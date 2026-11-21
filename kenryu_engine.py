import os, json, zipfile, asyncio, logging, re, io, base64, datetime, hmac, threading, time, subprocess
from contextlib import asynccontextmanager
from pathlib import Path
from typing import List, Literal, Optional
from collections import Counter, deque
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib_venn import venn2, venn3
import httpx
from fastapi import FastAPI, File, HTTPException, Request, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, JSONResponse
from pydantic import BaseModel, Field
import uvicorn

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

class AnalysisRequest(BaseModel):
    mirnas: List[str] = Field(min_length=1, max_length=50)
    years: int = Field(default=10, ge=1, le=50)
    mode: Literal["strict", "n-1", "n-2"] = "strict"
    month: Optional[str] = None
    offline: bool = False
    species: str = "human"
    # Ronda 3 (v1.48): interruptor opcional, apagado por defecto. Cuando
    # está activo, el consenso solo conserva genes con evidencia
    # experimental en miRTarBase (fuerte o débil) para al menos un miRNA
    # del panel. No cambia el método por defecto (el consenso sigue
    # siendo TargetScan + respaldo miRTarBase); es un filtro opcional.
    require_experimental: bool = False

BASE_DIR = Path(__file__).resolve().parent

# Unreleased: hash corto del commit de GitHub que produjo este servicio.
# Lo escribe el paso de sincronización en local_db/build_info.json con
# `git rev-parse --short HEAD` (el hash de GitHub, nunca el del Space).
# Si el archivo no existe o está mal formado, load_commit_hash()
# devuelve "unknown": nunca se inventa un valor.
BUILD_INFO_PATH = BASE_DIR / "local_db" / "build_info.json"
_COMMIT_RE = re.compile(r"[0-9a-f]{7,40}")


def _running_in_space():
    """True when this process runs in a Hugging Face Space.

    Primary signal: the SPACE_ID environment variable, which Hugging
    Face sets at runtime for every Space (docs: Spaces Overview,
    "Built-in environment variables"). Fallback: a checkout whose
    git origin is huggingface.co — that checkout tracks the Space
    repo, so its HEAD is the Space's commit, not GitHub's. A fork
    on another service (GitLab, ...) matches neither, so its HEAD
    is the code that runs.
    """
    if os.environ.get("SPACE_ID"):
        return True
    try:
        remote = subprocess.run(["git", "config", "--get", "remote.origin.url"],
                                cwd=BASE_DIR, capture_output=True, text=True,
                                timeout=5)
        return "huggingface.co" in remote.stdout
    except Exception:
        return False


def _git_head():
    """Short HEAD of this checkout, or None when there is no usable git."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             cwd=BASE_DIR, capture_output=True, text=True,
                             timeout=5)
        if out.returncode == 0:
            commit = out.stdout.strip()
            if _COMMIT_RE.fullmatch(commit):
                return commit
    except Exception:
        pass
    return None


def load_commit_hash():
    """Short GitHub commit hash of the code running here, or "unknown".

    Rule: in a Space (SPACE_ID set, or a Hugging Face origin),
    build_info.json is the only record of the GitHub commit — the
    Space's own HEAD is a different history. Everywhere else the
    git HEAD is the code that runs and always wins, even against a
    same-version build_info.json. When git is unavailable (a Docker
    image without .git), the version-matched build_info.json is the
    fallback.
    """
    if not _running_in_space():
        head = _git_head()
        if head:
            return head
    try:
        data = json.loads(BUILD_INFO_PATH.read_text(encoding="utf-8"))
        if isinstance(data, dict) and data.get("version") == APP_VERSION:
            commit = data.get("commit")
            if isinstance(commit, str) and _COMMIT_RE.fullmatch(commit):
                return commit
    except Exception:
        pass
    return "unknown"
# Directorio de datos: contiene targetscan_full.txt.zip y los archivos
# hsa-miR-*.txt de ejemplo. Compatible hacia atrás: si data/ no existe
# se cae a BASE_DIR (estructura antigua).
DATA_DIR = BASE_DIR / "data" if (BASE_DIR / "data").exists() else BASE_DIR
RAW_DATA_DIR = DATA_DIR
TARGETSCAN_DB = {}
# Ratón (TargetScanMouse 8.0, Species ID 10090). Las capas de validación
# y enriquecimiento locales son humanas; para ratón el análisis corre
# solo con TargetScan local + el respaldo vivo de Harmonizome.
MOUSE_DB = {}
MOUSE_ZIP = "data/mouse/Predicted_Targets_Info.default_predictions.txt.zip"


def load_mouse_data():
    """Carga la base local de ratón (una sola vez)."""
    global MOUSE_DB
    if not MOUSE_DB:
        path = BASE_DIR / MOUSE_ZIP
        if path.exists():
            with zipfile.ZipFile(path) as zf:
                internal = [f for f in zf.namelist() if f.endswith(".txt")][0]
                with zf.open(internal) as fh:
                    MOUSE_DB = _parse_targetscan_txt(fh, species_id="10090")
            logger.info(f"🐭 TargetScanMouse local: {len(MOUSE_DB)} familias.")
        else:
            logger.warning(f"⚠️ Sin base de ratón ({path} no existe).")
    return MOUSE_DB
# miRTarBase validado (local: índice OFICIAL 9.0 si existe, si no el GMT de
# Harmonizome) y librerías GMT de enriquecimiento (Reactome/WikiPathways/GO,
# de Enrichr). KEGG NO está en el repo: Kanehisa Labs prohíbe redistribuir
# KEGG sin acuerdo (ver DATA_LICENSES.md); el actualizador lo descarga en
# ejecución si el usuario lo pide, y el enriquecimiento en vivo lo consulta
# vía Enrichr. Se cargan bajo demanda y se cachean como TARGETSCAN_DB.
MIRTARBASE_DB = {}
MIRTARBASE_FORMAT = None  # "official" | "gmt" | None
# Build metadata of the loaded official index (rows/pairs/miRNAs and
# whether it is the synthetic CI fixture); None when no official index
# is loaded. Drives the honest version string in loaded_data_versions().
MIRTARBASE_META = None
GMT_LIBRARIES = {}
GMT_FILES = {
    "KEGG": "local_db/gmt_KEGG_2021_Human.gmt",
    "Reactome": "local_db/gmt_Reactome_Pathways_2024.gmt",
    "WikiPathways": "local_db/gmt_WikiPathways_2019_Human.gmt",
    "GO": "local_db/gmt_GO_Biological_Process_2023.gmt",
}
MIRTARBASE_GMT = "local_db/mirtarbase_gene_set_library.gmt.gz"
MIRTARBASE_OFFICIAL = "local_db/mirtarbase_official.json.gz"


def load_mirtarbase_local():
    """Carga la librería local de miRTarBase (una sola vez). Prefiere el
    índice oficial 9.0 (con niveles fuerte/débil); si no existe, usa el
    GMT de Harmonizome (nivel "legacy", sin distinción de evidencia),
    descargado bajo demanda y no incluido en el repo desde v1.52."""
    global MIRTARBASE_DB, MIRTARBASE_FORMAT, MIRTARBASE_META
    if not MIRTARBASE_DB:
        official = BASE_DIR / MIRTARBASE_OFFICIAL
        if official.exists():
            MIRTARBASE_DB, MIRTARBASE_META = load_mirtarbase_official_with_meta(official)
            MIRTARBASE_FORMAT = "official"
            logger.info(f"📚 miRTarBase local: índice oficial 9.0, {len(MIRTARBASE_DB)} miRNAs.")
        else:
            path = BASE_DIR / MIRTARBASE_GMT
            if path.exists():
                MIRTARBASE_DB = load_mirtarbase_gmt(path)
                MIRTARBASE_FORMAT = "gmt"
                MIRTARBASE_META = None
                logger.info(f"📚 miRTarBase local: {len(MIRTARBASE_DB)} miRNAs validados (GMT Harmonizome).")
            else:
                MIRTARBASE_META = None
                logger.warning(f"⚠️ Sin miRTarBase local ({official} ni {path} existen).")
        # v1.59: las claves de gen se normalizan con el mapa HGNC al
        # cargar, igual que los conjuntos de TargetScan.
        MIRTARBASE_DB = _apply_symbol_map_to_mirtarbase(MIRTARBASE_DB)
    return MIRTARBASE_DB


async def _download_gmt_bytes(url):
    """GET the GMT file; returns raw bytes. Raises on any failure."""
    async with make_async_client(timeout=30.0) as client:
        r = await client.get(url)
        r.raise_for_status()
        return r.content


async def _ensure_legacy_gmt():
    """Best-effort startup download of the Harmonizome GMT (v1.53).

    The CC BY-NC-SA snapshot is not distributed in the repo (v1.52) nor
    in the deployed image: when neither the official miRTarBase index nor
    the GMT is on disk, fetch the GMT from Harmonizome so the legacy
    fallback still works. Returns True when a local miRTarBase source is
    present afterwards. Never raises: a failed download only logs a
    warning and the service starts without the fallback. Disabled with
    KENRYU_STARTUP_GMT_DOWNLOAD=0 (the test suite disables it)."""
    import gzip

    official = BASE_DIR / MIRTARBASE_OFFICIAL
    gmt = BASE_DIR / MIRTARBASE_GMT
    if official.exists() or gmt.exists():
        return True
    if os.environ.get("KENRYU_STARTUP_GMT_DOWNLOAD", "1") == "0":
        return False
    url = next((f["url"] for s in DATA_SOURCES if s.get("key") == "mirtarbase"
                for f in s.get("files", [])), None)
    if not url:
        logger.warning("⚠️ Sin URL de descarga para el GMT de Harmonizome.")
        return False
    try:
        data = await _download_gmt_bytes(url)
        if data[:2] != b"\x1f\x8b":
            raise ValueError("la descarga no es un gzip válido")
        gzip.decompress(data)  # valida la integridad antes de instalar
        tmp = gmt.with_suffix(".tmp")
        tmp.write_bytes(data)
        os.replace(tmp, gmt)
        logger.info(f"📥 GMT de Harmonizome descargado al arrancar ({len(data)} B).")
        return True
    except Exception as e:  # noqa: BLE001 - best-effort, nunca frena el arranque
        logger.warning(f"⚠️ No se pudo descargar el GMT de Harmonizome al arrancar: {e}")
        return False


def load_gmt_libraries():
    """Carga las librerías GMT locales de enriquecimiento (una sola vez)."""
    global GMT_LIBRARIES
    if not GMT_LIBRARIES:
        for label, rel in GMT_FILES.items():
            path = BASE_DIR / rel
            if path.exists():
                GMT_LIBRARIES[label] = load_gmt_library(path)
        logger.info(f"📚 GMT locales: {', '.join(f'{k} ({len(v)} términos)' for k, v in GMT_LIBRARIES.items()) or 'ninguna'}.")
    return GMT_LIBRARIES
pubmed_semaphore = asyncio.Semaphore(3)
PERSISTENT_CACHE = {"trans": {}, "pubmed": {}, "enrichr": {}, "gene_research": {}}
CACHE_FILE = Path("local_db/analysis_cache.json")
# H-15 (Ronda 3, v1.49): la caché en disco debe ser atómica (escribir a
# un temporal y renombrar), estar acotada (las secciones viejas crecen
# sin límite) y tener frescura (la investigación multi-fuente caduca;
# las traducciones y los metadatos de PubMed son hechos estables y no
# caducan). El candado es de proceso: el servidor es un solo proceso y
# las escrituras concurrentes vendrían de sus propias tareas.
_CACHE_LOCK = threading.Lock()
_CACHE_LIMITS = {"trans": 2000, "pubmed": 1000, "enrichr": 200, "gene_research": 500}
CACHE_TTL_DAYS = 90

def load_persistent_cache():
    global PERSISTENT_CACHE
    if CACHE_FILE.exists():
        try:
            with open(CACHE_FILE, "r", encoding="utf-8") as f:
                data = json.load(f)
            # Frescura: la investigación multi-fuente caduca a los
            # CACHE_TTL_DAYS; sin _saved_at (archivos viejos) se trata
            # como fresca para no borrar datos sin avisar.
            saved_at = data.get("_saved_at")
            expired = False
            if saved_at:
                try:
                    saved = datetime.datetime.fromisoformat(saved_at)
                    if saved.tzinfo is None:
                        saved = saved.replace(tzinfo=datetime.timezone.utc)
                    expired = (datetime.datetime.now(datetime.timezone.utc) - saved).days >= CACHE_TTL_DAYS
                except ValueError:
                    expired = False
            if expired:
                data.pop("gene_research", None)
                logger.info("💾 Caché: investigación multi-fuente caducada (>90 días), descartada.")
            # Fusionar de forma defensiva — versiones antiguas del cache
            # no tienen 'gene_research', no debe fallar la carga.
            for key in ("trans", "pubmed", "enrichr", "gene_research"):
                if key in data:
                    PERSISTENT_CACHE[key].update(data[key])
            logger.info(
                f"💾 Caché persistente cargada "
                f"({len(PERSISTENT_CACHE['trans'])} traducciones, "
                f"{len(PERSISTENT_CACHE['pubmed'])} PubMed, "
                f"{len(PERSISTENT_CACHE['gene_research'])} genes investigados)."
            )
        except Exception as e:
            logger.warning(f"No se pudo cargar la caché: {e}")

def save_persistent_cache():
    """Guarda la caché de forma atómica y acotada: escribe a un archivo
    temporal en el mismo directorio y lo renombra (os.replace), bajo un
    candado de proceso. Cada sección se recorta a _CACHE_LIMITS
    (sobreviven las entradas más recientes). Un fallo a mitad de
    escritura deja intacto el archivo anterior."""
    try:
        with _CACHE_LOCK:
            CACHE_FILE.parent.mkdir(exist_ok=True)
            trimmed = {}
            for key, limit in _CACHE_LIMITS.items():
                section = PERSISTENT_CACHE.get(key, {})
                if len(section) > limit:
                    section = dict(list(section.items())[-limit:])
                trimmed[key] = section
            payload = dict(trimmed)
            payload["_saved_at"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            tmp = CACHE_FILE.with_suffix(CACHE_FILE.suffix + ".tmp")
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(payload, f, ensure_ascii=False, indent=2)
            os.replace(tmp, CACHE_FILE)
    except Exception as e:
        logger.warning(f"No se pudo guardar la caché: {e}")

async def translate_to_spanish(text: str, client: httpx.AsyncClient) -> str:
    """Traduce texto usando la API MyMemory con persistencia en disco."""
    if not text or len(text) < 3: return text
    if text in PERSISTENT_CACHE["trans"]: return PERSISTENT_CACHE["trans"][text]
    
    try:
        url = "https://api.mymemory.translated.net/get"
        params = {"q": text, "langpair": "en|es"}
        r = await client.get(url, params=params, timeout=8.0)
        if r.status_code == 200:
            trans = r.json().get("responseData", {}).get("translatedText", text)
            PERSISTENT_CACHE["trans"][text] = trans
            # Guardar periódicamente o al final
            return trans
    except Exception as e:
        logger.warning(f"Error traducción: {e}")
    return text

# ── LIMPIEZA DE NOMBRES ─────────────────────────────────────────────────────
# ── PURE HELPERS (live in kenryu_core.py, re-exported here) ─────────────
from kenryu_core import (
    CONVERGENCE_LABEL,
    CONVERGENCE_LIMIT,
    build_export,
    build_synthesis,
    build_validation_report,
    clean_mirna_name_cloud,
    consensus_genes,
    create_heatmap,
    create_visuals,
    deep_clean_mirna,
    get_validated_evidence,
    get_validated_local,
    load_gmt_library,
    load_mirtarbase_gmt,
    load_mirtarbase_official,
    load_mirtarbase_official_with_meta,
    local_enrichment,
    mirna_seed,
    normalize_isomir,
    panel_convergence,
    parse_family_seeds,
    parse_mirbase_sequences,
    prioritize_genes,
    reverse_lookup,
    term_specificity,
    verify_arm_seed,
    group_redundant_terms,
)
from kenryu_updater import (
    DATA_SOURCES,
    check_updates,
    install_source_update,
    load_local_versions,
    make_async_client,
)

def load_local_data(mapped=True):
    """TargetScan DB (familia -> lista de genes). Con mapped=True (por
    defecto) los símbolos se normalizan con el mapa HGNC al cargar; con
    mapped=False se devuelve el contenido crudo (lo usan la auditoría
    de símbolos y el constructor del mapa, que miden sobre datos sin
    mapear). La caché .pkl guarda siempre el contenido crudo."""
    cache_path = Path("local_db/targetscan_db.pkl")
    if cache_path.exists():
        logger.info("⚡ Cargando caché .pkl...")
        try:
            raw = pd.read_pickle(cache_path)
            return (_apply_symbol_map_to_targetscan(raw) if mapped
                    else raw)
        except:
            pass

    candidates = [
        BASE_DIR / "Predicted_Targets_Info.default_predictions.txt",
        DATA_DIR / "targetscan_full.txt.zip",
        BASE_DIR / "targetscan_full.txt.zip",
        # H-19 (v1.49): the zip was renamed from targetscan_full.json.zip
        # (it holds a .txt, not JSON). The old name stays as a fallback
        # so the Hugging Face Space keeps working until its next sync.
        DATA_DIR / "targetscan_full.json.zip",
        BASE_DIR / "targetscan_full.json.zip",
        RAW_DATA_DIR / "Predicted_Targets_Info.default_predictions.txt",
    ]

    for src in candidates:
        if not src.exists():
            continue
        logger.info(f"📄 Procesando: {src.name}")
        try:
            if src.suffix == '.zip':
                with zipfile.ZipFile(src, 'r') as z:
                    internal = [f for f in z.namelist() if f.endswith('.txt') or f.endswith('.json')][0]
                    with z.open(internal) as f:
                        if internal.endswith('.json'):
                            data = json.load(f)
                            temp_db = {k.lower(): set(v) for k, v in data.items()}
                        else:
                            temp_db = _parse_targetscan_txt(f)
            else:
                with open(src, 'rb') as f:
                    temp_db = _parse_targetscan_txt(f)

            final_db = {k: list(v) for k, v in temp_db.items()}
            Path("local_db").mkdir(exist_ok=True)
            pd.to_pickle(final_db, cache_path)
            logger.info(f"✅ {len(final_db)} familias de miRNA cargadas.")
            return (_apply_symbol_map_to_targetscan(final_db) if mapped
                    else final_db)
        except Exception as e:
            logger.error(f"❌ Error con {src.name}: {e}")

    # Último recurso: archivos individuales por miRNA en repo_tesis/data/raw/
    if RAW_DATA_DIR.exists():
        temp_db = {}
        for txt in RAW_DATA_DIR.glob("hsa-*.txt"):
            key = txt.stem.lower()
            with open(txt, 'r', encoding='utf-8-sig') as f:
                genes = {l.strip() for l in f if l.strip()}
            if genes:
                temp_db[key] = list(genes)
        if temp_db:
            logger.info(f"📁 DB desde archivos individuales: {len(temp_db)} miRNAs.")
            return temp_db

    logger.error("❌ No se encontró ninguna fuente de datos.")
    return {}

def _parse_targetscan_txt(fileobj, species_id="9606"):
    """
    Parsea el archivo global TargetScan (multi-especie).

    CORRECCIONES vs versión original:
    1. Filtra SOLO la Species ID pedida (9606 humano, 10090 ratón) antes
       de indexar.
    2. Indexa por familia completa (clave exacta, lowercase).
    3. Indexa también por sub-partes del nombre compuesto separadas por '/'
       para que variantes como 'mir-106b-5p' sean encontrables dentro de
       'mir-17-5p/20-5p/93-5p/106-5p/519-3p'.
    4. NO mezcla genes de otras especies en los conjuntos humanos.
    """
    temp_db = {}
    chunks = pd.read_csv(
        fileobj, sep='\t',
        usecols=['miR Family', 'Gene Symbol', 'Species ID'],
        dtype=str, chunksize=400_000
    )
    for chunk in chunks:
        # ── BUG 1 CORREGIDO: filtrar solo la especie pedida ANTES de indexar ──
        rows = chunk[chunk['Species ID'].str.strip() == species_id].copy()
        if rows.empty:
            continue

        for mir, gene in zip(rows['miR Family'], rows['Gene Symbol']):
            gene = str(gene).strip()
            if not gene or gene == 'nan':
                continue
            fam = str(mir).strip().lower()

            # Índice por familia completa exacta
            if fam not in temp_db:
                temp_db[fam] = set()
            temp_db[fam].add(gene)

            # Índice por cada sub-nombre dentro de la familia compuesta
            # Ej: "mir-17-5p/20-5p/93-5p/106-5p/519-3p"
            #     → ["mir-17-5p", "mir-20-5p", "mir-93-5p", "mir-106-5p", "mir-519-3p"]
            # ── BUG 2 CORREGIDO: reconstruir 'mir-' correctamente ──────────
            parts = fam.split('/')
            base_prefix = None
            for idx, part in enumerate(parts):
                part = part.strip()
                if idx == 0:
                    # Primera parte ya tiene "mir-" completo
                    key = part
                    # Extrae prefijo: "mir-"
                    base_prefix = re.match(r'^(mir-)', part)
                    base_prefix = base_prefix.group(1) if base_prefix else 'mir-'
                else:
                    # Partes siguientes pueden ser solo el número: "20-5p", "93-5p"
                    # Necesitan el prefijo "mir-"
                    if not part.startswith('mir-'):
                        key = base_prefix + part
                    else:
                        key = part

                if key not in temp_db:
                    temp_db[key] = set()
                temp_db[key].add(gene)

    return temp_db

# ── VERIFICACIÓN DE BRAZOS (-5p/-3p) POR SEMILLA ───────────────────────────
# Audit 2026-11-21 (v1.56): 33 miRNAs con brazo se resolvían en silencio a
# una familia TargetScan sin brazo; en 14 horquillas ambos brazos
# devolvían el mismo conjunto. La verificación compara la semilla de la
# familia (miR_Family_Info.txt) con la secuencia madura del miRNA
# (miRBase 23.1 mature.fa): la semilla debe aparecer al inicio de la
# secuencia (offset 0-3, ver kenryu_core.verify_arm_seed). Sin los dos
# archivos la verificación no está disponible: se devuelve el resultado
# de antes con un aviso (fail-open ruidoso), y la respuesta de /analyze
# lo dice en "arm_verification". Los dos archivos son datos con versión
# (local_db/, el FASTA ya se distribuía; el de familias lo descarga el
# actualizador junto con TargetScan 8.0).
MIRBASE_MATURE_PATH = "local_db/mirbase_mature.fa"
FAMILY_INFO_PATH = "local_db/mir_family_info.txt"
FAMILY_INFO_ZIP = "data/mir_family_info.txt.zip"
_ARM_MAPS = None  # (mirbase_sequences, family_seeds) | (None, None)


def _load_arm_seed_maps():
    """Lee los dos archivos de semillas (una sola vez). Devuelve
    (secuencias miRBase, semillas de familia) o (None, None) si falta
    alguno: la verificación queda no disponible, nunca inventada."""
    global _ARM_MAPS
    if _ARM_MAPS is not None:
        return _ARM_MAPS
    seqs = None
    fasta = BASE_DIR / MIRBASE_MATURE_PATH
    if fasta.exists():
        try:
            seqs = parse_mirbase_sequences(fasta)
        except Exception as e:
            logger.warning(f"⚠️ No se pudo leer {fasta}: {e}")
    fam = None
    info = BASE_DIR / FAMILY_INFO_PATH
    if info.exists():
        try:
            fam, conflicts = parse_family_seeds(info)
            if conflicts:
                logger.warning(
                    f"⚠️ miR_Family_Info: {len(conflicts)} claves con dos "
                    f"semillas (se acepta cualquiera): {conflicts[:4]}")
        except Exception as e:
            logger.warning(f"⚠️ No se pudo leer {info}: {e}")
    elif (BASE_DIR / FAMILY_INFO_ZIP).exists():
        try:
            with zipfile.ZipFile(BASE_DIR / FAMILY_INFO_ZIP) as z:
                inner = [n for n in z.namelist() if n.endswith(".txt")][0]
                with z.open(inner) as fh:
                    text = fh.read().decode("utf-8")
            import io as _io
            fam, conflicts = parse_family_seeds(_io.StringIO(text))
        except Exception as e:
            logger.warning(f"⚠️ No se pudo leer {FAMILY_INFO_ZIP}: {e}")
    _ARM_MAPS = (seqs, fam)
    if seqs is None or fam is None:
        logger.warning(
            "⚠️ Verificación de brazos NO DISPONIBLE (falta "
            "local_db/mirbase_mature.fa o local_db/mir_family_info.txt): "
            "las dianas se devuelven como antes, sin comprobar el brazo.")
    return _ARM_MAPS


def arm_seed_maps():
    """Mapas de semillas para la verificación de brazos (o (None, None)
    si los archivos no están). Pensado para tests y para la auditoría."""
    return _load_arm_seed_maps()


def arm_verification_status() -> str:
    """'active' cuando los dos archivos de semillas están cargados,
    'unavailable' cuando falta alguno."""
    seqs, fam = _load_arm_seed_maps()
    return "active" if (seqs is not None and fam is not None) else "unavailable"


# ── NORMALIZACIÓN DE SÍMBOLOS DE GEN (HGNC, v1.59) ──────────────────────────
# 2.496 de 13.035 genes de TargetScan no aparecen en ninguna biblioteca
# local de enriquecimiento; HGNC recupera 527 con un símbolo actual
# único en las bibliotecas (+11 solo vía KEGG). El mapa
# (local_db/hgnc_symbol_map.json, 538 entradas, CC0) se aplica al
# cargar: a los conjuntos de TargetScan y a las claves del índice
# miRTarBase (y al fetch vivo de Harmonizome), para que las dos capas
# hablen con símbolos actuales. Los genes ambiguos nunca se mapean.
# Sin el archivo, los símbolos pasan tal cual, con un aviso (fail-open
# ruidoso). El mapa es un dato con versión; se reconstruye con
# scripts/build_hgnc_symbol_map.py.
HGNC_MAP_PATH = "local_db/hgnc_symbol_map.json"
_HGNC_MAP = None  # dict old -> current | {} (archivo ausente)


def _load_hgnc_map():
    """Lee el mapa de símbolos HGNC (una sola vez). Devuelve {} si el
    archivo falta o no se puede leer: los símbolos pasan sin mapear."""
    global _HGNC_MAP
    if _HGNC_MAP is not None:
        return _HGNC_MAP
    path = BASE_DIR / HGNC_MAP_PATH
    _HGNC_MAP = {}
    if path.exists():
        try:
            doc = json.loads(path.read_text(encoding="utf-8"))
            _HGNC_MAP = dict(doc.get("mapping") or {})
        except Exception as e:
            logger.warning(f"⚠️ No se pudo leer {path}: {e}")
    if not _HGNC_MAP:
        logger.warning(
            "⚠️ Mapa de símbolos HGNC no disponible "
            f"({HGNC_MAP_PATH}): los símbolos antiguos no se normalizan.")
    return _HGNC_MAP


def map_gene_symbol(symbol: str) -> str:
    """Símbolo actual HGNC para `symbol`, o el mismo símbolo si no
    está en el mapa (ambiguo, sin mapear, o ya actual)."""
    return _load_hgnc_map().get(symbol, symbol)


def _apply_symbol_map_to_targetscan(db):
    """Un conjunto de dianas por familia con los símbolos normalizados
    (sin duplicados; los consumidores usan sets)."""
    return {k: list(dict.fromkeys(map_gene_symbol(g) for g in v))
            for k, v in db.items()}


_LEVEL_RANK = {"strong": 3, "weak": 2, "legacy": 1, "none": 0}


def _apply_symbol_map_to_mirtarbase(db):
    """Claves de gen del índice validado con los símbolos normalizados.
    Si dos símbolos antiguos caen en el mismo actual dentro de un
    miRNA, se fusiona la evidencia quedándose el mejor nivel (con los
    PMIDs y experimentos unidos)."""
    out = {}
    for mirna, genes in db.items():
        if isinstance(genes, dict):
            merged = {}
            for gene, ev in genes.items():
                new = map_gene_symbol(gene)
                if new not in merged:
                    merged[new] = dict(ev)
                else:
                    cur = merged[new]
                    if (_LEVEL_RANK.get(ev.get("level"), 0)
                            > _LEVEL_RANK.get(cur.get("level"), 0)):
                        cur["level"] = ev["level"]
                    cur["pmids"] = sorted(
                        set(cur.get("pmids") or [])
                        | set(ev.get("pmids") or []))
                    cur["experiments"] = sorted(
                        set(cur.get("experiments") or [])
                        | set(ev.get("experiments") or []))
            out[mirna] = merged
        else:
            out[mirna] = {map_gene_symbol(g) for g in genes}
    return out


def check_arm_resolution(mirna: str, db=None) -> dict:
    """¿La familia TargetScan que responde a `mirna` modela de verdad
    el brazo pedido? Devuelve {status, family_key, reason}:
    - verified: la semilla de la familia está al inicio de la secuencia
      madura del miRNA (offset 0-3).
    - unverified: la semilla no está (el brazo pedido NO es el
      modelado) o no se pudo comprobar (sin secuencia o sin semilla).
    - not_applicable: el nombre no trae brazo (resuelve al brazo
      modelado, la hebra guía histórica), el miRNA es de ratón, o no
      hay datos locales.
    - unavailable: faltan los archivos de semillas."""
    m_clean = deep_clean_mirna(mirna)
    if m_clean.startswith("mmu-"):
        return {"status": "not_applicable", "family_key": None,
                "reason": "mouse: no arm data"}
    arm = re.search(r"-(5p|3p)$", m_clean)
    if not arm:
        return {"status": "not_applicable", "family_key": None,
                "reason": "bare name resolves to the modelled arm"}
    genes, key = _resolve_targetscan(mirna, db)
    if not genes or key is None:
        return {"status": "not_applicable", "family_key": key,
                "reason": "no local targets to verify"}
    if key == "__individual_file__":
        return {"status": "verified", "family_key": key,
                "reason": "individual TargetScan file named for this miRNA"}
    seqs, fam = _load_arm_seed_maps()
    if seqs is None or fam is None:
        return {"status": "unavailable", "family_key": key,
                "reason": "seed files missing"}
    query_seq = seqs.get(m_clean)
    if not query_seq:
        return {"status": "unverified", "family_key": key,
                "reason": f"{mirna} has no sequence in miRBase 23.1"}
    family_seeds = fam.get(key)
    if not family_seeds:
        return {"status": "unverified", "family_key": key,
                "reason": f"family '{key}' not in miR_Family_Info"}
    if verify_arm_seed(query_seq, family_seeds):
        return {"status": "verified", "family_key": key,
                "reason": f"family seed {sorted(family_seeds)} in sequence"}
    return {"status": "unverified", "family_key": key,
            "reason": (f"family seed {sorted(family_seeds)} not in the "
                       f"{arm.group(1)} sequence (offset 0-3)")}


def _resolve_targetscan(mirna: str, db=None):
    """Resuelve un miRNA a sus genes diana locales y a la clave de
    familia TargetScan que respondió. Devuelve (genes, family_key);
    family_key es "__individual_file__" cuando los genes vienen de un
    archivo individual (el archivo lleva el nombre del miRNA, así que
    el brazo está en el nombre del archivo) y None cuando no hay datos."""
    m_clean = deep_clean_mirna(mirna)            # ej. "mir-33a-5p"
    m_no_letter = re.sub(r'(\d+)[a-z]+(-\d+)', r'\1\2', m_clean)  # "mir-33-5p"
    m_no_suffix = re.sub(r'-(5p|3p)$', '', m_clean)         # "mir-33a"
    m_no_letter_no_suffix = re.sub(r'-(5p|3p)$', '', m_no_letter)  # "mir-33"

    # 1. Archivos individuales (repo_tesis/data/raw/)  ← MÁXIMA PRIORIDAD
    #    Estos son los archivos descargados directamente de TargetScan,
    #    ya filtrados y limpios. Si existen, usarlos siempre.
    if RAW_DATA_DIR.exists():
        for variant in [mirna, f"hsa-{m_clean}", m_clean, m_no_letter, m_no_suffix]:
            fp = RAW_DATA_DIR / f"{variant}.txt"
            if fp.exists():
                with open(fp, 'r', encoding='utf-8-sig') as f:
                    genes = {l.strip() for l in f if l.strip()}
                logger.info(f"📂 {mirna}: {len(genes)} genes (archivo individual)")
                return genes, "__individual_file__"

    # 2. DB global — variantes de nombre exactas
    db = TARGETSCAN_DB if db is None else db
    for key in [m_clean, m_no_letter, m_no_suffix, m_no_letter_no_suffix]:
        res = db.get(key)
        if res:
            logger.info(f"🗄 {mirna} → '{key}': {len(res)} genes")
            return set(res), key

    # 3. BUG 3 CORREGIDO: búsqueda en familias compuestas
    #    Busca el nombre exacto del miRNA dentro de las claves de la DB
    #    (que ya están indexadas por sub-parte en _parse_targetscan_txt)
    #    Esto cubre el caso de hsa-miR-106b-5p → familia compuesta.
    for key_variant in [m_clean, m_no_letter]:
        for db_key, db_genes in db.items():
            # Match exacto de sub-parte: la clave en DB debe SER el miRNA
            if db_key == key_variant:
                logger.info(f"🔍 {mirna} → match exacto sub-familia '{db_key}': {len(db_genes)} genes")
                return set(db_genes), db_key

    # 4. Último recurso: búsqueda por número base SOLO si es suficientemente
    #    específico (incluye el sufijo 5p/3p para evitar falsos positivos).
    #    CORREGIDO (auditoría externa, 2026-11-19): la subcadena simple
    #    `specific_num in k` asignaba dianas de familias ajenas — p. ej.
    #    "1-5p" aparece dentro de "mir-211-5p" y devolvía 793 genes de
    #    miR-204/211 para un hsa-miR-1-5p inexistente. Ahora el número solo
    #    vale en un LÍMITE DE TOKEN (tras mir-/let-/ o /, y antes de / o .
    #    o fin de clave): así "873-5p" sigue encontrando "mir-873-5p.2"
    #    (variante real sin clave exacta), pero "1-5p" ya no casa con
    #    "mir-211-5p". Un miRNA sin familia real devuelve vacío con aviso.
    num_match = re.search(r'(\d+[ab]?(?:-[35]p)?)', m_clean)
    if num_match:
        specific_num = num_match.group(1)
        # Solo hacer difuso si el número incluye sufijo de strand (5p/3p)
        if '5p' in specific_num or '3p' in specific_num:
            token_pat = re.compile(
                rf'(?:^|mir-|let-|/){re.escape(specific_num)}(?:$|[/.])')
            for k, v in db.items():
                if token_pat.search(k):
                    logger.warning(f"⚠️ {mirna} → difuso (último recurso) '{k}': {len(v)} genes")
                    return set(v), k

    logger.warning(f"⚠️ {mirna}: sin datos en DB local")
    return set(), None


# ── BÚSQUEDA DE GENES POR miRNA ───────────────────────────────────────────────
def get_targets_local(mirna: str, db=None) -> set:
    """
    Busca los genes diana de un miRNA.

    CORRECCIONES vs versión original:
    - BUG 3 CORREGIDO: la búsqueda "difusa por número" devolvía el primer
      match que contenía el número, sin importar si era la familia correcta.
      Ej: buscar "106b" matcheaba "mir-106-5p/17-5p/..." que es la familia
      INCORRECTA, devolviendo miles de genes mezclados.
      Ahora la búsqueda difusa es más específica y solo se usa como último
      recurso con log de advertencia.
    - BRAZO CORREGIDO (v1.58): si el nombre trae brazo (-5p/-3p) y la
      familia que responde no lo indica, la semilla de la familia debe
      aparecer al inicio de la secuencia madura del miRNA; si no, el
      brazo pedido NO es el modelado y se devuelve vacío con aviso
      (check_arm_resolution). Un nombre sin brazo no se comprueba.
    """
    genes, key = _resolve_targetscan(mirna, db)
    if not genes or key is None:
        return genes
    check = check_arm_resolution(mirna, db)
    if check["status"] == "unverified":
        logger.warning(
            f"⚠️ {mirna}: brazo no verificado (familia '{key}'): "
            f"{check['reason']} — devolviendo vacío, no las dianas del otro brazo")
        return set()
    return genes


def _validated_fallback_arm_ok(mirna: str, fallback_key: str) -> bool:
    """The miRTarBase lookup fell back from an armed name to the bare
    name (pre-23 naming: the bare key is the guide arm). Verify by seed,
    as for TargetScan families: the query's mature sequence must carry
    the bare key's family seed at offset 0-3. When the key is not a
    TargetScan family at all (268 bare keys in the v9.0 index — viral,
    removed or renamed miRNAs) the seed cannot be checked: the validated
    pairs are kept, with a warning, because experimental evidence is not
    dropped for an unanswerable question. Same when the seed files are
    missing: keep, loudly."""
    seqs, fam = _load_arm_seed_maps()
    if seqs is None or fam is None:
        logger.warning(
            f"⚠️ {mirna}: no se pudo comprobar el brazo del respaldo "
            f"validado a {fallback_key} (faltan los archivos de semillas)")
        return True
    qseq = seqs.get(deep_clean_mirna(mirna))
    fseeds = fam.get(deep_clean_mirna(fallback_key))
    if not qseq or not fseeds:
        logger.warning(
            f"⚠️ {mirna}: el respaldo validado a {fallback_key} no se "
            f"puede comprobar por semilla (sin secuencia o sin familia) — "
            f"se conserva")
        return True
    ok = verify_arm_seed(qseq, fseeds)
    if not ok:
        logger.warning(
            f"⚠️ {mirna}: el respaldo validado a {fallback_key} es del "
            f"otro brazo (semilla {sorted(fseeds)} ausente) — se descarta")
    return ok


def get_validated_evidence_checked(mirna, db, species="hsa"):
    """get_validated_evidence + arm guard on the bare-name fallback:
    an armed query whose exact name is not in the library falls back to
    the bare name; if that bare key is the other arm's, the evidence is
    dropped (it belongs to the other arm). An exact-name hit is the
    miRNA's own evidence and is returned as-is."""
    name = clean_mirna_name_cloud(mirna, species)
    if name in db:
        return get_validated_evidence(mirna, db, species)
    ev = get_validated_evidence(mirna, db, species)
    if not ev:
        return {}
    base = re.sub(r"-(5p|3p)$", "", name, flags=re.IGNORECASE)
    return ev if _validated_fallback_arm_ok(mirna, base) else {}


def get_validated_local_checked(mirna, db, species="hsa"):
    """get_validated_local + the same arm guard (only levels
    strong/weak/legacy count)."""
    return {g for g, ev in get_validated_evidence_checked(mirna, db, species).items()
            if ev.get("level") in ("strong", "weak", "legacy")}

def family_labels_for_targets(mirnas, gene_sets, db):
    """TargetScan family key per miRNA (same order as gene_sets), for the
    convergence seed and for excluding observed families from the null
    pools. The label is canonical PER SET: miRNAs with identical target
    sets (miR-33a and miR-33b) share one label, so the seed depends on
    the families, not on which name came first. A set that matches no
    family key (e.g. individual raw files) falls back to the first
    sorted cleaned miRNA name sharing that set — the label only needs
    to be deterministic, not a true family name."""
    groups = {}
    for m, genes in zip(mirnas, gene_sets):
        groups.setdefault(frozenset(genes), []).append(m)
    canonical = {}
    for target, names in groups.items():
        match = sorted(k for k, v in db.items() if v and set(v) == set(target))
        canonical[target] = (match[0] if match
                             else sorted(deep_clean_mirna(n) for n in names)[0])
    return [canonical[frozenset(g)] for g in gene_sets]

async def fetch_mirtarbase(mirna: str, species: str = "hsa") -> set:
    m_name = clean_mirna_name_cloud(mirna, species)
    variants = [m_name, re.sub(r'-(5p|3p)$', '', m_name, flags=re.IGNORECASE)]
    try:
        async with make_async_client(timeout=12.0) as client:
            for v in variants:
                url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{v}/miRTarBase+microRNA+Targets"
                r = await client.get(url)
                if r.status_code == 200:
                    genes = {a['gene']['symbol'] for a in r.json().get('associations', [])}
                    if genes:
                        if v != variants[0] and not _validated_fallback_arm_ok(mirna, v):
                            continue
                        logger.info(f"🌐 miRTarBase {mirna}: {len(genes)} genes")
                        return {map_gene_symbol(g) for g in genes}
    except:
        pass
    return set()

# ── UTILIDADES DE RED RESILIENTE ─────────────────────────────────────────────
async def safe_pubmed_request(client: httpx.AsyncClient, url: str, params: dict, retries: int = 3):
    """Realiza peticiones a NCBI con reintentos y manejo de error 429."""
    for attempt in range(retries):
        try:
            async with pubmed_semaphore:
                r = await client.get(url, params=params, timeout=15.0)
                if r.status_code == 200:
                    return r.json()
                elif r.status_code == 429:
                    wait = (attempt + 1) * 2
                    logger.warning(f"⚠️ PubMed 429 (Too Many Requests). Reintentando en {wait}s...")
                    await asyncio.sleep(wait)
                else:
                    logger.error(f"❌ PubMed Error {r.status_code}: {r.text[:100]}")
        except Exception as e:
            if attempt == retries - 1: logger.error(f"❌ Fallo persistente en PubMed: {e}")
            await asyncio.sleep(1)
    return None

def annotate_enrichment_specificity(enrichment, gene_list):
    """Add a 'specificity' field to each enrichment term (v1.54): the
    exact probability (hypergeometric tail, unreleased 2026-11-19) that a
    random gene list of the same size, drawn from the TargetScan
    universe, overlaps the term at least as much as observed. 1.0 = the term appears by chance as often as observed
    (generic, e.g. "regulation of transcription"); 0.0 = it never
    appears by chance (specific). Terms not found in the local GMT
    libraries get None. Mutates the items in place and returns them."""
    if not enrichment or not gene_list:
        return enrichment
    ts_universe = (set().union(*TARGETSCAN_DB.values())
                   if TARGETSCAN_DB else set())
    if not ts_universe:
        return enrichment
    gene_set = set(gene_list)
    for item in enrichment:
        lib = GMT_LIBRARIES.get(item.get("Source"))
        # The live path translates Term to Spanish; the English name
        # survives in RawTerm (local path: RawTerm is the GMT key).
        term_name = item.get("RawTerm") or item.get("Term")
        term_genes = lib.get(term_name) if lib else None
        if not term_genes:
            item["specificity"] = None
            continue
        bg = ts_universe & set().union(*lib.values())
        observed = len(gene_set & term_genes)
        item["specificity"] = term_specificity(
            term_genes, bg, len(gene_list), observed)
    return enrichment


def annotate_enrichment_groups(enrichment, gene_list):
    """Mark redundant enrichment terms (unreleased, 2026-11-19). Terms
    whose panel hit sets overlap (overlap coefficient >= 0.75) collapse
    into the most significant one. Adds to each item `is_group_representative`
    (bool); representatives also get `group_size` and `group_members`
    (names); members get `group_of`. No p-value changes and every item
    stays in the response. Terms not found in the local GMT libraries stay
    alone. Mutates the items in place and returns them."""
    if not enrichment or not gene_list:
        return enrichment
    gene_set = set(gene_list)
    hits = {}
    for item in enrichment:
        lib = GMT_LIBRARIES.get(item.get("Source"))
        name = item.get("RawTerm") or item.get("Term")
        genes = lib.get(name) if lib else None
        if genes:
            hits[(item.get("Source"), name)] = gene_set & set(genes)
    for g in group_redundant_terms(enrichment, hits):
        rep = g["representative"]
        rep["is_group_representative"] = True
        rep["group_size"] = 1 + len(g["members"])
        rep["group_members"] = [m.get("Term") for m in g["members"]]
        for m in g["members"]:
            m["is_group_representative"] = False
            m["group_of"] = rep.get("Term")
    return enrichment


# ── ENRIQUECIMIENTO DIRECTO (BYPASS GSEApy) ──────────────────────────────────
async def compute_enrichment_routes(gene_list: list, offline: bool, client: httpx.AsyncClient):
    """Split the enrichment into its two routes (Unreleased, 2026-11-19).

    The LOCAL Fisher test with the TargetScan-universe background is the
    method; the live Enrichr route is an exploratory cross-check whose
    adjusted p-values are not a standard FDR at small list sizes (see
    docs/validation/calibration-enrichr-route.md: single-gene overlaps
    already reach adjusted p < 0.05, and the values are incompatible
    with a Benjamini-Hochberg correction over the whole library).
    Offline mode returns the local list only.
    Returns (local_terms, enrichr_terms)."""
    load_gmt_libraries()
    # Fondo del enriquecimiento local (v1.54): el universo de TargetScan
    # — todos los genes con algún sitio predicho (13.035 humanos), no la
    # unión del panel (v1.49, fondo dependiente del panel y sin potencia)
    # ni la unión de cada GMT (~20k genes que el análisis nunca probó,
    # que infla la significancia). La ruta viva de Enrichr usa el fondo
    # propio de Enrichr (su API no admite un fondo personalizado); el
    # local es el que controlamos.
    targetscan_universe = set().union(*TARGETSCAN_DB.values()) if TARGETSCAN_DB else None
    local_terms = local_enrichment(gene_list, GMT_LIBRARIES, background=targetscan_universe)
    if offline:
        return local_terms, []
    enrichr_terms = await enrich_genes_direct(gene_list, client)
    return local_terms, enrichr_terms


async def enrich_genes_direct(gene_list: list, client: httpx.AsyncClient):
    """Llamada directa a la API de Enrichr con caché persistente."""
    if len(gene_list) < 3: return []
    
    # Crear una clave única basada en los genes (ordenados). Versionada
    # (enr-v2): la v1 guardaba términos con el corte viejo de p-adj < 0,3
    # y se seguían sirviendo como significativos (v1.46, H-05).
    cache_key = "enr-v2|" + ",".join(sorted(gene_list[:300]))
    if cache_key in PERSISTENT_CACHE["enrichr"]:
        logger.info("⚡ Usando resultados de enriquecimiento desde caché.")
        return PERSISTENT_CACHE["enrichr"][cache_key]
    
    results = []
    try:
        # 1. Subir la lista de genes (un reintento: la primera llamada
        # en frío a veces agota el timeout de 20 s y la segunda responde
        # en ~1 s; verificado en vivo 2026-11-18).
        payload = {"list": (None, "\n".join(gene_list[:300])), "description": (None, "KENRYU Analysis")}
        r_add = None
        for attempt in range(2):
            try:
                r_add = await client.post("https://maayanlab.cloud/Enrichr/addList", files=payload, timeout=20.0)
                break
            except Exception as exc:
                if attempt == 0:
                    logger.warning(f"⏳ Enrichr addList reintentando tras fallo: {exc}")
                    await asyncio.sleep(2.0)
                else:
                    raise
        if r_add is None or r_add.status_code != 200: return []
        user_list_id = r_add.json().get("userListId")
        if not user_list_id: return []

        # 2. Consultar bases de datos clave
        sources = {
            "GO_Biological_Process_2023": "GO",
            "Reactome_Pathways_2024": "Reactome",
            "WikiPathways_2019_Human": "WikiPathways",
            "KEGG_2021_Human": "KEGG"
        }
        
        for gs_full, label in sources.items():
            try:
                r_enr = await client.get(
                    "https://maayanlab.cloud/Enrichr/enrich",
                    params={"userListId": user_list_id, "backgroundType": gs_full},
                    timeout=20.0
                )
                if r_enr.status_code == 200:
                    data = r_enr.json().get(gs_full, [])
                    # Enrichr: [Rank, Term, P-val, Z-score, Combined, Genes, AdjP, ...]
                    for row in data[:8]:
                        if row[6] < 0.05: # p-adj < 0.05 (auditoría externa, 2026-11-19: 0.3 no es significancia)
                            results.append({
                                "Term": row[1],
                                "Source": label,
                                "Pval": row[6],
                                "RawTerm": row[1]
                            })
            except: continue
        
        results.sort(key=lambda x: x["Pval"])
        final_results = results[:24]
        PERSISTENT_CACHE["enrichr"][cache_key] = final_results
        return final_results
    except Exception as e:
        logger.error(f"❌ Fallo Enrichr Directo: {e}")
        return []

# ── PUBMED: BÚSQUEDA DE METADATOS REALES DE ARTÍCULO ──────────────────────────
async def fetch_real_article(pmid: str, client: httpx.AsyncClient) -> dict:
    """Recupera título, autores, año y revista reales desde PubMed esummary.
    Devuelve {} si el PMID no resuelve o si la API falla."""
    if not pmid:
        return {}
    cache_key = f"esum_{pmid}"
    if cache_key in PERSISTENT_CACHE.get("pubmed", {}):
        return PERSISTENT_CACHE["pubmed"][cache_key]
    try:
        async with pubmed_semaphore:
            await asyncio.sleep(0.15)
            r = await client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
                params={"db": "pubmed", "id": pmid, "retmode": "json"},
                timeout=10.0,
            )
        if r.status_code != 200:
            return {}
        doc = r.json().get("result", {}).get(pmid, {})
        if not doc or doc.get("error"):
            return {}
        authors_list = doc.get("authors", []) or []
        author_names = [a.get("name") for a in authors_list if a.get("name")]
        if len(author_names) == 0:
            authors_str = ""
        elif len(author_names) == 1:
            authors_str = author_names[0]
        elif len(author_names) <= 3:
            authors_str = ", ".join(author_names)
        else:
            authors_str = f"{author_names[0]} et al."
        pubdate = doc.get("pubdate", "") or doc.get("epubdate", "")
        year = ""
        m = re.search(r"\b(19|20)\d{2}\b", pubdate)
        if m:
            year = m.group(0)
        meta = {
            "title": doc.get("title", "").strip().rstrip("."),
            "authors": authors_str,
            "year": year,
            "journal": doc.get("source", "") or doc.get("fulljournalname", ""),
        }
        PERSISTENT_CACHE["pubmed"][cache_key] = meta
        return meta
    except Exception as e:
        logger.warning(f"esummary falló para PMID {pmid}: {e}")
        return {}

# ── INVESTIGACIÓN MULTI-FUENTE PARA GENES CORE ────────────────────────────────
# Para cada gen core del Venn, buscamos evidencia REAL y FRESCA en:
#   - PubMed (artículos del gen + contexto biológico, con filtro de años)
#   - OMIM (catálogo de enfermedades mendelianas asociadas)
#   - ClinVar (variantes patogénicas clasificadas)
#   - ClinicalTrials.gov (ensayos clínicos activos)
# Resultados se cachean en PERSISTENT_CACHE["gene_research"] indexados por
# clave compuesta {gene}|{years}, para acelerar análisis futuros.

async def search_pubmed_for_gene(gene: str, term: str, years: int, client: httpx.AsyncClient) -> tuple:
    """Busca el PMID más relevante para un gen + contexto biológico, respetando
    el filtro de años. Cascada robusta con reintentos NCBI (429-aware).
    Devuelve (pmid, year_window_used) o (None, None).

    NOTA CRÍTICA: PubMed indexa SÓLO en inglés. Si `term` viene en español
    (frecuente en PRESETS), traducimos a una clave en inglés con un mapeo
    estático. Si no hay mapeo, usamos el gen + 'miRNA' como contexto natural
    al dominio de la app.
    """
    start_year = datetime.datetime.now().year - years
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # Mapeo ES → EN para los `term` típicos del PRESET. Sólo necesitamos
    # las keywords científicas — PubMed entiende ese subset perfectamente.
    ES_TO_EN = {
        "metabolismo del colesterol": "cholesterol metabolism",
        "ensamblaje de hdl": "HDL assembly",
        "srebf y mir33": "SREBF miR-33",
        "importación de proteínas al núcleo": "nuclear protein import",
        "efectos mediados por ns1": "NS1 protein effects",
        "sistema neuronal": "neuronal system",
        "percepción sensorial del dolor": "sensory pain perception",
        "respuesta al estrés osmótico": "osmotic stress response",
        "actividad de factores de transcripción": "transcription factor activity",
        "complejo asociado a distrofina": "dystrophin-associated complex",
        "organización del citoesqueleto": "cytoskeleton organization",
        "señalización mapk/erk": "MAPK ERK signaling",
        "regulación del complejo arp2/3": "Arp2/3 complex regulation",
        "señalización wnt/beta-catenina": "Wnt beta-catenin signaling",
        "evasión inmune pd-l1/pd-1": "PD-L1 PD-1 immune evasion",
    }

    # Limpiar: quitar IDs entre paréntesis, GO/WP/R-HSA codes
    clean_term_raw = re.sub(r'\([^)]*\)|GO:\d+|hsa\d+|WP\d+|R-HSA-\d+', '', term).strip()

    # Traducción a inglés (clave del mapping o fallback)
    en_term = ES_TO_EN.get(clean_term_raw.lower(), "")
    if not en_term:
        # Sin mapeo conocido: usar contexto natural del dominio
        en_term = "microRNA"

    # Cascada: 3 estrategias de búsqueda con reintentos NCBI
    strategies = [
        # 1. Estricto: gen + contexto inglés + años + humano
        (f'"{gene}"[Gene/Protein Name] AND "{en_term}"[Title/Abstract] AND ("{start_year}"[PDAT] : "3000"[PDAT]) AND humans[MeSH Terms]', f"{years}y"),
        # 2. Ampliado: gen + contexto inglés + humano (sin filtro de años)
        (f'"{gene}"[Gene/Protein Name] AND "{en_term}"[Title/Abstract] AND humans[MeSH Terms]', "any"),
        # 3. Mínimo: gen + miRNA + humano (sin contexto específico)
        (f'"{gene}"[Gene/Protein Name] AND "microRNA"[Title/Abstract] AND humans[MeSH Terms]', "any"),
        # 4. Último recurso: solo gen + humano (literatura general)
        (f'"{gene}"[Gene/Protein Name] AND humans[MeSH Terms]', "any"),
    ]

    for query, window in strategies:
        data = await safe_pubmed_request(client, url, {
            "db": "pubmed", "term": query, "retmax": 1,
            "retmode": "json", "sort": "relevance",
        }, retries=3)
        if not data:
            continue
        ids = data.get("esearchresult", {}).get("idlist", [])
        if ids:
            return ids[0], window
    return None, None


async def fetch_omim_for_gene(gene: str, client: httpx.AsyncClient) -> list:
    """Busca entradas OMIM (enfermedades hereditarias) asociadas al gen vía
    NCBI eUtils con reintentos 429-aware. Devuelve lista de
    [{omim_id, title, url}] (máximo 3 entradas)."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    search_data = await safe_pubmed_request(client, url, {
        "db": "omim", "term": f"{gene}[Gene Name]", "retmax": 3, "retmode": "json",
    }, retries=3)
    if not search_data:
        return []
    ids = search_data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return []

    # esummary con reintentos
    s_data = await safe_pubmed_request(
        client,
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
        {"db": "omim", "id": ",".join(ids), "retmode": "json"},
        retries=3,
    )
    if not s_data:
        # No pudimos resolver el título — devolver con título vacío
        # para que el caller pueda decidir si descartar o no.
        return [{"omim_id": i, "title": "", "url": f"https://www.omim.org/entry/{i}"} for i in ids]

    result = s_data.get("result", {})
    out = []
    for oid in ids:
        doc = result.get(oid, {})
        title = (doc.get("title", "") or doc.get("alttitles", "") or "").strip()
        out.append({
            "omim_id": oid,
            "title": title,
            "url": f"https://www.omim.org/entry/{oid}",
        })
    return out


async def fetch_clinvar_for_gene(gene: str, client: httpx.AsyncClient) -> dict:
    """Busca variantes patogénicas del gen en ClinVar vía NCBI eUtils con
    reintentos. Devuelve {count, url}."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    data = await safe_pubmed_request(client, url, {
        "db": "clinvar",
        "term": f"{gene}[Gene Name] AND (pathogenic[Clinical_Significance] OR likely_pathogenic[Clinical_Significance])",
        "retmax": 0, "retmode": "json",
    }, retries=3)
    if not data:
        return {"count": 0, "url": ""}
    try:
        count = int(data.get("esearchresult", {}).get("count", 0))
    except (ValueError, TypeError):
        count = 0
    return {
        "count": count,
        "url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}%5BGene+Name%5D",
    }


async def fetch_pubmed_count(gene: str, client: httpx.AsyncClient) -> int:
    """Conteo ligero de artículos en PubMed para un gen (esearch retmax=0).
    Solo el número, sin descargar artículos: alimenta el volcán real
    (eje X = frecuencia en PubMed, adaptación propia de la tesis Fig. 15)."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    data = await safe_pubmed_request(client, url, {
        "db": "pubmed",
        "term": f'"{gene}"[Gene/Protein Name] AND humans[MeSH Terms]',
        "retmax": 0, "retmode": "json",
    }, retries=2)
    if not data:
        return 0
    try:
        return int(data.get("esearchresult", {}).get("count", 0))
    except (ValueError, TypeError):
        return 0


async def fetch_string_network(genes: list, client: httpx.AsyncClient,
                               add_nodes: int = 0) -> list:
    """Interacciones reales entre `genes` vía STRING-DB (API pública, sin clave).
    Devuelve [(genA, genB, score 0-1)] con confianza >= 0,4
    (required_score=400, especie humana 9606). Solo aristas entre los genes
    pedidos; si STRING no devuelve ninguna, la lista queda vacía y el
    gráfico PPI no se dibuja (no se inventan aristas).
    Con add_nodes=N, STRING añade N proteínas vecinas de primera corona
    (las que más interactúan con los genes pedidos): las aristas hacia
    ellas también son reales y el dibujo las marca como vecinas."""
    genes = [g for g in dict.fromkeys((genes or [])) if g]
    if len(genes) < 2:
        return []
    try:
        params = {"identifiers": "%0d".join(genes),
                  "species": 9606, "required_score": 400}
        if add_nodes:
            params["add_nodes"] = add_nodes
        r = await client.get("https://string-db.org/api/json/network",
                             params=params, timeout=20.0)
        if r.status_code != 200:
            return []
        wanted = {g.upper() for g in genes}
        out = []
        for row in r.json():
            a = (row.get("preferredName_A") or "").upper()
            b = (row.get("preferredName_B") or "").upper()
            if a == b:
                continue
            if add_nodes:
                if a not in wanted and b not in wanted:
                    continue
            elif a not in wanted or b not in wanted:
                continue
            try:
                score = float(row.get("score", 0))
            except (TypeError, ValueError):
                continue
            out.append((a, b, score))
        seen, dedup = set(), []
        for a, b, s in out:
            key = tuple(sorted((a, b)))
            if key in seen:
                continue
            seen.add(key)
            dedup.append((key[0], key[1], s))
        return dedup
    except Exception:
        return []


def parse_uniprot_entry(entry: dict, gene: str) -> dict | None:
    """Extrae la relevancia de una proteína de una entrada UniProtKB
    (revisada). Puro: recibe el dict tal como lo devuelve la API.
    Devuelve {uniprot_id, protein_name, function, locations, diseases, url}."""
    if not entry:
        return None
    name = (entry.get("proteinDescription", {}).get("recommendedName", {})
            .get("fullName", {}).get("value", gene))
    function, locations, diseases = "", [], []
    for c in entry.get("comments", []):
        t = c.get("commentType")
        if t == "FUNCTION":
            function = " ".join(x.get("value", "") for x in c.get("texts", []))
        elif t == "SUBCELLULAR LOCATION":
            for loc in c.get("subcellularLocations", []):
                v = (loc.get("location") or {}).get("value", "")
                if v:
                    locations.append(v)
        elif t == "DISEASE":
            d = c.get("disease") or {}
            label = d.get("diseaseId") or d.get("acronym") or ""
            desc = d.get("description", "")
            if label:
                diseases.append(f"{label}: {desc}" if desc else label)
    acc = entry.get("primaryAccession", "")
    return {
        "uniprot_id": acc,
        "protein_name": name,
        "function": function[:600],
        "locations": locations[:3],
        "diseases": diseases[:3],
        "url": f"https://rest.uniprot.org/uniprotkb/{acc}" if acc else "",
    }


async def fetch_uniprot_for_gene(gene: str, client: httpx.AsyncClient) -> dict | None:
    """Relevancia de la proteína del gen vía UniProt (API pública, sin clave).
    Devuelve el dict de parse_uniprot_entry o None si UniProt no tiene
    entrada revisada para el gen (o si la red falla)."""
    try:
        r = await client.get("https://rest.uniprot.org/uniprotkb/search", params={
            "query": f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true",
            "fields": "accession,protein_name,cc_function,cc_subcellular_location,cc_disease",
            "format": "json", "size": 1,
        }, timeout=15.0)
        if r.status_code != 200:
            return None
        results = r.json().get("results", [])
        if not results:
            return None
        return parse_uniprot_entry(results[0], gene)
    except Exception:
        return None


async def fetch_clinicaltrials_for_gene(gene: str, client: httpx.AsyncClient) -> list:
    """Busca ensayos clínicos relacionados al gen vía ClinicalTrials.gov API v2.
    Filtra estrictamente por título que mencione el gen para descartar
    estudios donde el gen aparece sólo como mención tangencial.
    Devuelve lista de [{nct_id, title, status, url}] (máximo 3 ensayos).
    """
    try:
        url = "https://clinicaltrials.gov/api/v2/studies"
        # Pedimos más estudios y filtramos client-side por título
        async with pubmed_semaphore:
            await asyncio.sleep(0.15)
            r = await client.get(url, params={
                "query.term": gene,
                "filter.overallStatus": "RECRUITING,ACTIVE_NOT_RECRUITING,COMPLETED",
                "pageSize": 20,
                "fields": "NCTId,BriefTitle,OverallStatus",
            }, timeout=12.0)
        if r.status_code != 200:
            return []
        studies = r.json().get("studies", [])
        out = []
        gene_upper = gene.upper()
        for s in studies:
            ident = s.get("protocolSection", {}).get("identificationModule", {})
            status_mod = s.get("protocolSection", {}).get("statusModule", {})
            nct = ident.get("nctId", "")
            title = ident.get("briefTitle", "").strip()
            if not nct or not title:
                continue
            # FILTRO ESTRICTO: el título DEBE mencionar el gen (case-insensitive,
            # word-boundary para evitar matches parciales).
            if not re.search(rf'\b{re.escape(gene_upper)}\b', title.upper()):
                continue
            out.append({
                "nct_id": nct,
                "title": title,
                "status": status_mod.get("overallStatus", ""),
                "url": f"https://clinicaltrials.gov/study/{nct}",
            })
            if len(out) >= 3:
                break
        return out
    except Exception as e:
        logger.warning(f"ClinicalTrials.gov falló para {gene}: {e}")
        return []


async def fetch_disgenet(gene: str, client: httpx.AsyncClient, api_key=None):
    """Gene-disease associations from DisGeNET (disgenet.org).
    Requires a FREE API key (DISGENET_API_KEY env var): the public API
    answers without one only with the website HTML. Returns
    {"diseases": [{"name": ..., "score": ...}]} or {"diseases": []}.
    The response shape is parsed defensively (list, {"data": [...]},
    or {"results": [...]}) because the API version is unpinned."""
    if not api_key:
        return {"diseases": []}
    try:
        r = await client.get(
            f"https://www.disgenet.org/api/gda/gene/{gene}",
            headers={"Authorization": api_key}, timeout=12.0)
        if r.status_code != 200:
            return {"diseases": []}
        data = r.json()
        rows = data if isinstance(data, list) else data.get("data") or data.get("results") or []
        diseases = []
        for row in rows[:10]:
            name = row.get("disease_name") or row.get("disease") or ""
            if name:
                diseases.append({"name": name, "score": row.get("score")})
        return {"diseases": diseases}
    except Exception:
        return {"diseases": []}


async def fetch_opentargets_diseases(gene: str, client: httpx.AsyncClient) -> dict:
    """Disease associations from Open Targets (api.platform.opentargets.org,
    free GraphQL, no key): the gene symbol is resolved to an Ensembl ID
    via mygene.info, then the target's associated diseases are fetched
    with their association scores. Returns {"diseases": [{"name", "score"}],
    "ensembl_id": ...} (top 3 by score) or {"diseases": []}. Never raises:
    any failure (no Ensembl ID, API error, timeout) is an empty result,
    not an error — the disease tier is informative, never blocking."""
    try:
        r = await client.get(
            "https://mygene.info/v3/query",
            params={"q": f"symbol:{gene}", "fields": "ensembl.gene",
                    "species": "human"},
            timeout=15.0)
        if r.status_code != 200:
            return {"diseases": []}
        ensembl_id = None
        for hit in r.json().get("hits") or []:
            ens = hit.get("ensembl")
            if isinstance(ens, dict) and ens.get("gene"):
                ensembl_id = ens["gene"]
                break
            if isinstance(ens, list):
                for entry in ens:
                    if isinstance(entry, dict) and entry.get("gene"):
                        ensembl_id = entry["gene"]
                        break
                if ensembl_id:
                    break
        if not ensembl_id:
            return {"diseases": []}
        query = ("query TargetDiseases($ensemblId: String!) { "
                 "target(ensemblId: $ensemblId) { associatedDiseases { "
                 "rows { disease { name } score } } } }")
        r2 = await client.post(
            "https://api.platform.opentargets.org/api/v4/graphql",
            json={"query": query, "variables": {"ensemblId": ensembl_id}},
            timeout=20.0)
        if r2.status_code != 200:
            return {"diseases": []}
        target = (r2.json().get("data") or {}).get("target") or {}
        rows = (target.get("associatedDiseases") or {}).get("rows") or []
        diseases = []
        for row in rows:
            name = (row.get("disease") or {}).get("name")
            score = row.get("score")
            if name and score is not None:
                diseases.append({"name": name, "score": round(float(score), 3)})
        diseases.sort(key=lambda d: -d["score"])
        return {"diseases": diseases[:3], "ensembl_id": ensembl_id}
    except Exception:
        return {"diseases": []}


def dedupe_pubmed_refs(refs):
    """Elimina referencias duplicadas por PMID, conservando la primera.
    Pura, sin red. Cada término de evidencia clínica busca por separado en
    PubMed y dos términos pueden devolver el mismo artículo; sin esto, un
    mismo PMID sumaba 5 puntos dos veces en el puntaje de validación."""
    seen = set()
    out = []
    for r in refs or []:
        pmid = r.get("pmid") if isinstance(r, dict) else None
        if pmid and pmid in seen:
            continue
        if pmid:
            seen.add(pmid)
        out.append(r)
    return out


async def research_gene_complete(gene: str, preset: dict, years: int, client: httpx.AsyncClient) -> dict:
    """Orquesta TODA la investigación de un gen core: PubMed (con cascada de años)
    + OMIM + ClinVar + ClinicalTrials.gov. Cachea el resultado por {gene}|{years}.

    Devuelve dict con:
      - pubmed: [{pmid, term, source, editorial_note, year_window}]
      - omim: [{omim_id, title, url}]
      - clinvar: {count, url}
      - trials: [{nct_id, title, status, url}]
      - year_window_used: ventana real usada en pubmed (puede diferir del input)
    """
    # Cache key versionado: si subimos la versión, invalida resultados viejos
    # con queries en español/PMIDs incorrectos. v2 = búsquedas con términos
    # en inglés + filtro de título en Trials. v3 = v2 + deduplicación de
    # PubMed por PMID (sin ella, un mismo artículo sumaba 5 puntos dos veces).
    # v4 = v3 + H-06 (orden de prioridad por validación; los genes no
    # investigados ya no se puntúan). v5 = v4 + Open Targets: el resultado
    # gana "diseases" (asociaciones de enfermedad) y "diseases_source".
    cache_key = f"v5|{gene}|{years}"
    if cache_key in PERSISTENT_CACHE["gene_research"]:
        logger.info(f"  💾 Cache hit: {gene} (años={years})")
        return PERSISTENT_CACHE["gene_research"][cache_key]

    logger.info(f"  🔬 Investigando {gene} en multi-fuente (años={years})...")

    # 1. PubMed: una búsqueda por cada term del PRESET (contexto biológico real)
    pubmed_refs = []
    evidence_list = preset.get("clinical_evidence", [])
    if evidence_list:
        # Buscar PMID fresco para CADA evidencia clínica registrada en PRESETS
        pubmed_tasks = [
            search_pubmed_for_gene(gene, ev["term"], years, client)
            for ev in evidence_list
        ]
        pubmed_results = await asyncio.gather(*pubmed_tasks, return_exceptions=True)
        for ev, res in zip(evidence_list, pubmed_results):
            if isinstance(res, Exception) or not res or not res[0]:
                continue
            pmid, window = res
            # La descripción del preset es una NOTA EDITORIAL (contexto del
            # equipo, sin cita): va en su propio campo y NUNCA se presenta
            # como respaldada por este PMID. El PMID es solo el primer
            # resultado de la búsqueda en PubMed por gen + término.
            pubmed_refs.append({
                "pmid": pmid,
                "term": ev["term"],
                "source": ev["source"],
                "editorial_note": ev.get("desc", ""),
                "year_window": window,
            })
    else:
        # Sin clinical_evidence en PRESETS: una búsqueda genérica
        res = await search_pubmed_for_gene(gene, "function", years, client)
        if res and res[0]:
            pmid, window = res
            pubmed_refs.append({
                "pmid": pmid,
                "term": "Función biológica",
                "source": "PubMed",
                "editorial_note": "",
                "year_window": window,
            })

    # Deduplicar por PMID: cada término de evidencia busca por separado y
    # dos términos pueden devolver el mismo artículo; sin esto, un mismo
    # PMID sumaba 5 puntos dos veces en el puntaje de validación.
    pubmed_refs = dedupe_pubmed_refs(pubmed_refs)

    # 2-4. OMIM + ClinVar + ClinicalTrials + enfermedades en paralelo.
    # Enfermedades: DisGeNET cuando tiene clave; si no, Open Targets
    # (gratis, sin clave). El origen queda registrado en diseases_source.
    omim, clinvar, trials, disgenet, opentargets = await asyncio.gather(
        fetch_omim_for_gene(gene, client),
        fetch_clinvar_for_gene(gene, client),
        fetch_clinicaltrials_for_gene(gene, client),
        fetch_disgenet(gene, client, os.environ.get("DISGENET_API_KEY")),
        fetch_opentargets_diseases(gene, client),
        return_exceptions=True,
    )
    if isinstance(omim, Exception): omim = []
    if isinstance(clinvar, Exception): clinvar = {"count": 0, "url": ""}
    if isinstance(trials, Exception): trials = []
    if isinstance(disgenet, Exception): disgenet = {"diseases": []}
    if isinstance(opentargets, Exception): opentargets = {"diseases": []}
    disgenet_diseases = (disgenet or {}).get("diseases") or []
    opentargets_diseases = (opentargets or {}).get("diseases") or []
    if disgenet_diseases:
        diseases, diseases_source = disgenet_diseases, "disgenet"
    elif opentargets_diseases:
        diseases, diseases_source = opentargets_diseases, "opentargets"
    else:
        diseases, diseases_source = [], "none"

    # Ventana real usada: si alguna ref usó "any" (fallback), el campo lo refleja
    year_windows = [r["year_window"] for r in pubmed_refs]
    if not year_windows:
        year_window_used = "none"
    elif all(w == f"{years}y" for w in year_windows):
        year_window_used = f"{years}y"
    elif any(w == "any" for w in year_windows):
        year_window_used = f"{years}y_ampliado"
    else:
        year_window_used = f"{years}y"

    result = {
        "pubmed": pubmed_refs,
        "omim": omim,
        "clinvar": clinvar,
        "trials": trials,
        "diseases": diseases,
        "diseases_source": diseases_source,
        "year_window_used": year_window_used,
    }
    PERSISTENT_CACHE["gene_research"][cache_key] = result
    logger.info(
        f"  ✅ {gene}: {len(pubmed_refs)} PubMed, {len(omim)} OMIM, "
        f"{clinvar.get('count', 0)} variantes ClinVar, {len(trials)} ensayos clínicos"
    )
    return result


# ── PUBMED EXPERTO ────────────────────────────────────────────────────────────
async def get_pubmed_for_term(term: str, years: int, client: httpx.AsyncClient):
    """Query PubMed con caché persistente y filtros oficiales."""
    cache_key = f"{term}_{years}"
    if cache_key in PERSISTENT_CACHE["pubmed"]:
        return PERSISTENT_CACHE["pubmed"][cache_key]
    
    # Limpiar el término: quitar paréntesis, IDs, texto extra de Enrichr
    clean = re.sub(r'\s*-\s*(Homo sapiens|KEGG|Reactome|WikiPathways).*$', '', term, flags=re.IGNORECASE)
    clean = re.sub(r'hsa\d+|GO:\d+|\(.*?\)|R-HSA-\d+', '', clean).strip()
    
    # Tomar solo las palabras clave para evitar queries vacías por exceso de especificidad
    words = clean.split()[:5]
    simple_term = " ".join(words)
    if len(simple_term) < 3: return []
    
    current_year = datetime.datetime.now().year
    start_year = current_year - years
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    # Estrategia de búsqueda en cascada
    strategies = [
        f'({simple_term}[Title/Abstract]) AND ("{start_year}"[PDAT] : "3000"[PDAT]) AND human[MH]',
        f'({simple_term}) AND human[Organism]',
        f'{simple_term}'
    ]
    
    for query in strategies:
        data = await safe_pubmed_request(client, url, {"db": "pubmed", "term": query, "retmax": 1, "retmode": "json"})
        if data:
            ids = data.get("esearchresult", {}).get("idlist", [])
            if ids:
                res = [{"title": "", "id": ids[0]}]
                PERSISTENT_CACHE["pubmed"][cache_key] = res
                return res
    return []
# ── DETALLES DE GEN ───────────────────────────────────────────────────────────
SYSTEM_MAP = {
    "Cardiovascular": ["heart", "cardio", "vessel", "blood", "artery", "atherosclerosis", "cholesterol", "hdl", "ldl", "lipoprotein"],
    "Neurológico": ["brain", "neuron", "synapse", "nervous", "epilepsy", "pain", "sensory", "dravet", "channel", "axon"],
    "Metabólico": ["metabolism", "lipid", "glucose", "insulin", "diabetes", "obesity", "fat", "eflujo", "transport"],
    "Oncológico": ["cancer", "tumor", "oncogene", "proliferation", "apoptosis", "metastasis", "carcinoma", "growth"],
    "Inmunológico": ["immune", "inflammation", "cytokine", "leukocyte", "t-cell", "b-cell", "infection", "interferon"],
    "Celular": ["transport", "nucleus", "importin", "cytoskeleton", "mitochondria", "cell cycle", "autophagy"],
    "Renal": ["nephropathy", "kidney", "renal", "osmotic", "urea", "filtration"],
}

async def get_gene_details(gene: str, client: httpx.AsyncClient) -> dict:
    info = {"full_name": gene, "system": "Multisistémico",
            "pathology": "Nodo regulador crítico identificado en convergencia de miRNAs.",
            "associated_routes": [], "pmid": None,
            "protein": None}
    
    # Descripciones de ÉLITE con EVIDENCIA ESTRUCTURADA
    # NOTA (2026-11-19): los PMIDs escritos a mano de clinical_evidence se
    # eliminaron porque los 11 resolvían en PubMed a artículos SIN relación
    # con el gen ni con la frase (p. ej. veterinaria, salmón, escoliosis).
    # La evidencia PubMed verificada la aporta la capa viva de investigación
    # (research_gene_complete), no estos presets. No volver a escribir PMIDs
    # a mano: todo PMID que se muestre debe venir de una consulta en vivo.
    PRESETS = {
        "ABCA1": {
            "full_name": "ATP binding cassette subfamily A member 1",
            "system": "Cardiovascular / Metabólico",
            "pathology": "Regulador maestro del eflujo de colesterol. Su deficiencia causa el Síndrome de Tangier y aumenta drásticamente el riesgo de aterosclerosis coronaria.",
            "clinical_evidence": [
                {"term": "Metabolismo del colesterol", "source": "KEGG", "pmid": None, "desc": "miR-33 suprime ABCA1, reduciendo el eflujo de colesterol en macrófagos."},
                {"term": "SREBF y miR33", "source": "WikiPathways", "pmid": None, "desc": "La inhibición de miR-33 restaura ABCA1 y aumenta HDL en modelos animales."},
                {"term": "Ensamblaje de HDL", "source": "Reactome", "pmid": None, "desc": "Ruta crítica para la formación de lipoproteínas mediada por ABCA1."}
            ],
            "full_routes": [
                "Ensamblaje de HDL", "SREBF y miR33 en la homeostasis de lípidos WP2011",
                "Ruta metabólica de LDL, HDL y TG WP4522", "Ensamblaje de lipoproteínas plasmáticas",
                "Ruta de las estatinas WP430", "Receptores nucleares en metabolismo lipídico WP299",
                "Digestión y absorción de grasas", "Regulación de expresión génica por NR1H3 y NR1H2",
                "Transportadores ABC", "Metabolismo del colesterol", "Metabolismo de Vitamina B12 WP1533",
                "Señalización mediada por NR1H2 y NR1H3", "Trastornos de transportadores ABC",
                "Metabolismo del Folato WP176"
            ],
            "conclusion": "La evidencia sugiere que antagonistas de miR-33 podrían ser una estrategia terapéutica en enfermedades cardiovasculares asociadas a ABCA1."
        },
        "KPNA3": {
            "full_name": "Karyopherin subunit alpha 3",
            "system": "Celular / Inmunológico",
            "pathology": "Importina clave en el transporte nucleocitoplasmático de factores de transcripción. Implicada en la respuesta al estrés celular y la progresión oncológica.",
            "clinical_evidence": [
                {"term": "Importación de proteínas al núcleo", "source": "GO", "pmid": None, "desc": "Modulación de la importación de proteínas críticas para la respuesta al estrés."},
                {"term": "Efectos mediados por NS1", "source": "WikiPathways", "pmid": None, "desc": "Interacción con rutas virales que afectan la estabilidad genómica."}
            ],
            "full_routes": [
                "Importación de proteínas con NLS al núcleo (GO:0006607)", "Efectos mediados por NS1 en rutas del huésped",
                "Importación de proteínas al núcleo", "Propagación de ARNm viral", "Interacción huésped-patógeno"
            ],
            "conclusion": "Su regulación es vital para el mantenimiento de la homeostasis proteica bajo estrés celular."
        },
        "SCN1A": {
            "full_name": "Sodium voltage-gated channel alpha subunit 1",
            "system": "Neurológico",
            "pathology": "Canal de sodio esencial para la excitabilidad neuronal. Mutaciones están vinculadas al Síndrome de Dravet y diversas formas de epilepsia severa.",
            "clinical_evidence": [
                {"term": "Sistema Neuronal", "source": "Reactome", "pmid": None, "desc": "Disfunción de interneuronas GABAérgicas por alteración del canal de sodio."},
                {"term": "Percepción sensorial del dolor", "source": "GO", "pmid": None, "desc": "Implicado en la detección y procesamiento de estímulos nociceptivos."}
            ],
            "full_routes": [
                "Detección de estímulos en percepción del dolor (GO:0062149)",
                "Detección de estímulos mecánicos en percepción sensorial (GO:0050974)",
                "Canales de sodio dependientes de voltaje", "Sistema Neuronal", "Canalopatías"
            ],
            "conclusion": "La modulación de SCN1A por miRNAs representa un eje terapéutico en epilepsias refractarias."
        },
        "TSC22D2": {
            "full_name": "TSC22 domain family member 2",
            "system": "Metabólico / Renal",
            "pathology": "Regulador de la transcripción inducido por estrés. Vinculado a la respuesta osmótica en nefropatías.",
            "clinical_evidence": [
                {"term": "Respuesta al estrés osmótico", "source": "GO", "pmid": None, "desc": "Regulación de la homeostasis de aminoácidos bajo estrés osmótico renal."},
                {"term": "Actividad de factores de transcripción", "source": "KEGG", "pmid": None, "desc": "Control transcripcional en rutas de supervivencia celular."}
            ],
            "full_routes": ["Regulación por miRNA", "Actividad de factores de transcripción", "Respuesta al estrés osmótico", "Homeostasis renal"],
            "pmid": None,
            "conclusion": "Su papel como regulador de estrés lo posiciona como biomarcador de progresión en nefropatías."
        },
        "SNTB2": {
            "full_name": "Syntrophin beta 2",
            "system": "Cardiovascular / Estructural",
            "pathology": "Proteína adaptadora del complejo distrofina. Su alteración se asocia a cardiomiopatías y disfunción del citoesqueleto.",
            "clinical_evidence": [
                {"term": "Complejo asociado a distrofina", "source": "Reactome", "pmid": None, "desc": "Afecta la integridad de la membrana celular en cardiomiocitos."},
                {"term": "Organización del citoesqueleto", "source": "GO", "pmid": None, "desc": "Regula el anclaje de proteínas de membrana al citoesqueleto de actina."}
            ],
            "full_routes": ["Regulación por miRNA", "Complejo de glicoproteína asociado a distrofina", "Organización del citoesqueleto", "Rutas de cardiomiopatía"],
            "conclusion": "La pérdida de integridad estructural por desregulación de SNTB2 es un factor clave en cardiomiopatías estructurales."
        },
        "GMFB": {
            "full_name": "Glia maturation factor beta",
            "system": "Neurológico",
            "pathology": "Factor neurotrófico esencial para la maduración glial y la regeneración neural. Su sobreexpresión está vinculada a neuroinflamación en Alzheimer y Parkinson.",
            "clinical_evidence": [
                {"term": "Señalización MAPK/ERK", "source": "KEGG", "pmid": None, "desc": "Regula la diferenciación y supervivencia celular en el sistema nervioso central."},
                {"term": "Regulación del complejo Arp2/3", "source": "GO", "pmid": None, "desc": "Modula la dinámica de actina, vital para la plasticidad sináptica y movilidad glial."}
            ],
            "full_routes": [
                "Señalización MAPK/ERK", "Diferenciación de astrocitos", "Regulación del complejo Arp2/3",
                "Respuesta a la neuroinflamación", "Mantenimiento de la red de actina"
            ],
            "conclusion": "GMFB actúa como un mediador clave en la respuesta neuroinflamatoria y el desarrollo glial."
        },
        "KLF12": {
            "full_name": "Kruppel-like factor 12",
            "system": "Oncológico / Inmunológico",
            "pathology": "Factor de transcripción tipo dedos de zinc que actúa como represor génico (AP-2rep). Vinculado a susceptibilidad en Artritis Reumatoide y progresión tumoral.",
            "clinical_evidence": [
                {"term": "Señalización Wnt/beta-catenina", "source": "WikiPathways", "pmid": None, "desc": "Regula la proliferación y metástasis en diversos tipos de cáncer y endometriosis."},
                {"term": "Evasión Inmune PD-L1/PD-1", "source": "Reactome", "pmid": None, "desc": "Participa en los mecanismos de escape inmunitario del microambiente tumoral."}
            ],
            "full_routes": [
                "Señalización Wnt/beta-catenina", "Represión de la transcripción por ARN Pol II",
                "Evasión inmunitaria PD-L1/PD-1", "Regulación del ciclo celular G1/S", "Susceptibilidad a Artritis Reumatoide"
            ],
            "conclusion": "KLF12 es un nodo regulador crítico en la interfaz de la respuesta inmunitaria y la progresión neoplásica."
        }
    }
    
    if gene in PRESETS:
        p = PRESETS[gene]
        # Combinar rutas de evidencia y rutas completas para el Sidebar
        all_routes = sorted(list(set([f"{ev['source']}: {ev['term']}" for ev in p["clinical_evidence"]] + p["full_routes"])))
        info.update({
            "full_name": p["full_name"],
            "system": p["system"],
            "pathology": p["pathology"],
            "associated_routes": all_routes,
            "clinical_evidence": p["clinical_evidence"],
            "conclusion": p.get("conclusion", ""),
            "pmid": p["clinical_evidence"][0]["pmid"] if p["clinical_evidence"] else None,
        })
        # Relevancia de la proteína (UniProt): qué hace, dónde trabaja,
        # con qué enfermedad se asocia.
        info["protein"] = await fetch_uniprot_for_gene(gene, client)
        return info

    try:
        r = await client.get(f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&fields=name,summary,pathway,go", timeout=12.0)
        if r.status_code == 200:
            hits = r.json().get('hits', [])
            if hits:
                h = hits[0]
                if h.get('name') and gene not in PRESETS: info["full_name"] = h['name']
                
                # Determinar sistema si no es preset
                summary_text = (h.get('summary') or "").lower()
                if gene not in PRESETS:
                    for sys_name, keywords in SYSTEM_MAP.items():
                        if any(k in summary_text for k in keywords):
                            info["system"] = sys_name
                            break
                
                if h.get('summary') and gene not in PRESETS:
                    # Traducción/Adaptación profesional al español con más patrones
                    raw_sum = h['summary']
                    if "involved in" in raw_sum.lower():
                        clean_desc = raw_sum.split("involved in")[-1].split(".")[0].strip()
                        trans_desc = await translate_to_spanish(clean_desc, client)
                        info["pathology"] = f"Proteína crítica implicada en {trans_desc}. Actúa como un nodo de regulación esencial en la homeostasis celular y la estabilidad tisular."
                    elif "is a member of" in raw_sum.lower():
                        clean_desc = raw_sum.split("is a member of")[-1].split(".")[0].strip()
                        trans_desc = await translate_to_spanish(clean_desc, client)
                        info["pathology"] = f"Miembro clave de la familia {trans_desc}. Desempeña un papel fundamental en la señalización celular y la respuesta adaptativa."
                    elif "encodes a protein that" in raw_sum.lower():
                        clean_desc = raw_sum.split("encodes a protein that")[-1].split(".")[0].strip()
                        trans_desc = await translate_to_spanish(clean_desc, client)
                        info["pathology"] = f"Codifica una proteína que {trans_desc}. Actúa como componente clave en rutas biológicas fundamentales identificadas mediante análisis genómico."
                    else:
                        trans_full = await translate_to_spanish(raw_sum[:250], client)
                        info["pathology"] = f"Nodo regulador identificado mediante análisis genómico. {trans_full}... integrando señales para el mantenimiento del equilibrio fisiológico."

                pways = h.get('pathway', {})
                routes_with_sources = []
                # 1. Rutas clásicas
                for lib_key, lib_name in [('kegg', 'KEGG'), ('reactome', 'Reactome'), ('wikipathways', 'WikiPathways')]:
                    if lib_key in pways:
                        d = pways[lib_key]; items = d if isinstance(d, list) else [d]
                        for x in items[:3]:
                            trans_route = await translate_to_spanish(x['name'], client)
                            routes_with_sources.append(f"{lib_name}: {trans_route}")
                
                # 2. Fallback a GO (Gene Ontology) si hay pocas rutas
                if len(routes_with_sources) < 4:
                    go_terms = h.get('go', {}).get('BP', [])
                    if isinstance(go_terms, dict): go_terms = [go_terms]
                    for x in go_terms[:5]:
                        if x.get('term'):
                            trans_go = await translate_to_spanish(x['term'], client)
                            routes_with_sources.append(f"GO: {trans_go}")

                info["associated_routes"] = routes_with_sources
                
                # Crear conclusión genérica para nuevos genes: lenguaje de
                # hipótesis computacional (H-12, v1.47) — candidato para
                # validación experimental, nunca "punto crítico de intervención".
                info["conclusion"] = (f"La convergencia de miRNAs sobre {gene} lo señala "
                    f"como candidato para validación experimental en procesos de {info['system'].lower()}.")
    except:
        pass
        
    try:
        async with pubmed_semaphore:
            await asyncio.sleep(0.2)
            pm = await client.get(
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
                params={"db": "pubmed", "term": f"{gene}[Gene/Protein Name] AND microRNA[Title/Abstract] AND human[MH]", "retmax": 1, "retmode": "json"}, timeout=10.0)
        ids = pm.json().get('esearchresult', {}).get('idlist', [])
        if ids: info["pmid"] = ids[0]
    except:
        pass
    info["protein"] = await fetch_uniprot_for_gene(gene, client)
    return info

# ── GRÁFICOS ──────────────────────────────────────────────────────────────────
# H-21 (v1.47): una sola constante de versión — la usan FastAPI, /health y
# la respuesta de /analyze (engine_version), para que la exportación no
# dependa de un literal fijo en el frontend.
APP_VERSION = "1.59"


@asynccontextmanager
async def lifespan(app):
    # v1.52: @app.on_event("startup") is deprecated (it raised a
    # DeprecationWarning on every run); the lifespan handler replaces
    # it. Same work as the old startup_event: load the cache and the
    # local databases before serving the first request.
    global TARGETSCAN_DB
    load_persistent_cache()
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()
    # v1.53: el GMT de Harmonizome ya no se distribuye; si no hay índice
    # oficial ni GMT en disco, se intenta descargar al arrancar (best-effort).
    await _ensure_legacy_gmt()
    load_mirtarbase_local()
    load_gmt_libraries()
    logger.info(f"🚀 KENRYU v{APP_VERSION} listo — {len(TARGETSCAN_DB)} familias miRNA cargadas.")
    yield


app = FastAPI(title="KENRYU Bioinformatics Engine", version=APP_VERSION,
              lifespan=lifespan)

# Seguridad CORS (H-14, v1.47): solo el origen exacto del Space. Antes el
# regex https://.*\.hf\.space con allow_credentials=True confiaba en
# cualquier Space de terceros. La app no usa cookies de sesión, así que
# allow_credentials=False.
SPACE_ORIGIN = "https://kenryu007-bioinformatica.hf.space"
origins = [
    "http://localhost",
    "http://localhost:3001",
    "http://127.0.0.1",
    "http://127.0.0.1:3001",
    "https://huggingface.co",
    SPACE_ORIGIN,
    "https://abrangel.github.io",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=False,
    allow_methods=["GET", "POST", "HEAD", "OPTIONS"],
    allow_headers=["*"],
)


@app.middleware("http")
async def security_headers(request, call_next):
    """Cabeceras de seguridad básicas (H-14, v1.47).
    La CSP permite los CDN que usa la interfaz (cdnjs, Google Fonts) y
    'unsafe-inline' en scripts y estilos porque la UI usa manejadores
    onclick y estilos en línea; la defensa contra XSS es escapeHtml en
    el frontend, no el bloqueo de scripts. frame-ancestors permite el
    iframe de Hugging Face Spaces (la vista del Space) y nada más."""
    response = await call_next(request)
    response.headers["X-Content-Type-Options"] = "nosniff"
    response.headers["Referrer-Policy"] = "strict-origin-when-cross-origin"
    response.headers["Content-Security-Policy"] = (
        "default-src 'self'; "
        "script-src 'self' https://cdnjs.cloudflare.com 'unsafe-inline'; "
        "style-src 'self' 'unsafe-inline' https://cdnjs.cloudflare.com https://fonts.googleapis.com; "
        "font-src 'self' https://fonts.gstatic.com https://cdnjs.cloudflare.com; "
        "img-src 'self' data:; connect-src 'self'; object-src 'none'; base-uri 'self'; "
        "frame-ancestors 'self' https://huggingface.co https://*.hf.space;"
    )
    return response


# ── RATE LIMIT (H-10, v1.50) ─────────────────────────────────────────────────
# La auditoría encontró que ningún endpoint tenía límite de tasa: un solo
# cliente podía inundar /api/v1/analyze (cada llamada quema 30-90 s de
# peticiones vivas a NCBI/Enrichr) y agotar la cuota de esas APIs para
# todos. Ventana deslizante por IP, en memoria, sin dependencias nuevas.
# /health queda fuera: los health checks de Hugging Face no deben fallar.
# Detrás del proxy de HF todos los clientes llegan con la IP del proxy,
# así que la clave es X-Forwarded-For (primer salto) cuando está
# presente. La cabecera es falsificable, pero el objetivo es frenar
# inundaciones, no un sistema anti-abuso perfecto.
# Configuración por entorno: KENRYU_RATE_LIMIT=0 lo desactiva,
# KENRYU_RATE_LIMIT_PER_MINUTE cambia el presupuesto (60 por defecto).
RATE_LIMIT_PER_MINUTE = int(os.environ.get("KENRYU_RATE_LIMIT_PER_MINUTE", "60"))
RATE_LIMIT_ENABLED = os.environ.get("KENRYU_RATE_LIMIT", "1") != "0"


class _RateLimiter:
    def __init__(self, per_minute=60, window=60.0):
        self.per_minute = per_minute
        self.window = window
        self._hits = {}  # ip -> deque[timestamps]
        self._lock = threading.Lock()

    def _client_ip(self, request):
        xff = request.headers.get("x-forwarded-for")
        if xff:
            return xff.split(",")[0].strip()
        return request.client.host if request.client else "unknown"

    def allowed(self, request):
        if not RATE_LIMIT_ENABLED:
            return True
        ip = self._client_ip(request)
        now = time.monotonic()
        with self._lock:
            q = self._hits.setdefault(ip, deque())
            while q and q[0] <= now - self.window:
                q.popleft()
            if len(q) >= self.per_minute:
                return False
            q.append(now)
            return True

    def reset(self):
        with self._lock:
            self._hits.clear()


RATE_LIMITER = _RateLimiter(RATE_LIMIT_PER_MINUTE)


@app.middleware("http")
async def rate_limit(request, call_next):
    """429 con Retry-After cuando una IP supera el presupuesto en
    /api/v1/*. /health y la UI estática quedan fuera, y también los
    preflight OPTIONS (CORS, sin coste de cuota)."""
    path = request.url.path
    if (path.startswith("/api/v1/") and path != "/health"
            and request.method != "OPTIONS"):
        if not RATE_LIMITER.allowed(request):
            return JSONResponse(
                {"detail": "rate limit exceeded, try again later"},
                status_code=429,
                headers={"Retry-After": "60"},
            )
    return await call_next(request)


app.mount("/static", StaticFiles(directory="static"), name="static")

@app.get("/")
@app.head("/")
async def root():
    return FileResponse('static/index.html')

@app.get("/script.js")
async def serve_script():
    return FileResponse('static/script.js')

@app.get("/style.css")
async def serve_style():
    return FileResponse('static/style.css')

def mirtarbase_mode():
    """How miRTarBase evidence is sourced in this process (v1.54):
    "official" (index 9.0, verified levels), "legacy" (Harmonizome GMT
    snapshot, no verified level) or "none" (no local miRTarBase)."""
    return {"official": "official", "gmt": "legacy"}.get(MIRTARBASE_FORMAT, "none")


LEGACY_MODE_WARNING = (
    "Modo legacy: miRTarBase procede del snapshot de Harmonizome (sin nivel "
    "de evidencia verificado); los genes con evidencia legacy quedan en "
    "'media' como máximo y no pueden llegar a 'alta'. Para reproducir los "
    "resultados del artículo, construye el índice oficial: descarga "
    "hsa_MTI.xlsx desde miRTarBase y ejecuta `python "
    "build_mirtarbase_official.py` (ver DATA_LICENSES.md).")


@app.get("/health")
async def health():
    return {"status": "OK", "version": APP_VERSION,
            "commit": load_commit_hash(),
            "db_families": len(TARGETSCAN_DB),
            "raw_files": len(list(RAW_DATA_DIR.glob("*.txt"))) if RAW_DATA_DIR.exists() else 0,
            "data_sources": loaded_data_versions(),
            "mirtarbase_mode": mirtarbase_mode()}


@app.get("/api/v1/updates")
def updates():
    """Self-update check: which reference databases are outdated, and
    where to download the new release. Never fails the request: with no
    network every source reports 'unknown'."""
    try:
        return check_updates(timeout=8.0)
    except Exception as exc:  # pragma: no cover - belt and braces
        logger.error(f"❌ Fallo chequeo de actualizaciones: {exc}")
        return {"checked_at": None, "updates_available": 0,
                "unpinned": 0, "unknown": len(DATA_SOURCES),
                "sources": [{**s, "latest_checked": None, "status": "unknown"}
                            for s in DATA_SOURCES]}


@app.post("/api/v1/updates/{source}/install")
def install_update_endpoint(source: str, request: Request):
    """Install the new release for one source. Only the registry's own
    whitelisted files are downloaded; the old file is kept as a backup and
    the new version is recorded. Returns the install report (the notice).

    PROTEGIDO (auditoría externa, 2026-11-19): antes cualquier cliente
    anónimo podía llamar a este endpoint y sobrescribir archivos del
    servidor. Ahora exige el secreto KENRYU_ADMIN_TOKEN en la cabecera
    X-Admin-Token; si la variable no está configurada, las
    actualizaciones están desactivadas (403). Tras instalar, las DB en
    memoria se invalidan para recargar los archivos nuevos en la
    siguiente petición."""
    global TARGETSCAN_DB, MIRTARBASE_DB
    admin_token = os.environ.get("KENRYU_ADMIN_TOKEN")
    if not admin_token:
        raise HTTPException(403, "updates are disabled: KENRYU_ADMIN_TOKEN is not configured")
    if not hmac.compare_digest(request.headers.get("X-Admin-Token", ""), admin_token):
        raise HTTPException(403, "invalid admin token")
    try:
        report = install_source_update(source)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    except Exception as exc:
        logger.error(f"❌ Fallo instalando actualización {source}: {exc}")
        raise HTTPException(502, f"install failed: {exc}")
    TARGETSCAN_DB = {}
    MIRTARBASE_DB = {}
    return report


@app.get("/api/v1/reverse")
def reverse(gene: str, species: str = "human"):
    """Búsqueda inversa: qué miRNAs (familias TargetScan locales) tienen
    a `gene` como diana predicha."""
    global TARGETSCAN_DB
    if species == "mouse":
        db = load_mouse_data()
    else:
        if not TARGETSCAN_DB:
            TARGETSCAN_DB = load_local_data()
        db = TARGETSCAN_DB
    hits = reverse_lookup(gene, db)
    return {"gene": (gene or "").strip().upper(), "species": species,
            "mirnas": hits, "n": len(hits)}


def panel_core(mirnas, mode="strict"):
    """Local-only core of one panel: TargetScan sets + consensus.
    No network, no live calls — the fast path behind /compare and /batch.
    Returns (found_mirnas, common_genes)."""
    global TARGETSCAN_DB
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()
    load_mirtarbase_local()
    names, sets = [], []
    for m in [normalize_isomir(x.strip()) for x in mirnas if x.strip()]:
        local = get_targets_local(m)
        remote = get_validated_local_checked(m, MIRTARBASE_DB)
        combined = local if local else remote
        if combined:
            names.append(m)
            sets.append(combined)
    if not sets:
        return names, []
    genes = consensus_genes(sets, mode)
    return names, sorted(genes)


class PanelItem(BaseModel):
    name: str = "panel"
    mirnas: List[str] = Field(min_length=1, max_length=30)


class CompareRequest(BaseModel):
    panels: List[PanelItem] = Field(min_length=1, max_length=20)
    mode: Literal["strict", "n-1", "n-2"] = "strict"


@app.post("/api/v1/compare")
def compare(req: CompareRequest):
    """Compare several panels (local data only): each panel's core genes,
    the union, the intersection of all panels, and every pairwise overlap."""
    out_panels = []
    for p in req.panels:
        name = str(p.name or "panel")
        found, genes = panel_core(p.mirnas or [], req.mode)
        out_panels.append({"name": name, "mirnas": found,
                           "genes": genes, "n": len(genes)})
    union = sorted(set().union(*[set(p["genes"]) for p in out_panels]) or set())
    inter = sorted(set(out_panels[0]["genes"]).intersection(
        *[set(p["genes"]) for p in out_panels[1:]]) if out_panels else [])
    pairwise = []
    for i, a in enumerate(out_panels):
        for b in out_panels[i + 1:]:
            shared = sorted(set(a["genes"]) & set(b["genes"]))
            pairwise.append({"a": a["name"], "b": b["name"],
                             "shared": shared, "n": len(shared)})
    return {"panels": out_panels, "union": union, "n_union": len(union),
            "intersection_all": inter, "n_intersection": len(inter),
            "pairwise": pairwise, "mode": req.mode}


@app.post("/api/v1/batch")
async def batch(file: UploadFile = File(...)):
    """Batch of panels from a file: CSV with a `panel,mirna` header (one
    row per miRNA), or plain text with one panel per line
    (`Name: mir1, mir2` or just `mir1, mir2`). Local data only."""
    raw_bytes = await file.read()
    if len(raw_bytes) > 256 * 1024:
        raise HTTPException(400, "file too large (max 256 KB)")
    raw = raw_bytes.decode("utf-8", errors="replace")
    panels = []
    lines = [ln for ln in raw.splitlines() if ln.strip()]
    if len(lines) > 500:
        raise HTTPException(400, "too many rows (max 500)")
    if lines and lines[0].strip().lower().startswith("panel"):
        import csv, io as _io
        reader = csv.DictReader(_io.StringIO(raw))
        grouped = {}
        for row in reader:
            name = (row.get("panel") or "panel").strip()
            mirna = (row.get("mirna") or "").strip()
            if mirna:
                grouped.setdefault(name, []).append(mirna)
        panels = [{"name": n, "mirnas": ms} for n, ms in grouped.items()]
    else:
        for ln in lines:
            if ":" in ln:
                name, rest = ln.split(":", 1)
                panels.append({"name": name.strip(),
                               "mirnas": [x.strip() for x in rest.split(",") if x.strip()]})
            else:
                panels.append({"name": f"panel-{len(panels) + 1}",
                               "mirnas": [x.strip() for x in ln.split(",") if x.strip()]})
    if not panels:
        raise HTTPException(400, "no panels found in the uploaded file")
    if len(panels) > 20:
        raise HTTPException(400, "too many panels (max 20)")
    for p in panels:
        if not p["mirnas"]:
            raise HTTPException(400, f"panel '{p['name']}' has no miRNAs")
        if len(p["mirnas"]) > 30:
            raise HTTPException(400, f"panel '{p['name']}' has too many miRNAs (max 30)")
    return compare(CompareRequest(panels=panels))
# Tope de investigación multi-fuente por análisis (v1.46, H-06): cada gen
# investigado cuesta varias llamadas a PubMed/OMIM/ClinVar/ClinicalTrials.
# Se investigan solo los RESEARCH_CAP genes priorizados; los demás figuran
# como "no evaluado" y la respuesta lo dice.
RESEARCH_CAP = 5


def loaded_data_versions():
    """Flat {source: version} of the reference data actually loaded in
    this process (v1.51, idea 2): the static record in data_versions.json
    says what SHOULD be there; this says what IS there. miRTarBase and
    KEGG are the two that can be absent in a fresh clone, so they report
    "not loaded" instead of the declared version; TargetScan reports the
    loaded family count. Same keys as load_local_versions(), so /health,
    the analysis JSON and the export all agree."""
    versions = dict(load_local_versions())
    if MIRTARBASE_FORMAT == "official":
        # The counts come from the loaded index's own metadata, so a
        # synthetic fixture can never claim the real file's numbers.
        meta = MIRTARBASE_META or {}
        if meta.get("synthetic"):
            versions["mirtarbase"] = (
                "miRTarBase official index (synthetic CI fixture, "
                f"{meta.get('n_mirnas', 0)} miRNAs, {meta.get('pairs', 0)} pairs)")
        elif meta.get("rows_read") is not None:
            versions["mirtarbase"] = (
                "miRTarBase 9.0 official index "
                f"({meta.get('built_from')}, {meta.get('rows_read'):,} rows, "
                f"{meta.get('n_mirnas'):,} miRNAs, {meta.get('pairs'):,} pairs)")
        else:  # index without build metadata (defensive fallback)
            versions["mirtarbase"] = ("miRTarBase 9.0 official index "
                                      "(hsa_MTI.xlsx, 957,040 rows, 2,947 miRNAs, "
                                      "396,758 pairs)")
    elif MIRTARBASE_FORMAT == "gmt":
        versions["mirtarbase"] = "Harmonizome GMT snapshot (legacy)"
    else:
        versions["mirtarbase"] = "not loaded"
    if "KEGG" not in GMT_LIBRARIES:
        versions["kegg"] = "not loaded"
    versions["targetscan"] = (f"8.0 ({len(TARGETSCAN_DB)} families)"
                              if TARGETSCAN_DB else "not loaded")
    return versions


@app.post("/api/v1/analyze")
async def analyze(req: AnalysisRequest):
    global TARGETSCAN_DB
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()

    # ... resto del procesamiento ...


    mirnas = [normalize_isomir(m.strip(), "mmu" if req.species == "mouse" else "hsa")
              for m in req.mirnas if m.strip()]
    # H-17 (v1.47): deduplicar preservando el orden — los duplicados
    # sesgaban el umbral de N-1/N-2 (contaban como miRNAs distintos).
    mirnas = list(dict.fromkeys(mirnas))
    if not mirnas:
        raise HTTPException(400, "Debe ingresar al menos un miRNA.")
    if req.species not in ("human", "mouse"):
        raise HTTPException(400, "species must be human or mouse")

    # ── RECOLECCIÓN DE GENES POR miRNA ───────────────────────────────────────────
    gene_sets, found_names = [], []
    missing_mirnas = []  # H-17 (v1.47): miRNAs sin datos en ninguna fuente.
    arm_unverified = []  # v1.58: brazos que TargetScan no modela (ver
    # check_arm_resolution): sus dianas predichas se bloquean, nunca se
    # sustituyen por las del otro brazo.
    arm_status = arm_verification_status()
    local_by_mirna, remote_by_mirna = {}, {}
    evidence_by_mirna = {}  # Ronda 3 (v1.48): {miRNA: {gen: evidencia}}
    source_by_mirna = {}  # por miRNA, qué fuente se usó (v1.46): un consenso
    # strict puede mezclar predicción de TargetScan con validación
    # experimental de miRTarBase; la respuesta lo dice por miRNA.
    if req.species == "mouse":
        db = load_mouse_data()
        species_prefix = "mmu"
    else:
        if not TARGETSCAN_DB:
            TARGETSCAN_DB = load_local_data()
        db = TARGETSCAN_DB
        species_prefix = "hsa"
    load_mirtarbase_local()
    for m in mirnas:
        local  = get_targets_local(m, db)
        if req.species != "mouse":
            check = check_arm_resolution(m, db)
            if check["status"] == "unverified":
                arm_unverified.append({"mirna": m, "reason": check["reason"],
                                       "family_key": check["family_key"]})
        # Validado: primero la librería local de miRTarBase; la API viva
        # de Harmonizome solo como respaldo (y nunca en modo offline).
        remote = get_validated_local_checked(m, MIRTARBASE_DB, species_prefix)
        evidence = get_validated_evidence_checked(m, MIRTARBASE_DB, species_prefix)
        if not remote and not req.offline:
            remote = await fetch_mirtarbase(m, species_prefix)
            # La API viva no distingue niveles: sus genes son "legacy".
            evidence = {g: {"level": "legacy", "pmids": [], "experiments": []}
                        for g in remote}
        local_by_mirna[m] = local
        remote_by_mirna[m] = remote
        evidence_by_mirna[m] = evidence

        # ── BUG 4 CORREGIDO: si hay datos locales, NO mezclar con remoto ────
        # miRTarBase contiene genes validados experimentalmente (no predichos),
        # mezclarlos con TargetScan diluyó la intersección en la versión anterior.
        # Ahora: usamos local (TargetScan) como fuente primaria.
        # Solo si local está vacío usamos miRTarBase como fallback.
        if local:
            combined = local
            source_by_mirna[m] = "TargetScan (predicción local)"
            if remote:
                logger.info(f"🔬 {m}: {len(local)} local (TargetScan) | {len(remote)} miRTarBase — usando solo TargetScan")
        else:
            combined = remote
            source_by_mirna[m] = (
                "miRTarBase (validación experimental; respaldo porque TargetScan no tenía este miRNA)"
                if remote else "sin datos en ninguna fuente"
            )
            logger.info(f"🔬 {m}: sin datos locales, usando miRTarBase: {len(remote)} genes")

        if combined:
            gene_sets.append(combined)
            found_names.append(m)
        else:
            missing_mirnas.append(m)
            logger.warning(f"⚠️ {m}: sin datos en ninguna fuente")

    if not gene_sets:
        raise HTTPException(404, "No se encontraron datos para los miRNAs ingresados.")

    # ── INTERSECCIÓN ───────────────────────────────────────────────────────────
    # BUG 5 CORREGIDO: el fallback automático a N-1 cuando la intersección
    # estricta era vacía hacía que con 5 miRNAs se devolvieran ~100 genes
    # en lugar de los 5 correctos.
    # Ahora: modo strict = solo intersección estricta, sin fallback silencioso.
    mode = req.mode
    common = consensus_genes(gene_sets, mode)
    logger.info(f"🎯 Modo {mode.upper()}: {len(common)} genes")

    # ── PRUEBA DE CONVERGENCIA (Unreleased) ──────────────────────────────────
    # p empírico del número de genes en común contra paneles aleatorios de
    # familias emparejadas por tamaño (ver panel_convergence). Solo tiene
    # sentido sobre dianas TargetScan: si algún miRNA del panel no tiene
    # datos locales, se omite con su motivo en lugar de comparar peras
    # con manzanas.
    if all(local_by_mirna[m] for m in found_names):
        convergence = panel_convergence(
            gene_sets, family_labels_for_targets(found_names, gene_sets, db),
            db, mode=mode)
    else:
        convergence = {
            "k_names": len(found_names), "k_effective": None,
            "same_family_warning": False, "genes_observed": len(common),
            "perms": 0, "seed": None, "resolution": None, "p_by_tolerance": {},
            "p_min": None, "p_max": None,
            "label": CONVERGENCE_LABEL, "limit": CONVERGENCE_LIMIT,
            "skipped": "needs_local_targets",
        }

    # ── MANIFIESTO DE EJECUCIÓN (Unreleased) ─────────────────────────────────
    # Todo lo necesario para repetir este análisis con un comando:
    # versión del motor, commit, versiones de los datos, parámetros y
    # semillas (incluida la de convergencia). La respuesta, el JSON del
    # CLI y el informe Markdown llevan el mismo manifiesto.
    manifest = {
        "engine_version": APP_VERSION,
        "commit": load_commit_hash(),
        "data_versions": loaded_data_versions(),
        "parameters": {
            "mirnas": found_names,
            "missing_mirnas": missing_mirnas,
            "arm_unverified": [e["mirna"] for e in arm_unverified],
            "mode": mode,
            "years": req.years,
            "offline": req.offline,
            "species": req.species,
            "require_experimental": req.require_experimental,
        },
        "seeds": {"convergence": convergence.get("seed")},
        "reproduce": (
            "curl -s -X POST http://localhost:7860/api/v1/analyze "
            "-H 'Content-Type: application/json' -d '"
            + json.dumps({"mirnas": found_names, "mode": mode,
                          "years": req.years, "offline": req.offline,
                          "species": req.species,
                          "require_experimental": req.require_experimental})
            + "'"
        ),
    }

    # Ronda 3 (v1.48): interruptor opcional "exigir evidencia experimental"
    # (apagado por defecto). Cuando está activo, el consenso solo conserva
    # genes con evidencia en miRTarBase (fuerte o débil) para al menos un
    # miRNA del panel. Es un filtro sobre el consenso, no un cambio de método.
    if req.require_experimental:
        supported = {g for ev in evidence_by_mirna.values() for g, e in ev.items()
                     if e.get("level") in ("strong", "weak", "legacy")}
        before = len(common)
        common = [g for g in common if g in supported]
        logger.info(f"🔬 require_experimental: {before} → {len(common)} genes")

    # ── BUG 5 CORREGIDO: eliminado el fallback automático silencioso ──────────
    # Si el usuario eligió "strict" y la intersección es vacía, informar
    # correctamente en lugar de cambiar el modo sin avisar.
    if not common and mode == "strict":
        logger.warning("⚠️ Intersección estricta vacía. Considere usar modo N-1 o N-2.")
        # Devolver respuesta vacía pero informativa en lugar de cambiar modo
        empty_validation = build_validation_report(
            found_names, [], {}, local_by_mirna, remote_by_mirna,
            loaded_data_versions(), evidence_by_mirna)
        return {
            "engine_version": APP_VERSION,
            "commit": load_commit_hash(),
            "data_versions": loaded_data_versions(),
            "mirtarbase_mode": mirtarbase_mode(),
            "warnings": ([LEGACY_MODE_WARNING] if mirtarbase_mode() == "legacy" else []),
            "common_genes": [],
            "found_mirnas": found_names,
            "missing_mirnas": missing_mirnas,
            "arm_unverified": arm_unverified,
            "arm_verification": arm_status,
            "convergence": convergence,
            "execution_manifest": manifest,
            "single_mirna_warning": (
                "El panel tiene un solo miRNA con datos: no hay convergencia "
                "entre miRNAs; el resultado son todas sus dianas."
                if len(found_names) == 1 else None
            ),
            "source_by_mirna": source_by_mirna,
            "gene_details": {},
            "validation": empty_validation,
            "export": build_export(found_names, [], {}, [], empty_validation, gene_sets),
            "heatmap_plot": None,
            "offline": req.offline,
            "species": req.species,
            "enrichment_truncated": False,
            "enrichment_note": None,
            "research_scope": {
                "cap": RESEARCH_CAP,
                "researched_genes": [],
                "not_evaluated_genes": [],
                "note": "Sin genes comunes: no hubo investigación multi-fuente.",
            },
            "scientific_synthesis": (
                f"No se encontraron genes diana comunes a todos los miRNAs analizados "
                f"({', '.join(found_names)}) en modo de intersección estricta. "
                f"Considere usar el modo N-1 para obtener genes presentes en al menos "
                f"{len(gene_sets)-1} de los {len(gene_sets)} miRNAs."
            ),
            "enrichment": [],
            "enrichment_enrichr": [],
            "venn_plot": None,
            "volcano_plot": None,
            "ppi_plot": None,
            "report_references": []
        }

    # Orden de prioridad (v1.46, H-06): por nº de miRNAs del panel que
    # regulan el gen y por validación en miRTarBase — nunca alfabético.
    # Todas las listas recortadas (investigación, enriquecimiento,
    # gráficos, detalles) usan este orden.
    prioritized_common = prioritize_genes(list(common), gene_sets, remote_by_mirna)
    logger.info(f"✅ {len(prioritized_common)} genes comunes identificados: {prioritized_common[:10]}")

    # ── ENRIQUECIMIENTO FUNCIONAL (REESCRITO PARA ESTABILIDAD) ─────────────────
    enrichment_results = []
    enrichment_enrichr = []
    gene_list_enr = prioritized_common[:300]
    enrichment_truncated = len(prioritized_common) > 300
    
    if len(gene_list_enr) >= 3:
        enrichment_enrichr = []
        async with make_async_client(timeout=45.0) as client:
            try:
                logger.info(f"🧪 Iniciando Enriquecimiento para {len(gene_list_enr)} genes...")
                # Unreleased (2026-11-19): la ruta local es el método y la
                # viva de Enrichr es solo una comprobación exploratoria.
                local_terms, enrichr_terms = await compute_enrichment_routes(
                    gene_list_enr, req.offline, client)
                logger.info(f"📚 Enriquecimiento local: {len(local_terms)} rutas; "
                            f"Enrichr (exploratorio): {len(enrichr_terms)} rutas")
                if local_terms:
                    # Procesamiento en PARALELO para PubMed y Traducción
                    async def process_row(term_data):
                        term = term_data['Term']
                        pval = term_data['Pval']
                        gs_label = term_data['Source']
                        if req.offline:
                            # Sin red: término sin traducir y sin PubMed.
                            pubmed, trans_term = [], term
                        else:
                            try:
                                # Motor experto con reintentos y 429-handling
                                pubmed_task = get_pubmed_for_term(term, req.years, client)
                                trans_task = translate_to_spanish(term, client)
                                pubmed, trans_term = await asyncio.gather(pubmed_task, trans_task)
                            except:
                                pubmed, trans_term = [], term

                        return {
                            "Term": trans_term, "Source": gs_label, "Pval": float(pval),
                            "RawTerm": term,
                            "ScientificDesc": f"Ruta identificada con p-adj={pval:.4f}",
                            "Evidence": pubmed[0] if pubmed else None,
                            "Evidences": [a['id'] for a in pubmed]
                        }

                    tasks = [process_row(t) for t in local_terms]
                    enrichment_results = await asyncio.gather(*tasks)
                    logger.info(f"✅ Enriquecimiento local: {len(enrichment_results)} rutas, {sum(1 for e in enrichment_results if e['Evidence'])} con PubMed")
                else:
                    logger.warning("🔬 No se encontraron rutas significativas en la ruta local.")
                # La ruta viva de Enrichr es una comprobación exploratoria:
                # se lista aparte, sin las llamadas por término (PubMed y
                # traducción), y su p ajustado se rotula como exploratorio.
                enrichment_enrichr = [{
                    "Term": t["Term"], "Source": t["Source"],
                    "Pval": float(t["Pval"]), "RawTerm": t["RawTerm"],
                    "ScientificDesc": "p ajustado de Enrichr (exploratorio)",
                    "Evidence": None, "Evidences": [],
                } for t in enrichr_terms]
            except Exception as e:
                logger.error(f"❌ Error crítico enriquecimiento: {e}")

    # v1.54: columna de especificidad — en qué fracción de paneles
    # aleatorios del mismo tamaño aparece cada término.
    if enrichment_results:
        annotate_enrichment_specificity(enrichment_results, gene_list_enr)
    if enrichment_enrichr:
        annotate_enrichment_specificity(enrichment_enrichr, gene_list_enr)
    # Unreleased: agrupar términos redundantes (mismos genes del panel).
    if enrichment_results:
        annotate_enrichment_groups(enrichment_results, gene_list_enr)
    if enrichment_enrichr:
        annotate_enrichment_groups(enrichment_enrichr, gene_list_enr)

    # ── DATOS REALES PARA LOS GRÁFICOS ────────────────────────────────────────
    # Volcán: conteos de PubMed (top 25 genes de la unión por nº de miRNAs).
    # PPI: interacciones reales STRING-DB entre los genes core.
    # Sin red (offline) ambos quedan en None y la página lo dice:
    # no se inventa ningún gráfico.
    pubmed_counts, ppi_edges = None, None
    if not req.offline:
        union_counts = {}
        for s in gene_sets:
            for g in s:
                union_counts[g] = union_counts.get(g, 0) + 1
        top_genes = sorted(union_counts,
                           key=lambda g: (-union_counts[g], g))[:25]
        async with make_async_client(timeout=20.0) as _viz_client:
            sem = asyncio.Semaphore(5)

            async def _count(g):
                async with sem:
                    return g, await fetch_pubmed_count(g, _viz_client)

            pairs = await asyncio.gather(*[_count(g) for g in top_genes])
            pubmed_counts = dict(pairs)
            # Primera corona: 5 proteínas vecinas que más interactúan con
            # los genes core (la red que la tesis describía como dos módulos).
            ppi_edges = await fetch_string_network(
                prioritized_common[:12], _viz_client, add_nodes=5)

    # ── GRÁFICOS ──────────────────────────────────────────────────────────────
    v_p, volc_p, ppi_p = create_visuals(gene_sets, found_names, prioritized_common,
                                        pubmed_counts=pubmed_counts,
                                        ppi_edges=ppi_edges)

    # ── DETALLES DE GENES ──────────────────────────────────────────────────────
    if req.offline:
        gene_details = {}
    else:
        async with make_async_client(timeout=20.0) as client:
            # Aumentar a 40 para cubrir toda la tabla del informe
            tasks = [get_gene_details(g, client) for g in prioritized_common[:40]]
            details = await asyncio.gather(*tasks)
        gene_details = {g: d for g, d in zip(prioritized_common[:40], details)}

    # ── REFERENCIAS — INVESTIGACIÓN MULTI-FUENTE ──────────────────────────────
    # Para CADA gen core del Venn (prioritized_common[:RESEARCH_CAP]):
    #   - PubMed: búsqueda fresca (gen + contexto biológico de PRESETS + filtro
    #     de años req.years). Cascada: estricto → ampliado → fallback.
    #   - OMIM: enfermedades mendelianas asociadas (NCBI eUtils).
    #   - ClinVar: variantes patogénicas (NCBI eUtils).
    #   - ClinicalTrials.gov: ensayos clínicos activos (API v2, gratis).
    # Resultados cacheados en PERSISTENT_CACHE["gene_research"] por {gene}|{years}
    # para que análisis repetidos sean instantáneos.
    # Para las rutas globales: PMIDs del enrichment_results (ya filtrados por
    # req.years en get_pubmed_for_term).
    report_references, ref_id = [], 1
    gene_research_data = {}  # {gene: {pubmed: [...], omim: [...], clinvar: {...}, trials: [...]}}
    pubmed_candidates = []   # [{pmid, source, fallback_title}] — para metadata fetch en batch
    seen_pmids = set()

    # 1. Investigación multi-fuente PARALELA de los RESEARCH_CAP genes core del Venn
    if req.offline:
        # Sin red: investigación vacía, marcada como offline (no inventada).
        research_results = [
            {"pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""},
             "trials": [], "diseases": [], "diseases_source": "none", "year_window_used": "offline"}
            for g in prioritized_common[:RESEARCH_CAP]
        ]
    else:
        async with make_async_client(timeout=25.0) as client:
            research_tasks = [
                research_gene_complete(g, gene_details.get(g, {}), req.years, client)
                for g in prioritized_common[:RESEARCH_CAP]
            ]
            research_results = await asyncio.gather(*research_tasks, return_exceptions=True)

    for g, res in zip(prioritized_common[:RESEARCH_CAP], research_results):
        if isinstance(res, Exception):
            logger.warning(f"⚠️  Investigación falló para {g}: {res}")
            gene_research_data[g] = {"pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""}, "trials": [], "diseases": [], "diseases_source": "none", "year_window_used": "error"}
            continue
        gene_research_data[g] = res
        # Acumular PMIDs frescos como candidatos para el bloque de Bibliografía
        for ref in res.get("pubmed", []):
            pmid = ref["pmid"]
            if pmid in seen_pmids: continue
            seen_pmids.add(pmid)
            pubmed_candidates.append({
                "pmid": pmid,
                "source": ref["source"],
                "fallback_title": f"{ref['term']} — evidencia clínica para {g}",
                "gene": g,
            })

    # 2. PMIDs de Rutas Globales (Enrichr + PubMed, ya filtrados por req.years
    #    dentro de get_pubmed_for_term — vienen de enrichment_results).
    for item in enrichment_results:
        if len(pubmed_candidates) >= 40: break
        if item.get("Evidence"):
            pmid = item["Evidence"]["id"]
            if pmid in seen_pmids: continue
            seen_pmids.add(pmid)
            pubmed_candidates.append({
                "pmid": pmid,
                "source": item["Source"],
                "fallback_title": f"Ruta funcional: {item['Term']}",
                "gene": None,
            })

    logger.info(f"📚 PMIDs frescos: {len(pubmed_candidates)} candidatos (genes core + rutas, filtro={req.years}y).")

    # 3. Batch fetch de metadatos PubMed reales (título, autores, año, revista)
    resolved_count = 0
    fallback_count = 0
    if pubmed_candidates:
        async with make_async_client(timeout=20.0) as client:
            meta_tasks = [fetch_real_article(c["pmid"], client) for c in pubmed_candidates]
            metas = await asyncio.gather(*meta_tasks, return_exceptions=True)
        for c, meta in zip(pubmed_candidates, metas):
            if isinstance(meta, Exception) or not meta or not meta.get("title"):
                # Fallback: esummary falló por rate-limit; conservamos PMID real.
                report_references.append({
                    "id": ref_id,
                    "title": c["fallback_title"],
                    "authors": "",
                    "year": "",
                    "journal": "",
                    "pmid": c["pmid"],
                    "source": c["source"],
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{c['pmid']}",
                    "metadata_resolved": False,
                    "ref_type": "pubmed",
                })
                fallback_count += 1
            else:
                report_references.append({
                    "id": ref_id,
                    "title": meta["title"],
                    "authors": meta.get("authors", ""),
                    "year": meta.get("year", ""),
                    "journal": meta.get("journal", ""),
                    "pmid": c["pmid"],
                    "source": c["source"],
                    "url": f"https://pubmed.ncbi.nlm.nih.gov/{c['pmid']}",
                    "metadata_resolved": True,
                    "ref_type": "pubmed",
                })
                resolved_count += 1
            ref_id += 1

    # 4. Añadir OMIM, ClinVar, ClinicalTrials como referencias separadas
    seen_omim = set()
    seen_nct = set()
    seen_clinvar_gene = set()
    extras_count = 0
    for g, res in gene_research_data.items():
        # OMIM
        for o in res.get("omim", []):
            if o["omim_id"] in seen_omim or not o.get("title"): continue
            seen_omim.add(o["omim_id"])
            report_references.append({
                "id": ref_id,
                "title": o["title"],
                "authors": "",
                "year": "",
                "journal": "OMIM — Online Mendelian Inheritance in Man",
                "pmid": "",
                "source": f"OMIM:{o['omim_id']} ({g})",
                "url": o["url"],
                "metadata_resolved": True,
                "ref_type": "omim",
            })
            ref_id += 1
            extras_count += 1
        # ClinVar: una entrada por gen (resumen del conteo)
        cv = res.get("clinvar", {})
        if cv.get("count", 0) > 0 and g not in seen_clinvar_gene:
            seen_clinvar_gene.add(g)
            report_references.append({
                "id": ref_id,
                "title": f"Variantes patogénicas reportadas en {g}: {cv['count']} entradas en ClinVar",
                "authors": "",
                "year": "",
                "journal": "ClinVar — NCBI",
                "pmid": "",
                "source": f"ClinVar ({g})",
                "url": cv["url"],
                "metadata_resolved": True,
                "ref_type": "clinvar",
            })
            ref_id += 1
            extras_count += 1
        # ClinicalTrials.gov
        for t in res.get("trials", []):
            if t["nct_id"] in seen_nct or not t.get("title"): continue
            seen_nct.add(t["nct_id"])
            report_references.append({
                "id": ref_id,
                "title": t["title"],
                "authors": "",
                "year": "",
                "journal": f"ClinicalTrials.gov ({t.get('status', '')})",
                "pmid": "",
                "source": f"NCT:{t['nct_id']} ({g})",
                "url": t["url"],
                "metadata_resolved": True,
                "ref_type": "clinicaltrials",
            })
            ref_id += 1
            extras_count += 1

    logger.info(
        f"📚 Referencias finales: {len(report_references)} totales — "
        f"PubMed: {resolved_count} resueltas + {fallback_count} fallback, "
        f"Extras (OMIM/ClinVar/Trials): {extras_count}."
    )

    synthesis_obj = build_synthesis(found_names, prioritized_common, gene_details, enrichment_results, report_references, gene_research_data)

    validation_report = build_validation_report(
        found_names, prioritized_common, gene_research_data,
        local_by_mirna, remote_by_mirna, loaded_data_versions(),
        evidence_by_mirna)

    export_data = build_export(
        found_names, prioritized_common, gene_details, enrichment_results,
        validation_report, gene_sets)

    # Persistir descubrimientos (Traducciones, PubMed, Enriquecimiento)
    save_persistent_cache()

    return {
        "engine_version": APP_VERSION,
        "commit": load_commit_hash(),
        "data_versions": loaded_data_versions(),
        "mirtarbase_mode": mirtarbase_mode(),
        "warnings": ([LEGACY_MODE_WARNING] if mirtarbase_mode() == "legacy" else []),
        "common_genes": prioritized_common,
        "found_mirnas": found_names,
        "missing_mirnas": missing_mirnas,
        "arm_unverified": arm_unverified,
        "arm_verification": arm_status,
        "convergence": convergence,
        "execution_manifest": manifest,
        "single_mirna_warning": (
            "El panel tiene un solo miRNA con datos: no hay convergencia "
            "entre miRNAs; el resultado son todas sus dianas."
            if len(found_names) == 1 else None
        ),
        "source_by_mirna": source_by_mirna,
        "gene_details": gene_details,
        "validation": validation_report,
        "export": export_data,
        "heatmap_plot": create_heatmap(found_names, prioritized_common, gene_sets),
        "offline": req.offline,
        "species": req.species,
        "scientific_synthesis": synthesis_obj["academic"],
        "functional_context": synthesis_obj["functional"],
        "references_text": synthesis_obj["references"],
        "enrichment": enrichment_results,
        "enrichment_enrichr": enrichment_enrichr,
        "enrichment_truncated": enrichment_truncated,
        "enrichment_note": (
            f"Enriquecimiento calculado sobre los 300 genes priorizados de "
            f"{len(prioritized_common)} (orden por nº de miRNAs reguladores "
            f"y validación en miRTarBase)."
            if enrichment_truncated else None
        ),
        "research_scope": {
            "cap": RESEARCH_CAP,
            "researched_genes": prioritized_common[:RESEARCH_CAP],
            "not_evaluated_genes": prioritized_common[RESEARCH_CAP:],
            "note": (
                f"Investigación multi-fuente (PubMed, OMIM, ClinVar, "
                f"ClinicalTrials) limitada a los {RESEARCH_CAP} genes "
                f"priorizados por costo de API — la prioridad es por número "
                f"de miRNAs del panel que regulan el gen y por validación "
                f"en el snapshot local de miRTarBase (Harmonizome, 2026-11-18); "
                f"los {len(prioritized_common) - RESEARCH_CAP} restantes figuran "
                f"como no evaluados, sin puntaje."
                if len(prioritized_common) > RESEARCH_CAP else
                "Todos los genes del panel fueron investigados (dentro del tope por costo de API)."
            ),
        },
        "venn_plot": v_p,
        "volcano_plot": volc_p,
        "ppi_plot": ppi_p,
        "report_references": report_references,
        "gene_research": gene_research_data
    }

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 7860))
    uvicorn.run(app, host="0.0.0.0", port=port)

























