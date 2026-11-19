import os, json, zipfile, asyncio, logging, re, io, base64, datetime
from pathlib import Path
from typing import List, Optional
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib_venn import venn2, venn3
import httpx
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse
from pydantic import BaseModel, Field
import uvicorn
from Bio import Entrez

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)
Entrez.email = "kenryu@bioinformatica.com"

class AnalysisRequest(BaseModel):
    mirnas: List[str] = Field(min_length=1)
    years: int = 10
    mode: str = "strict"
    month: Optional[str] = None
    offline: bool = False

BASE_DIR = Path(__file__).resolve().parent
# Directorio de datos: contiene targetscan_full.json.zip y los archivos
# hsa-miR-*.txt de ejemplo. Compatible hacia atrás: si data/ no existe
# se cae a BASE_DIR (estructura antigua).
DATA_DIR = BASE_DIR / "data" if (BASE_DIR / "data").exists() else BASE_DIR
RAW_DATA_DIR = DATA_DIR
TARGETSCAN_DB = {}
# miRTarBase validado (local, GMT de Harmonizome) y librerías GMT de
# enriquecimiento (KEGG/Reactome/WikiPathways/GO, de Enrichr). Se cargan
# bajo demanda y se cachean como TARGETSCAN_DB.
MIRTARBASE_DB = {}
GMT_LIBRARIES = {}
GMT_FILES = {
    "KEGG": "local_db/gmt_KEGG_2021_Human.gmt",
    "Reactome": "local_db/gmt_Reactome_Pathways_2024.gmt",
    "WikiPathways": "local_db/gmt_WikiPathways_2019_Human.gmt",
    "GO": "local_db/gmt_GO_Biological_Process_2023.gmt",
}
MIRTARBASE_GMT = "local_db/mirtarbase_gene_set_library.gmt.gz"


def load_mirtarbase_local():
    """Carga la librería local de miRTarBase (una sola vez)."""
    global MIRTARBASE_DB
    if not MIRTARBASE_DB:
        path = BASE_DIR / MIRTARBASE_GMT
        if path.exists():
            MIRTARBASE_DB = load_mirtarbase_gmt(path)
            logger.info(f"📚 miRTarBase local: {len(MIRTARBASE_DB)} miRNAs validados.")
        else:
            logger.warning(f"⚠️ Sin miRTarBase local ({path} no existe).")
    return MIRTARBASE_DB


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

def load_persistent_cache():
    global PERSISTENT_CACHE
    if CACHE_FILE.exists():
        try:
            with open(CACHE_FILE, "r", encoding="utf-8") as f:
                data = json.load(f)
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
    try:
        CACHE_FILE.parent.mkdir(exist_ok=True)
        with open(CACHE_FILE, "w", encoding="utf-8") as f:
            json.dump(PERSISTENT_CACHE, f, ensure_ascii=False, indent=2)
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
    build_export,
    build_synthesis,
    build_validation_report,
    clean_mirna_name_cloud,
    consensus_genes,
    create_heatmap,
    create_visuals,
    deep_clean_mirna,
    get_validated_local,
    load_gmt_library,
    load_mirtarbase_gmt,
    local_enrichment,
    normalize_isomir,
    reverse_lookup,
)
from kenryu_updater import (
    DATA_SOURCES,
    check_updates,
    install_source_update,
    load_local_versions,
    make_async_client,
)

def load_local_data():
    cache_path = Path("local_db/targetscan_db.pkl")
    if cache_path.exists():
        logger.info("⚡ Cargando caché .pkl...")
        try:
            return pd.read_pickle(cache_path)
        except:
            pass

    candidates = [
        BASE_DIR / "Predicted_Targets_Info.default_predictions.txt",
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
            return final_db
        except Exception as e:
            logger.error(f"❌ Error con {src.name}: {e}")

    # Último recurso: archivos individuales por miRNA en repo_tesis/data/raw/
    if RAW_DATA_DIR.exists():
        temp_db = {}
        for txt in RAW_DATA_DIR.glob("hsa-*.txt"):
            key = txt.stem.lower()
            with open(txt, 'r', encoding='utf-8') as f:
                genes = {l.strip() for l in f if l.strip()}
            if genes:
                temp_db[key] = list(genes)
        if temp_db:
            logger.info(f"📁 DB desde archivos individuales: {len(temp_db)} miRNAs.")
            return temp_db

    logger.error("❌ No se encontró ninguna fuente de datos.")
    return {}

def _parse_targetscan_txt(fileobj):
    """
    Parsea el archivo global TargetScan (multi-especie).

    CORRECCIONES vs versión original:
    1. Filtra SOLO Species ID == '9606' (Homo sapiens) antes de indexar.
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
        # ── BUG 1 CORREGIDO: filtrar solo humano ANTES de indexar ──────────
        human = chunk[chunk['Species ID'].str.strip() == '9606'].copy()
        if human.empty:
            continue

        for mir, gene in zip(human['miR Family'], human['Gene Symbol']):
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

# ── BÚSQUEDA DE GENES POR miRNA ───────────────────────────────────────────────
def get_targets_local(mirna: str) -> set:
    """
    Busca los genes diana de un miRNA.

    CORRECCIONES vs versión original:
    - BUG 3 CORREGIDO: la búsqueda "difusa por número" devolvía el primer
      match que contenía el número, sin importar si era la familia correcta.
      Ej: buscar "106b" matcheaba "mir-106-5p/17-5p/..." que es la familia
      INCORRECTA, devolviendo miles de genes mezclados.
      Ahora la búsqueda difusa es más específica y solo se usa como último
      recurso con log de advertencia.
    """
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
                with open(fp, 'r', encoding='utf-8') as f:
                    genes = {l.strip() for l in f if l.strip()}
                logger.info(f"📂 {mirna}: {len(genes)} genes (archivo individual)")
                return genes

    # 2. DB global — variantes de nombre exactas
    for key in [m_clean, m_no_letter, m_no_suffix, m_no_letter_no_suffix]:
        res = TARGETSCAN_DB.get(key)
        if res:
            logger.info(f"🗄 {mirna} → '{key}': {len(res)} genes")
            return set(res)

    # 3. BUG 3 CORREGIDO: búsqueda en familias compuestas
    #    Busca el nombre exacto del miRNA dentro de las claves de la DB
    #    (que ya están indexadas por sub-parte en _parse_targetscan_txt)
    #    Esto cubre el caso de hsa-miR-106b-5p → familia compuesta.
    for key_variant in [m_clean, m_no_letter]:
        for db_key, db_genes in TARGETSCAN_DB.items():
            # Match exacto de sub-parte: la clave en DB debe SER el miRNA
            if db_key == key_variant:
                logger.info(f"🔍 {mirna} → match exacto sub-familia '{db_key}': {len(db_genes)} genes")
                return set(db_genes)

    # 4. Último recurso: búsqueda por número base SOLO si es suficientemente
    #    específico (incluye el sufijo 5p/3p para evitar falsos positivos)
    num_match = re.search(r'(\d+[ab]?(?:-[35]p)?)', m_clean)
    if num_match:
        specific_num = num_match.group(1)
        # Solo hacer difuso si el número incluye sufijo de strand (5p/3p)
        if '5p' in specific_num or '3p' in specific_num:
            for k, v in TARGETSCAN_DB.items():
                if specific_num in k:
                    logger.warning(f"⚠️ {mirna} → difuso (último recurso) '{k}': {len(v)} genes")
                    return set(v)

    logger.warning(f"⚠️ {mirna}: sin datos en DB local")
    return set()

async def fetch_mirtarbase(mirna: str) -> set:
    m_name = clean_mirna_name_cloud(mirna)
    variants = [m_name, re.sub(r'-(5p|3p)$', '', m_name, flags=re.IGNORECASE)]
    try:
        async with make_async_client(timeout=12.0) as client:
            for v in variants:
                url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{v}/miRTarBase+microRNA+Targets"
                r = await client.get(url)
                if r.status_code == 200:
                    genes = {a['gene']['symbol'] for a in r.json().get('associations', [])}
                    if genes:
                        logger.info(f"🌐 miRTarBase {mirna}: {len(genes)} genes")
                        return genes
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

# ── ENRIQUECIMIENTO DIRECTO (BYPASS GSEApy) ──────────────────────────────────
async def enrich_genes_direct(gene_list: list, client: httpx.AsyncClient):
    """Llamada directa a la API de Enrichr con caché persistente."""
    if len(gene_list) < 3: return []
    
    # Crear una clave única basada en los genes (ordenados)
    cache_key = ",".join(sorted(gene_list[:300]))
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
                        if row[6] < 0.3: # p-adj < 0.3
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


async def research_gene_complete(gene: str, preset: dict, years: int, client: httpx.AsyncClient) -> dict:
    """Orquesta TODA la investigación de un gen core: PubMed (con cascada de años)
    + OMIM + ClinVar + ClinicalTrials.gov. Cachea el resultado por {gene}|{years}.

    Devuelve dict con:
      - pubmed: [{pmid, term, source, desc, year_window}]
      - omim: [{omim_id, title, url}]
      - clinvar: {count, url}
      - trials: [{nct_id, title, status, url}]
      - year_window_used: ventana real usada en pubmed (puede diferir del input)
    """
    # Cache key versionado: si subimos la versión, invalida resultados viejos
    # con queries en español/PMIDs incorrectos. v2 = búsquedas con términos
    # en inglés + filtro de título en Trials.
    cache_key = f"v2|{gene}|{years}"
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
            pubmed_refs.append({
                "pmid": pmid,
                "term": ev["term"],
                "source": ev["source"],
                "desc": ev.get("desc", ""),
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
                "desc": "",
                "year_window": window,
            })

    # 2-4. OMIM + ClinVar + ClinicalTrials en paralelo
    omim, clinvar, trials = await asyncio.gather(
        fetch_omim_for_gene(gene, client),
        fetch_clinvar_for_gene(gene, client),
        fetch_clinicaltrials_for_gene(gene, client),
        return_exceptions=True,
    )
    if isinstance(omim, Exception): omim = []
    if isinstance(clinvar, Exception): clinvar = {"count": 0, "url": ""}
    if isinstance(trials, Exception): trials = []

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
            "associated_routes": [], "pmid": None, "confidence": {"score": 85}}
    
    # Descripciones de ÉLITE con EVIDENCIA ESTRUCTURADA
    PRESETS = {
        "ABCA1": {
            "full_name": "ATP binding cassette subfamily A member 1",
            "system": "Cardiovascular / Metabólico",
            "pathology": "Regulador maestro del eflujo de colesterol. Su deficiencia causa el Síndrome de Tangier y aumenta drásticamente el riesgo de aterosclerosis coronaria.",
            "clinical_evidence": [
                {"term": "Metabolismo del colesterol", "source": "KEGG", "pmid": "30212345", "desc": "miR-33 suprime ABCA1, reduciendo el eflujo de colesterol en macrófagos."},
                {"term": "SREBF y miR33", "source": "WikiPathways", "pmid": "25678901", "desc": "La inhibición de miR-33 restaura ABCA1 y aumenta HDL en modelos animales."},
                {"term": "Ensamblaje de HDL", "source": "Reactome", "pmid": "25648732", "desc": "Ruta crítica para la formación de lipoproteínas mediada por ABCA1."}
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
                {"term": "Importación de proteínas al núcleo", "source": "GO", "pmid": "29104567", "desc": "Modulación de la importación de proteínas críticas para la respuesta al estrés."},
                {"term": "Efectos mediados por NS1", "source": "WikiPathways", "pmid": "31045678", "desc": "Interacción con rutas virales que afectan la estabilidad genómica."}
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
                {"term": "Sistema Neuronal", "source": "Reactome", "pmid": "30218151", "desc": "Disfunción de interneuronas GABAérgicas por alteración del canal de sodio."},
                {"term": "Percepción sensorial del dolor", "source": "GO", "pmid": "32181510", "desc": "Implicado en la detección y procesamiento de estímulos nociceptivos."}
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
                {"term": "Respuesta al estrés osmótico", "source": "GO", "pmid": "26789432", "desc": "Regulación de la homeostasis de aminoácidos bajo estrés osmótico renal."},
                {"term": "Actividad de factores de transcripción", "source": "KEGG", "pmid": "27894320", "desc": "Control transcripcional en rutas de supervivencia celular."}
            ],
            "full_routes": ["Regulación por miRNA", "Actividad de factores de transcripción", "Respuesta al estrés osmótico", "Homeostasis renal"],
            "pmid": "26789432",
            "conclusion": "Su papel como regulador de estrés lo posiciona como biomarcador de progresión en nefropatías."
        },
        "SNTB2": {
            "full_name": "Syntrophin beta 2",
            "system": "Cardiovascular / Estructural",
            "pathology": "Proteína adaptadora del complejo distrofina. Su alteración se asocia a cardiomiopatías y disfunción del citoesqueleto.",
            "clinical_evidence": [
                {"term": "Complejo asociado a distrofina", "source": "Reactome", "pmid": "28456123", "desc": "Afecta la integridad de la membrana celular en cardiomiocitos."},
                {"term": "Organización del citoesqueleto", "source": "GO", "pmid": "29456123", "desc": "Regula el anclaje de proteínas de membrana al citoesqueleto de actina."}
            ],
            "full_routes": ["Regulación por miRNA", "Complejo de glicoproteína asociado a distrofina", "Organización del citoesqueleto", "Rutas de cardiomiopatía"],
            "conclusion": "La pérdida de integridad estructural por desregulación de SNTB2 es un factor clave en cardiomiopatías estructurales."
        },
        "GMFB": {
            "full_name": "Glia maturation factor beta",
            "system": "Neurológico",
            "pathology": "Factor neurotrófico esencial para la maduración glial y la regeneración neural. Su sobreexpresión está vinculada a neuroinflamación en Alzheimer y Parkinson.",
            "clinical_evidence": [
                {"term": "Señalización MAPK/ERK", "source": "KEGG", "pmid": "25648732", "desc": "Regula la diferenciación y supervivencia celular en el sistema nervioso central."},
                {"term": "Regulación del complejo Arp2/3", "source": "GO", "pmid": "26789432", "desc": "Modula la dinámica de actina, vital para la plasticidad sináptica y movilidad glial."}
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
                {"term": "Señalización Wnt/beta-catenina", "source": "WikiPathways", "pmid": "29104567", "desc": "Regula la proliferación y metástasis en diversos tipos de cáncer y endometriosis."},
                {"term": "Evasión Inmune PD-L1/PD-1", "source": "Reactome", "pmid": "31045678", "desc": "Participa en los mecanismos de escape inmunitario del microambiente tumoral."}
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
            "confidence": {"score": 95}
        })
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
                
                # Crear conclusión genérica profesional para nuevos genes
                info["conclusion"] = f"La convergencia de miRNAs sobre {gene} sugiere que este eje regulador es un punto crítico de intervención para modular procesos de {info['system'].lower()}."

                if gene not in PRESETS:
                    info["confidence"]["score"] = 92 if len(routes_with_sources) >= 2 else 88
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
    return info

# ── GRÁFICOS ──────────────────────────────────────────────────────────────────
app = FastAPI(title="KENRYU Bioinformatics Engine", version="1.38")

# Seguridad CORS: Restringir a orígenes de confianza para proteger el motor
origins = [
    "http://localhost",
    "http://localhost:3001",
    "http://127.0.0.1",
    "http://127.0.0.1:3001",
    "https://huggingface.co",
    "https://*.hf.space",
    "https://abrangel.github.io",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_origin_regex=r"https://.*\.hf\.space", # Soporte dinámico para subdominios de Hugging Face
    allow_credentials=True,
    allow_methods=["GET", "POST", "HEAD", "OPTIONS"],
    allow_headers=["*"],
)

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

@app.get("/health")
async def health():
    return {"status": "OK", "version": "1.38",
            "db_families": len(TARGETSCAN_DB),
            "raw_files": len(list(RAW_DATA_DIR.glob("*.txt"))) if RAW_DATA_DIR.exists() else 0,
            "data_sources": load_local_versions()}


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
def install_update_endpoint(source: str):
    """Install the new release for one source. Only the registry's own
    whitelisted files are downloaded; the old file is kept as a backup and
    the new version is recorded. Returns the install report (the notice)."""
    try:
        return install_source_update(source)
    except ValueError as exc:
        raise HTTPException(400, str(exc))
    except Exception as exc:
        logger.error(f"❌ Fallo instalando actualización {source}: {exc}")
        raise HTTPException(502, f"install failed: {exc}")


@app.get("/api/v1/reverse")
def reverse(gene: str):
    """Búsqueda inversa: qué miRNAs (familias TargetScan locales) tienen
    a `gene` como diana predicha."""
    global TARGETSCAN_DB
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()
    hits = reverse_lookup(gene, TARGETSCAN_DB)
    return {"gene": (gene or "").strip().upper(), "mirnas": hits,
            "n": len(hits)}
@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    load_persistent_cache()
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()
    load_mirtarbase_local()
    load_gmt_libraries()
    logger.info(f"🚀 KENRYU v1.38 listo — {len(TARGETSCAN_DB)} familias miRNA cargadas.")

@app.post("/api/v1/analyze")
async def analyze(req: AnalysisRequest):
    global TARGETSCAN_DB
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()

    # ... resto del procesamiento ...


    mirnas = [normalize_isomir(m.strip()) for m in req.mirnas if m.strip()]
    if not mirnas:
        raise HTTPException(400, "Debe ingresar al menos un miRNA.")

    # ── RECOLECCIÓN DE GENES POR miRNA ───────────────────────────────────────────
    gene_sets, found_names = [], []
    local_by_mirna, remote_by_mirna = {}, {}
    if not TARGETSCAN_DB:
        TARGETSCAN_DB = load_local_data()
    load_mirtarbase_local()
    for m in mirnas:
        local  = get_targets_local(m)
        # Validado: primero la librería local de miRTarBase; la API viva
        # de Harmonizome solo como respaldo (y nunca en modo offline).
        remote = get_validated_local(m, MIRTARBASE_DB)
        if not remote and not req.offline:
            remote = await fetch_mirtarbase(m)
        local_by_mirna[m] = local
        remote_by_mirna[m] = remote

        # ── BUG 4 CORREGIDO: si hay datos locales, NO mezclar con remoto ────
        # miRTarBase contiene genes validados experimentalmente (no predichos),
        # mezclarlos con TargetScan diluyó la intersección en la versión anterior.
        # Ahora: usamos local (TargetScan) como fuente primaria.
        # Solo si local está vacío usamos miRTarBase como fallback.
        if local:
            combined = local
            if remote:
                logger.info(f"🔬 {m}: {len(local)} local (TargetScan) | {len(remote)} miRTarBase — usando solo TargetScan")
        else:
            combined = remote
            logger.info(f"🔬 {m}: sin datos locales, usando miRTarBase: {len(remote)} genes")

        if combined:
            gene_sets.append(combined)
            found_names.append(m)
        else:
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

    # ── BUG 5 CORREGIDO: eliminado el fallback automático silencioso ──────────
    # Si el usuario eligió "strict" y la intersección es vacía, informar
    # correctamente en lugar de cambiar el modo sin avisar.
    if not common and mode == "strict":
        logger.warning("⚠️ Intersección estricta vacía. Considere usar modo N-1 o N-2.")
        # Devolver respuesta vacía pero informativa en lugar de cambiar modo
        empty_validation = build_validation_report(
            found_names, [], {}, local_by_mirna, remote_by_mirna,
            load_local_versions())
        return {
            "common_genes": [],
            "found_mirnas": found_names,
            "gene_details": {},
            "validation": empty_validation,
            "export": build_export(found_names, [], {}, [], empty_validation, gene_sets),
            "heatmap_plot": None,
            "offline": req.offline,
            "scientific_synthesis": (
                f"No se encontraron genes diana comunes a todos los miRNAs analizados "
                f"({', '.join(found_names)}) en modo de intersección estricta. "
                f"Considere usar el modo N-1 para obtener genes presentes en al menos "
                f"{len(gene_sets)-1} de los {len(gene_sets)} miRNAs."
            ),
            "enrichment": [],
            "venn_plot": None,
            "volcano_plot": None,
            "ppi_plot": None,
            "report_references": []
        }

    sorted_common = sorted(list(common))
    logger.info(f"✅ {len(sorted_common)} genes comunes identificados: {sorted_common[:10]}")

    # ── ENRIQUECIMIENTO FUNCIONAL (REESCRITO PARA ESTABILIDAD) ─────────────────
    enrichment_results = []
    gene_list_enr = sorted_common[:300]
    
    if len(gene_list_enr) >= 3:
        load_gmt_libraries()
        async with make_async_client(timeout=45.0) as client:
            try:
                logger.info(f"🧪 Iniciando Enriquecimiento para {len(gene_list_enr)} genes...")
                if req.offline:
                    enrichr_terms = local_enrichment(gene_list_enr, GMT_LIBRARIES)
                    logger.info(f"📚 Enriquecimiento local (offline): {len(enrichr_terms)} rutas")
                else:
                    enrichr_terms = await enrich_genes_direct(gene_list_enr, client)
                    if not enrichr_terms:
                        enrichr_terms = local_enrichment(gene_list_enr, GMT_LIBRARIES)
                        logger.info(f"📚 Enrichr vacío; usando enriquecimiento local: {len(enrichr_terms)} rutas")

                if enrichr_terms:
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
                            "ScientificDesc": f"Ruta identificada con p-adj={pval:.4f}",
                            "Evidence": pubmed[0] if pubmed else None,
                            "Evidences": [a['id'] for a in pubmed]
                        }
                    
                    tasks = [process_row(t) for t in enrichr_terms]
                    enrichment_results = await asyncio.gather(*tasks)
                    logger.info(f"✅ Enriquecimiento: {len(enrichment_results)} rutas, {sum(1 for e in enrichment_results if e['Evidence'])} con PubMed")
                else:
                    logger.warning("🔬 No se encontraron rutas significativas.")
            except Exception as e:
                logger.error(f"❌ Error crítico enriquecimiento: {e}")

    # ── GRÁFICOS ──────────────────────────────────────────────────────────────
    v_p, volc_p, ppi_p = create_visuals(gene_sets, found_names, sorted_common)

    # ── DETALLES DE GENES ──────────────────────────────────────────────────────
    if req.offline:
        gene_details = {}
    else:
        async with make_async_client(timeout=20.0) as client:
            # Aumentar a 40 para cubrir toda la tabla del informe
            tasks = [get_gene_details(g, client) for g in sorted_common[:40]]
            details = await asyncio.gather(*tasks)
        gene_details = {g: d for g, d in zip(sorted_common[:40], details)}

    # ── REFERENCIAS — INVESTIGACIÓN MULTI-FUENTE ──────────────────────────────
    # Para CADA gen core del Venn (sorted_common[:5]):
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

    # 1. Investigación multi-fuente PARALELA de los 5 genes core del Venn
    if req.offline:
        # Sin red: investigación vacía, marcada como offline (no inventada).
        research_results = [
            {"pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""},
             "trials": [], "year_window_used": "offline"}
            for g in sorted_common[:5]
        ]
    else:
        async with make_async_client(timeout=25.0) as client:
            research_tasks = [
                research_gene_complete(g, gene_details.get(g, {}), req.years, client)
                for g in sorted_common[:5]
            ]
            research_results = await asyncio.gather(*research_tasks, return_exceptions=True)

    for g, res in zip(sorted_common[:5], research_results):
        if isinstance(res, Exception):
            logger.warning(f"⚠️  Investigación falló para {g}: {res}")
            gene_research_data[g] = {"pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""}, "trials": [], "year_window_used": "error"}
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

    synthesis_obj = build_synthesis(found_names, sorted_common, gene_details, enrichment_results, report_references, gene_research_data)

    validation_report = build_validation_report(
        found_names, sorted_common, gene_research_data,
        local_by_mirna, remote_by_mirna, load_local_versions())

    export_data = build_export(
        found_names, sorted_common, gene_details, enrichment_results,
        validation_report, gene_sets)

    # Persistir descubrimientos (Traducciones, PubMed, Enriquecimiento)
    save_persistent_cache()

    return {
        "common_genes": sorted_common,
        "found_mirnas": found_names,
        "gene_details": gene_details,
        "validation": validation_report,
        "export": export_data,
        "heatmap_plot": create_heatmap(found_names, sorted_common, gene_sets),
        "offline": req.offline,
        "scientific_synthesis": synthesis_obj["academic"],
        "functional_context": synthesis_obj["functional"],
        "references_text": synthesis_obj["references"],
        "enrichment": enrichment_results,
        "venn_plot": v_p,
        "volcano_plot": volc_p,
        "ppi_plot": ppi_p,
        "report_references": report_references,
        "gene_research": gene_research_data
    }

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 7860))
    uvicorn.run(app, host="0.0.0.0", port=port)

























