import os
import json
import zipfile
import asyncio
import logging
import re
import io
import base64
from pathlib import Path
from typing import List, Optional, Dict, Set
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib_venn import venn2, venn3
import httpx
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn
import gseapy as gp
from Bio import Entrez
from collections import Counter

# Configuración de Logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

Entrez.email = "kenryu@bioinformatica.com"

# --- CONFIGURACIÓN DE RUTAS (COMO EN EL ESCRITORIO) ---
BASE_DIR = Path(__file__).resolve().parent
RAW_DATA_DIR = BASE_DIR / "repo_tesis" / "data" / "raw"
TARGETSCAN_DB = {}
pubmed_semaphore = asyncio.Semaphore(2)

# --- DEFINICIÓN DE MODELOS ---
class PipelineRequest(BaseModel):
    mirnas: List[str]
    years: int = 10
    month: int = 3
    mode: str = "strict"

# --- LÓGICA DE NEGOCIO IDENTICA AL RESPALDO ---

def load_local_data():
    """Carga masiva idéntica al código perfeccionado."""
    cache_path = Path("local_db/targetscan_db.pkl")
    if cache_path.exists():
        logger.info("⚡ Cargando caché .pkl...")
        try: return pd.read_pickle(cache_path)
        except: pass

    txt_file = Path("Predicted_Targets_Info.default_predictions.txt")
    if txt_file.exists():
        logger.info(f"📄 Procesando archivo masivo: {txt_file.name}")
        try:
            temp_db = {}
            for chunk in pd.read_csv(txt_file, sep='\t', usecols=['miR Family', 'Gene Symbol', 'Species ID'], dtype=str, chunksize=300000):
                human = chunk[chunk['Species ID'] == '9606']
                for mir, gene in zip(human['miR Family'], human['Gene Symbol']):
                    m_low = str(mir).lower().strip()
                    if m_low not in temp_db: temp_db[m_low] = set()
                    temp_db[m_low].add(str(gene))
            final_db = {k: list(v) for k, v in temp_db.items()}
            Path("local_db").mkdir(exist_ok=True)
            pd.to_pickle(final_db, cache_path)
            return final_db
        except Exception as e:
            logger.error(f"❌ Error carga masiva: {e}")
    return {}

def clean_mirna_name(m: str) -> str:
    """Limpieza de nombres idéntica al escritorio."""
    m = m.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', '-').strip()
    if not m.lower().startswith('hsa-'): m = 'hsa-' + m
    parts = m.split('-')
    if len(parts) >= 3:
        # Forzar formato hsa-miR-XXX...
        parts[1] = "miR"
        return "-".join(parts)
    return m

def get_targetscan_genes_local(mirna: str) -> set:
    """Búsqueda en archivos RAW (Dianas reales de la tesis)."""
    # 1. Intentar cargar desde los archivos TXT específicos de la tesis
    file_path = RAW_DATA_DIR / f"{mirna}.txt"
    if file_path.exists():
        logger.info(f"📂 Cargando dianas reales desde archivo: {mirna}.txt")
        with open(file_path, 'r', encoding='utf-8') as f:
            return {line.strip() for line in f if line.strip()}
    
    # 2. Fallback a la base de datos masiva por familia
    query = mirna.lower().replace('hsa-', '')
    num_match = re.search(r'mir-(\d+)', query)
    if not num_match: return set()
    num = num_match.group(1)
    results = set()
    for family, genes in TARGETSCAN_DB.items():
        if f"mir-{num}" in family.lower():
            if "-5p" in query and "-3p" in family.lower() and "-5p" not in family.lower(): continue
            if "-3p" in query and "-5p" in family.lower() and "-3p" not in family.lower(): continue
            results.update(genes)
    return results

async def fetch_mirtarbase_cloud(mirna: str):
    """miRTarBase via Harmonizome."""
    url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{mirna}/miRTarBase+microRNA+Targets"
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=12.0)
            if resp.status_code == 200:
                return {assoc['gene']['symbol'] for assoc in resp.json().get('associations', [])}
    except: pass
    return set()

async def get_detailed_gene_info(gene: str, client: httpx.AsyncClient):
    info = {"full_name": gene, "system": "Multisistémico", "pathology": "Nodo regulador crítico.", "associated_routes": [], "pmid": None}
    try:
        url = f"https://mygene.info/v3/query?q=symbol:{gene}&species=human&fields=name,summary,pathway"
        resp = await client.get(url, timeout=10.0)
        if resp.status_code == 200:
            hits = resp.json().get('hits', [])
            if hits:
                h = hits[0]
                info["full_name"] = h.get('name', gene)
                if h.get('summary'): info["pathology"] = f"Proteína con función en {h.get('summary')[:300]}..."
                pways = h.get('pathway', {})
                routes = []
                for lib in ['kegg', 'reactome']:
                    if lib in pways:
                        d = pways[lib]
                        routes.extend([x['name'] for x in (d if isinstance(d, list) else [d])[:2]])
                info["associated_routes"] = routes

        async with pubmed_semaphore:
            await asyncio.sleep(0.3)
            pm_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={gene}+AND+microRNA&retmax=1&retmode=json"
            pm_resp = await client.get(pm_url)
            if pm_resp.status_code == 200:
                ids = pm_resp.json().get('esearchresult', {}).get('idlist', [])
                if ids: info["pmid"] = ids[0]
    except: pass
    return info

def create_visuals(gene_sets, mirnas, common_all):
    """Motor de visualización idéntico al perfeccionado."""
    plt.style.use('dark_background')
    
    # 1. Venn / Pétalos
    fig1 = plt.figure(figsize=(10, 10))
    if len(mirnas) <= 3:
        if len(mirnas) == 2: venn2(gene_sets, set_labels=mirnas)
        elif len(mirnas) == 3: venn3(gene_sets, set_labels=mirnas)
    else:
        ax = plt.gca(); ax.axis('off')
        colors = sns.color_palette("husl", len(mirnas))
        for i in range(len(mirnas)):
            angle = i * (360 / len(mirnas))
            ax.add_patch(patches.Circle((3 * np.cos(np.radians(angle)), 3 * np.sin(np.radians(angle))), 4, color=colors[i], alpha=0.3))
            plt.text(6 * np.cos(np.radians(angle)), 6 * np.sin(np.radians(angle)), f"{mirnas[i]}", ha='center', color='white', fontweight='bold', fontsize=10)
        plt.text(0, 0, f"DIANAS\n{len(common_all)}", ha='center', va='center', color='#ffcc00', fontweight='bold', fontsize=16)
        ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
    
    buf_v = io.BytesIO()
    plt.savefig(buf_v, format='png', transparent=True, bbox_inches='tight', dpi=120)
    plt.close(fig1)

    # 2. Interactoma Mock de alta calidad
    fig3, ax3 = plt.subplots(figsize=(8, 8))
    ax3.axis('off')
    for i in range(min(len(common_all), 12)):
        angle = 2 * np.pi * i / 12
        plt.plot([0.5, 0.5 + 0.3*np.cos(angle)], [0.5, 0.5 + 0.3*np.sin(angle)], c='#58a6ff', alpha=0.3)
        plt.text(0.5 + 0.35*np.cos(angle), 0.5 + 0.35*np.sin(angle), common_all[i], color='white', ha='center', fontsize=9)
    plt.scatter([0.5], [0.5], c='#ffcc00', s=150, zorder=5)
    
    buf_p = io.BytesIO()
    plt.savefig(buf_p, format='png', transparent=True, bbox_inches='tight', dpi=120)
    plt.close(fig3)

    return base64.b64encode(buf_v.getvalue()).decode(), base64.b64encode(buf_p.getvalue()).decode()

# --- APLICACIÓN ---

app = FastAPI(title="KENRYU SaaS - Cloud Edition")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    logger.info("🚀 Iniciando Motor Kenryu...")
    TARGETSCAN_DB = load_local_data()
    logger.info("✅ Sistema Online.")

@app.post("/api/v1/analyze")
async def analyze(request: PipelineRequest):
    global TARGETSCAN_DB
    mirnas = [clean_mirna_name(m) for m in request.mirnas if m.strip()]
    
    gene_sets = []
    for m in mirnas:
        genes_pred = get_targetscan_genes_local(m)
        genes_exp = await fetch_mirtarbase_cloud(m)
        combined = genes_pred.union(genes_exp)
        if combined:
            gene_sets.append(combined)
    
    if not gene_sets: raise HTTPException(status_code=404, detail="No genomic data found.")
            
    # Lógica de Consenso
    if request.mode == "strict":
        common_all = list(set.intersection(*gene_sets))
    else:
        from collections import Counter
        counts = Counter([g for s in gene_sets for g in s])
        min_m = len(gene_sets) - (1 if request.mode == "n-1" else 2)
        common_all = [g for g, c in counts.items() if c >= max(1, min_m)]
    
    # BASE DE CONOCIMIENTO EXPERTA
    KNOWLEDGE_BASE = {
        "TSC22D2": ["Sistémico / Oncológico", "Regula el ciclo celular; posible implicación en oncogénesis."],
        "KPNA3": ["Inmunológico / Celular", "Alteraciones en transporte nucleocitoplasmático."],
        "ABCA1": ["Cardiovascular / Metabólico", "Riesgo elevado de aterosclerosis y dislipidemia."],
        "SNTB2": ["Cardiovascular / Muscular", "Vinculado a cardiomiopatías y defectos en sarcómero."],
        "SCN1A": ["Neurológico", "Asociado a epilepsia y síndrome de Dravet."],
        "LDLR": ["Cardiovascular / Metabólico", "Regulación crítica de receptores de LDL."],
        "HMGCR": ["Metabólico / Cardiovascular", "Diana farmacológica clave en dislipidemias."],
        "PCSK9": ["Cardiovascular / Sistémico", "Regulación de la degradación del receptor de LDL."],
        "INS": ["Metabólico / Diabetes", "Implicación directa en diabetes mellitus tipo 2."]
    }

    detected_systems = {}
    for g in common_all:
        detected_systems[g] = KNOWLEDGE_BASE.get(g, ["Sistémico", "Alteración en procesos biológicos generales."])

    # Enriquecimiento Real (GSEAPY)
    enrichment_data = []
    try:
        enr = gp.enrichr(gene_list=common_all[:300], gene_sets='KEGG_2021_Human', organism='human')
        res = enr.results[enr.results['Adjusted P-value'] < 0.05]
        for _, row in res.head(10).iterrows():
            enrichment_data.append({
                "Term": row['Term'], "Source": "KEGG", "Pval": row['Adjusted P-value'], 
                "ScientificDesc": f"Ruta biológica significativa.",
                "Evidences": ["PubMed"], "OddsRatio": row['Odds Ratio'], "GeneCount": int(row['Overlap'].split('/')[0])
            })
    except: pass

    # Síntesis Académica
    expert_synthesis = f"El análisis de convergencia molecular ha revelado {len(common_all)} biomarcadores core. " \
                       f"Se observa una afectación coordinada en sistemas metabólicos y celulares críticos, " \
                       f"vinculando la regulación de genes como {', '.join(common_all[:4])} con la homeostasis sistémica."

    venn_b64, ppi_b64 = create_visuals(gene_sets, mirnas, common_all)

    return {
        "common_genes": common_all,
        "found_mirnas": mirnas,
        "executive_summary": expert_synthesis,
        "detected_systems": detected_systems,
        "enrichment": enrichment_data,
        "venn_plot": venn_b64,
        "ppi_plot": ppi_b64,
        "volcano_plot": None
    }

@app.get("/")
async def root():
    return {"status": "online", "author": "Cesar Manzo"}

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=7860)
