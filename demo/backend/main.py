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

# --- LÓGICA DE NEGOCIO IDENTICA AL RESPALDO ---

def load_local_data():
    """Carga masiva idéntica al código perfeccionado de escritorio."""
    cache_path = Path("local_db/targetscan_db.pkl")
    if cache_path.exists():
        logger.info("⚡ Cargando caché .pkl...")
        try: return pd.read_pickle(cache_path)
        except: pass

    txt_file = Path("targetscan_full.json.zip")
    if not txt_file.exists():
        txt_file = Path("Predicted_Targets_Info.default_predictions.txt")

    if txt_file.exists():
        logger.info(f"📄 Procesando archivo masivo: {txt_file.name}")
        try:
            temp_db = {}
            if txt_file.suffix == '.zip':
                with zipfile.ZipFile(txt_file, 'r') as z:
                    file_name = [f for f in z.namelist() if f.endswith('.txt') or f.endswith('.json')][0]
                    with z.open(file_name) as f:
                        # Leer como CSV si es .txt, o JSON si es .json
                        if file_name.endswith('.txt'):
                            chunks = pd.read_csv(f, sep='\t', usecols=['miR Family', 'Gene Symbol', 'Species ID'], dtype=str, chunksize=300000)
                            for chunk in chunks:
                                human = chunk[chunk['Species ID'] == '9606']
                                for mir, gene in zip(human['miR Family'], human['Gene Symbol']):
                                    m_low = str(mir).lower().strip()
                                    if m_low not in temp_db: temp_db[m_low] = set()
                                    temp_db[m_low].add(str(gene))
                        else:
                            data = json.load(f)
                            temp_db = {k.lower(): set(v) for k, v in data.items()}
            else:
                chunks = pd.read_csv(txt_file, sep='\t', usecols=['miR Family', 'Gene Symbol', 'Species ID'], dtype=str, chunksize=300000)
                for chunk in chunks:
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
    m = m.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', '-').strip()
    if not m.lower().startswith('hsa-'): m = 'hsa-' + m
    parts = m.split('-')
    if len(parts) >= 3:
        parts[1] = "miR"
        return "-".join(parts)
    return m

def get_targetscan_genes_local(mirna: str) -> set:
    global TARGETSCAN_DB
    # 1. Archivos RAW específicos de la tesis
    file_path = RAW_DATA_DIR / f"{mirna}.txt"
    if file_path.exists():
        with open(file_path, 'r', encoding='utf-8') as f:
            return {line.strip() for line in f if line.strip()}
    
    # 2. Base de datos masiva
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
    url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{mirna}/miRTarBase+microRNA+Targets"
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=12.0)
            if resp.status_code == 200:
                return {assoc['gene']['symbol'] for assoc in resp.json().get('associations', [])}
    except: pass
    return set()

async def get_kegg_scientific_summary(pathway_name: str) -> str:
    p_low = pathway_name.lower()
    if "lipid" in p_low or "atherosclerosis" in p_low or "cholesterol" in p_low:
        return "Resumen Experto: Alteración en el transporte de lípidos y eflujo de colesterol. Riesgo cardiovascular/metabólico."
    if "insulin" in p_low or "diabetes" in p_low:
        return "Resumen Experto: Regulación de la señalización metabólica de glucosa y sensibilidad a la insulina."
    if "epilep" in p_low or "neuron" in p_low:
        return "Resumen Experto: Disregulación de la excitabilidad neuronal mediada por canales iónicos."
    return f"Proceso biológico significativo identificado en {pathway_name}."

async def get_pubmed_details(term: str, years: int, client: httpx.AsyncClient):
    clean_term = re.sub(r'hsa\d+|GO:\d+|\(.*?\)', '', term).strip()
    query = f'("{clean_term}"[Title/Abstract]) AND ("{2026-years}"[Date - Publication] : "3000"[Date - Publication]) AND human[Organism]'
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        r = await client.get(url, params={"db": "pubmed", "term": query, "retmax": 10, "retmode": "json"})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return ids if ids else []
    except: return []

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=120, bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')

app = FastAPI(title="KENRYU SaaS - Cesar Manzo")
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"])

@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    logger.info("🚀 Iniciando Motor Kenryu Perfeccionado...")
    TARGETSCAN_DB = load_local_data()
    logger.info("✅ Sistema Online.")

@app.get("/")
async def root():
    return {"status": "online", "author": "Cesar Manzo", "db_size": len(TARGETSCAN_DB)}

@app.post("/api/v1/analyze")
async def analyze(request: PipelineRequest):
    global TARGETSCAN_DB
    mirnas = [clean_mirna_name(m) for m in request.mirnas if m.strip()]
    
    gene_sets = []
    logs = ["Iniciando orquestación de dianas consenso..."]
    for m in mirnas:
        genes_pred = get_targetscan_genes_local(m)
        genes_exp = await fetch_mirtarbase_cloud(m)
        combined = genes_pred.union(genes_exp)
        if combined:
            gene_sets.append(combined)
            logs.append(f"Sincronizado {m}: {len(combined)} dianas analizadas.")
        else:
            logs.append(f"Fallo en datos para {m}.")

    if not gene_sets: raise HTTPException(status_code=404, detail="No genomic data found.")
            
    # Lógica de Consenso
    common_all = list(set.intersection(*gene_sets))
    if len(common_all) == 0 and len(gene_sets) >= 3:
        counts = Counter([g for s in gene_sets for g in s])
        common_all = [g for g, c in counts.items() if c >= len(gene_sets) - 1]
    
    # BASE DE CONOCIMIENTO EXPERTA
    KNOWLEDGE_BASE = {
        "TSC22D2": ["Sistémico / Oncológico", "Regula el ciclo celular; posible implicación en oncogénesis."],
        "KPNA3": ["Inmunológico / Celular", "Alteraciones en transporte nucleocitoplasmático."],
        "ABCA1": ["Cardiovascular / Metabólico", "Riesgo elevado de aterosclerosis y dislipidemia."],
        "SNTB2": ["Cardiovascular / Muscular", "Vinculado a cardiomiopatías y defectos en sarcómero."],
        "SCN1A": ["Neurológico", "Asociado a epilepsia y síndrome de Dravet."]
    }

    enrichment_data = []
    detected_systems = {}
    for g in common_all:
        detected_systems[g] = KNOWLEDGE_BASE.get(g, ["Sistémico", "Alteración en procesos biológicos generales."])

    async with httpx.AsyncClient() as client:
        try:
            enr = gp.enrichr(gene_list=common_all[:300], gene_sets='KEGG_2021_Human', organism='human')
            res = enr.results[enr.results['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
            for _, row in res.head(10).iterrows():
                evidences = await get_pubmed_details(row['Term'], request.years, client)
                scientific_desc = await get_kegg_scientific_summary(row['Term'])
                enrichment_data.append({
                    "Term": row['Term'], "Source": "KEGG", "Pval": row['Adjusted P-value'], 
                    "ScientificDesc": scientific_desc, "Evidences": evidences,
                    "OddsRatio": row['Odds Ratio'], "GeneCount": int(row['Overlap'].split('/')[0])
                })
        except: pass

    # Síntesis
    expert_synthesis = f"La orquestación de dianas consenso ha revelado {len(common_all)} biomarcadores clave. " \
                       f"Se observa una convergencia funcional sobre las rutas metabólicas y celulares, " \
                       f"vinculando la regulación de genes como {', '.join(common_all[:4])} con la homeostasis sistémica."

    # Visualizaciones
    plt.figure(figsize=(10, 10)); plt.gcf().set_facecolor('white')
    if len(mirnas) <= 3:
        if len(mirnas) == 2: venn2(gene_sets, set_labels=mirnas)
        elif len(mirnas) == 3: venn3(gene_sets, set_labels=mirnas)
    else:
        ax = plt.gca(); ax.axis('off'); colors = sns.color_palette("husl", len(mirnas))
        for i in range(len(mirnas)):
            angle = i * (360 / len(mirnas))
            ax.add_patch(patches.Circle((3 * np.cos(np.radians(angle)), 3 * np.sin(np.radians(angle))), 4, color=colors[i], alpha=0.3))
            plt.text(6 * np.cos(np.radians(angle)), 6 * np.sin(np.radians(angle)), f"{mirnas[i]}", ha='center', fontweight='bold', fontsize=10)
        plt.text(0, 0, f"DIANAS\n{len(common_all)}", ha='center', va='center', fontweight='bold')
        ax.set_xlim(-10, 10); ax.set_ylim(-10, 10)
    plt.title("Identificación de Dianas Consenso", fontsize=14, fontweight='bold')
    venn_b64 = fig_to_base64(plt.gcf())

    volcano_b64 = None
    if enrichment_data:
        df_v = pd.DataFrame(enrichment_data)
        df_v['-log10P'] = -np.log10(df_v['Pval'].replace(0, 1e-10))
        plt.figure(figsize=(12, 7)); sns.set_style("whitegrid")
        sns.scatterplot(data=df_v, x='OddsRatio', y='-log10P', size='GeneCount', sizes=(100, 500), hue='Source', palette='bright')
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
        plt.title("Significancia de Rutas Biológicas")
        volcano_b64 = fig_to_base64(plt.gcf())

    ppi_b64 = None
    try:
        url = f"https://string-db.org/api/image/network?identifiers={'%0d'.join(common_all[:12])}&species=9606"
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=12.0)
            if resp.status_code == 200: ppi_b64 = base64.b64encode(resp.content).decode('utf-8')
    except: pass

    return {
        "common_genes": common_all, "venn_plot": venn_b64, "ppi_plot": ppi_b64,
        "volcano_plot": volcano_b64, "enrichment": enrichment_data,
        "executive_summary": expert_synthesis, "logs": logs,
        "detected_systems": detected_systems
    }

if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=7860)
