from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Set, Optional
import httpx
import asyncio
import pandas as pd
import numpy as np
import io
import zipfile
from pathlib import Path
import base64
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from Bio import Entrez
import gseapy as gp
import re
import datetime
from matplotlib_venn import venn2, venn3

Entrez.email = "carva@tesis.com"

app = FastAPI(title="KENRYU SaaS")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

TARGETSCAN_DB = {}
LOCAL_ZIP_PATH = Path(__file__).parent / "targetscan_data.zip"
if not LOCAL_ZIP_PATH.exists():
    LOCAL_ZIP_PATH = Path(__file__).parent.parent.parent / "targetscan_data.zip"

RAW_DATA_DIR = Path(__file__).parent.parent.parent / "repo_tesis" / "data" / "raw"

def load_local_data():
    if not LOCAL_ZIP_PATH.exists(): return {}
    try:
        with zipfile.ZipFile(LOCAL_ZIP_PATH, 'r') as z:
            file_name = [f for f in z.namelist() if f.endswith('.txt')][0]
            with z.open(file_name) as f:
                chunks = pd.read_csv(f, sep='\t', usecols=['miR Family', 'Gene Symbol', 'Species ID'], 
                                   chunksize=500000, nrows=5000000)
                temp_db = {}
                for chunk in chunks:
                    human = chunk[chunk['Species ID'] == 9606]
                    for _, row in human.iterrows():
                        m, g = str(row['miR Family']), str(row['Gene Symbol'])
                        if m not in temp_db: temp_db[m] = set()
                        temp_db[m].add(g)
                return temp_db
    except: return {}

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
    file_path = RAW_DATA_DIR / f"{mirna}.txt"
    if file_path.exists():
        with open(file_path, 'r', encoding='utf-8') as f:
            return {line.strip() for line in f if line.strip()}
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
            resp = await client.get(url, timeout=15.0)
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
        r = await client.get(url, params={"db": "pubmed", "term": query, "retmax": 50, "retmode": "json"})
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return ids if ids else []
    except: return []

def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')

@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    print("Iniciando KENRYU Engine...")
    TARGETSCAN_DB = load_local_data()
    print("Sistema Online.")

class PipelineRequest(BaseModel):
    mirnas: List[str]
    years: int = 10
    month: int = 3

@app.get("/")
async def root():
    return {"status": "online", "name": "KENRYU"}

@app.post("/api/v1/analyze")
async def run_pipeline(request: PipelineRequest):
    global TARGETSCAN_DB
    logs = ["Iniciando orquestación de dianas consenso..."]
    mirnas = [clean_mirna_name(m) for m in request.mirnas if m.strip()]
    
    gene_sets = []
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
            
    common_all = list(set.intersection(*gene_sets))
    if len(common_all) == 0 and len(gene_sets) >= 3:
        from collections import Counter
        counts = Counter([g for s in gene_sets for g in s])
        common_all = [g for g, c in counts.items() if c >= len(gene_sets) - 1]
    
    logs.append(f"Identificados {len(common_all)} biomarcadores consenso.")
    
    # BASE DE CONOCIMIENTO EXPERTA (PARA PRECISIÓN CLÍNICA TOTAL)
    KNOWLEDGE_BASE = {
        "TSC22D2": ["Sistémico / Oncológico", "Regula el ciclo celular; posible implicación en oncogénesis; asociado a alteraciones de homeostasis y transporte molecular."],
        "KPNA3": ["Inmunológico / Celular / Sistémico", "Asociado a susceptibilidad a infecciones; alteraciones en transporte nucleocitoplasmático; disrupción de señalización celular."],
        "ABCA1": ["Cardiovascular / Metabólico", "Riesgo elevado de aterosclerosis; dislipidemia; síndrome de Tangier; alteraciones del perfil lipídico."],
        "SNTB2": ["Cardiovascular / Muscular / Sistémico", "Vinculado a cardiomiopatías, distrofias musculares y defectos en arquitectura del sarcómero."],
        "SCN1A": ["Neurológico", "Asociado a epilepsia, síndrome de Dravet, convulsiones febriles; alteraciones en canales de sodio neuronales."],
        "LDLR": ["Cardiovascular / Metabólico", "Regulación crítica de receptores de LDL; riesgo de hipercolesterolemia familiar y placas ateromatosas."],
        "HMGCR": ["Metabólico / Cardiovascular", "Enzima limitante en la biosíntesis de colesterol; diana farmacológica clave en dislipidemias."],
        "PCSK9": ["Cardiovascular / Sistémico", "Modulador de la homeostasis lipídica; regulación de la degradación del receptor de LDL hepático."],
        "INS": ["Metabólico / Diabetes", "Regulación hormonal de la glucemia; implicación directa en diabetes mellitus tipo 2 y síndrome metabólico."]
    }

    enrichment_data = []
    detected_systems = {} # {Gene: [System, Relevance]}
    
    # Inicializar con Base de Conocimiento
    for g in common_all:
        if g in KNOWLEDGE_BASE:
            detected_systems[g] = KNOWLEDGE_BASE[g]

    async with httpx.AsyncClient() as client:
        try:
            # EXPANSIÓN MÁXIMA DE BASES DE DATOS
            libs = {
                'KEGG': 'KEGG_2021_Human', 
                'Reactome': 'Reactome_Pathways_2024', 
                'GO_BP': 'GO_Biological_Process_2023',
                'GO_MF': 'GO_Molecular_Function_2023',
                'GO_CC': 'GO_Cellular_Component_2023',
                'WikiPathways': 'WikiPathways_2024_Human'
            }
            for label, lib in libs.items():
                enr = gp.enrichr(gene_list=common_all[:500], gene_sets=lib, organism='human')
                res = enr.results
                if not res.empty:
                    # RECOGIDA EXHAUSTIVA: Tomamos todos los significativos (p < 0.05)
                    top = res[res['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
                    for _, row in top.iterrows():
                        evidences = await get_pubmed_details(row['Term'], request.years, client)
                        scientific_desc = await get_kegg_scientific_summary(row['Term'])
                        
                        # Lógica de detección de sistemas EXHAUSTIVA con PRIORIDAD
                        t_low = row['Term'].lower()
                        system = None
                        relevance = None
                        
                        if any(x in t_low for x in ["lipid", "athero", "cholesterol", "cardiac", "heart", "vascular", "artery", "coronary", "blood pressure", "apolipoprotein", "lipoprotein"]):
                            system = "Cardiovascular / Metabólico"
                            relevance = "Riesgo elevado de aterosclerosis, dislipidemia y eventos coronarios agudos."
                        elif any(x in t_low for x in ["insulin", "diabetes", "glucose", "metabolic", "metabolism", "glycogen", "ppar", "adipose", "fatty acid", "hyperglycemia"]):
                            system = "Metabólico / Diabetes"
                            relevance = "Disregulación de la homeostasis de glucosa y resistencia a la insulina periférica."
                        elif any(x in t_low for x in ["cancer", "tumor", "oncogene", "cell cycle", "proliferation", "p53", "apoptosis", "growth factor", "carcinoma", "metastasis", "wnt", "pi3k", "akt"]):
                            system = "Sistémico / Oncológico"
                            relevance = "Activación de rutas proliferativas y evasión de mecanismos de muerte celular programada."
                        elif any(x in t_low for x in ["neuron", "synapse", "brain", "nervous", "epilep", "axon", "neuro", "glial", "parkinson", "alzheimer", "glutamate", "gaba"]):
                            system = "Neurológico"
                            relevance = "Alteración en la excitabilidad neuronal y riesgo de trastornos del espectro epiléptico."
                        elif any(x in t_low for x in ["immune", "cytokine", "inflammation", "t-cell", "b-cell", "interleukin", "tnf", "nf-kappa", "autoimmune", "leukocyte", "nucleocytoplasmic"]):
                            system = "Inmunológico / Celular"
                            relevance = "Respuesta inflamatoria exacerbada y posibles defectos en el transporte molecular."

                        # Si se detectó algo clínico, lo guardamos. Si no, queda como Sistémico por defecto.
                        if not system:
                            system = "Sistémico"
                            relevance = "Alteración en procesos biológicos generales y señalización celular."

                        # Mapear genes a estos sistemas (Damos prioridad a lo clínico sobre lo Sistémico)
                        row_genes = row['Genes'].split(';')
                        for g in row_genes:
                            if g in common_all:
                                # Si el gen NO está en la base de conocimiento y no tiene sistema clínico aún
                                if g not in KNOWLEDGE_BASE:
                                    if g not in detected_systems or (detected_systems[g][0] == "Sistémico" and system != "Sistémico"):
                                        detected_systems[g] = [system, relevance]

                        enrichment_data.append({
                            "Term": row['Term'], "Source": label, "Pval": row['Adjusted P-value'], 
                            "Genes": row['Genes'], "ScientificDesc": scientific_desc,
                            "Evidences": evidences,
                            "OddsRatio": row['Odds Ratio'], "GeneCount": int(row['Overlap'].split('/')[0])
                        })
        except: pass

    # Ordenar por p-val para la síntesis
    enrichment_data.sort(key=lambda x: x['Pval'])
    
    # SÍNTESIS CIENTÍFICA EXTENSA Y DINÁMICA
    top_terms = [e['Term'] for e in enrichment_data[:10]]
    systems_found = list(set([s[0] for s in detected_systems.values()]))
    systems_summary = ", ".join(systems_found[:5])
    
    expert_synthesis = f"La orquestación de dianas consenso ha revelado una arquitectura regulatoria compleja que afecta a {len(common_all)} biomarcadores clave. " \
                       f"Se observa una convergencia funcional significativa sobre las rutas de {', '.join(top_terms[:5])}. " \
                       f"Este perfil molecular sugiere una afectación coordinada en los sistemas {systems_summary}, " \
                       f"lo que implica un riesgo clínico sustancial para el desarrollo de patologías {systems_found[0].split('/')[0].strip().lower()}s. " \
                       f"La represión sinérgica de genes como {', '.join(common_all[:4])} por parte del panel de miRNAs analizado, " \
                       f"representa un nodo crítico en la fisiopatología de enfermedades complejas, " \
                       f"vinculando la regulación epigenética con la manifestación de fenotipos relacionados con {enrichment_data[0]['ScientificDesc'] if enrichment_data else 'la homeostasis sistémica'}. " \
                       f"Además, la identificación de dianas en rutas secundarias de {', '.join(top_terms[5:10])} " \
                       f"refuerza la hipótesis de un control post-transcripcional multi-nivel que predispone a la cronicidad en estados {systems_summary.split(',')[0].lower()}s."

    # Visualizaciones
    plt.figure(figsize=(10, 10)); plt.gcf().set_facecolor('white')
    if len(mirnas) <= 3:
        if len(mirnas) == 2: venn2([gene_sets[0], gene_sets[1]], set_labels=mirnas)
        else: venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=mirnas)
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
        plt.figure(figsize=(14, 9)); sns.set_style("whitegrid")
        sns.scatterplot(data=df_v, x='OddsRatio', y='-log10P', size='GeneCount', sizes=(100, 800), hue='Source', palette='bright', alpha=0.8, edgecolors="black")
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.6)
        plt.xscale('log'); plt.title("Significancia de Rutas Biológicas", fontsize=18, fontweight='bold')
        volcano_b64 = fig_to_base64(plt.gcf())

    ppi_b64 = None
    try:
        url = f"https://string-db.org/api/image/network?identifiers={'%0d'.join(common_all[:12])}&species=9606&add_nodes=10"
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=15.0)
            if resp.status_code == 200: ppi_b64 = base64.b64encode(resp.content).decode('utf-8')
    except: pass

    return {
        "common_genes": common_all, "venn_plot": venn_b64, "ppi_plot": ppi_b64,
        "volcano_plot": volcano_b64, "enrichment": enrichment_data,
        "executive_summary": expert_synthesis, "logs": logs,
        "detected_systems": detected_systems
    }

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)
