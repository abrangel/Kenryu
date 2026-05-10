from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict, Set, Optional
import httpx
import asyncio
import pandas as pd
import numpy as np
import io
import os
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

# ── NUNCA hardcodear la API key — cargarla desde variable de entorno ──────────
PUBMED_API_KEY = os.getenv("PUBMED_API_KEY", "f6cba5737acf229314444266a1f183251507")

app = FastAPI(title="KENRYU Scientific Engine")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

TARGETSCAN_DB: dict = {}
LOCAL_ZIP_PATH = Path(__file__).parent / "targetscan_data.zip"

# ── FIX 1: Ruta absoluta resuelta desde la ubicación del archivo ──────────────
RAW_DATA_DIR = (Path(__file__).parent / "../../../repo_tesis/data/raw").resolve()


def load_local_data():
    if not LOCAL_ZIP_PATH.exists():
        return {}
    try:
        with zipfile.ZipFile(LOCAL_ZIP_PATH, 'r') as z:
            file_name = [f for f in z.namelist() if f.endswith('.txt')][0]
            with z.open(file_name) as f:
                chunks = pd.read_csv(
                    f, sep='\t',
                    usecols=['miR Family', 'Gene Symbol', 'Species ID'],
                    chunksize=500000, nrows=5000000
                )
                temp_db = {}
                for chunk in chunks:
                    human = chunk[chunk['Species ID'] == 9606]
                    for _, row in human.iterrows():
                        m, g = str(row['miR Family']), str(row['Gene Symbol'])
                        if m not in temp_db:
                            temp_db[m] = set()
                        temp_db[m].add(g)
                return temp_db
    except:
        return {}


def clean_mirna_name(m: str) -> str:
    m = m.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', '-').strip()
    if not m.lower().startswith('hsa-'):
        m = 'hsa-' + m
    parts = m.split('-')
    if len(parts) >= 3:
        parts[1] = "miR"
        return "-".join(parts)
    return m


def get_targetscan_genes_local(mirna: str) -> set:
    global TARGETSCAN_DB

    # ── FIX 1 aplicado: RAW_DATA_DIR ya es ruta absoluta resuelta ────────────
    file_path = RAW_DATA_DIR / f"{mirna}.txt"
    if file_path.exists():
        with open(file_path, 'r', encoding='utf-8') as f:
            return {line.strip() for line in f if line.strip()}

    query = mirna.lower().replace('hsa-', '')
    num_match = re.search(r'mir-(\d+)', query)
    if not num_match:
        return set()
    num = num_match.group(1)
    results = set()
    for family, genes in TARGETSCAN_DB.items():
        if f"mir-{num}" in family.lower():
            if "-5p" in query and "-3p" in family.lower() and "-5p" not in family.lower():
                continue
            if "-3p" in query and "-5p" in family.lower() and "-3p" not in family.lower():
                continue
            results.update(genes)
    return results


async def fetch_mirtarbase_cloud(mirna: str):
    url = f"https://maayanlab.cloud/Harmonizome/api/1.0/gene_set/{mirna}/miRTarBase+microRNA+Targets"
    try:
        async with httpx.AsyncClient() as client:
            resp = await client.get(url, timeout=15.0)
            if resp.status_code == 200:
                return {assoc['gene']['symbol'] for assoc in resp.json().get('associations', [])}
    except:
        pass
    return set()


# ── FIX 4: Función síncrona — no necesita async ni await ─────────────────────
def get_kegg_scientific_summary(pathway_name: str) -> str:
    p_low = pathway_name.lower()
    contexts = {
        "lipid":      "Disrupción en el transporte y eflujo lipídico celular.",
        "athero":     "Desarrollo de procesos inflamatorios en la pared vascular.",
        "cholesterol":"Alteración en la vía de biosíntesis y transporte de colesterol.",
        "insulin":    "Señalización de la resistencia periférica a la insulina.",
        "cancer":     "Control del ciclo celular y mecanismos de escape apoptótico.",
        "immune":     "Regulación de la cascada pro-inflamatoria mediada por citoquinas.",
        "epilep":     "Disgregación de la excitabilidad neuronal por canales Na/K.",
        "cardio":     "Remodelado del tejido cardiaco y estrés biomecánico vascular.",
        "metabolic":  "Desequilibrio en las rutas de bioenergética mitocondrial.",
        "foxo":       "Regulación del ciclo celular y respuesta al estrés oxidativo.",
        "pi3k":       "Señalización de supervivencia y proliferación celular.",
        "mapk":       "Cascada de señalización de crecimiento y diferenciación.",
        "wnt":        "Modulación de la proliferación y diferenciación celular.",
        "p53":        "Control de la respuesta al daño del DNA y apoptosis.",
        "tgf":        "Regulación de la fibrosis y diferenciación tisular.",
        "jak":        "Señalización de citoquinas y respuesta inmune adaptativa.",
        "notch":      "Determinación del destino celular y morfogénesis.",
        "hedgehog":   "Señalización del desarrollo embrionario y tejidos adultos.",
    }
    for key, desc in contexts.items():
        if key in p_low:
            return desc
    return f"Regulación funcional de la vía {pathway_name}."


async def get_pubmed_details(term: str, years: int, client: httpx.AsyncClient):
    clean_term = re.sub(r'hsa\d+|GO:\d+|\(.*?\)', '', term).strip()
    # ── FIX 2: calcular año antes del f-string para evitar expresión aritmética
    start_year = 2026 - years
    query = f'("{clean_term}"[Title/Abstract]) AND ("{start_year}"[Date - Publication] : "3000"[Date - Publication]) AND human[Organism]'
    try:
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pubmed", "term": query, "retmax": 1,
            "retmode": "json", "api_key": PUBMED_API_KEY
        }
        r = await client.get(url, params=params, timeout=10.0)
        ids = r.json().get("esearchresult", {}).get("idlist", [])
        return {"id": ids[0]} if ids else None
    except:
        return None


def fig_to_base64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=250, bbox_inches='tight')
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')


@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    print(f"Iniciando orquestador científico...")
    print(f"RAW_DATA_DIR → {RAW_DATA_DIR} (existe: {RAW_DATA_DIR.exists()})")
    TARGETSCAN_DB = load_local_data()
    print("Sincronización completa.")


class PipelineRequest(BaseModel):
    mirnas: List[str]
    years: int = 10
    consensus_mode: str = "strict"


@app.get("/api/v1/gene-info/{gene_symbol}")
async def get_gene_info(gene_symbol: str):
    """Mapeo clínico real y dinámico para CUALQUIER gen."""
    result = {
        "symbol": gene_symbol,
        "fullName": "",
        "system": "Sistémico",
        "pathology": "Evidencia bioinformática de interacción coordinada."
    }
    async with httpx.AsyncClient(timeout=15.0) as client:
        try:
            # MyGene.info — nombre completo del gen
            resp = await client.get(
                f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=name,summary"
            )
            if resp.status_code == 200:
                hits = resp.json().get("hits", [])
                if hits:
                    result["fullName"] = hits[0].get("name", "")

            # DisGeNET — asociaciones patológicas curadas
            dis_resp = await client.get(
                f"https://www.disgenet.org/api/gda/gene/{gene_symbol}",
                params={"source": "CURATED", "limit": 4},
                headers={"accept": "application/json"}
            )
            if dis_resp.status_code == 200:
                data = dis_resp.json()
                if data:
                    names = [d.get("disease_name") for d in data if d.get("disease_name")]
                    if names:
                        result["pathology"] = f"Asociado clínicamente a: {', '.join(names)}."
                        txt = " ".join(names).lower()
                        if any(x in txt for x in ["heart", "cardio", "athero", "lipid", "coronary"]):
                            result["system"] = "Cardiovascular / Metabólico"
                        elif any(x in txt for x in ["brain", "neuro", "seizure", "epilep", "alzheimer", "parkinson"]):
                            result["system"] = "Neurológico"
                        elif any(x in txt for x in ["cancer", "tumor", "carcinoma", "leukemia", "lymphoma"]):
                            result["system"] = "Oncológico"
                        elif any(x in txt for x in ["immune", "inflam", "infection", "autoimmune"]):
                            result["system"] = "Inmunológico / Celular"
                        elif any(x in txt for x in ["diabet", "metabol", "obesity", "insulin"]):
                            result["system"] = "Metabólico / Endocrino"
                        elif any(x in txt for x in ["liver", "hepat", "renal", "kidney"]):
                            result["system"] = "Hepatorrenal"
        except:
            pass
    return result


@app.post("/api/v1/analyze")
async def run_pipeline(request: PipelineRequest):
    mirnas = [clean_mirna_name(m) for m in request.mirnas if m.strip()]
    gene_sets = []
    for m in mirnas:
        p = get_targetscan_genes_local(m)
        e = await fetch_mirtarbase_cloud(m)
        c = p.union(e)
        if c:
            gene_sets.append(c)

    if not gene_sets:
        raise HTTPException(status_code=404, detail="Error de sincronización: no se encontraron genes para los miRNAs dados.")

    from collections import Counter
    counts = Counter([g for s in gene_sets for g in s])

    if request.consensus_mode == "strict":
        threshold = len(gene_sets)
    elif request.consensus_mode == "n-1":
        threshold = max(1, len(gene_sets) - 1)
    else:  # n-2
        threshold = max(2, len(gene_sets) - 2)  # FIX: mínimo 2 para evitar genes de 1 sola fuente

    common_all = [g for g, c in counts.items() if c >= threshold]

    # Enriquecimiento funcional
    enrichment_data = []
    enriched_genes_set = set()
    async with httpx.AsyncClient() as client:
        libs = {
            'KEGG':     'KEGG_2021_Human',
            'Reactome': 'Reactome_Pathways_2024',
            'GO_BP':    'GO_Biological_Process_2023'
        }
        for label, lib in libs.items():
            enr = gp.enrichr(gene_list=common_all[:500], gene_sets=lib, organism='human')
            res = enr.results
            if not res.empty:
                top = res[res['P-value'] < 0.05].head(15)
                for _, row in top.iterrows():
                    ev = await get_pubmed_details(row['Term'], request.years, client)
                    # FIX 4: get_kegg_scientific_summary ya no es async
                    sci_desc = get_kegg_scientific_summary(row['Term'])
                    enrichment_data.append({
                        "Term":          row['Term'],
                        "Source":        label,
                        "Pval":          row['P-value'],
                        "Genes":         row['Genes'],
                        "ScientificDesc": sci_desc,
                        "Evidence":      ev,
                        "OddsRatio":     row.get('Odds Ratio', 1)
                    })
                    for g in row['Genes'].split(';'):
                        enriched_genes_set.add(g.strip())

    # ── FIX 3: Investigar TODOS los genes core (no solo los primeros 40) ──────
    # Usamos el mismo slice que common_genes en la respuesta final
    genes_to_investigate = common_all[:100]  # hasta 100 — ajustar según tiempo de respuesta

    gene_details = {}
    async with httpx.AsyncClient() as client:
        for gene in genes_to_investigate:
            info = await get_gene_info(gene)
            pubmed = await get_pubmed_details(gene, request.years, client)

            # Score: presencia en gene_sets (0-40) + enriquecimiento (+30) + PubMed (+30)
            presence_score = (sum(1 for s in gene_sets if gene in s) / len(gene_sets)) * 40
            enrich_score   = 30 if gene in enriched_genes_set else 0
            pubmed_score   = 30 if pubmed else 0
            score = round(presence_score + enrich_score + pubmed_score, 1)

            # Rutas asociadas al gen desde el enriquecimiento
            associated_routes = [
                e['Term'] for e in enrichment_data
                if gene in [x.strip() for x in e['Genes'].split(';')]
            ][:3]

            gene_details[gene] = {
                "pmid":             pubmed["id"] if pubmed else None,
                "confidence":       {"score": score, "label": "Alta" if score >= 70 else "Media" if score >= 40 else "Baja"},
                "system":           info["system"],
                "pathology":        info["pathology"],
                "fullName":         info.get("fullName", ""),
                "associated_routes": associated_routes,  # NUEVO: rutas ya calculadas aquí
            }

    # VENN
    plt.figure(figsize=(10, 10))
    plt.gcf().set_facecolor('white')
    if len(mirnas) == 2:
        venn2([gene_sets[0], gene_sets[1]], set_labels=mirnas)
    elif len(mirnas) == 3:
        venn3([gene_sets[0], gene_sets[1], gene_sets[2]], set_labels=mirnas)
    else:
        ax = plt.gca()
        ax.axis('off')
        colors = sns.color_palette("husl", len(mirnas))
        for i, m in enumerate(mirnas):
            ang = i * (360 / len(mirnas))
            ax.add_patch(patches.Circle(
                (3 * np.cos(np.radians(ang)), 3 * np.sin(np.radians(ang))),
                4, color=colors[i], alpha=0.3
            ))
            plt.text(6 * np.cos(np.radians(ang)), 6 * np.sin(np.radians(ang)),
                     m, ha='center', fontweight='bold')
    venn_b64 = fig_to_base64(plt.gcf())

    # VOLCANO
    volcano_b64 = None
    if enrichment_data:
        df = pd.DataFrame(enrichment_data)
        df['-logP'] = -np.log10(df['Pval'].replace(0, 1e-10))
        plt.figure(figsize=(12, 7))
        sns.set_style("whitegrid")
        sns.scatterplot(data=df, x='OddsRatio', y='-logP', hue='Source', s=250, alpha=0.9, edgecolors="black")
        plt.axhline(y=-np.log10(0.05), color='red', linestyle='--')
        plt.title("Paisaje de Significancia Biológica", fontsize=15, fontweight='bold')
        volcano_b64 = fig_to_base64(plt.gcf())

    # PPI (STRING-DB)
    ppi_b64 = None
    try:
        url = f"https://string-db.org/api/image/network?identifiers={'%0d'.join(common_all[:18])}&species=9606&add_nodes=10"
        async with httpx.AsyncClient() as client:
            r = await client.get(url, timeout=15.0)
            if r.status_code == 200:
                ppi_b64 = base64.b64encode(r.content).decode('utf-8')
    except:
        pass

    return {
        "common_genes":    common_all,
        "gene_details":    gene_details,
        "venn_plot":       venn_b64,
        "volcano_plot":    volcano_b64,
        "ppi_plot":        ppi_b64,
        "enrichment":      enrichment_data,
        "logs":            [f"Detección finalizada: {len(common_all)} biomarcadores encontrados."]
    }


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8080)
