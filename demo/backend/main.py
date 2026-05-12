import os
import json
import zipfile
import asyncio
import logging
from pathlib import Path
from typing import List, Optional
import pandas as pd
import httpx
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import uvicorn

# Configuración de Logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = FastAPI(title="KENRYU API - Hugging Face Edition")

# CORS abierto para que tu GitHub Pages se conecte sin problemas
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "local_db"
DATA_DIR.mkdir(exist_ok=True)

LOCAL_ZIP_PATH = BASE_DIR / "targetscan_full.json.zip"
CACHE_PATH = DATA_DIR / "targetscan_db.pkl"

TARGETSCAN_DB = {}

class AnalysisRequest(BaseModel):
    mirnas: List[str]
    years: int = 25
    mode: str = "strict"

def load_local_data():
    if CACHE_PATH.exists():
        logger.info("⚡ Cargando caché ultra-rápida .pkl...")
        return pd.read_pickle(CACHE_PATH).to_dict('index')

    if not LOCAL_ZIP_PATH.exists():
        logger.error(f"❌ ARCHIVO NO ENCONTRADO: {LOCAL_ZIP_PATH}")
        return {}

    logger.info("📦 Procesando archivo pesado por primera vez...")
    with zipfile.ZipFile(LOCAL_ZIP_PATH, 'r') as z:
        filename = z.namelist()[0]
        with z.open(filename) as f:
            data = json.load(f)
            df = pd.DataFrame.from_dict(data, orient='index')
            df.to_pickle(CACHE_PATH)
            return data

@app.get("/")
async def root():
    return {"status": "online", "message": "Kenryu Backend de Cesar Manzo activo en Hugging Face"}

@app.post("/api/v1/analyze")
async def analyze(req: AnalysisRequest):
    logger.info(f"Analizando: {req.mirnas}")
    gene_sets = [set(TARGETSCAN_DB.get(m, [])) for m in req.mirnas]
    if not any(gene_sets): return {"common_genes": [], "logs": "No hay coincidencias."}
    
    common = set.intersection(*[s for s in gene_sets if s])
    sorted_common = sorted(list(common))[:100]
    
    return {
        "common_genes": sorted_common,
        "scientific_synthesis": f"Análisis completado exitosamente para {len(req.mirnas)} miRNAs.",
        "report_references": []
    }

@app.on_event("startup")
async def startup_event():
    global TARGETSCAN_DB
    TARGETSCAN_DB = load_local_data()
    logger.info("✅ Sistema listo para consultas.")

if __name__ == "__main__":
    # Hugging Face requiere el puerto 7860
    uvicorn.run(app, host="0.0.0.0", port=7860)
