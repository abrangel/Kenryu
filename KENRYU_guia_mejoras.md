# KENRYU — Guía Completa de Mejoras (v1.0 → v2.0)

> Diagnóstico exhaustivo y código corregido para cada problema identificado en el sistema bioinformático KENRYU.
> Cubre: correcciones críticas, mejoras de backend, nuevas funciones y exportación Markdown.

---

## Índice

1. [Correcciones críticas](#1-correcciones-críticas)
   - [1.1 PMIDs inventados en el reporte](#11-pmids-inventados-en-el-reporte)
   - [1.2 Mapeo de enfermedades hardcodeado](#12-mapeo-de-enfermedades-hardcodeado)
2. [Mejoras al backend (main.py)](#2-mejoras-al-backend-mainpy)
   - [2.1 Umbral de enriquecimiento relajado](#21-umbral-de-enriquecimiento-relajado)
   - [2.2 Más bases de datos de predicción](#22-más-bases-de-datos-de-predicción)
   - [2.3 Conexión a OMIM y DisGeNET](#23-conexión-a-omim-y-disgenet)
   - [2.4 Score de confianza por gen](#24-score-de-confianza-por-gen)
   - [2.5 Tests mejorados](#25-tests-mejorados)
3. [Mejoras al frontend (script.js)](#3-mejoras-al-frontend-scriptjs)
   - [3.1 Gene pills interactivos](#31-gene-pills-interactivos)
   - [3.2 Modo consenso configurable](#32-modo-consenso-configurable)
4. [Exportación Markdown en el Editor de Reporte](#4-exportación-markdown-en-el-editor-de-reporte)
   - [4.1 Función exportMarkdown()](#41-función-exportmarkdown)
   - [4.2 Botón en la barra de herramientas](#42-botón-en-la-barra-de-herramientas)
   - [4.3 Helpers de conversión HTML→Markdown](#43-helpers-de-conversión-htmlmarkdown)
5. [Resumen de archivos modificados](#5-resumen-de-archivos-modificados)

---

## 1. Correcciones críticas

### 1.1 PMIDs inventados en el reporte

**Problema:** En `script.js`, la función `populateReport()` genera PMIDs con `Math.random()`. Esto produce citas bibliográficas falsas que apuntan a artículos inexistentes o de otros autores. En una tesis esto constituye fabricación de datos.

```javascript
// ❌ CÓDIGO ACTUAL — genera números arbitrarios, NO son PMIDs reales
let pmid1 = Math.floor(Math.random() * 10000000) + 30000000;
let pmid2 = pmid1 + 1234;
detailHtml += `<p>Evidencia: PMID ${pmid1} | PMID ${pmid2}</p>`;
```

**Causa raíz:** El backend ya tiene la función `get_pubmed_details()` que busca PMIDs reales en eutils de NCBI. Sin embargo, el frontend no usa los IDs devueltos en `data.enrichment[i].Evidence` para los bloques de gen, sino que genera números al azar.

**Solución:** Reemplazar completamente la generación aleatoria por el uso de los datos reales devueltos por el backend.

```javascript
// ✅ SOLUCIÓN — usa los PMIDs reales del backend

/**
 * Construye los bloques de "Traducción Patológica" usando evidencia real de PubMed.
 * El backend ya cargó los PMIDs en data.enrichment[i].Evidence.id
 *
 * @param {Object} data - Respuesta completa del endpoint /api/v1/analyze
 */
function buildPathologicalBlocks(data) {
  const container = document.getElementById('rep-pathological-translation');
  if (!container) return;

  // Construimos un mapa gen → lista de PMIDs reales a partir del enriquecimiento
  // Cada término de enriquecimiento tiene un campo "Genes" (ej: "ABCA1;SCN1A")
  // y un campo "Evidence" con el PMID buscado para ese término.
  const genePmidMap = {};

  data.enrichment.forEach(item => {
    if (!item.Evidence?.id) return;          // solo si hay PMID verificado
    const genesInTerm = item.Genes
      .split(/[;,]/)
      .map(g => g.trim())
      .filter(Boolean);

    genesInTerm.forEach(gene => {
      if (!genePmidMap[gene]) genePmidMap[gene] = new Set();
      genePmidMap[gene].add(item.Evidence.id);
    });
  });

  // Mapeo de información base para los genes conocidos
  const GENE_INFO = {
    ABCA1:   { system: "Cardiovascular / Metabólico",
               desc:   "Transportador de eflujo de colesterol. Su represión post-transcripcional reduce la biogénesis de HDL y aumenta el riesgo aterosclerótico." },
    SCN1A:   { system: "Neurológico",
               desc:   "Subunidad alfa del canal Nav1.1. Su disregulación reduce la excitabilidad de interneuronas GABAérgicas, favoreciendo estados epilépticos." },
    SNTB2:   { system: "Cardiovascular / Muscular",
               desc:   "Componente del complejo distrofina-glicoproteína. Mantiene la integridad estructural del sarcómero cardíaco." },
    KPNA3:   { system: "Inmunológico / Celular",
               desc:   "Importina alfa-4. Media el transporte nuclear de proteínas de señalización inmune, incluyendo factores de transcripción NF-κB." },
    TSC22D2: { system: "Sistémico / Oncológico",
               desc:   "Factor de transcripción regulador del ciclo celular. Actúa como supresor tumoral potencial en varios tipos de cáncer." },
  };

  let html = '';
  const genesToShow = data.common_genes.slice(0, 8);  // máx 8 en el reporte

  genesToShow.forEach(gene => {
    const info   = GENE_INFO[gene] ?? { system: "Sistémico", desc: "Análisis funcional pendiente de validación experimental." };
    const pmids  = [...(genePmidMap[gene] ?? [])].slice(0, 3);  // máx 3 PMIDs por gen

    // Si no hay PMIDs reales, indicarlo explícitamente (NO inventar)
    const evidenceHtml = pmids.length > 0
      ? pmids.map(id =>
          `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}" target="_blank"
             style="color:#1f6feb;text-decoration:none;">
             📄 PMID: ${id}
          </a>`
        ).join(' · ')
      : `<span style="color:#888;font-style:italic;">
           Sin evidencia PubMed recuperada para este análisis — verificar manualmente.
         </span>`;

    html += `
      <div style="margin-bottom:20px; border-left:4px solid #1f6feb;
                  padding-left:15px; page-break-inside:avoid; color:black;">
        <p style="margin:0; font-size:13px; color:#1f6feb; font-weight:bold;">
          Marcador identificado:
        </p>
        <p style="margin:0; font-size:15px;"><b>mRNA: ${gene}</b></p>
        <p style="margin:5px 0; font-size:11px;">
          <b>Sistema:</b> ${info.system}<br>
          ${info.desc}
        </p>
        <p style="margin:0; font-size:11px;">
          <b>Evidencia bibliográfica:</b> ${evidenceHtml}
        </p>
      </div>`;
  });

  container.innerHTML = html;
}
```

**Dónde llamarlo:** En `populateReport()`, reemplaza la sección que construye `detailHtml`:

```javascript
// En populateReport(), al final, en lugar del bloque con Math.random():
buildPathologicalBlocks(data);
```

---

### 1.2 Mapeo de enfermedades hardcodeado

**Problema:** La tabla de biomarcadores solo reconoce 5 genes (ABCA1, SCN1A, SNTB2, KPNA3, TSC22D2). Cualquier otro gen que devuelva el análisis aparece como `"Requiere correlación diagnóstica específica."`, sin información real.

```javascript
// ❌ CÓDIGO ACTUAL — tabla fija de 5 genes
let mapping = {
    "ABCA1": ["Cardiovascular / Metabólico", "Enfermedad cardiovascular prematura..."],
    "SCN1A": ["Neurológico", "Trastornos del espectro epiléptico..."],
    // Solo 3 más... todos los demás quedan vacíos
};
```

**Solución en dos capas:**

**Capa A — Backend:** Agregar un endpoint que consulte OMIM/DisGeNET en tiempo real.

```python
# main.py — nuevo endpoint para información de genes
import httpx

@app.get("/api/v1/gene-info/{gene_symbol}")
async def get_gene_info(gene_symbol: str):
    """
    Obtiene información clínica de un gen desde múltiples fuentes:
    - MyGene.info (función, nombre completo, resumen)
    - OpenTargets (enfermedades asociadas, score de evidencia)
    
    Nota: OMIM requiere clave de API. Usamos alternativas open-access.
    """
    result = {
        "symbol": gene_symbol,
        "fullName": "",
        "summary": "",
        "diseases": [],
        "source": []
    }

    async with httpx.AsyncClient(timeout=10.0) as client:

        # ── 1. MyGene.info (sin clave de API, muy completo) ──
        try:
            url = f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=name,summary,entrezgene"
            resp = await client.get(url)
            if resp.status_code == 200:
                hits = resp.json().get("hits", [])
                if hits:
                    result["fullName"] = hits[0].get("name", "")
                    result["summary"]  = hits[0].get("summary", "")[:400]  # primeros 400 chars
                    result["source"].append("MyGene.info")
        except Exception:
            pass

        # ── 2. OpenTargets GraphQL (enfermedades con score > 0.1) ──
        # OpenTargets es open-access, no requiere API key
        try:
            query = """
            query GeneDisease($geneId: String!) {
              target(ensemblId: $geneId) {
                associatedDiseases(page: {index: 0, size: 5}) {
                  rows {
                    disease { name id }
                    score
                  }
                }
              }
            }
            """
            # Primero obtenemos el Ensembl ID desde MyGene
            ensembl_resp = await client.get(
                f"https://mygene.info/v3/query?q=symbol:{gene_symbol}&species=human&fields=ensembl.gene"
            )
            ensembl_id = None
            if ensembl_resp.status_code == 200:
                hits = ensembl_resp.json().get("hits", [])
                if hits:
                    ensembl_data = hits[0].get("ensembl", {})
                    ensembl_id = ensembl_data.get("gene") if isinstance(ensembl_data, dict) else None

            if ensembl_id:
                ot_resp = await client.post(
                    "https://api.platform.opentargets.org/api/v4/graphql",
                    json={"query": query, "variables": {"geneId": ensembl_id}}
                )
                if ot_resp.status_code == 200:
                    rows = (ot_resp.json()
                            .get("data", {})
                            .get("target", {})
                            .get("associatedDiseases", {})
                            .get("rows", []))
                    result["diseases"] = [
                        {"name": r["disease"]["name"], "score": round(r["score"], 3)}
                        for r in rows if r["score"] > 0.1
                    ]
                    result["source"].append("OpenTargets")
        except Exception:
            pass

    return result
```

**Capa B — Frontend:** Consultar el nuevo endpoint al construir la tabla.

```javascript
// script.js — tabla de biomarcadores dinámica

/**
 * Construye la tabla de biomarcadores consultando el backend para cada gen.
 * Muestra un spinner mientras carga y llena la tabla progresivamente.
 *
 * @param {string[]} genes - Lista de genes comunes del análisis
 */
async function buildBiomarkerTable(genes) {
  const container = document.getElementById('rep-biomarker-table-container');
  if (!container) return;

  // Tabla esqueleto mientras carga
  let html = `
    <table style="width:100%; border-collapse:collapse; font-size:11px; color:black; border:1px solid #ddd;">
      <thead style="background:#f8f9fa;">
        <tr>
          <th style="border:1px solid #ddd; padding:8px; text-align:left;">Gen Core</th>
          <th style="border:1px solid #ddd; padding:8px; text-align:left;">Nombre completo</th>
          <th style="border:1px solid #ddd; padding:8px; text-align:left;">Enfermedades asociadas (OpenTargets)</th>
        </tr>
      </thead>
      <tbody id="dynamic-bio-tbody">
        <tr><td colspan="3" style="padding:12px; text-align:center; color:#888;">
          ⏳ Consultando bases de datos clínicas...
        </td></tr>
      </tbody>
    </table>`;
  container.innerHTML = html;

  const tbody = document.getElementById('dynamic-bio-tbody');
  tbody.innerHTML = '';

  // Consultamos en paralelo (máx 10 genes para no saturar)
  const promises = genes.slice(0, 10).map(gene =>
    fetch(`${API_URL}/api/v1/gene-info/${gene}`)
      .then(r => r.ok ? r.json() : null)
      .catch(() => null)
  );

  const results = await Promise.all(promises);

  results.forEach((info, idx) => {
    const gene = genes[idx];
    const fullName = info?.fullName || '—';
    const diseases = info?.diseases?.length > 0
      ? info.diseases
          .map(d => `<span style="background:#e8eef7;color:#1a3a6b;padding:2px 6px;
                                   border-radius:10px;margin:2px;display:inline-block;
                                   font-size:10px;">${d.name}</span>`)
          .join(' ')
      : '<span style="color:#888;font-size:10px;">Sin asociaciones significativas en OpenTargets</span>';

    const row = document.createElement('tr');
    row.innerHTML = `
      <td style="border:1px solid #ddd; padding:8px; font-weight:bold;">${gene}</td>
      <td style="border:1px solid #ddd; padding:8px; font-size:10px; color:#555;">${fullName}</td>
      <td style="border:1px solid #ddd; padding:8px;">${diseases}</td>`;
    tbody.appendChild(row);
  });
}
```

**Dónde llamarlo:** En `populateReport()`:

```javascript
// En populateReport(), reemplaza el bloque de mapping estático:
await buildBiomarkerTable(data.common_genes);
// (necesitas hacer populateReport() async también)
```

---

## 2. Mejoras al backend (main.py)

### 2.1 Umbral de enriquecimiento relajado

**Problema:** Con solo 5 genes de entrada, el p-valor ajustado (corrección FDR/Benjamini-Hochberg) rara vez cae por debajo de 0.05, por lo que el análisis devuelve solo 4 términos en el CSV. Falta mostrar resultados con p-valor nominal y ordenar por Combined Score.

```python
# ❌ CÓDIGO ACTUAL — filtra solo p-adj < 0.05
top = res[res['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
```

**Solución:** Estrategia de captura en tres niveles, tomando siempre los más significativos.

```python
# ✅ SOLUCIÓN — umbral adaptativo y múltiples criterios de ordenamiento

def get_top_enrichment_results(res: pd.DataFrame, label: str, max_results: int = 15) -> list:
    """
    Extrae los términos de enriquecimiento más relevantes con estrategia adaptativa.
    
    Nivel 1 (ideal):    p-adj < 0.05  — estadísticamente sólido con corrección múltiple
    Nivel 2 (aceptable): p-adj < 0.25  — sugestivo, informativo para exploración
    Nivel 3 (exploratorio): p-valor nominal < 0.05 — para conjuntos pequeños de genes
    
    Siempre ordena por Combined Score (= Odds Ratio × -log10(P-valor)),
    que es la métrica más balanceada para conjuntos pequeños.
    
    Args:
        res: DataFrame devuelto por gseapy.enrichr
        label: Nombre de la base de datos (ej: 'KEGG', 'GO_BP')
        max_results: Máximo de términos a devolver por base de datos
        
    Returns:
        Lista de diccionarios con los términos más significativos
    """
    if res.empty:
        return []

    # Normalizar nombres de columnas (gseapy puede variar entre versiones)
    res = res.copy()
    if 'P-value' not in res.columns and 'P_value' in res.columns:
        res.rename(columns={'P_value': 'P-value'}, inplace=True)
    if 'Adjusted P-value' not in res.columns and 'Adjusted_P-value' in res.columns:
        res.rename(columns={'Adjusted_P-value': 'Adjusted P-value'}, inplace=True)

    # Estrategia adaptativa: usar el nivel más estricto que tenga resultados
    for p_col, threshold, level in [
        ('Adjusted P-value', 0.05,  'strict'),
        ('Adjusted P-value', 0.25,  'suggestive'),
        ('P-value',          0.05,  'nominal'),
        ('P-value',          0.20,  'exploratory'),
    ]:
        subset = res[res[p_col] < threshold].copy()
        if not subset.empty:
            # Calcular Combined Score si no existe
            if 'Combined Score' not in subset.columns:
                import numpy as np
                subset['Combined Score'] = (
                    subset['Odds Ratio'] * (-np.log10(subset['P-value'].replace(0, 1e-10)))
                )
            top = subset.sort_values('Combined Score', ascending=False).head(max_results)
            
            results = []
            for _, row in top.iterrows():
                # Extraer número de genes del overlap "3/152" → 3
                overlap_str = str(row.get('Overlap', '0/0'))
                gene_count = int(overlap_str.split('/')[0]) if '/' in overlap_str else 0
                
                results.append({
                    "Term":          row['Term'],
                    "Source":        label,
                    "Pval":          float(row.get('P-value', 1.0)),
                    "PvalAdj":       float(row.get('Adjusted P-value', 1.0)),
                    "EvidenceLevel": level,          # nuevo campo: qué tan sólida es la evidencia
                    "Genes":         row.get('Genes', ''),
                    "OddsRatio":     float(row.get('Odds Ratio', 0)),
                    "CombinedScore": float(row.get('Combined Score', 0)),
                    "GeneCount":     gene_count,
                    "Overlap":       overlap_str,
                })
            return results
    
    return []


# En el endpoint run_pipeline, reemplazar el bloque de enriquecimiento:
# ❌ ANTES:
#   top = res[res['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
#   for _, row in top.iterrows():
#       enrichment_data.append({...})

# ✅ DESPUÉS:
for label, lib in libs.items():
    try:
        enr = gp.enrichr(gene_list=common_all[:500], gene_sets=lib, organism='human')
        term_results = get_top_enrichment_results(enr.results, label, max_results=15)
        
        # Añadir evidencia PubMed a los top 5 por base de datos (para no sobrecargar la API)
        for i, term_item in enumerate(term_results[:5]):
            evidence = await get_pubmed_details(term_item["Term"], request.years, client)
            scientific_desc = await get_kegg_scientific_summary(term_item["Term"])
            term_item["ScientificDesc"] = scientific_desc
            term_item["Evidence"] = evidence
        
        enrichment_data.extend(term_results)
    except Exception as e:
        logs.append(f"Enriquecimiento {label} falló: {str(e)[:100]}")
```

---

### 2.2 Más bases de datos de predicción

**Problema:** Solo se usan TargetScan + miRTarBase. Agregar Miranda, DIANA-microT y TarBase aumenta la confianza estadística al requerir que los genes aparezcan en múltiples predictores independientes.

```python
# main.py — funciones adicionales de predicción

async def fetch_diana_microt(mirna: str, client: httpx.AsyncClient) -> set:
    """
    Consulta DIANA-microT-CDS vía su API REST.
    Devuelve genes con threshold > 0.7 (alta confianza).
    
    Documentación: http://diana.imis.athena-innovation.gr/DianaTools/index.php?r=microT_CDS/index
    """
    # Normalizar: hsa-miR-33a-5p → hsa-miR-33a-5p (ya está bien)
    url = "http://diana.imis.athena-innovation.gr/DianaTools/index.php"
    params = {
        "r":        "microT_CDS/getTargets",
        "mirna":    mirna,
        "threshold": 0.7,
        "format":   "json",
        "genes_only": 1
    }
    try:
        resp = await client.get(url, params=params, timeout=12.0)
        if resp.status_code == 200:
            data = resp.json()
            # La API devuelve {"targets": [{"gene_name": "ABCA1", ...}, ...]}
            return {t["gene_name"] for t in data.get("targets", []) if t.get("gene_name")}
    except Exception:
        pass
    return set()


async def fetch_tarbase(mirna: str, client: httpx.AsyncClient) -> set:
    """
    Consulta TarBase v8 — solo interacciones con evidencia experimental (CLASH, PAR-CLIP, etc.).
    Estas son las dianas con mayor validez biológica.
    
    API: https://dianalab.e-ce.uth.gr/tarbasev8/api
    """
    url = "https://dianalab.e-ce.uth.gr/tarbasev8/api/interactions/"
    params = {
        "mirna":    mirna,
        "organism": "Homo sapiens",
        "type":     "direct",   # solo interacciones directas validadas
        "page_size": 500
    }
    try:
        resp = await client.get(url, params=params, timeout=12.0)
        if resp.status_code == 200:
            results = resp.json().get("results", [])
            return {r["gene_name"] for r in results if r.get("gene_name")}
    except Exception:
        pass
    return set()


async def get_all_predictions(mirna: str, client: httpx.AsyncClient) -> dict:
    """
    Reúne predicciones de todas las fuentes disponibles para un miRNA.
    
    Devuelve un diccionario con los genes por fuente, para calcular
    el score de consenso (cuántas fuentes concuerdan en un gen).
    
    Returns:
        {
          "TargetScan":  {"ABCA1", "SCN1A", ...},
          "miRTarBase":  {"ABCA1", ...},
          "DIANA":       {"SCN1A", ...},
          "TarBase":     {"ABCA1", "SCN1A", ...}
        }
    """
    targetscan = get_targetscan_genes_local(mirna)
    
    # Ejecutar en paralelo las consultas remotas
    mirtarbase_task = fetch_mirtarbase_cloud(mirna)
    diana_task      = fetch_diana_microt(mirna, client)
    tarbase_task    = fetch_tarbase(mirna, client)
    
    mirtarbase, diana, tarbase = await asyncio.gather(
        mirtarbase_task, diana_task, tarbase_task,
        return_exceptions=True  # no falla si una fuente falla
    )
    
    return {
        "TargetScan": targetscan,
        "miRTarBase": mirtarbase if isinstance(mirtarbase, set) else set(),
        "DIANA":      diana      if isinstance(diana, set)      else set(),
        "TarBase":    tarbase    if isinstance(tarbase, set)    else set(),
    }
```

**Integración en el endpoint principal:**

```python
# En run_pipeline, reemplazar el bloque de gene_sets:

gene_sets        = []
source_breakdown = {}   # para mostrar en el reporte cuántas fuentes aportaron cada gen

async with httpx.AsyncClient() as prediction_client:
    for m in mirnas:
        source_data = await get_all_predictions(m, prediction_client)
        
        # Unión de todas las fuentes para este miRNA
        all_genes = set.union(*source_data.values()) if source_data else set()
        
        if all_genes:
            gene_sets.append(all_genes)
            source_breakdown[m] = {src: len(genes) for src, genes in source_data.items()}
            
            total = len(all_genes)
            srcs  = ", ".join(f"{src}:{n}" for src, n in source_breakdown[m].items() if n > 0)
            logs.append(f"Sincronizado {m}: {total} dianas totales ({srcs})")
        else:
            logs.append(f"Sin datos para {m} en ninguna fuente.")
```

---

### 2.3 Conexión a OMIM y DisGeNET

**Problema:** La información de enfermedades está escrita fija en el JS. Integrar OMIM y DisGeNET provee datos reales y actualizados.

> **Nota sobre APIs:** OMIM requiere clave de API gratuita (registro en omim.org). DisGeNET tiene una API abierta y una versión autenticada.

```python
# main.py — nuevos endpoints de información clínica

OMIM_API_KEY = "TU_CLAVE_OMIM_AQUI"  # Obtener en: https://www.omim.org/api

@app.get("/api/v1/diseases/{gene_symbol}")
async def get_diseases_for_gene(gene_symbol: str):
    """
    Obtiene enfermedades asociadas al gen desde DisGeNET (open-access) y OMIM.
    
    DisGeNET no requiere autenticación para consultas básicas de gen → enfermedad.
    OMIM requiere API key gratuita.
    """
    diseases = []

    async with httpx.AsyncClient(timeout=12.0) as client:

        # ── DisGeNET (GDAs — Gene-Disease Associations) ──
        # Documentación: https://www.disgenet.org/api/
        try:
            resp = await client.get(
                f"https://www.disgenet.org/api/gda/gene/{gene_symbol}",
                params={
                    "source": "CURATED",   # solo asociaciones curadas, no inferidas
                    "min_score": 0.3,       # score mínimo de evidencia
                    "limit": 10
                },
                headers={"accept": "application/json"}
            )
            if resp.status_code == 200:
                for gda in resp.json():
                    diseases.append({
                        "name":   gda.get("disease_name", ""),
                        "id":     gda.get("diseaseid", ""),
                        "score":  gda.get("score", 0),
                        "source": "DisGeNET"
                    })
        except Exception:
            pass

        # ── OMIM (si tienes clave de API) ──
        if OMIM_API_KEY and OMIM_API_KEY != "TU_CLAVE_OMIM_AQUI":
            try:
                resp = await client.get(
                    "https://api.omim.org/api/geneMap/search",
                    params={
                        "search":       gene_symbol,
                        "retrieve":     "geneMap",
                        "format":       "json",
                        "apiKey":       OMIM_API_KEY,
                        "start":        0,
                        "limit":        5
                    }
                )
                if resp.status_code == 200:
                    entries = (resp.json()
                               .get("omim", {})
                               .get("searchResponse", {})
                               .get("geneMapList", []))
                    for entry in entries:
                        gm = entry.get("geneMap", {})
                        for pheno in gm.get("phenotypeMapList", []):
                            pm = pheno.get("phenotypeMap", {})
                            diseases.append({
                                "name":   pm.get("phenotype", ""),
                                "id":     f"OMIM:{pm.get('mimNumber', '')}",
                                "score":  1.0,          # OMIM es curado manualmente
                                "source": "OMIM"
                            })
            except Exception:
                pass

    return {
        "gene":     gene_symbol,
        "diseases": diseases,
        "count":    len(diseases)
    }
```

---

### 2.4 Score de confianza por gen

**Problema:** No existe un criterio cuantitativo para priorizar genes. Los usuarios ven todos los genes con el mismo peso visual.

```python
# main.py — función para calcular score de confianza 0-100

def compute_confidence_score(
    gene: str,
    source_counts: dict,       # {"TargetScan": 1, "miRTarBase": 1, "DIANA": 0, "TarBase": 1}
    is_in_enrichment: bool,    # si el gen aparece en algún término enriquecido
    has_pubmed: bool,          # si tiene PMID verificado
    n_mirnas_total: int        # total de miRNAs analizados
) -> dict:
    """
    Calcula un score de confianza para un gen diana.

    Criterios y pesos:
    - Número de fuentes que predicen el gen:  hasta 40 pts (10 por fuente)
    - Aparece en análisis de enriquecimiento: 30 pts
    - Tiene evidencia en PubMed:              20 pts
    - Es intersección de todos los miRNAs:    10 pts extra

    Returns:
        {"score": int, "label": str, "breakdown": dict}
    """
    score = 0
    breakdown = {}

    # Fuentes de predicción (máx 40 pts)
    n_sources   = sum(1 for n in source_counts.values() if n > 0)
    source_pts  = min(n_sources * 10, 40)
    score      += source_pts
    breakdown["fuentes"] = f"{source_pts}/40 ({n_sources} de {len(source_counts)} fuentes)"

    # Enriquecimiento funcional (30 pts)
    enrich_pts  = 30 if is_in_enrichment else 0
    score      += enrich_pts
    breakdown["enriquecimiento"] = f"{enrich_pts}/30"

    # Evidencia PubMed (20 pts)
    pubmed_pts  = 20 if has_pubmed else 0
    score      += pubmed_pts
    breakdown["pubmed"] = f"{pubmed_pts}/20"

    # Intersección completa entre todos los miRNAs (bonus 10 pts)
    # (este bonus lo asigna el caller si el gen está en la intersección estricta)
    breakdown["interseccion_total"] = "pendiente (asignar externamente)"

    # Etiqueta cualitativa
    if score >= 80:
        label = "Alta confianza"
    elif score >= 50:
        label = "Confianza media"
    else:
        label = "Exploratoria"

    return {
        "score":     min(score, 100),
        "label":     label,
        "breakdown": breakdown
    }
```

---

### 2.5 Tests mejorados

**Problema:** El archivo `test_main.py` solo tiene 3 tests, uno de ellos vacío (`test_pipeline_request_validation`).

```python
# demo/backend/tests/test_main.py — versión mejorada

import pytest
import asyncio
from unittest.mock import AsyncMock, patch, MagicMock
from demo.backend.main import (
    clean_mirna_name,
    get_top_enrichment_results,
    compute_confidence_score,
    PipelineRequest
)
import pandas as pd


# ─── Tests de clean_mirna_name ───────────────────────────────────────────────

class TestCleanMirnaName:
    def test_basic_normalization(self):
        """Agrega prefijo hsa- y capitaliza miR."""
        assert clean_mirna_name("mir-33a") == "hsa-miR-33a"
        assert clean_mirna_name("hsa-mir-33a-5p") == "hsa-miR-33a-5p"

    def test_unicode_hyphens(self):
        """Convierte guiones Unicode a ASCII estándar."""
        assert clean_mirna_name(" mir\u201133a ") == "hsa-miR-33a"  # guión no-ruptura
        assert clean_mirna_name("mir\u201333a") == "hsa-miR-33a"    # guión em

    def test_already_correct(self):
        """No modifica nombres ya bien formateados."""
        result = clean_mirna_name("hsa-miR-33a-5p")
        assert result == "hsa-miR-33a-5p"

    def test_whitespace_stripped(self):
        """Elimina espacios al inicio y final."""
        assert clean_mirna_name("  hsa-mir-106b-5p  ") == "hsa-miR-106b-5p"


# ─── Tests de get_top_enrichment_results ─────────────────────────────────────

class TestGetTopEnrichmentResults:
    def _make_df(self, n: int, p_adj: float, p_nom: float = None) -> pd.DataFrame:
        """Crea un DataFrame de resultados de enriquecimiento de prueba."""
        if p_nom is None:
            p_nom = p_adj * 0.5
        return pd.DataFrame({
            "Term":              [f"Term_{i}" for i in range(n)],
            "P-value":           [p_nom] * n,
            "Adjusted P-value":  [p_adj] * n,
            "Odds Ratio":        [10.0] * n,
            "Combined Score":    [50.0] * n,
            "Genes":             ["ABCA1;SCN1A"] * n,
            "Overlap":           ["2/100"] * n,
        })

    def test_returns_empty_for_empty_df(self):
        result = get_top_enrichment_results(pd.DataFrame(), "KEGG")
        assert result == []

    def test_strict_threshold_used_when_available(self):
        """Usa p-adj < 0.05 cuando hay resultados."""
        df = self._make_df(5, p_adj=0.01)
        result = get_top_enrichment_results(df, "KEGG", max_results=10)
        assert len(result) == 5
        assert all(r["EvidenceLevel"] == "strict" for r in result)

    def test_falls_back_to_nominal_pvalue(self):
        """Cuando p-adj > 0.05 pero p-nominal < 0.05, usa nivel nominal."""
        df = self._make_df(3, p_adj=0.8, p_nom=0.02)
        result = get_top_enrichment_results(df, "GO_BP")
        assert len(result) > 0
        assert result[0]["EvidenceLevel"] == "nominal"

    def test_respects_max_results(self):
        """No devuelve más del límite especificado."""
        df = self._make_df(20, p_adj=0.01)
        result = get_top_enrichment_results(df, "Reactome", max_results=5)
        assert len(result) <= 5


# ─── Tests de compute_confidence_score ───────────────────────────────────────

class TestComputeConfidenceScore:
    def test_max_score(self):
        """Gen con todas las fuentes + enriquecimiento + PubMed = ≥ 90."""
        result = compute_confidence_score(
            gene="ABCA1",
            source_counts={"TargetScan": 50, "miRTarBase": 10, "DIANA": 5, "TarBase": 8},
            is_in_enrichment=True,
            has_pubmed=True,
            n_mirnas_total=5
        )
        assert result["score"] >= 90
        assert result["label"] == "Alta confianza"

    def test_low_score_no_evidence(self):
        """Gen sin evidencia tiene score bajo."""
        result = compute_confidence_score(
            gene="UNKNOWN_GENE",
            source_counts={"TargetScan": 0, "miRTarBase": 0, "DIANA": 0, "TarBase": 0},
            is_in_enrichment=False,
            has_pubmed=False,
            n_mirnas_total=5
        )
        assert result["score"] < 20
        assert result["label"] == "Exploratoria"

    def test_score_never_exceeds_100(self):
        """El score está capeado en 100."""
        result = compute_confidence_score(
            gene="X",
            source_counts={"TargetScan": 999, "miRTarBase": 999, "DIANA": 999, "TarBase": 999},
            is_in_enrichment=True,
            has_pubmed=True,
            n_mirnas_total=5
        )
        assert result["score"] <= 100


# ─── Tests de integración del endpoint ───────────────────────────────────────

class TestPipelineEndpoint:
    """
    Tests de integración que mockean las dependencias externas.
    Requieren instalar: pip install pytest-asyncio httpx
    """

    @pytest.mark.asyncio
    async def test_pipeline_validates_empty_mirna_list(self):
        """El endpoint rechaza listas vacías de miRNAs."""
        from fastapi.testclient import TestClient
        from demo.backend.main import app

        client = TestClient(app)
        response = client.post("/api/v1/analyze", json={"mirnas": [], "years": 10})
        # Debe fallar antes de procesar (no hay genes)
        assert response.status_code in [404, 422]

    @pytest.mark.asyncio
    async def test_pipeline_request_model_validation(self):
        """Valida que el modelo PipelineRequest rechaza tipos incorrectos."""
        with pytest.raises(Exception):
            # years debe ser int, no string
            PipelineRequest(mirnas=["hsa-miR-33a-5p"], years="diez")
```

---

## 3. Mejoras al frontend (script.js)

### 3.1 Gene pills interactivos

**Problema:** Los "gene pills" (chips de genes) en el dashboard solo muestran el nombre. No son interactivos.

```javascript
// script.js — gene pills con panel de detalles al clic

/**
 * Renderiza los gene pills con interactividad.
 * Al hacer clic en un gen, se muestra un panel lateral con:
 * - Nombre completo (desde MyGene.info vía backend)
 * - Enfermedades asociadas (OpenTargets)
 * - Score de confianza
 */
function renderCorGenes(genes, enrichmentData) {
  const container = document.getElementById('core-genes-list');
  if (!container) return;

  // Crear mapa de genes que aparecen en enriquecimiento
  const enrichedGenes = new Set();
  enrichmentData?.forEach(item => {
    item.Genes.split(/[;,]/).forEach(g => enrichedGenes.add(g.trim()));
  });

  container.innerHTML = genes.map(gene => {
    const isEnriched = enrichedGenes.has(gene);
    return `
      <span class="gene-pill ${isEnriched ? 'gene-pill--enriched' : ''}"
            onclick="showGeneSidebar('${gene}')"
            title="${isEnriched ? 'Aparece en términos enriquecidos' : 'Gen diana'}">
        ${gene}
        ${isEnriched ? '<span class="pill-dot"></span>' : ''}
      </span>`;
  }).join('');
}

/**
 * Muestra un panel lateral con información detallada del gen.
 * Consulta el backend de forma asíncrona.
 */
async function showGeneSidebar(geneSymbol) {
  // Crear o reutilizar el panel
  let panel = document.getElementById('gene-detail-panel');
  if (!panel) {
    panel = document.createElement('div');
    panel.id = 'gene-detail-panel';
    panel.style.cssText = `
      position: fixed; right: 0; top: 0; height: 100vh; width: 320px;
      background: #1e1e1e; border-left: 1px solid #333; z-index: 200;
      overflow-y: auto; padding: 20px; font-family: 'DM Sans', sans-serif;
      transform: translateX(100%); transition: transform 0.25s ease;
      color: #d4d4d4;`;
    document.body.appendChild(panel);
  }

  // Mostrar panel con estado de carga
  panel.innerHTML = `
    <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
      <h3 style="margin:0; color:#58a6ff; font-size:18px;">${geneSymbol}</h3>
      <button onclick="document.getElementById('gene-detail-panel').style.transform='translateX(100%)'"
              style="background:none; border:none; color:#666; font-size:18px; cursor:pointer;">✕</button>
    </div>
    <div style="color:#666; font-size:12px;">⏳ Consultando bases de datos...</div>`;
  panel.style.transform = 'translateX(0)';

  // Consultar info
  try {
    const [infoResp, diseaseResp] = await Promise.all([
      fetch(`${API_URL}/api/v1/gene-info/${geneSymbol}`),
      fetch(`${API_URL}/api/v1/diseases/${geneSymbol}`)
    ]);

    const info     = infoResp.ok     ? await infoResp.json()     : {};
    const diseases = diseaseResp.ok  ? await diseaseResp.json()  : { diseases: [] };

    const diseasesHtml = diseases.diseases?.length > 0
      ? diseases.diseases.slice(0, 6).map(d => `
          <div style="padding:6px 8px; background:#252525; border-radius:4px;
                      margin-bottom:4px; font-size:11px;">
            <span style="color:#e8eef7;">${d.name}</span>
            <span style="float:right; color:#555; font-size:10px;">${d.source}</span>
          </div>`).join('')
      : '<p style="color:#555; font-size:11px;">Sin asociaciones curadas disponibles</p>';

    panel.innerHTML = `
      <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
        <h3 style="margin:0; color:#58a6ff; font-size:18px;">${geneSymbol}</h3>
        <button onclick="document.getElementById('gene-detail-panel').style.transform='translateX(100%)'"
                style="background:none; border:none; color:#666; font-size:18px; cursor:pointer;">✕</button>
      </div>

      ${info.fullName ? `<p style="color:#aaa; font-size:12px; margin:0 0 12px;">${info.fullName}</p>` : ''}

      ${info.summary ? `
        <div style="margin-bottom:16px;">
          <div style="font-size:10px; color:#555; text-transform:uppercase; margin-bottom:6px;">Función</div>
          <p style="font-size:11px; color:#ccc; line-height:1.5; margin:0;">${info.summary}</p>
        </div>` : ''}

      <div style="margin-bottom:16px;">
        <div style="font-size:10px; color:#555; text-transform:uppercase; margin-bottom:8px;">
          Enfermedades asociadas
        </div>
        ${diseasesHtml}
      </div>

      <div style="font-size:10px; color:#444; margin-top:16px; padding-top:12px; border-top:1px solid #2a2a2a;">
        Fuentes: MyGene.info · DisGeNET · OpenTargets
      </div>`;
  } catch (err) {
    panel.innerHTML += `<p style="color:#ff7b72; font-size:11px;">Error al cargar: ${err.message}</p>`;
  }
}

// CSS adicional para gene pills (agregar al <style> del index.html)
const genePillStyle = `
  .gene-pill {
    display: inline-flex; align-items: center; gap: 4px;
    padding: 4px 10px; background: #1a3a6b22; color: #58a6ff;
    border: 1px solid #1a3a6b55; border-radius: 12px;
    font-size: 12px; cursor: pointer; transition: all 0.15s;
    position: relative;
  }
  .gene-pill:hover { background: #1a3a6b44; transform: translateY(-1px); }
  .gene-pill--enriched { border-color: #0f5e3a55; background: #0f5e3a22; color: #56d364; }
  .pill-dot {
    width: 5px; height: 5px; background: #56d364;
    border-radius: 50%; display: inline-block;
  }
`;
```

---

### 3.2 Modo consenso configurable

**Problema:** La intersección actual es estricta: un gen debe aparecer en TODOS los miRNAs. Con 5 miRNAs esto reduce mucho la lista. Agregar modo "≥N de M".

```javascript
// Agregar al HTML en la sección de configuración:
/*
<div>
  <label style="font-size:0.7rem; color:var(--text-dim);">MODO CONSENSO</label>
  <select id="consensus-mode" style="width:100%; background:var(--bg-panel); color:white; border:1px solid var(--border); padding:5px;">
    <option value="strict">Estricto (todos los miRNAs)</option>
    <option value="4of5" selected>Alto (≥4 de 5 miRNAs)</option>
    <option value="3of5">Moderado (≥3 de 5 miRNAs)</option>
    <option value="2of5">Amplio (≥2 de 5 miRNAs)</option>
  </select>
</div>
*/

// En analyzePro():
const consensusMode = document.getElementById('consensus-mode')?.value || 'strict';
const body = {
  mirnas,
  years:          parseInt(years),
  month:          parseInt(month),
  consensus_mode: consensusMode   // nuevo campo
};
```

```python
# main.py — soporte para consensus_mode en PipelineRequest

class PipelineRequest(BaseModel):
    mirnas:         List[str]
    years:          int = 10
    month:          int = 3
    consensus_mode: str = "strict"   # "strict" | "4of5" | "3of5" | "2of5"

# En run_pipeline, después de obtener gene_sets:
def apply_consensus(gene_sets: list, mode: str) -> list:
    """
    Aplica el modo de consenso seleccionado.
    
    "strict" → intersección de todos
    "4of5"   → aparece en al menos 4 conjuntos
    "3of5"   → aparece en al menos 3
    "2of5"   → aparece en al menos 2
    """
    from collections import Counter
    
    if mode == "strict" or len(gene_sets) <= 2:
        return list(set.intersection(*gene_sets)) if gene_sets else []
    
    threshold_map = {"4of5": 4, "3of5": 3, "2of5": 2}
    threshold = threshold_map.get(mode, len(gene_sets))
    
    counts = Counter(gene for s in gene_sets for gene in s)
    return [gene for gene, count in counts.items() if count >= threshold]

# Reemplaza la lógica actual de common_all:
common_all = apply_consensus(gene_sets, request.consensus_mode)
logs.append(f"Consenso ({request.consensus_mode}): {len(common_all)} biomarcadores.")
```

---

## 4. Exportación Markdown en el Editor de Reporte

Esta es la mejora principal solicitada: el editor actualmente exporta a PDF y TXT. Se agrega exportación a **Markdown** (.md) con:

- Frontmatter YAML (metadatos del informe)
- Secciones como encabezados `##`
- Bloques de genes con sub-encabezados `###`
- Tablas convertidas a sintaxis Markdown
- Referencias PubMed como links Markdown
- Imágenes como notas de referencia

### 4.1 Función exportMarkdown()

```javascript
// ─────────────────────────────────────────────────────────────────────────────
// EXPORTACIÓN MARKDOWN — Agregar al bloque <script> de kenryu_report_editor.html
// Reemplaza / complementa la función exportText() existente
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Exporta el informe activo como archivo Markdown (.md).
 *
 * Estructura del archivo generado:
 *   ── frontmatter YAML (metadatos)
 *   ── encabezado del informe
 *   ── una sección ## por cada .report-section
 *      ── bloques de genes como sub-secciones ###
 *      ── tablas convertidas a Markdown GFM
 *      ── citas PubMed como links [PMID: XXXXX](url)
 *   ── pie de página
 */
function exportMarkdown() {
  // ── 1. Recolectar metadatos del encabezado ──────────────────────────────
  const reportId   = document.getElementById('meta-id')?.textContent?.trim()   || 'KR-DRAFT';
  const inst       = document.getElementById('meta-inst')?.textContent?.trim() || 'KENRYU';
  const version    = document.getElementById('meta-ver')?.textContent?.trim()  || 'v1.0';
  const dateLabel  = document.getElementById('meta-date')?.textContent?.trim() || new Date().toLocaleDateString('es-ES');
  const isoDate    = new Date().toISOString().split('T')[0];

  // ── 2. Frontmatter YAML ─────────────────────────────────────────────────
  let md = `---
title: "Informe Bioinformático KENRYU"
id: "${reportId}"
institucion: "${inst}"
version: "${version}"
fecha: "${dateLabel}"
generado: "${isoDate}"
sistema: "KENRYU Bioinformatics Engine"
formato: "Markdown / GitHub Flavored Markdown (GFM)"
---

`;

  // ── 3. Encabezado principal ─────────────────────────────────────────────
  const titleEl    = document.querySelector('.report-title');
  const subtitleEl = document.querySelector('.report-subtitle');

  md += `# ${titleEl?.textContent?.trim() || 'Informe Bioinformático'}\n\n`;
  if (subtitleEl) {
    md += `_${subtitleEl.textContent.trim()}_\n\n`;
  }
  md += `**ID:** ${reportId} · **Institución:** ${inst} · **Versión:** ${version} · **Fecha:** ${dateLabel}\n\n`;
  md += `---\n\n`;

  // ── 4. Procesar cada sección del documento ──────────────────────────────
  const sections = document.querySelectorAll('.report-section');

  sections.forEach((sec) => {
    // 4a. Encabezado de sección
    const headingEl = sec.querySelector('.section-heading');
    if (headingEl) {
      // Limpiar el número romano del span.s-num si existe
      const sNum = headingEl.querySelector('.s-num');
      const headingText = sNum
        ? headingEl.textContent.replace(sNum.textContent, '').trim()
        : headingEl.textContent.trim();
      md += `## ${headingText}\n\n`;
    }

    // 4b. Bloques de gen (elemento .gene-block)
    const geneBlocks = sec.querySelectorAll('.gene-block');
    geneBlocks.forEach(gb => {
      const geneName  = gb.querySelector('.gene-block-name')?.textContent?.trim()  || '';
      const geneType  = gb.querySelector('.gene-block-type')?.textContent?.trim()  || '';
      const geneText  = gb.querySelector('.gene-block-text');
      const genePmids = gb.querySelector('.gene-pmid');

      if (geneName) md += `### ${geneName}\n\n`;
      if (geneType) md += `**Función:** ${geneType}\n\n`;

      if (geneText) {
        md += innerHtmlToMarkdown(geneText.innerHTML) + '\n\n';
      }

      if (genePmids) {
        // Convertir los <a href="..."> de PMIDs a links Markdown
        const pmidLinks = genePmids.querySelectorAll('a');
        if (pmidLinks.length > 0) {
          const linksStr = Array.from(pmidLinks)
            .map(a => `[${a.textContent.trim()}](${a.href})`)
            .join(' · ');
          md += `> **Evidencia bibliográfica:** ${linksStr}\n\n`;
        } else {
          md += `> ${genePmids.textContent.trim()}\n\n`;
        }
      }
    });

    // 4c. Bloques de texto editable genérico
    const editableBlocks = sec.querySelectorAll(
      '.editable-block, .quote-block'
    );
    editableBlocks.forEach(block => {
      const content = innerHtmlToMarkdown(block.innerHTML).trim();
      if (content) {
        const isQuote = block.classList.contains('quote-block');
        md += isQuote
          ? `> ${content}\n\n`
          : `${content}\n\n`;
      }
    });

    // 4d. Tablas HTML → Markdown GFM
    const tables = sec.querySelectorAll('table');
    tables.forEach(table => {
      md += tableToMarkdown(table) + '\n\n';
    });

    // 4e. Notas (note-block)
    const noteBlock = sec.querySelector('.note-block');
    if (noteBlock) {
      const noteText = noteBlock.querySelector('.editable-block')?.innerText?.trim();
      if (noteText) {
        md += `> ✏️ **Nota del investigador:** ${noteText}\n\n`;
      }
    }

    // 4f. Divisores
    if (sec.querySelector('hr')) {
      md += `---\n\n`;
    }
  });

  // ── 5. Nota sobre imágenes (base64 no se incluyen en MD) ────────────────
  md += `---\n\n`;
  md += `## Notas sobre este documento\n\n`;
  md += `- Este archivo fue generado automáticamente por **KENRYU Bioinformatics Engine**.\n`;
  md += `- Las visualizaciones (diagrama de Venn, volcano plot, red PPI) se exportan por separado en formato PDF.\n`;
  md += `- Para la versión completa con gráficos, use la opción **Exportar PDF**.\n`;
  md += `- Los PMIDs enlazados apuntan directamente a PubMed.\n\n`;
  md += `_Documento generado el ${new Date().toLocaleString('es-ES')} por KENRYU v2.0_\n`;

  // ── 6. Descargar el archivo ─────────────────────────────────────────────
  const blob = new Blob([md], { type: 'text/markdown;charset=utf-8' });
  const a    = document.createElement('a');
  a.href     = URL.createObjectURL(blob);
  a.download = `Informe_KENRYU_${reportId}.md`;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(a.href);
}
```

### 4.2 Botón en la barra de herramientas

```html
<!-- En kenryu_report_editor.html, dentro de #toolbar,
     agregar DESPUÉS del botón de PDF y ANTES del botón DOCX -->

<!-- Botón Exportar Markdown -->
<button class="tb-export tb-export-md" onclick="exportMarkdown()" title="Ctrl+M">
  <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
    <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
    <polyline points="14 2 14 8 20 8"/>
    <line x1="8" y1="13" x2="16" y2="13"/>
    <line x1="8" y1="17" x2="16" y2="17"/>
  </svg>
  ↓ MD
</button>
```

```css
/* Agregar a la sección <style> de kenryu_report_editor.html */

.tb-export-md {
  background: #4a1a6b;        /* púrpura — color tradicional de Markdown */
}
.tb-export-md:hover {
  background: #3a0f55;
}
```

```javascript
// Agregar al listener de keydown en el DOMContentLoaded:
if ((e.ctrlKey || e.metaKey) && e.key === 'm') {
  e.preventDefault();
  exportMarkdown();
}
```

### 4.3 Helpers de conversión HTML→Markdown

```javascript
// ─────────────────────────────────────────────────────────────────────────────
// HELPERS — Agregar también al bloque <script> de kenryu_report_editor.html
// ─────────────────────────────────────────────────────────────────────────────

/**
 * Convierte innerHTML de un elemento a texto Markdown.
 *
 * Maneja los casos más comunes en el editor:
 *   <b>/<strong> → **texto**
 *   <i>/<em>     → _texto_
 *   <a href>     → [texto](url)
 *   <br>         → salto de línea
 *   <p>          → párrafo
 *   <code>       → `código`
 *   entidades    → caracteres reales (&lt; → <, etc.)
 *
 * @param {string} html - innerHTML del elemento
 * @returns {string} Texto en formato Markdown
 */
function innerHtmlToMarkdown(html) {
  if (!html) return '';

  return html
    // Párrafos: <p>texto</p> → texto\n\n
    .replace(/<p[^>]*>([\s\S]*?)<\/p>/gi, (_, content) =>
      stripTags(content).trim() + '\n\n'
    )
    // Negritas
    .replace(/<(b|strong)[^>]*>([\s\S]*?)<\/(b|strong)>/gi, (_, _t, c) =>
      `**${stripTags(c).trim()}**`
    )
    // Cursivas
    .replace(/<(i|em)[^>]*>([\s\S]*?)<\/(i|em)>/gi, (_, _t, c) =>
      `_${stripTags(c).trim()}_`
    )
    // Código inline
    .replace(/<code[^>]*>([\s\S]*?)<\/code>/gi, (_, c) =>
      `\`${stripTags(c).trim()}\``
    )
    // Links
    .replace(/<a[^>]+href=["']([^"']+)["'][^>]*>([\s\S]*?)<\/a>/gi, (_, url, text) =>
      `[${stripTags(text).trim()}](${url})`
    )
    // Saltos de línea
    .replace(/<br\s*\/?>/gi, '\n')
    // Quitar etiquetas restantes
    .replace(/<[^>]+>/g, '')
    // Decodificar entidades HTML
    .replace(/&lt;/g,   '<')
    .replace(/&gt;/g,   '>')
    .replace(/&amp;/g,  '&')
    .replace(/&nbsp;/g, ' ')
    .replace(/&quot;/g, '"')
    .replace(/&#39;/g,  "'")
    // Normalizar espacios múltiples (mantener saltos de línea)
    .replace(/[ \t]{2,}/g, ' ')
    // Normalizar más de 2 saltos de línea consecutivos
    .replace(/\n{3,}/g, '\n\n')
    .trim();
}

/**
 * Quita todas las etiquetas HTML de una cadena (helper interno).
 */
function stripTags(html) {
  return html.replace(/<[^>]+>/g, '');
}

/**
 * Convierte una tabla HTML a Markdown GFM (GitHub Flavored Markdown).
 *
 * Entrada:   <table><thead><tr><th>A</th><th>B</th>...</tr></thead>
 *            <tbody><tr><td>1</td><td>2</td></tr>...</tbody></table>
 *
 * Salida:    | A | B |
 *            |---|---|
 *            | 1 | 2 |
 *
 * @param {HTMLTableElement} table - Elemento <table> del DOM
 * @returns {string} Tabla en formato Markdown
 */
function tableToMarkdown(table) {
  const rows = Array.from(table.querySelectorAll('tr'));
  if (rows.length === 0) return '';

  let md = '';

  rows.forEach((row, rowIndex) => {
    const cells = Array.from(row.querySelectorAll('th, td'));
    if (cells.length === 0) return;

    // Limpiar contenido de cada celda
    const cellTexts = cells.map(cell => {
      // Reemplazar saltos de línea internos por espacio para que quede en una línea
      const text = (cell.innerText || cell.textContent || '')
        .trim()
        .replace(/\n+/g, ' ')
        .replace(/\|/g, '\\|');  // escapar pipes que romperían la tabla
      return text || ' ';       // celda vacía → espacio para que GFM la renderice
    });

    // Línea de la tabla
    md += '| ' + cellTexts.join(' | ') + ' |\n';

    // Separador después del encabezado (primera fila o fila con <th>)
    const isHeader = row.querySelector('th') !== null || rowIndex === 0;
    if (isHeader) {
      md += '| ' + cells.map(() => '---').join(' | ') + ' |\n';
    }
  });

  return md;
}
```

---

## 5. Resumen de archivos modificados

| Archivo | Cambio | Sección |
|---|---|---|
| `demo/frontend/script.js` | Eliminar `Math.random()` para PMIDs | 1.1 |
| `demo/frontend/script.js` | Tabla de biomarcadores dinámica | 1.2 |
| `demo/frontend/script.js` | Gene pills interactivos + sidebar | 3.1 |
| `demo/frontend/script.js` | Modo consenso configurable | 3.2 |
| `demo/frontend/index.html` | Select de modo consenso | 3.2 |
| `demo/backend/main.py` | `get_top_enrichment_results()` | 2.1 |
| `demo/backend/main.py` | `fetch_diana_microt()`, `fetch_tarbase()` | 2.2 |
| `demo/backend/main.py` | Endpoint `/api/v1/diseases/{gene}` | 2.3 |
| `demo/backend/main.py` | Endpoint `/api/v1/gene-info/{gene}` | 1.2 |
| `demo/backend/main.py` | `compute_confidence_score()` | 2.4 |
| `demo/backend/main.py` | `PipelineRequest.consensus_mode` | 3.2 |
| `demo/backend/tests/test_main.py` | Tests completos con pytest | 2.5 |
| `kenryu_report_editor.html` | `exportMarkdown()` + helpers | 4.1 |
| `kenryu_report_editor.html` | Botón MD en toolbar | 4.2 |
| `kenryu_report_editor.html` | `innerHtmlToMarkdown()`, `tableToMarkdown()` | 4.3 |

---

_Guía generada para el sistema KENRYU · Tesis de Maestría · Mayo 2026_
