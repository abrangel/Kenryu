# KENRYU — Motor Bioinformático de microARNs

<div align="center">

[![Hugging Face Spaces](https://img.shields.io/badge/🤗_Hugging_Face-Spaces-FFD21E?style=flat-square)](https://huggingface.co/spaces/Kenryu007/Bioinformatica)
[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-181717?style=flat-square&logo=github)](https://abrangel.github.io/Kenryu)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-3776AB?style=flat-square&logo=python&logoColor=white)](https://www.python.org/)
[![FastAPI](https://img.shields.io/badge/FastAPI-0.115-009688?style=flat-square&logo=fastapi&logoColor=white)](https://fastapi.tiangolo.com/)
[![License](https://img.shields.io/badge/License-Academic-blue?style=flat-square)](#licencia)
[![Status](https://img.shields.io/badge/Status-Production-green?style=flat-square)](https://huggingface.co/spaces/Kenryu007/Bioinformatica)

**Plataforma de análisis bioinformático para la identificación de dianas convergentes de microARNs e integración de evidencia clínica multi-fuente.**

[Demo en vivo](https://huggingface.co/spaces/Kenryu007/Bioinformatica) · [Reportar bug](https://github.com/abrangel/Kenryu/issues) · [Solicitar feature](https://github.com/abrangel/Kenryu/issues)

</div>

---

## Índice

- [Sobre el proyecto](#sobre-el-proyecto)
- [Características principales](#características-principales)
- [Capturas](#capturas)
- [Stack tecnológico](#stack-tecnológico)
- [Arquitectura](#arquitectura)
- [Flujo de análisis](#flujo-de-análisis)
- [APIs y bases de datos integradas](#apis-y-bases-de-datos-integradas)
- [Instalación](#instalación)
- [Uso](#uso)
- [Exportación de resultados](#exportación-de-resultados)
- [Estructura del repositorio](#estructura-del-repositorio)
- [Roadmap](#roadmap)
- [Contribuir](#contribuir)
- [Licencia](#licencia)
- [Autor](#autor)

---

## Sobre el proyecto

Los microARNs (miRNAs) son reguladores post-transcripcionales que controlan la expresión de redes génicas complejas mediante silenciamiento dirigido. Identificar qué genes son regulados **simultáneamente por múltiples miRNAs** es fundamental para descubrir nodos críticos en patologías como cáncer, enfermedades cardiovasculares, neurodegeneración y desórdenes metabólicos.

**KENRYU** automatiza este flujo:

1. **Integra** predicciones termodinámicas (TargetScan) con validación experimental (miRTarBase).
2. **Identifica** la intersección Venn dinámica de los genes core regulados por el panel ingresado.
3. **Enriquece** funcionalmente las dianas (KEGG · Reactome · WikiPathways · GO).
4. **Investiga** cada gen core en cuatro fuentes en paralelo: PubMed, OMIM, ClinVar y ClinicalTrials.gov.
5. **Genera** un informe académico paginado con bibliografía estilo Vancouver, exportable a PDF y Markdown.

A diferencia de pipelines tradicionales, KENRYU usa **cascada de búsqueda con fallback** (estricto → ampliado → mínimo) y **traducción de términos al inglés** para superar limitaciones de PubMed con queries en español. Todo el flujo está optimizado para **Hugging Face Spaces**, incluyendo manejo robusto de rate-limits de NCBI eUtils.

---

## Características principales

### Análisis bioinformático

- **Convergencia dinámica** — intersección Venn de N miRNAs (sin límite hardcoded en el número de genes core).
- **Tres modos de consenso** — Estricto, N-1, N-2 para detectar redes regulatorias robustas.
- **Enriquecimiento funcional balanceado** — algoritmo de equidad entre KEGG, Reactome, WikiPathways y GO.
- **Volcano plot, diagrama de Venn e interactoma STRING-DB** generados al vuelo.

### Investigación multi-fuente (nuevo en v1.38)

| Fuente | Qué aporta |
|---|---|
| **PubMed** | Artículos científicos relevantes por gen + contexto biológico, filtrados por año |
| **OMIM** | Catálogo Mendelian Inheritance in Man — enfermedades hereditarias asociadas |
| **ClinVar** | Conteo de variantes patogénicas/probablemente patogénicas reportadas |
| **ClinicalTrials.gov** | Ensayos clínicos activos donde el gen es protagonista (filtrado por título) |

### Robustez y rendimiento

- **Reintentos 429-aware** — manejo automático de rate-limits de NCBI eUtils.
- **Cache persistente** — `local_db/analysis_cache.json` almacena traducciones, búsquedas PubMed y resultados de enriquecimiento entre sesiones.
- **Mapeo ES→EN** — traduce automáticamente términos biológicos del PRESET al inglés para queries efectivas en PubMed.
- **Cascada de filtros** — si la búsqueda estricta no devuelve resultados, amplía progresivamente la ventana de años.

### Generación de informes

- **Editor A4 profesional** — vista previa WYSIWYG con paginación que se puede editar y mejorar el contenido con las
- propias palabras del investigador.
- **Bibliografía estilo Vancouver** — `[N] Autores. Título. Revista. Año. PMID: X. URL`.
- **Citas multi-tipo** — PMID / OMIM / ClinVar / NCT con identificador correcto por fuente.
- **Exportación PDF** — vía `window.print()` con CSS de impresión optimizado, footer anclado al fondo de cada hoja A4.
- **Exportación Markdown profesional** — ZIP con `Reporte.md` + carpeta `/assets/` (imágenes PNG separadas) + `README.md` explicativo. Compatible con VSCode, Obsidian, Pandoc y GitHub.

---

## Capturas

> **Sugerencia:** captura screenshots del editor, del Venn y del informe exportado y agrégalos en `docs/screenshots/` para mostrar el resultado real.

```
docs/
└── screenshots/
    ├── editor.png
    ├── venn-plot.png
    ├── report-page.png
    └── narrative-evidence.png
```

---

## Stack tecnológico

**Backend:**
- Python 3.11
- FastAPI + Uvicorn
- httpx (HTTP async)
- Matplotlib + Seaborn + matplotlib-venn (visualización)
- Pandas + NumPy (manipulación de datos)

**Frontend:**
- HTML5 + CSS3 (sin framework, diseño dark mode con acentos gold/teal)
- JavaScript ES2020+ (vanilla)
- JSZip 3.10 (empaquetado de exportación Markdown)
- html2pdf.js 0.10 (cargado pero no usado — exporta vía `window.print()`)
- FontAwesome 6 (iconografía)

**Despliegue:**
- Docker (Python slim)
- Hugging Face Spaces (CPU)
- GitHub Pages (presentación)

---

## Arquitectura

```
                    ┌─────────────────────────────────┐
                    │   Frontend (HTML/JS/CSS)        │
                    │  - Editor A4 paginado           │
                    │  - Exportación PDF y MD-ZIP     │
                    └────────────┬────────────────────┘
                                 │ POST /api/v1/analyze
                                 ▼
                    ┌─────────────────────────────────┐
                    │  FastAPI Backend (kenryu_engine)│
                    │  - Intersección Venn            │
                    │  - Enriquecimiento funcional    │
                    │  - Investigación multi-fuente   │
                    └────────────┬────────────────────┘
                                 │
        ┌───────────┬────────────┼────────────┬────────────┐
        ▼           ▼            ▼            ▼            ▼
   ┌─────────┐ ┌─────────┐ ┌──────────┐ ┌─────────┐ ┌──────────┐
   │TargetSc.│ │miRTarBs.│ │ Enrichr  │ │ PubMed  │ │ OMIM /   │
   │ (local) │ │ (remoto)│ │  (KEGG…) │ │(eUtils) │ │ClinVar / │
   └─────────┘ └─────────┘ └──────────┘ └─────────┘ │ CTrials  │
                                                    └──────────┘
                                 │
                                 ▼
                    ┌─────────────────────────────────┐
                    │  Cache persistente              │
                    │  local_db/analysis_cache.json   │
                    └─────────────────────────────────┘
```

---

## Flujo de análisis

### 1. Entrada de miRNAs

El usuario ingresa N miRNAs (formato `hsa-miR-XX-Yp`), selecciona año de corte y modo de consenso (Estricto / N-1 / N-2).

### 2. Recolección de dianas

Por cada miRNA, KENRYU:
- Consulta **TargetScanHuman v8.0** local para predicciones termodinámicas validadas.
- Consulta **miRTarBase** (vía Harmonizome) para interacciones experimentales (CLIP-seq, Luciferasa, Western Blot).

### 3. Intersección Venn

Aplica el algoritmo de consenso seleccionado:

- **Estricto:** genes regulados por TODOS los miRNAs del panel.
- **N-1:** genes regulados por N-1 o más miRNAs.
- **N-2:** genes regulados por N-2 o más miRNAs.

### 4. Enriquecimiento funcional

Para los genes core obtenidos:
- Consulta Enrichr con balanceo entre KEGG, Reactome, WikiPathways y GO Biological Process.
- Para cada ruta significativa, busca un PMID de respaldo en PubMed (con filtro de años).

### 5. Investigación multi-fuente por gen core

Cada gen del Venn dispara cuatro búsquedas paralelas:

| Búsqueda | API | Estrategia |
|---|---|---|
| PubMed | NCBI eUtils esearch + esummary | Cascada: gen + contexto EN + año → ampliado → mínimo |
| OMIM | NCBI eUtils (db=omim) | Búsqueda directa por símbolo de gen |
| ClinVar | NCBI eUtils (db=clinvar) | Filtro: `pathogenic[Clinical_Significance]` |
| ClinicalTrials.gov | API v2 | Filtro estricto: título DEBE contener el símbolo del gen |

Todos los resultados se cachean en `PERSISTENT_CACHE["gene_research"]` con clave `v2|{gene}|{years}`.

### 6. Generación del informe

`build_synthesis()` construye 3 secciones narrativas:

- **I. Síntesis académica** — descripción de cada gen core + rutas validadas + evidencia genómica complementaria (narrativa integrada con OMIM/ClinVar/Trials) + conclusión clínica.
- **2.1 Contexto funcional global** — cada ruta enriquecida con su p-adj real y PMID de respaldo.
- **I.1 Bibliografía** — cita estilo Vancouver adaptada al tipo de fuente.

---

## APIs y bases de datos integradas

| Fuente | Función | Tipo |
|---|---|---|
| **TargetScanHuman v8.0** | Predicción termodinámica de sitios de unión | Local (zip) |
| **miRTarBase** (vía Harmonizome) | Validación experimental de interacciones | API remota |
| **Enrichr** | Enriquecimiento funcional (KEGG/Reactome/WikiPathways/GO) | API remota |
| **NCBI PubMed eUtils** | Artículos científicos por gen + contexto | API remota (con reintentos) |
| **NCBI OMIM** | Enfermedades hereditarias asociadas | API remota |
| **NCBI ClinVar** | Variantes patogénicas reportadas | API remota |
| **ClinicalTrials.gov v2** | Ensayos clínicos activos | API remota |
| **MyGene.info** | Anotación clínica de genes | API remota |
| **MyMemory API** | Traducción ES↔EN en tiempo real | API remota |
| **STRING-DB** | Generación del interactoma proteico | API remota |

---

## Instalación

### Requisitos

- Python 3.11+
- pip
- (Opcional) Docker para despliegue en contenedor

### Instalación local

```bash
# 1. Clonar repositorio
git clone https://github.com/abrangel/Kenryu.git
cd Kenryu

# 2. Crear entorno virtual (recomendado)
python3.11 -m venv venv
source venv/bin/activate   # Linux/macOS
# venv\Scripts\activate    # Windows

# 3. Instalar dependencias
pip install -r requirements.txt

# 4. Iniciar el servidor
uvicorn kenryu_engine:app --host 0.0.0.0 --port 7860 --reload

# 5. Abrir en navegador
# http://localhost:7860
```

### Despliegue con Docker

```bash
docker build -t kenryu .
docker run -p 7860:7860 -v $(pwd)/local_db:/app/local_db kenryu
```

El volumen `local_db` persiste el cache entre reinicios del contenedor.

### Despliegue en Hugging Face Spaces

El repo está listo para HF Spaces — el Dockerfile expone el puerto 7860 y monta `local_db` para cache persistente. Solo se necesita conectar el Space al repo de GitHub.

---

## Uso

### Vía interfaz web

1. Abrir [el Space en producción](https://huggingface.co/spaces/Kenryu007/Bioinformatica) o `http://localhost:7860`.
2. Ingresar miRNAs separados por comas: `hsa-miR-33a-5p, hsa-miR-144-3p, hsa-miR-758-3p`.
3. Seleccionar año de corte (filtro de antigüedad de evidencia PubMed).
4. Seleccionar modo de consenso.
5. Click en **Ejecutar**.

### Vía API REST

```bash
curl -X POST "http://localhost:7860/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "mirnas": ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p"],
    "years": 10,
    "mode": "strict"
  }'
```

**Respuesta:**

```json
{
  "common_genes": ["ABCA1", "KPNA3", "SCN1A", "..."],
  "gene_details": { "ABCA1": { "system": "...", "pathology": "...", ... } },
  "gene_research": {
    "ABCA1": {
      "pubmed": [{"pmid": "...", "term": "...", "year_window": "10y"}],
      "omim": [{"omim_id": "600046", "title": "..."}],
      "clinvar": {"count": 201, "url": "..."},
      "trials": [{"nct_id": "NCT01456650", "title": "..."}]
    }
  },
  "enrichment": [...],
  "scientific_synthesis": "...",
  "functional_context": "...",
  "references_text": "...",
  "report_references": [...],
  "venn_plot": "data:image/png;base64,...",
  "volcano_plot": "data:image/png;base64,...",
  "ppi_plot": "data:image/png;base64,..."
}
```

---

## Exportación de resultados

KENRYU exporta el informe en tres formatos desde el editor:

### PDF

Usa `window.print()` del navegador → "Guardar como PDF". El CSS `@media print` configura:
- Hoja A4 fija (210 × 297 mm)
- Footer "KENRYU Bioinformatics Engine" anclado al fondo de cada página
- Paginación inteligente por medición de altura para evitar solapamientos

### Markdown profesional (.zip)

Genera un paquete ZIP estándar para bioinformáticos:

```
Reporte-KENRYU.zip
├── Reporte.md              # Markdown completo (Pandoc-compatible)
│   ├── Frontmatter YAML
│   ├── 6 encabezados ##
│   ├── 2 tablas GFM
│   ├── 40+ enlaces [texto](url)
│   ├── Referencias a imágenes ![alt](assets/...)
│   └── Narrativa completa con OMIM/ClinVar/Trials
├── assets/
│   ├── figura-1.png        # Diagrama de Venn
│   ├── figura-2.png        # Volcano plot
│   └── figura-3.png        # Interactoma STRING-DB
└── README.md               # Instrucciones de uso
```

**Casos de uso:**
- Editar en VSCode/Obsidian/JupyterLab sin re-renderizar gráficos
- Versionar en Git/GitHub (texto + imágenes binarias separadas)
- Convertir a otros formatos: `pandoc Reporte.md -o Reporte.pdf`
- Publicar directamente en GitHub Pages (renderizado automático)

### TXT (texto plano)

Para análisis rápido sin estructura. Vuelca el `innerText` del editor.

---

## Estructura del repositorio

```
Kenryu/
├── kenryu_engine.py              # Motor FastAPI + lógica bioinformática (~1600 líneas)
├── requirements.txt              # Dependencias Python
├── Dockerfile                    # Configuración de despliegue
├── targetscan_full.json.zip      # Base TargetScan v8.0 indexada
├── hsa-miR-XXX-Yp.txt            # Bases pre-procesadas por miRNA (ejemplos)
├── index.html                    # Interfaz principal del editor
├── script.js                     # Lógica del frontend (paginación, exportación)
├── style.css                     # Diseño visual (dark mode + acentos gold/teal)
├── static/
│   ├── index.html                # Copia servida por FastAPI
│   ├── script.js                 # Copia servida por FastAPI
│   └── style.css                 # Copia servida por FastAPI
├── local_db/
│   └── analysis_cache.json       # Cache persistente (creado en runtime)
└── README.md                     # Este archivo
```

> **Nota:** El backend sirve los archivos del frontend desde `static/`. Mantener `static/*` sincronizado con la raíz del repo si se editan ambos. La copia en la raíz existe por compatibilidad con GitHub Pages.

---

## Roadmap

### Implementado

- [x] Convergencia Venn dinámica (modos Estricto / N-1 / N-2)
- [x] Enriquecimiento funcional balanceado (KEGG / Reactome / WikiPathways / GO)
- [x] Visualización: Venn, Volcano, Interactoma STRING-DB
- [x] Bibliografía estilo Vancouver
- [x] Investigación multi-fuente (PubMed + OMIM + ClinVar + ClinicalTrials)
- [x] Cache persistente versionado
- [x] Reintentos 429-aware en NCBI eUtils
- [x] Mapeo ES→EN para queries PubMed efectivas
- [x] Exportación Markdown profesional con ZIP + assets
- [x] Footer A4 anclado correctamente al fondo en PDF
- [x] Narrativa de evidencia genómica integrada al cuerpo del informe

### Planificado

- [ ] Soporte para isómeros de miRNA (isomiRs)
- [ ] Integración con DisGeNET para asociación gen-enfermedad
- [ ] Exportación a JSON estructurado (intercambio con otras herramientas)
- [ ] Comparación side-by-side de paneles de miRNAs
- [ ] Heatmap de expresión cruzada (genes × miRNAs)
- [ ] Modo "Batch" para procesar múltiples paneles desde CSV
- [ ] Integración con TargetScan v9.0 cuando esté disponible

---

## Contribuir

Las contribuciones son bienvenidas. Para reportar bugs, sugerir features o enviar pull requests:

1. **Fork** el repositorio
2. Crea una rama: `git checkout -b feat/nombre-feature`
3. Commit con mensajes descriptivos
4. Push: `git push origin feat/nombre-feature`
5. Abre un Pull Request

### Convenciones de commit

- `feat:` nueva funcionalidad
- `fix:` corrección de bug
- `docs:` cambios en documentación
- `refactor:` refactorización sin cambio funcional
- `chore:` cambios en infraestructura/build

### Reportar bugs

Por favor incluye:
- Pasos para reproducir
- Resultado esperado vs resultado obtenido
- Logs del navegador (consola) si es un error de frontend
- Logs del Hugging Face Space si es un error de backend

---

## Licencia

Este proyecto se distribuye bajo una **licencia académica de uso libre para investigación y educación**. Para uso comercial, contactar al autor.

KENRYU integra datos de fuentes públicas (NCBI, Enrichr, OMIM, ClinVar, ClinicalTrials.gov, STRING-DB) sujetos a sus respectivas políticas de uso. El usuario es responsable de cumplir con los términos de servicio de cada base de datos.

---

## Autor

**Cesar Manzo**

- Bioinformática Clínica · Análisis Genómico · Medicina Traslacional
- Plataforma: [Kenryu en Hugging Face](https://huggingface.co/spaces/Kenryu007/Bioinformatica)
- Presentación: [Kenryu en GitHub Pages](https://abrangel.github.io/Kenryu)
- Otros proyectos: [@abrangel](https://github.com/abrangel)

---

## Reconocimientos

- **TargetScanHuman** — Lewis Lab, Whitehead Institute
- **miRTarBase** — Chou *et al.*, Nucleic Acids Research
- **NCBI** — National Center for Biotechnology Information
- **Enrichr** — Ma'ayan Lab, Mount Sinai
- **ClinicalTrials.gov** — U.S. National Library of Medicine
- **STRING-DB** — European Molecular Biology Laboratory

---

<div align="center">

**KENRYU Bioinformatics Engine** · 2026

*Hecho con rigor científico para la comunidad de investigación genómica.*

</div>

