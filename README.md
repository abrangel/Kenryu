<div align="center">

<p align="center">
     <a href="README.es.md">🇪🇸   Leer en Español</a> | 🇬🇧   <b>Read in English</b>
    </p>

<img src="https://capsule-render.vercel.app/api?type=waving&height=200&color=0:1a1d24,50:c8a96e,100:1a1d24&text=KENRYU&fontColor=ffffff&fontSize=80&fontAlignY=40&desc=Advanced%20miRNA%20Bioinformatics%20Engine&descAlignY=68&descSize=18&animation=fadeIn" alt="Kenryu banner" width="100%"/>

<p align="center">
  <strong>High-precision bioinformatics platform for microRNA research and biomarker discovery</strong><br/>
  <em>Genomic Convergence · PubMed Evidence · Functional Enrichment · Automated Scientific Reporting</em>
</p>

<p align="center">
  <a href="https://huggingface.co/spaces/Kenryu007/Bioinformatica"><img src="https://img.shields.io/badge/Demo-Hugging%20Face%20Space-c8a96e?style=for-the-badge&logo=huggingface&logoColor=white" alt="Demo"/></a>
  <a href="#"><img src="https://img.shields.io/badge/Python-3.11-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python"/></a>
  <a href="#"><img src="https://img.shields.io/badge/FastAPI-0.109-009688?style=for-the-badge&logo=fastapi&logoColor=white" alt="FastAPI"/></a>
  <a href="https://www.targetscan.org"><img src="https://img.shields.io/badge/Data-TargetScan%208.0-4fc3a1?style=for-the-badge&logo=googlesheets&logoColor=white" alt="TargetScan"/></a>
</p>

<p align="center">
  <a href="https://abrangel.github.io/Kenryu/"><img src="https://img.shields.io/badge/Live%20Demo-🚀-c8a96e?style=for-the-badge&labelColor=1a1d24" alt="Live Demo"/></a>
  <a href="https://github.com/abrangel/Kenryu/issues"><img src="https://img.shields.io/badge/Report%20Bug-🐛-e05c5c?style=for-the-badge&labelColor=1a1d24" alt="Report Bug"/></a>
  <a href="https://github.com/abrangel/Kenryu/issues"><img src="https://img.shields.io/badge/Request%20Feature-💡-4fc3a1?style=for-the-badge&labelColor=1a1d24" alt="Request Feature"/></a>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/miRNA%20families-285-c8a96e?style=flat-square&labelColor=1a1d24" alt="285 families"/>
  <img src="https://img.shields.io/badge/genomic%20targets-18k+-4fc3a1?style=flat-square&labelColor=1a1d24" alt="18k+ targets"/>
  <img src="https://img.shields.io/badge/validation-miRTarBase-6eaadc?style=flat-square&labelColor=1a1d24" alt="miRTarBase"/>
  <img src="https://img.shields.io/badge/evidence-PubMed%20Integrated-f5a623?style=flat-square&labelColor=1a1d24" alt="PubMed Integrated"/>
  <img src="https://img.shields.io/badge/analysis-Functional%20Enrichment-c8a96e?style=flat-square&labelColor=1a1d24" alt="Functional Enrichment"/>
  <img src="https://img.shields.io/badge/DOI-Ready-005b9f?style=flat-square&labelColor=1a1d24" alt="DOI Ready"/>
</p>

<p align="center">
  <a href="#-overview"><b>Overview</b></a> ·
  <a href="#-scientific-rigor"><b>Scientific Rigor</b></a> ·
  <a href="#%EF%B8%8F-engine-architecture"><b>Architecture</b></a> ·
  <a href="#-local-installation"><b>Installation</b></a> ·
  <a href="#-api-documentation"><b>API</b></a> ·
  <a href="#-datasets"><b>Databases</b></a> ·
  <a href="#-citation--doi"><b>Citation & DOI</b></a>
</p>

</div>

---

## About the Project

MicroRNAs (miRNAs) are post-transcriptional regulators that control the expression of complex gene networks via targeted silencing. Identifying which genes are **simultaneously regulated by multiple miRNAs** is fundamental to discovering critical nodes in pathologies such as cancer, cardiovascular diseases, neurodegeneration, and metabolic disorders.

**KENRYU** automates this workflow:

1. **Integrates** thermodynamic predictions (TargetScan) with experimental validation (miRTarBase).
2. **Identifies** the dynamic Venn intersection of core genes regulated by the input panel.
3. **Enriches** the targets functionally (KEGG · Reactome · WikiPathways · GO).
4. **Investigates** each core gene across four sources in parallel: PubMed, OMIM, ClinVar, and ClinicalTrials.gov.
5. **Generates** an editable academic report with Vancouver-style bibliography, exportable to PDF, Markdown, and plain text (.txt).

Unlike traditional pipelines, KENRYU uses a **search cascade with fallback** (strict → expanded → minimal) and **term translation into English** to overcome PubMed's limitations with non-English queries. The entire flow is optimized for **Hugging Face Spaces**, including robust handling of NCBI eUtils rate-limits.

---
### ⚠️ Medical Disclaimer & Safety Note

KENRYU is designed for clinical education and decision support. It is not a substitute for professional medical judgment, individualized patient evaluation, local clinical guidelines, or review by a qualified healthcare professional. The evidence presented by the system comes from public databases (NCBI, miRTarBase, TargetScan, STRING-DB) and must be rigorously verified by the clinician before any therapeutic or research intervention.

*   **Regulatory Status:** KENRYU is for **Research Use Only (RUO)**. It is not intended for diagnostic or therapeutic procedures.
*   **Data Privacy:** The platform works as a decision support tool and does not store or process personally identifiable information (PII) or protected health information (PHI).
*   **Infrastructure Note:** Implemented and hosted on Hugging Face Spaces.
  
## Key Features

### Bioinformatics Analysis

- **Dynamic Convergence** — Venn intersection of N miRNAs (no hardcoded limit on core genes count).
- **Three Consensus Modes** — Strict, N-1, N-2 to detect robust regulatory networks.
- **Balanced Functional Enrichment** — fairness algorithm between KEGG, Reactome, WikiPathways, and GO.
- **Volcano plot, Venn diagram, and STRING-DB interactome** generated on-the-fly.

### Multi-source Research

| Source | Contribution |
|---|---|
| **PubMed** | Relevant scientific articles per gene + biological context, filtered by year |
| **OMIM** | Mendelian Inheritance in Man catalog — associated hereditary diseases |
| **ClinVar** | Reported pathogenic/likely pathogenic variants count |
| **ClinicalTrials.gov** | Active clinical trials where the gene is a primary focus (title filtered) |

### Robustness and Performance

- **429-aware Retries** — automatic handling of NCBI eUtils rate-limits.
- **Persistent Cache** — `local_db/analysis_cache.json` stores translations, PubMed searches, and enrichment results between sessions.
- **ES→EN Mapping** — automatically translates biological terms from PRESETS to English for effective PubMed queries.
- **Filter Cascade** — if the strict search yields no results, it progressively expands the year window.

### Report Generation

- **Professional A4 Editor** — WYSIWYG preview with smart pagination via height measurement.
- **Vancouver Style Bibliography** — `[N] Authors. Title. Journal. Year. PMID: X. URL`.
- **Multi-type Citations** — PMID / OMIM / ClinVar / NCT with correct identifier per source.
- **PDF Export** — via browser `window.print()` with optimized print CSS and footer anchored to each A4 page.
- **Professional Markdown Export** — ZIP containing `Report.md` + `/assets/` folder (separate PNG images) + explanatory `README.md`. Compatible with VSCode, Obsidian, Pandoc, and GitHub.

---

## Screenshots

> **Demo:** Real result examples.

```
docs/
└── screenshots/
    ├── editor.png
    ├── venn-plot.png
    ├── report-page.png
    └── narrative-evidence.png
```
<img width="1364" height="679" alt="Kenryu Interface" src="https://github.com/user-attachments/assets/28c4df8b-8f05-4913-8205-9430c7327af6" />

*Kenryu Graphics*
<img width="1364" height="676" alt="Graphics" src="https://github.com/user-attachments/assets/68353d73-c970-4af9-a467-d0984efba1e0" />

*Gene information in sidebars*
<img width="1365" height="680" alt="Gene Info" src="https://github.com/user-attachments/assets/903f93c8-50c6-486d-8cc2-92130dea8d16" />

*Research Results*
<img width="1364" height="679" alt="Research" src="https://github.com/user-attachments/assets/b4306c1d-029d-442f-be2b-3b0cf9f2e054" />

*Narrative Report Example*
<img width="1357" height="641" alt="Report" src="https://github.com/user-attachments/assets/dfad6fd3-54df-4811-9243-fc516f4a5168" />

---

## Tech Stack

**Backend:**
- Python 3.11
- FastAPI + Uvicorn
- httpx (Async HTTP)
- Matplotlib + Seaborn + matplotlib-venn (Visualization)
- Pandas + NumPy (Data manipulation)

**Frontend:**
- HTML5 + CSS3 (No framework, dark mode design with gold/teal accents)
- JavaScript ES2020+ (Vanilla)
- JSZip 3.10 (Markdown export packaging)
- html2pdf.js 0.10 (Loaded but unused — exports via `window.print()`)
- FontAwesome 6 (Iconography)

**Deployment:**
- Docker (Python slim)
- Hugging Face Spaces (CPU)
- GitHub Pages (Presentation)

---

## Architecture

```
                    ┌─────────────────────────────────┐
                    │   Frontend (HTML/JS/CSS)        │
                    │  - Paged A4 Editor              │
                    │  - PDF and MD-ZIP Export        │
                    └────────────┬────────────────────┘
                                 │ POST /api/v1/analyze
                                 ▼
                    ┌─────────────────────────────────┐
                    │  FastAPI Backend (kenryu_engine)│
                    │  - Venn Intersection            │
                    │  - Functional Enrichment        │
                    │  - Multi-source Research        │
                    └────────────┬────────────────────┘
                                 │
        ┌───────────┬────────────┼────────────┬────────────┐
        ▼           ▼            ▼            ▼            ▼
   ┌─────────┐ ┌─────────┐ ┌──────────┐ ┌─────────┐ ┌──────────┐
   │TargetSc.│ │miRTarBs.│ │ Enrichr  │ │ PubMed  │ │ OMIM /   │
   │ (local) │ │ (remote)│ │  (KEGG…) │ │(eUtils) │ │ClinVar / │
   └─────────┘ └─────────┘ └──────────┘ └─────────┘ │ CTrials  │
                                                    └──────────┘
                                 │
                                 ▼
                    ┌─────────────────────────────────┐
                    │  Persistent Cache               │
                    │  local_db/analysis_cache.json   │
                    └─────────────────────────────────┘
```

---

## Analysis Flow

### 1. miRNA Input
User enters N miRNAs (format `hsa-miR-XX-Yp`), selects year cutoff and consensus mode (Strict / N-1 / N-2).

### 2. Target Collection
For each miRNA, KENRYU:
- Queries local **TargetScanHuman v8.0** for validated thermodynamic predictions.
- Queries **miRTarBase** (via Harmonizome) for experimental interactions (CLIP-seq, Luciferase, Western Blot).

### 3. Venn Intersection
Applies the selected consensus algorithm:
- **Strict:** genes regulated by ALL miRNAs in the panel.
- **N-1:** genes regulated by N-1 or more miRNAs.
- **N-2:** genes regulated by N-2 or more miRNAs.

### 4. Functional Enrichment
For the obtained core genes:
- Queries Enrichr with balancing across KEGG, Reactome, WikiPathways, and GO Biological Process.
- For each significant pathway, searches for a supporting PMID in PubMed (with year filtering).

### 5. Multi-source Research per Core Gen
Each Venn gene triggers four parallel searches:

| Search | API | Strategy |
|---|---|---|
| PubMed | NCBI eUtils esearch + esummary | Cascade: gene + EN context + year → expanded → minimal |
| OMIM | NCBI eUtils (db=omim) | Direct search by gene symbol |
| ClinVar | NCBI eUtils (db=clinvar) | Filter: `pathogenic[Clinical_Significance]` |
| ClinicalTrials.gov | API v2 | Strict filter: title MUST contain the gene symbol |

All results are cached in `PERSISTENT_CACHE["gene_research"]` with key `v2|{gene}|{years}`.

### 6. Report Generation
`build_synthesis()` constructs 3 narrative sections:
- **I. Academic Synthesis** — description of each core gene + validated pathways + complementary genomic evidence (narrative integrated with OMIM/ClinVar/Trials) + clinical conclusion.
- **2.1 Global Functional Context** — each pathway enriched with its real p-adj and supporting PMID.
- **I.1 Bibliography** — Vancouver style citation adapted to source type.

---

## Integrated APIs and Databases

| Source | Function | Type |
|---|---|---|
| **TargetScanHuman v8.0** | Thermodynamic prediction of binding sites | Local (zip) |
| **miRTarBase** (via Harmonizome) | Experimental validation of interactions | Remote API |
| **Enrichr** | Functional enrichment (KEGG/Reactome/WikiPathways/GO) | Remote API |
| **NCBI PubMed eUtils** | Scientific articles per gene + context | Remote API (with retries) |
| **NCBI OMIM** | Associated hereditary diseases | Remote API |
| **NCBI ClinVar** | Reported pathogenic variants | Remote API |
| **ClinicalTrials.gov v2** | Active clinical trials | Remote API |
| **MyGene.info** | Clinical gene annotation | Remote API |
| **MyMemory API** | Real-time ES↔EN translation | Remote API |
| **STRING-DB** | Protein interactome generation | Remote API |

---

## Installation

### Prerequisites
- Python 3.11+
- pip
- (Optional) Docker for containerized deployment

### Local Installation
```bash
# 1. Clone repository
git clone https://github.com/abrangel/Kenryu.git
cd Kenryu

# 2. Create virtual environment (recommended)
python3.11 -m venv venv
source venv/bin/activate   # Linux/macOS
# venv\Scripts\activate    # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Start the server
uvicorn kenryu_engine:app --host 0.0.0.0 --port 7860 --reload

# 5. Open in browser
# http://localhost:7860
```

### Docker Deployment
```bash
docker build -t kenryu .
docker run -p 7860:7860 -v $(pwd)/local_db:/app/local_db kenryu
```
The `local_db` volume persists the cache across container restarts.

---

## Usage

### Via Web Interface
1. Open [Production Space](https://huggingface.co/spaces/Kenryu007/Bioinformatica) or `http://localhost:7860`.
2. Enter miRNAs separated by commas: `hsa-miR-33a-5p, hsa-miR-144-3p, hsa-miR-758-3p`.
3. Select cutoff year (PubMed evidence age filter).
4. Select consensus mode.
5. Click **Execute**.

### Via REST API
```bash
curl -X POST "http://localhost:7860/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{
    "mirnas": ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p"],
    "years": 10,
    "mode": "strict"
  }'
```

---

## Roadmap

### Implemented
- [x] Dynamic Venn convergence (Strict / N-1 / N-2 modes)
- [x] Balanced functional enrichment (KEGG / Reactome / WikiPathways / GO)
- [x] Visualization: Venn, Volcano, STRING-DB interactome
- [x] Vancouver style bibliography
- [x] Multi-source research (PubMed + OMIM + ClinVar + ClinicalTrials)
- [x] Versioned persistent cache
- [x] 429-aware retries in NCBI eUtils
- [x] ES→EN mapping for effective PubMed queries
- [x] Professional Markdown export with ZIP + assets

### Planned
- [ ] Support for miRNA isomers (isomiRs)
- [ ] DisGeNET integration for gene-disease association
- [ ] Structured JSON export (exchange with other tools)
- [ ] Side-by-side miRNA panel comparison
- [ ] TargetScan v9.0 integration when available

---

## Contributing
Contributions are welcome. To report bugs, suggest features, or submit pull requests:
1. **Fork** the repository
2. Create a branch: `git checkout -b feat/feature-name`
3. Commit with descriptive messages
4. Push: `git push origin feat/feature-name`
5. Open a Pull Request

---

## License
This project is distributed under a **Free Academic License for Research and Education**. For commercial use, contact the author.

---

## Author
**Cesar Manzo**
- Clinical Bioinformatics · Genomic Analysis · Translational Medicine
- Platform: [Kenryu on Hugging Face](https://huggingface.co/spaces/Kenryu007/Bioinformatica)

---

<div align="center">

**KENRYU Bioinformatics Engine** · 2026

*Made with scientific rigor for the genomic research community.*

</div>
