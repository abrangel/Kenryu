<div align="center">

<p align="center">
  <a href="README.es.md">🇪🇸 Leer en Español</a> &nbsp;|&nbsp; 🇬🇧 <b>Read in English</b>
</p>

<img src="https://capsule-render.vercel.app/api?type=waving&height=220&color=0:0d0f14,40:1a1423,70:c8a96e,100:0d0f14&text=KENRYU&fontColor=ffffff&fontSize=90&fontAlignY=42&desc=Advanced%20miRNA%20Bioinformatics%20Engine&descAlignY=65&descSize=20&animation=fadeIn&stroke=c8a96e&strokeWidth=1" width="100%" alt="KENRYU"/>

<br/>

<p align="center">
  <b>High-precision bioinformatics platform for microRNA research and biomarker discovery</b><br/>
  <sub>Genomic Convergence &nbsp;·&nbsp; PubMed Evidence &nbsp;·&nbsp; Functional Enrichment &nbsp;·&nbsp; Automated Scientific Reporting</sub>
</p>

<br/>

<!-- PRIMARY ACTIONS -->
<p align="center">
  <a href="https://huggingface.co/spaces/Kenryu007/Bioinformatica">
    <img src="https://img.shields.io/badge/%F0%9F%A4%97%20Live%20Demo-Hugging%20Face%20Space-c8a96e?style=for-the-badge&labelColor=1a1d24" alt="Live Demo"/>
  </a>
  &nbsp;
  <a href="https://abrangel.github.io/Kenryu/">
    <img src="https://img.shields.io/badge/%F0%9F%8C%90%20Website-GitHub%20Pages-4fc3a1?style=for-the-badge&labelColor=1a1d24" alt="Website"/>
  </a>
  &nbsp;
  <a href="https://github.com/abrangel/Kenryu/issues">
    <img src="https://img.shields.io/badge/%F0%9F%90%9B%20Report-Bug-e05c5c?style=for-the-badge&labelColor=1a1d24" alt="Report Bug"/>
  </a>
  &nbsp;
  <a href="https://github.com/abrangel/Kenryu/issues">
    <img src="https://img.shields.io/badge/%F0%9F%92%A1%20Request-Feature-6eaadc?style=for-the-badge&labelColor=1a1d24" alt="Request Feature"/>
  </a>
</p>

<br/>

<!-- TECH BADGES -->
<p align="center">
  <img src="https://img.shields.io/badge/Python-3.11-3776AB?style=flat-square&logo=python&logoColor=white&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/FastAPI-0.109-009688?style=flat-square&logo=fastapi&logoColor=white&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/Docker-ready-2496ED?style=flat-square&logo=docker&logoColor=white&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/HF%20Spaces-deployed-FFD21E?style=flat-square&logo=huggingface&logoColor=black&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/License-Academic%20RUO-c8a96e?style=flat-square&labelColor=1a1d24"/>
</p>

<!-- STATS -->
<p align="center">
  <img src="https://img.shields.io/badge/miRNA%20families-285-c8a96e?style=flat-square&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/genomic%20targets-18k%2B-4fc3a1?style=flat-square&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/validation-miRTarBase-6eaadc?style=flat-square&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/evidence-PubMed%20Integrated-f5a623?style=flat-square&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/enrichment-KEGG%20%C2%B7%20GO%20%C2%B7%20Reactome-c8a96e?style=flat-square&labelColor=1a1d24"/>
  <img src="https://img.shields.io/badge/status-RUO%20Only-e05c5c?style=flat-square&labelColor=1a1d24"/>
</p>

<br/>

<!-- NAV -->
<p align="center">
  <a href="#-about-the-project"><b>Overview</b></a> &nbsp;·&nbsp;
  <a href="#-analysis-flow"><b>Analysis Flow</b></a> &nbsp;·&nbsp;
  <a href="#-key-features"><b>Features</b></a> &nbsp;·&nbsp;
  <a href="#-tech-stack"><b>Tech Stack</b></a> &nbsp;·&nbsp;
  <a href="#-integrated-apis--databases"><b>Databases</b></a> &nbsp;·&nbsp;
  <a href="#-installation"><b>Installation</b></a> &nbsp;·&nbsp;
  <a href="#-roadmap"><b>Roadmap</b></a>
</p>

</div>

<br/>

---

## ⚠️ Medical Disclaimer

> **KENRYU is for Research Use Only (RUO).** It is not a substitute for professional medical judgment, individualized patient evaluation, local clinical guidelines, or review by a qualified healthcare professional. Evidence comes from public databases (NCBI, miRTarBase, TargetScan, STRING-DB) and must be rigorously verified before any therapeutic or research intervention. This platform does not store or process PII or PHI.

---

## 🔬 About the Project

MicroRNAs (miRNAs) are post-transcriptional regulators that control the expression of complex gene networks via targeted silencing. Identifying which genes are **simultaneously regulated by multiple miRNAs** is fundamental to discovering critical nodes in pathologies such as cancer, cardiovascular diseases, neurodegeneration, and metabolic disorders.

**KENRYU** automates the entire discovery workflow in one platform:

| Step | Action | Sources |
|:----:|--------|---------|
| **1** | 🔗 **Integrate** thermodynamic predictions + experimental validation | TargetScan · miRTarBase |
| **2** | 🎯 **Identify** dynamic Venn intersection of core regulated genes | Local DB + Harmonizome |
| **3** | 🧬 **Enrich** targets functionally with balancing algorithm | KEGG · Reactome · WikiPathways · GO |
| **4** | 📡 **Investigate** each core gene across 4 sources in parallel | PubMed · OMIM · ClinVar · ClinicalTrials |
| **5** | 📄 **Generate** editable academic report with Vancouver bibliography | PDF · Markdown ZIP · TXT |

> Unlike traditional pipelines, KENRYU uses a **search cascade with fallback** (strict → expanded → minimal) and **real-time ES→EN translation** to overcome PubMed's limitations with non-English queries. Fully optimized for **Hugging Face Spaces** with robust NCBI eUtils rate-limit handling.

---

## ⚡ Analysis Flow

```
  INPUT                COLLECTION            INTERSECTION         ENRICHMENT
┌──────────┐         ┌────────────┐         ┌──────────────┐    ┌──────────────┐
│ N miRNAs │──────▶  │TargetScan  │ ──────▶ │  Venn Core   │──▶ │ KEGG/GO/     │
│hsa-miR-  │         │ v8.0 local │         │  Gene Set    │    │ Reactome/    │
│  XX-Yp   │         ├────────────┤  modes: │              │    │ WikiPathways │
└──────────┘         │ miRTarBase │  Strict  │  N-1 · N-2   │    └──────┬───────┘
                     │ (remote)   │  N-1     └──────────────┘           │
                     └────────────┘  N-2                                ▼
                                                               ┌──────────────────┐
  REPORT              RESEARCH (×4 parallel per gene)          │  PubMed PMID     │
┌──────────┐         ┌────────────────────────────────┐        │  per pathway     │
│ PDF/MD/  │◀──────  │ PubMed · OMIM · ClinVar ·      │        └──────────────────┘
│   TXT    │         │ ClinicalTrials.gov              │
└──────────┘         └────────────────────────────────┘
```

### Consensus Modes

| Mode | Description | Use Case |
|------|-------------|----------|
| **Strict** | Genes regulated by ALL miRNAs in the panel | High-confidence core targets |
| **N-1** | Genes regulated by at least N-1 miRNAs | Robust regulatory networks |
| **N-2** | Genes regulated by at least N-2 miRNAs | Broader pathway exploration |

---

## 🧩 Key Features

### 🔬 Bioinformatics Analysis

- **Dynamic Convergence** — Venn intersection of N miRNAs, no hardcoded limit on core gene count
- **Balanced Functional Enrichment** — fairness algorithm across KEGG, Reactome, WikiPathways, and GO Biological Process
- **Volcano plot · Venn diagram · STRING-DB interactome** — generated on-the-fly per analysis

### 📡 Multi-source Research

| Source | What it provides |
|--------|-----------------|
| **PubMed** | Scientific articles per gene + biological context, year-filtered with cascade fallback |
| **OMIM** | Mendelian Inheritance in Man — associated hereditary diseases |
| **ClinVar** | Reported pathogenic / likely pathogenic variant count |
| **ClinicalTrials.gov** | Active clinical trials where the gene symbol appears in the title |

### ⚡ Robustness & Performance

- **429-aware Retries** — automatic NCBI eUtils rate-limit handling with exponential backoff
- **Persistent Cache** — `local_db/analysis_cache.json` stores translations, PubMed results, and enrichment across sessions
- **ES→EN Mapping** — auto-translates biological terms from Spanish presets for effective PubMed queries
- **Filter Cascade** — strict → expanded → minimal year window fallback

### 📄 Report Generation

- **Professional A4 Editor** — WYSIWYG preview with smart height-based pagination
- **Vancouver Bibliography** — `[N] Authors. Title. Journal. Year. PMID: X. URL`
- **Multi-type Citations** — PMID / OMIM / ClinVar / NCT with correct identifier per source
- **PDF Export** — via `window.print()` with optimized `@media print` CSS, footer anchored to each A4 page
- **Markdown ZIP Export** — `Report.md` + `/assets/` PNGs + `README.md` — Pandoc/Obsidian/VSCode compatible

---

## 📸 Screenshots

<div align="center">

<img width="1364" alt="Kenryu main interface with analysis results" src="https://github.com/user-attachments/assets/28c4df8b-8f05-4913-8205-9430c7327af6"/>
<sub><b>Main interface</b> — Volcano plot, Venn diagram and STRING-DB interactome generated live</sub>

<br/><br/>

<img width="1364" alt="Gene charts and visualizations" src="https://github.com/user-attachments/assets/68353d73-c970-4af9-a467-d0984efba1e0"/>
<sub><b>Bioinformatics charts</b> — Functional enrichment visualization</sub>

<br/><br/>

<img width="1365" alt="Gene information in sidebars" src="https://github.com/user-attachments/assets/903f93c8-50c6-486d-8cc2-92130dea8d16"/>
<sub><b>Gene sidebars</b> — Quick access to OMIM, ClinVar, and trial data</sub>

<br/><br/>

<img width="1364" alt="Multi-source research results" src="https://github.com/user-attachments/assets/b4306c1d-029d-442f-be2b-3b0cf9f2e054"/>
<sub><b>Research results</b> — PubMed · OMIM · ClinVar · ClinicalTrials consolidated</sub>

<br/><br/>

<img width="1357" alt="Narrative academic report" src="https://github.com/user-attachments/assets/dfad6fd3-54df-4811-9243-fc516f4a5168"/>
<sub><b>Narrative report</b> — Academic synthesis with Vancouver bibliography (excerpt)</sub>

</div>

---

## 🏗️ Architecture

```
                    ┌─────────────────────────────────┐
                    │   Frontend (HTML / JS / CSS)    │
                    │  · Paged A4 WYSIWYG Editor      │
                    │  · PDF · Markdown ZIP · TXT     │
                    └────────────┬────────────────────┘
                                 │  POST /api/v1/analyze
                                 ▼
                    ┌─────────────────────────────────┐
                    │  FastAPI Backend — kenryu_engine│
                    │  · Venn Intersection Engine     │
                    │  · Functional Enrichment        │
                    │  · Multi-source Research        │
                    └────────────┬────────────────────┘
                                 │
       ┌──────────┬──────────────┼───────────────┬──────────────────┐
       ▼          ▼              ▼               ▼                  ▼
 ┌──────────┐ ┌──────────┐ ┌─────────┐ ┌─────────────┐ ┌────────────────┐
 │TargetScan│ │miRTarBase│ │ Enrichr │ │NCBI PubMed  │ │ OMIM · ClinVar │
 │ v8 local │ │(Harmoniz)│ │KEGG/GO… │ │  eUtils     │ │  ClinicalTrials│
 └──────────┘ └──────────┘ └─────────┘ └─────────────┘ └────────────────┘
                                 │
                                 ▼
                    ┌─────────────────────────────────┐
                    │  Persistent Cache               │
                    │  local_db/analysis_cache.json   │
                    └─────────────────────────────────┘
```

---

## 🛠️ Tech Stack

<table>
<tr>
<td valign="top" width="33%">

**Backend**

![Python](https://img.shields.io/badge/Python%203.11-3776AB?style=flat-square&logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-009688?style=flat-square&logo=fastapi&logoColor=white)
![httpx](https://img.shields.io/badge/httpx%20async-1a1d24?style=flat-square)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557c?style=flat-square)
![Pandas](https://img.shields.io/badge/Pandas%20+%20NumPy-150458?style=flat-square&logo=pandas&logoColor=white)

</td>
<td valign="top" width="33%">

**Frontend**

![HTML5](https://img.shields.io/badge/HTML5%20+%20CSS3-E34F26?style=flat-square&logo=html5&logoColor=white)
![JavaScript](https://img.shields.io/badge/Vanilla%20JS%20ES2020+-F7DF1E?style=flat-square&logo=javascript&logoColor=black)
![JSZip](https://img.shields.io/badge/JSZip%203.10-1a1d24?style=flat-square)
![FontAwesome](https://img.shields.io/badge/FontAwesome%206-528DD7?style=flat-square&logo=fontawesome&logoColor=white)

</td>
<td valign="top" width="33%">

**Deployment**

![Docker](https://img.shields.io/badge/Docker%20slim-2496ED?style=flat-square&logo=docker&logoColor=white)
![HuggingFace](https://img.shields.io/badge/Hugging%20Face%20Spaces-FFD21E?style=flat-square&logo=huggingface&logoColor=black)
![GitHub Pages](https://img.shields.io/badge/GitHub%20Pages-24292e?style=flat-square&logo=github&logoColor=white)

</td>
</tr>
</table>

---

## 🗄️ Integrated APIs & Databases

| Source | Function | Type |
|--------|----------|------|
| 🧬 **TargetScanHuman v8.0** | Thermodynamic prediction of binding sites | `Local (zip)` |
| 🔬 **miRTarBase** via Harmonizome | Experimental interaction validation (CLIP-seq, Luciferase, WB) | `Remote API` |
| 📊 **Enrichr** | Functional enrichment — KEGG / Reactome / WikiPathways / GO | `Remote API` |
| 📚 **NCBI PubMed eUtils** | Scientific articles per gene + context with cascade + retries | `Remote API` |
| 🏥 **NCBI OMIM** | Associated hereditary diseases by gene symbol | `Remote API` |
| 🧪 **NCBI ClinVar** | Reported pathogenic / likely pathogenic variants | `Remote API` |
| 💊 **ClinicalTrials.gov v2** | Active clinical trials filtered by gene symbol in title | `Remote API` |
| 🔗 **STRING-DB** | Protein-protein interaction network generation | `Remote API` |
| 🌐 **MyGene.info** | Clinical gene annotation | `Remote API` |
| 🔤 **MyMemory API** | Real-time ES↔EN biological term translation | `Remote API` |

---

## 📦 Installation

### Prerequisites

- Python 3.11+
- pip
- Docker *(optional — for containerized deployment)*

### Local Setup

```bash
# 1. Clone the repository
git clone https://github.com/abrangel/Kenryu.git
cd Kenryu

# 2. Create virtual environment (recommended)
python3.11 -m venv venv
source venv/bin/activate      # Linux / macOS
# venv\Scripts\activate       # Windows

# 3. Install dependencies
pip install -r requirements.txt

# 4. Start the server
uvicorn kenryu_engine:app --host 0.0.0.0 --port 7860 --reload

# 5. Open in browser → http://localhost:7860
```

### Docker Deployment

```bash
docker build -t kenryu .
docker run -p 7860:7860 -v $(pwd)/local_db:/app/local_db kenryu
```

> The `local_db` volume persists the analysis cache across container restarts.

### Hugging Face Spaces

The repo is fully ready for HF Spaces — the Dockerfile exposes port `7860` and mounts `local_db` for persistent caching. Simply connect the Space to this GitHub repository.

---

## 🚀 Usage

### Via Web Interface

1. Open [production Space](https://huggingface.co/spaces/Kenryu007/Bioinformatica) or `http://localhost:7860`
2. Enter miRNAs separated by commas: `hsa-miR-33a-5p, hsa-miR-144-3p, hsa-miR-758-3p`
3. Select year cutoff *(PubMed evidence age filter)*
4. Select consensus mode *(Strict / N-1 / N-2)*
5. Click **Execute**

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

<details>
<summary><b>📋 Example JSON Response</b></summary>

```json
{
  "common_genes": ["ABCA1", "KPNA3", "SCN1A"],
  "gene_research": {
    "ABCA1": {
      "pubmed": [{"pmid": "...", "term": "...", "year_window": "10y"}],
      "omim":   [{"omim_id": "600046", "title": "..."}],
      "clinvar": {"count": 201, "url": "..."},
      "trials":  [{"nct_id": "NCT01456650", "title": "..."}]
    }
  },
  "enrichment": ["..."],
  "scientific_synthesis": "...",
  "venn_plot":    "data:image/png;base64,...",
  "volcano_plot": "data:image/png;base64,...",
  "ppi_plot":     "data:image/png;base64,..."
}
```

</details>

---

## 📁 Repository Structure

```
Kenryu/
├── kenryu_engine.py              # FastAPI backend + bioinformatics logic (~1700 lines)
├── Dockerfile                    # HF Spaces / Docker configuration
├── requirements.txt
│
├── static/
│   ├── index.html                # A4 editor main interface
│   ├── script.js                 # Frontend logic (pagination, PDF/MD export)
│   └── style.css                 # Dark mode · gold/teal accents
│
├── data/
│   ├── targetscan_full.json.zip  # TargetScan v8.0 indexed DB (Git LFS, ~8 MB)
│   └── hsa-miR-*.txt             # Pre-processed per-miRNA example files
│
└── local_db/                     # Runtime cache (auto-created, not in repo)
    └── analysis_cache.json
```

---

## 🗺️ Roadmap

### ✅ Implemented

- [x] Dynamic Venn convergence — Strict / N-1 / N-2 modes
- [x] Balanced functional enrichment — KEGG · Reactome · WikiPathways · GO
- [x] Visualizations — Venn · Volcano · STRING-DB interactome
- [x] Vancouver-style bibliography with multi-source citation types
- [x] Multi-source research — PubMed · OMIM · ClinVar · ClinicalTrials
- [x] Versioned persistent cache across sessions
- [x] 429-aware retries with exponential backoff for NCBI eUtils
- [x] ES→EN mapping for effective PubMed queries
- [x] Professional Markdown ZIP export with separate image assets
- [x] A4 PDF footer correctly anchored across all pages
- [x] Integrated genomic evidence narrative in report body

### 🔜 Planned

- [ ] Support for miRNA isomers (isomiRs)
- [ ] DisGeNET integration for gene-disease association
- [ ] Structured JSON export for interoperability with other tools
- [ ] Side-by-side miRNA panel comparison
- [ ] Cross-expression heatmap (genes × miRNAs)
- [ ] Batch mode — process multiple panels from CSV
- [ ] TargetScan v9.0 integration when available

---

## 🤝 Contributing

Contributions are welcome. To report bugs, suggest features, or submit pull requests:

1. **Fork** the repository
2. Create a branch: `git checkout -b feat/feature-name`
3. Commit with descriptive messages following conventions below
4. Push: `git push origin feat/feature-name`
5. Open a **Pull Request**

**Commit conventions:** `feat:` · `fix:` · `docs:` · `refactor:` · `chore:`

When reporting bugs please include: steps to reproduce, expected vs actual behavior, browser console logs (frontend errors), and HF Space logs (backend errors).

---

## 📜 License

This project is distributed under a **Free Academic License for Research and Education**. For commercial use, contact the author.

KENRYU integrates data from public sources (NCBI, Enrichr, OMIM, ClinVar, ClinicalTrials.gov, STRING-DB) subject to their respective terms of service. Users are responsible for complying with each database's policies.

---

## 👤 Author

**Cesar Manzo** — Clinical Bioinformatics · Genomic Analysis · Translational Medicine

[![Hugging Face](https://img.shields.io/badge/🤗%20Kenryu-Hugging%20Face-c8a96e?style=flat-square&labelColor=1a1d24)](https://huggingface.co/spaces/Kenryu007/Bioinformatica)
[![GitHub Pages](https://img.shields.io/badge/🌐%20kenryu-GitHub%20Pages-4fc3a1?style=flat-square&labelColor=1a1d24)](https://abrangel.github.io/Kenryu)
[![GitHub](https://img.shields.io/badge/@abrangel-GitHub-6eaadc?style=flat-square&logo=github&labelColor=1a1d24)](https://github.com/abrangel)

---

## 🙏 Acknowledgements

- **TargetScanHuman** — Lewis Lab, Whitehead Institute
- **miRTarBase** — Chou *et al.*, Nucleic Acids Research
- **NCBI** — National Center for Biotechnology Information
- **Enrichr** — Ma'ayan Lab, Mount Sinai
- **ClinicalTrials.gov** — U.S. National Library of Medicine
- **STRING-DB** — European Molecular Biology Laboratory

---

<div align="center">

<img src="https://capsule-render.vercel.app/api?type=waving&height=100&color=0:0d0f14,50:c8a96e,100:0d0f14&section=footer" width="100%"/>

<sub><b>KENRYU Bioinformatics Engine</b> · 2026<br/>
Made with scientific rigor for the genomic research community.</sub>

</div>
