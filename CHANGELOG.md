# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.40] — 2026-11-18

### Added
- Every dashboard card now explains itself in plain language (how to read
  the Venn, volcano, heatmap, PPI, enrichment and validation table), so a
  newcomer understands each graphic without a manual.
- Reverse lookup on the page: type a gene (e.g. ABCA1) and see the 55 miRNA
  families that target it, with the panel's own families highlighted
  (uses the existing `GET /api/v1/reverse`).
- Heatmap now covers the UNION of targets (top 20 by miRNA count) instead
  of only the core genes — the old matrix was all 1s by construction.

### Fixed
- The volcano plot was a SIMULATED illustration (200 random points): it now
  plots real data — X = PubMed article count (log10), Y = number of panel
  miRNAs regulating the gene, top 25 union genes, core genes in red and
  labelled. Without live PubMed data (offline) it is honestly absent.
- The PPI was a decorative ring with invented edges: it now draws only real
  STRING-DB interactions (confidence ≥ 0,4) between the core genes, with
  edge width = confidence. All core genes are shown as nodes, including
  isolated ones. Finding: among the 5 thesis genes STRING records a single
  direct interaction (SNTB2–ABCA1, 0.836) — the thesis two-module picture
  came from a query that added neighbouring proteins.

## [1.39] — 2026-11-18

### Added
- Self-updater: `GET /api/v1/updates` (version check) + `POST /api/v1/updates/{source}/install`
  (download, verify, backup, install, record). miRBase 23.1 installed live.
- Validation dossier in every `/analyze` response: per-gene 0–100 score from
  independent layers (TargetScan, miRTarBase, PubMed, OMIM, ClinVar,
  ClinicalTrials) + the exact database versions used.
- miRTarBase served LOCALLY (Harmonizome GMT, 596 miRNAs); live API only a fallback.
- Local enrichment (KEGG/Reactome/WikiPathways/GO GMTs + Fisher + Benjamini-Hochberg).
- Offline mode (`offline=true`): zero network calls; same core genes as live.
- Structured JSON export (schema 1.0) + genes × miRNAs heatmap in every response.
- `GET /api/v1/reverse` (gene → miRNA families), isomiR normalization,
  `kenryu` CLI, `POST /api/v1/compare` + `/api/v1/batch`, DisGeNET layer
  (free API key), mouse species (`species=mouse`, TargetScanMouse 8.0).
- Web UI: validation table + heatmap card on the dashboard; the gene panel
  now shows the real validation score instead of a fixed number.
- `pyproject.toml` (`pip install -e .` → `kenryu` command); Dockerfile fixed
  (it copied only the engine before — the image could not import the core).

## [Unreleased]

### Added
- `CITATION.cff` so the software can be cited from GitHub.
- Cancer-panel case study (`docs/case-study-cancer-panel.md`): second
  computational validation (18-gene Strict core, live run 2026-11-11).
- `kenryu_core.py`: pure helpers extracted from the engine
  (`deep_clean_mirna`, `clean_mirna_name_cloud`, `consensus_genes`,
  `create_visuals`, `build_synthesis`), unit-testable in isolation.
- `consensus_genes()` pure function for the Strict / N-1 / N-2 modes,
  extracted verbatim from `analyze()`.
- API validation: `POST /api/v1/analyze` now rejects an empty miRNA panel
  with 422 instead of running a meaningless analysis.
- Comparison with multiMiR, miRWalk 3.0, DIANA-miRPath v3, and miRNet in
  `paper.md` (State of the field), with verified citations.

### Fixed
- `clean_mirna_name_cloud` kept the let-7 family name: `let-7a-5p` was
  converted to `hsa-miR-7a-5p`, so the remote miRTarBase lookup queried a
  name that does not exist. Now returns `hsa-let-7a-5p` (regression tests
  added).
- The TargetScan database (`data/targetscan_full.json.zip`, 8.3 MB) was
  only a Git LFS pointer in the repository: a fresh clone without
  `git lfs pull` silently fell back to the 5 example files instead of the
  285 families. The file is now stored directly in git (SHA-256 verified
  against the LFS object) and `.gitattributes` / the CI `lfs: true` step
  are gone, so any clone reproduces the full local analysis.

### Changed
- CI: removed the dead Rust toolchain step (gseapy is gone) and added
  coverage reporting (`pytest --cov=kenryu_engine --cov=kenryu_core`).
  Current coverage: 34% total (kenryu_core 75%, engine 17% — the
  network-dependent functions are not covered offline).
- Test suite grew from 22 to 31 tests, all offline.

## [1.38] — 2026-07-04

Current release. Deployed at https://huggingface.co/spaces/Kenryu007/Bioinformatica.

### Added
- Strict / N-1 / N-2 consensus modes for genomic convergence.
- Balanced functional enrichment across KEGG, Reactome, WikiPathways, and GO Biological Process.
- Multi-source gene research: PubMed, OMIM, ClinVar, ClinicalTrials.gov (parallel queries).
- Automated academic report with narrative synthesis and Vancouver bibliography (PDF, Markdown ZIP, TXT).
- Visualizations generated per analysis: Venn diagram, volcano plot, STRING-DB interactome.
- Persistent cache (`local_db/analysis_cache.json`) for translations, PubMed results, and enrichment.
- 429-aware retries with exponential backoff (2–32 s) for NCBI eUtils.
- Real-time ES→EN translation of biological terms for PubMed queries.
- A4 WYSIWYG report editor with smart pagination and anchored PDF footer.
- REST API: `POST /api/v1/analyze`, `GET /health`.
- Docker deployment and Hugging Face Spaces support (port 7860).

### Validated
- Lipid panel (hsa-miR-33a-5p, hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p, hsa-miR-33b-5p): 5-gene Strict core (ABCA1, KPNA3, SCN1A, SNTB2, TSC22D2), documented in the UNIR master's thesis (2026-03-29).
- Example panel (hsa-miR-375, hsa-miR-29a-3p, hsa-miR-21-5p): 6-gene Strict core (GMFB, KLF12, NFIA, PDE3A, SKI, SOX6), live test 2026-11-11.
- Cancer panel (hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p): 18-gene Strict core, 24 enrichment terms, 3 significant (p-adj < 0.05), live test 2026-11-11.

## [1.0] — 2026-05-14

Initial public release: FastAPI engine with TargetScan v8.0 (local), miRTarBase via Harmonizome, Venn convergence, functional enrichment, and automated reporting.

[1.38]: https://github.com/abrangel/Kenryu/releases/tag/v1.38
[1.0]: https://github.com/abrangel/Kenryu/releases/tag/v1.0
