# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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

### Changed
- CI: removed the dead Rust toolchain step (gseapy is gone) and added
  coverage reporting (`pytest --cov=kenryu_engine`).
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
