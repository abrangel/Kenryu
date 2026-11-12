# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

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
