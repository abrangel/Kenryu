# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.44] — 2026-11-19

### Fixed
- **La síntesis imprimía "(PMID: None)":** al quitar los PMIDs escritos a
  mano en v1.43, el informe quedó con 11 ocurrencias literales de
  "PMID: None" junto a cada descripción de ruta. Ya no se imprime ningún
  PMID junto a la descripción (ver el cambio de abajo).

### Changed
- **Texto y cita separados (auditoría externa, opción 1):** la descripción
  escrita a mano de cada ruta se pegaba al primer artículo que devolvía
  PubMed para gen + término, como si ese artículo la respaldara — buscar
  por término no garantiza que el primer resultado diga lo que dice la
  frase. Ahora van separados y se dice así: la descripción se muestra como
  "Notas editoriales sobre las rutas (contexto del equipo, sin cita)" y el
  artículo en vivo como "Primeros resultados en PubMed (búsqueda por gen +
  término; no se verificó que respalden la nota anterior)". En la API, el
  campo `desc` del registro de PubMed pasa a llamarse `editorial_note`. El
  texto del método del dossier lo explica. La opción 2 (eliminar la
  descripción y mostrar solo fuentes) queda para la ronda del miRTarBase
  oficial.

### Docs
- **Cifras de v1.42 matizadas:** el CHANGELOG y las notas de la release
  decían "SCN1A 60" como si fuera un valor fijo del parche. Es la corrida
  de verificación: en vivo los puntajes varían unos puntos entre corridas
  según PubMed (el Space dio SCN1A 65 en una corrida, dos artículos
  distintos sin duplicar). Lo que el parche garantiza es que no hay doble
  conteo, no un puntaje exacto.
- **Comentario del CI sin mención del incidente:** el escaneo de secretos
  sigue en el CI; el comentario ya no anuncia la clave purgada (revocada y
  fuera del historial; el mensaje del commit que añadió el escaneo conserva
  el motivo).

## [1.43] — 2026-11-19

### Fixed
- **PMIDs escritos a mano eliminados (auditoría externa):** los 7 genes con
  texto fijo (ABCA1, KPNA3, SCN1A, SNTB2, TSC22D2, GMFB, KLF12) llevaban 11
  PMIDs escritos a mano y TODOS resolvían en PubMed a artículos sin relación
  con el gen ni con la frase (p. ej. PMID 25648732 es un estudio sobre
  consejería de embarazo, citado para el ensamblaje de HDL de ABCA1 y para
  MAPK/ERK de GMFB; PMID 29104567 es microbioma de salmón, citado para
  KPNA3 y KLF12). Los PMIDs quedan en `None`: el panel del gen ya no muestra
  el enlace "VERIFICAR EVIDENCIA" con un artículo falso. La evidencia PubMed
  verificada la aporta la capa viva de investigación, no estos presets.
  Test de regresión: ningún gen preset devuelve un PMID escrito a mano.

### Docs
- **README.md y README.es.md actualizados:** la hoja de ruta marcaba como
  "planeado" lo ya existente (isomiRs, DisGeNET, exportación JSON, comparación
  de paneles, heatmap, batch) — ahora está en Implementado, con la hoja de
  ruta futura (miRTarBase oficial, fondo de enriquecimiento, Open Targets,
  miRDB, GTEx, red interactiva, GSEA, TargetScan v9.0). **Licencia
  corregida:** el README decía "Licencia Académica" pero el archivo LICENSE
  del repo siempre ha sido MIT y GitHub reporta MIT; el README ahora dice MIT.
  Sección nueva "Known Limitations / Limitaciones conocidas" con el hueco del
  snapshot local de miRTarBase (hsa-miR-758-3p sin dianas locales) y el fondo
  de enriquecimiento. Tabla de fuentes y estructura del repo actualizadas
  (miRTarBase y Enrichr ahora son GMT locales + API en vivo).

## [1.42] — 2026-11-18

### Fixed
- **Deduplicación por PMID (auditoría externa):** `research_gene_complete`
  buscaba en PubMed una vez por cada término de evidencia y se quedaba con el
  primer resultado de cada uno, así que el mismo artículo podía entrar dos
  veces y sumar 5 puntos dos veces. Verificado en vivo: KPNA3 (PMID 32492291),
  SCN1A (36153432), SNTB2 (37762054) y TSC22D2 (35100036) tenían el mismo PMID
  duplicado. Ahora `dedupe_pubmed_refs()` deduplica por PMID antes de puntuar
  (la deduplicación que existía solo afectaba a la bibliografía del informe,
  no al puntaje). Caché de investigación versionada v2 → v3 para invalidar
  los resultados viejos con duplicados. Panel de la tesis tras el arreglo
  (corrida de verificación; los puntajes en vivo varían unos puntos entre
  corridas según lo que devuelva PubMed ese día): ABCA1 100 (alta), SCN1A 60,
  KPNA3 50, SNTB2 50, TSC22D2 50 (media). El Space en vivo dio SCN1A 65 en una
  corrida — dos artículos distintos, sin duplicar —, así que el 60 no es un
  valor fijo del parche: lo que el parche garantiza es que no hay doble
  conteo, no un puntaje exacto.

### Changed
- **El nivel "alta" exige validación experimental de la interacción:** un gen
  con 70+ puntos por literatura y clínica (PubMed, OMIM, ClinVar, ensayos)
  pero sin miRTarBase ya no es "alta" — queda en "media" y la explicación lo
  dice. El puntaje numérico no cambia; solo el nivel. La literatura describe
  lo importante que es el gen, no que el microARN lo regule. El panel de la
  tesis no cambia de niveles (ABCA1 sigue alta porque SÍ tiene miRTarBase).

### Known limitations
- El GMT local de miRTarBase (snapshot Harmonizome) no contiene ninguna
  entrada para hsa-miR-758-3p, aunque la literatura reporta dianas validadas
  (ABCA1, MDM2, MTOR, NOTCH2). El panel de la tesis no se ve afectado
  (ABCA1 se valida vía miR-33a-5p), pero el snapshot subestima la evidencia
  experimental de 758-3p. Pendiente: refrescar el GMT desde el miRTarBase
  oficial cuando su descarga sea alcanzable.

## [1.41] — 2026-11-18

### Added
- Protein relevance (UniProt, public API, no key): every core gene now
  carries its protein — name, UniProt accession, what it does, where it
  works in the cell, associated disease — in a new dashboard card
  ("Relevancia de las proteínas") and in the gene side-panel. Examples:
  ABCA1 → O95477, Tangier disease; SCN1A → P35498, GEFS+; KPNA3 → O00505,
  Spastic paraplegia 88. Offline, the card honestly says there is no data.
- One-click examples on the input panel: Lípidos (tesis), Cáncer,
  Control — they fill the input and run the analysis immediately.
- PPI now includes the first shell of STRING neighbours (add_nodes=5):
  the thesis two-module picture reappears (dystrophin module SNTA1/SNTB1/
  SNTB2/UTRN/DTNB around ABCA1–SCN1A, and TSC22D2–FAM117B). Neighbour
  proteins are drawn grey; core genes blue.

### Fixed
- The reverse lookup moved BELOW the results grid so it no longer
  clutters the front of the page.
- PPI crash on an isolated core gene: KPNA3 has no STRING edges even with
  neighbours, and the size lookup raised KeyError (caught by the live run
  before publishing; regression test added).

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
