# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.47] — 2026-11-19

### Ronda 2, Commit A: afirmaciones y consistencia (H-12, H-18, H-11, H-17, H-16, H-21)
### H-04 NO está aquí: es decisión del usuario (corregir el texto o cambiar el código a evidencia dual), pendiente.

### Consenso (H-17)
- `consensus_genes` ya no lanza `TypeError` con lista vacía (devuelve `set()`).
- Los miRNAs duplicados en el panel no sesgan el umbral de N-1/N-2: se
  deduplican preservando el orden, en `analyze()` y en `consensus_genes`.
- La respuesta lleva `missing_mirnas` (miRNAs sin datos en ninguna fuente,
  antes descartados en silencio) y `single_mirna_warning` cuando el panel
  tiene un solo miRNA con datos (no hay convergencia multi-miRNA). Ambos
  se muestran en la interfaz.
- `mode` fuera de {strict, n-1, n-2} ya devolvía 422 desde v1.45 (test de
  regresión añadido).

### Puntaje de confianza fijo (H-18)
- Eliminado el campo `confidence` (85/88/92/95) de `get_gene_details`.
  La barra "Rigor Científico" usa solo el puntaje real del dossier de
  validación; si el gen es "no evaluado" muestra "Sin evaluar", sin barra.

### Síntesis académica (H-12)
- Sin `clinical_evidence`, la síntesis ya no afirma "consecuencias
  directas en la estabilidad tisular" ni cita con un índice `[i+1]` que
  no corresponde a nada: es una predicción computacional, pendiente de
  validación experimental.
- La etiqueta de significancia dice solo significancia estadística:
  nunca "significancia funcional" ni "equilibrio fisiológico".
- La conclusión genérica de un gen sin preset ya no lo llama "punto
  crítico de intervención": lo presenta como candidato para validación
  experimental.
- NO se tocó la frase de la introducción ("intersección de TargetScan y
  miRTarBase"): es parte de H-04, pendiente de la decisión del usuario.

### Bibliografía (H-11)
- `paper.bib` corregido contra Crossref (y PubMed para Sayers): Lewis
  2005 (título y autor — "Benjamin P.", no "Christopher B." — clave
  renombrada a `lewis2005conserved`), Huang 2022 (35 autores completos,
  clave `huang2022mirtarbase`), Sayers 2022 (páginas D20-D26, no
  D11-D18; 25 autores), Ru 2014 (autor "Yuanbin", no "Yi"; "Katerina",
  no "Katherina"; número 17, no 21), Sticht 2018 (De La Torre
  "Carolina", no "Cecilia"). Añadido McGeary et al. 2019 (TargetScan v8,
  Science 366:eaav1741) y citado en `paper.md`. Las claves renombradas
  se actualizaron en `paper.md`. Las entradas no marcadas por la
  auditoría no se tocaron (no verificadas en esta ronda).
- El DOI de Zenodo (10.5281/zenodo.20564689) SÍ resuelve (la auditoría
  solo consultó Crossref, que no indexa DOIs de DataCite).

### Símbolos génicos (H-16)
- `reverse_lookup` ya no distingue mayúsculas: `C8orf58` se encuentra
  aunque TargetScan lo guarde con minúsculas.
- Documentado en ambos README (Limitaciones conocidas): no hay mapeo de
  alias HGNC, y los FASTA locales de miRBase (14 MB) los descarga el
  actualizador pero no se usan para validar nombres.

### Exportación Markdown (H-21)
- Los miRNAs del frontmatter YAML se serializan con `JSON.stringify`
  (una comilla o un salto de línea ya no rompe el YAML).
- La versión ya no es un literal fijo `v1.38`: la respuesta de `/analyze`
  lleva `engine_version` (constante única `APP_VERSION`) y la
  exportación y la portada la leen de ahí.

### Regresión
- Panel de la tesis, sin cambio: ABCA1 100 alta, SCN1A 60 media,
  KPNA3/SNTB2/TSC22D2 50 media. Panel de cáncer, sin cambio:
  investigados AGO4, MEIS1, PURB, RPS6KA3, WNK3 (80 alta), 13 "no
  evaluado". La versión de la caché no sube: ningún resultado cambia.

### Ronda 2, Commit B: seguridad del Space (H-08, H-14)

### XSS (H-08)
- Nueva función `escapeHtml` en `static/script.js` y aplicada a todos
  los datos externos y de usuario antes de `innerHTML`: gen buscado en
  la búsqueda inversa (reflejaba la entrada del usuario), genes del
  consenso (ahora en `data-gene`, no en un `onclick` interpolado),
  tabla de validación, tabla de proteínas, términos de enriquecimiento,
  mensajes de error, ficha del gen (nombre, sistema, patología,
  proteína, rutas, PMID), esquema del informe, párrafos de síntesis y
  contexto funcional, referencias (se escapa antes de enlazar las URLs)
  y tabla de genes del informe.
- Los PMID solo enlazan si son numéricos; la URL de UniProt solo se
  renderiza si es `https://rest.uniprot.org/`; las imágenes solo se
  insertan si el base64 es válido.

### CORS y cabeceras (H-14)
- CORS: solo el origen exacto del Space
  (`https://kenryu007-bioinformatica.hf.space`); eliminado el regex
  `https://.*\.hf\.space` y `allow_credentials=True` (la app no usa
  cookies de sesión: `allow_credentials=False`).
- Middleware de cabeceras: `X-Content-Type-Options: nosniff`,
  `Referrer-Policy: strict-origin-when-cross-origin` y
  `Content-Security-Policy` (permite los CDN que usa la interfaz —
  cdnjs, Google Fonts — y `frame-ancestors` solo para Hugging Face, de
  modo que la vista del Space sigue funcionando). La CSP lleva
  `'unsafe-inline'` en scripts y estilos porque la UI usa manejadores
  `onclick` y estilos en línea; la defensa contra XSS es `escapeHtml`,
  no el bloqueo de scripts.

### Tests
- 120/120 (4 nuevos en Commit B: CORS rechaza otro Space, CORS acepta
  el origen exacto sin credenciales, cabeceras de seguridad presentes,
  escape en innerHTML).

## [1.46] — 2026-11-19

### Ronda 1, último hallazgo: H-06 (sesgo alfabético)

La investigación por gen y el enriquecimiento se recortaban por orden
alfabético: en el panel de cáncer (18 genes) solo se investigaban AAK1,
ACVR2A, AGO4, BNC2 y BRWD3, y los 13 restantes salían como "baja" con
25/100 por no haberse investigado, no por falta de evidencia.

### Fixed
- **Orden de prioridad en vez de alfabético (H-06, alto):** nueva
  función pura `prioritize_genes` — primero los genes que regulan más
  miRNAs del panel y, en empate, los validados en miRTarBase por más
  miRNAs (el nombre solo desempata). El mismo orden se usa en la
  investigación (top 5), el enriquecimiento (top 300), los gráficos
  (top 12) y los detalles de genes (top 40).
- **Genes no investigados: "no evaluado", nunca "baja":** el informe
  de validación marca `level: "not_evaluated"` con `score: null` a los
  genes fuera del tope; la media y los conteos de alta/media/baja los
  excluyen. La interfaz muestra "no evaluado" sin puntaje (en gris) y
  el resumen "N sin evaluar".
- **Tope por costo de API, declarado en la respuesta:** la investigación
  multi-fuente sigue limitada a 5 genes por corrida (constante
  `RESEARCH_CAP`), pero ahora la respuesta incluye `research_scope`
  (tope, genes investigados, genes no evaluados y la nota) y el
  enriquecimiento avisa con `enrichment_truncated` + `enrichment_note`
  cuando usa solo los 300 priorizados.
- **Caché de Enrichr versionada (hallazgo propio en la regresión):**
  la caché guardaba términos con el corte viejo de p-adj < 0,3 y se
  seguían sirviendo (el panel de la tesis devolvía 24 términos con
  p-adj máximo 0,0584). La clave pasa a `enr-v2|` y las entradas viejas
  se ignoran. Caché de investigación v3 → v4.

### Resultados, antes y después (panel de cáncer, 18 genes, en vivo)
- Antes (v1.45): investigados AAK1, ACVR2A, AGO4, BNC2, BRWD3
  (alfabéticos); alta 1, media 8, baja 9, media 40,3. Ocho genes
  (CHD7, DCUN1D3, FOXP1, IL6R, KPNA4, RDX, SESN3, UBN2, XKR4) salían
  "baja" sin haberse investigado.
- Después (v1.46): investigados AGO4, MEIS1, PURB, RPS6KA3, WNK3
  (los 5 con validación en miRTarBase) — los 5 con 80/100 "alta";
  13 "no evaluado" sin puntaje; media 80,0 sobre los 5 evaluados;
  baja 0. Enriquecimiento sin cambio: 3 términos significativos
  (0,0176 / 0,0424 / 0,0424).
- Panel de la tesis, sin cambio (regresión): ABCA1 100 alta, SCN1A 60
  media, KPNA3/SNTB2/TSC22D2 50 media, media 62,0; enriquecimiento
  19 términos, p-adj máximo 0,0349 (todos < 0,05).

### Ajustes pre-publicación (revisión externa)
- **La media es solo de los evaluados y la etiqueta lo dice:** el
  resumen lleva `mean_score_label` — "Media de los genes evaluados
  (5 de 18)" — y la interfaz la muestra con esa etiqueta. La subida
  de 40,3 a 80,0 no es una mejora de evidencia: los 5 investigados se
  eligieron por tener validación en miRTarBase, así que son "alta" por
  construcción. El orden de prioridad (nº de miRNAs del panel que
  regulan el gen y validación en el snapshot local de miRTarBase,
  Harmonizome 2026-11-18) se indica en la nota del tope, en el texto
  del método del dossier y aquí.
- **Fuente por miRNA en la respuesta:** `source_by_mirna` dice, por
  cada miRNA, si se usó TargetScan (predicción local) o miRTarBase
  como respaldo (validación experimental, porque TargetScan no tenía
  ese miRNA). Un consenso strict puede mezclar ambos tipos de
  evidencia y ahora se ve. La interfaz lo muestra bajo el dossier.
  (Conecta con H-04 y H-12, ronda 2.)

### Tests
- 103/103 (9 nuevos en v1.46: orden de prioridad, "no evaluado" en el
  informe, tope en `/analyze` offline con el panel de control y el de
  cáncer, truncado del enriquecimiento en N-2, caché de Enrichr
  versionada, etiqueta de la media, nota del tope con la base de
  prioridad y el snapshot, fuente por miRNA).

## [1.45] — 2026-11-19

### Ronda 1 de la auditoría externa (23 hallazgos; 8 corregidos aquí)

Se verificó cada hallazgo contra el código y en vivo antes de tocarlo.
Dos afirmaciones de la auditoría NO se confirmaron: el DOI de Zenodo
`10.5281/zenodo.20564689` sí resuelve (200 en doi.org y registro activo
en la API de Zenodo; el 404 era solo de Crossref, que no indexa DOIs de
DataCite) y la primera prueba del matching difuso dio 0 genes porque la
DB no estaba cargada en ese proceso — con la DB cargada, el bug sí se
reprodujo exactamente como decía la auditoría.

### Security
- **El endpoint `POST /api/v1/updates/{source}/install` ya no es
  anónimo (H-01, crítico):** cualquiera podía descargar y sobrescribir
  archivos del servidor. Ahora exige el secreto `KENRYU_ADMIN_TOKEN` en
  la cabecera `X-Admin-Token` (comparación en tiempo constante); sin la
  variable configurada, las actualizaciones están desactivadas (403). La
  UI no llamaba a este endpoint, así que nada cambia en la web. Tras
  instalar, las DB en memoria se invalidan y se recargan en la siguiente
  petición.

### Fixed
- **Matching difuso por subcadena (H-03, crítico):** `hsa-miR-1-5p`
  devolvía 793 genes de miR-204/211, `hsa-miR-9-3p` 1.385 de
  miR-17/20/93/106/519, etc. El número solo vale ahora en límite de
  token (tras `mir-`/`let-`/`/` y antes de `/`/`.`/fin): un miRNA sin
  familia real devuelve vacío con aviso, y la variante real
  `hsa-miR-873-5p` sigue encontrando `mir-873-5p.2`.
- **Fisher + Benjamini-Hochberg (H-05, alto):** BH se calculaba solo
  sobre términos con solape (m ~28× menor que el total en GO) y `genes`
  no se recortaba al universo, con falsos positivos en 40/40 listas
  nulas a p_adj < 0,3. Ahora todos los términos no vacíos de la
  biblioteca entran en BH (los de solape 0 con p = 1,0), los genes se
  recortan al universo y el corte de significancia es **p_adj < 0,05**
  (también en la ruta viva de Enrichr, que usaba 0,3).
- **La UI enviaba el año calendario como `years` (H-07, alto):** el
  selector mandaba `2026` y el backend calculaba 2026−2026 = 0, anulando
  el filtro temporal de PubMed. Ahora la UI convierte el año de inicio a
  edad en años (por defecto, ventana de 10 años) y el backend valida
  `years` entre 1 y 50. El selector de mes se eliminó de la UI: el
  backend nunca lo usaba.
- **Límites y validación de entradas (H-10, alto):** `/analyze` aceptaba
  2.000 miRNAs (22,5 s de CPU); `/batch` leía archivos sin límite
  (200.000 filas, 22,3 s); `/compare` devolvía 500 con `mirnas: [1, 2]`
  y aceptaba 300 paneles. Ahora: panel ≤ 50 miRNAs, `mode` es un
  literal (`strict`/`n-1`/`n-2`), `/compare` tipa cada panel (≤ 30
  miRNAs, ≤ 20 paneles) y `/batch` rechaza archivos > 256 KB o > 500
  filas.
- **BOM UTF-8 en 4 archivos de datos (H-13, medio):** `\ufeffANO4`,
  `\ufeffCTDSPL`, `\ufeffZNF772` corrompían símbolos génicos. Los
  archivos se limpiaron y la lectura usa `utf-8-sig`.
- **La caché de los tests ya no toca el archivo real (H-15, medio):**
  un test de regresión había sobrescrito `local_db/analysis_cache.json`
  con una sola entrada TESTGENE. `tests/conftest.py` redirige la caché
  a un directorio temporal y el archivo queda fuera de git (es estado de
  runtime, no contenido del repo).

### Changed
- **KEGG fuera del repositorio (H-02, crítico):** Kanehisa Labs prohíbe
  redistribuir KEGG sin acuerdo. Se eliminó `local_db/gmt_KEGG_2021_Human.gmt`
  de git y del cargador de GMTs; el enriquecimiento en vivo sigue
  consultando KEGG vía Enrichr y el modo offline usa Reactome,
  WikiPathways y GO. Nuevo `DATA_LICENSES.md`: la licencia MIT cubre
  solo el código; cada dataset conserva la suya.
- **Dependencias (H-09, alto):** `fastapi` 0.109.0 → 0.165.0 y
  `python-multipart` 0.0.6 → 0.0.34 (CVEs parcheados; `pydantic`
  2.5.3 → 2.14.3 por compatibilidad). `biopython` y `jinja2` eliminados:
  no se usaban (Biopython solo ponía `Entrez.email`, que nunca se
  llamaba).
- **Licencia en `pyproject.toml` alineada a MIT (H-23):** seguía
  diciendo "Free Academic License for Research and Education".

### Tests
- 94/94 (8 nuevos: matching difuso, Fisher/BH ×3, validación de
  `/analyze`/`/compare`/`/batch`, token de instalación).

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
