# Changelog

All notable changes to KENRYU are documented here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Fixed
- **Version reconciliation (paper + READMEs + controls).** The paper
  and READMEs described the data layer as the miRTarBase 10.0 index,
  but the public service runs the 9.0 index and the 10.0 index exists
  only on the author's machine (its redistribution terms are
  unverifiable — the site returns 404). The paper now reports both
  control sets — 6 of 8 verdicts on both, TargetScan sensitivity 405
  of 754 (54%) on 9.0 and 358 of 643 (56%) on 10.0 — and says a
  stranger can verify the 9.0 figures on the public service. The 9.0
  controls were re-run with engine 2.0
  (`docs/validation/strong-controls.json`; one target moved, 404 →
  405, via CEP83 recovered by the HGNC map). The READMEs now state
  the public service runs 9.0 and 10.0 is local-only; the stale
  "13,035 human target genes" is corrected to 13,016 (post-HGNC-map
  universe) in the paper and both READMEs. Paper word count 1,663
  after trimming detail into `docs/` (limit 750–1,750).

### Added
- **CI on Windows and macOS.** The `test` and `smoke` jobs now run
  on a matrix of ubuntu-latest, windows-latest and macos-latest
  (H-20 was only ever verified by a simulated-Windows test). The
  smoke steps run under bash on all three OSes.
- **Dependency audit in CI.** New `pip-audit` job (H-09): fails the
  build on any known vulnerability in the pinned requirements.
- **Repo-hygiene regression tests** (`tests/test_repo_hygiene.py`):
  no KEGG/Harmonizome GMT or derived miRTarBase index may be tracked
  in git (H-02), only the three allowlisted GMTs may be tracked, and
  no fixed confidence literal (85/88/92/95) may sit next to
  "confidence" in the engine or the UI (H-18). Each guard was
  verified fail-first.
- **Reviewer path** (`examples/`): the three reference panels with
  their expected results and a single verification command
  (`python examples/verify_examples.py`). It prints the miRTarBase
  version it used per panel (9.0 / 10.0 / legacy GMT / not loaded)
  and, without an index, verifies the consensus genes in reduced mode
  instead of comparing evidence levels. The command also runs in CI
  (smoke job, before the synthetic index is copied in).

### Fixed
- **Windows: tests read text files as UTF-8.** The first Windows CI
  run failed in 8 tests: `Path.read_text()` without an encoding
  defaults to cp1252 on Windows, which cannot decode the Spanish UI
  files (`static/index.html`, `static/script.js`) and silently
  mojibakes `report.md` and `docs/openapi.json`. All text reads in
  the tests now pass `encoding="utf-8"`; the engine and CLI already
  did. The smoke job's output directory is now relative instead of
  `/tmp`, which Python on Windows does not resolve.

## [2.0] - 2026-11-21

### Added
- **miRTarBase 10.0 support (major version).** The official index is
  now built from `hsa_MTI.csv` (miRTarBase 10.0, Cui et al., *Nucleic
  Acids Res* 2025;53(D1):D147–D156, doi 10.1093/nar/gkae1072):
  4,004,967 rows, 1,730,805 human pairs, 3,048 miRNAs (vs 396,758
  pairs / 2,947 miRNAs in 9.0). The engine loads the v10.0 index
  (`local_db/mirtarbase_v10.json.gz`) when its file is present and
  falls back to the 9.0 index otherwise; `KENRYU_MIRTARBASE_VERSION`
  forces `9.0` or `10.0` (reproducing old results needs 9.0). The
  loaded version is fixed in the execution manifest, in /health and
  in the dossier. The v10.0 adoption is documented step by step in
  `docs/validation/mirtarbase-v10-adoption.md` (licence, level-change
  matrix, base rate, levels, panels, performance).
- **New evidence level "clip_only".** miRTarBase 10.0 adds 1,411,229
  weak pairs with NO PMID — all of them CLIP-only (measured:
  82.0% of all weak pairs, 81.5% of the whole index; in v9.0 that
  figure was 0). They are shown in the dossier ("Solo CLIP sin PMID
  … visible, sin puntos, no desbloquea «alta»") with +0 points and
  never unlock "alta". Strong evidence (with PMID) still scores +30
  and weak evidence with PMID +15. The constructor is fail-first:
  an unknown Support Type, a Functional MTI row without PMID, or a
  weak row without PMID whose experiments are not all CLIP all raise
  instead of being silently mapped; a PMID on a weak row never
  rescues a strong claim (per-support PMID sets are tracked
  separately). 11 tests in `tests/test_mirtarbase_levels.py`.
- **Bare miRNA names are marked.** Evidence that comes from a
  miRTarBase key without a -5p/-3p arm is labelled "brazo no
  especificado" in the dossier table (whether the query fell back
  to the bare key or the query itself was bare, e.g. the control
  panel's hsa-miR-375). v10.0 has 1,212 bare-name miRNAs (59.9% of
  its pairs, vs 26.7% in v9.0).

### Changed
- **Reference panels, before/after (each change explained gene by
  gene in the adoption doc).** Lipids: identical (5 genes; the 8
  benchmark pairs and the ABCA1 five-miRNA check unchanged — the
  thesis panel is untouched). Cancer: 19 genes; AGO4 best stays
  strong (the miR-155-5p pair downgraded strong → weak, same PMID
  20584899, but let-7a-5p → AGO4 is still strong); **IL6R drops
  from weak (40/media) to none (25/baja)** — both its pairs are
  absent from the v10.0 raw file — so the 5 researched genes are
  now {AAK1, AGO4, CHD7, MEIS1, PURB} (PURB replaces IL6R on the
  alphabetical tiebreak). Control: 6 genes; NFIA best stays strong
  (miR-21-5p → NFIA strong → weak, PMID 24980967) with a new
  miR-375 → NFIA clip_only; **SOX6 rises from none (25/baja) to
  strong (55/media)** — miR-21-5p → SOX6 went from no Support Type
  to "Functional MTI" (PMID 32757994) — so the prioritised order is
  now [NFIA, SOX6, GMFB, KLF12, PDE3A, SKI]. No gene reaches "high"
  on clip_only evidence. Positive controls re-run on v10.0: still
  **6 of 8 recover their expected biology** (same two partial,
  miR-155 and miR-34a); pooled TargetScan sensitivity on strong
  targets 358/643 (0.557) vs 404/754 (0.536) on v9.0; miR-155-5p
  drops 173 → 101 strong targets. miR-107-3p keeps 289 validated
  pairs (was 298); the -5p arm still keeps 0 (the arm guard holds).
  The AARS pair for hsa-miR-124-3p is absent from v10.0 entirely
  (the raw CSV has 0 AARS rows).
- **Licence (step A).** The v10.0 terms could not be verified: the
  live site returns 404 and the Wayback Machine holds no readable
  copy of the 2025 download page or LICENSE (the v9.0 LICENSE,
  academic use free, "AS IS", is archived and readable; the v10.0
  article is CC BY-NC 4.0, which covers the article, not
  necessarily the data). So the v10.0 index is built on the user's
  machine from the user's own download and is NEVER committed to
  the repo and NEVER shipped in the Space sync package; the public
  Space keeps serving v9.0 by default. The alternative — the Space
  downloading the official CSV itself at runtime and building the
  index there — is left open for the user to decide.
- **Performance (step F, this VM, cold start).** With the v10.0
  index the app takes 20.7 s to first /health (was 6.2 s) and peaks
  at 1,397 MB RSS (was 477 MB); the index file is 8.3 MB gz (was
  2.3 MB). Stable for local/CLI use; the public Space is unaffected
  (it serves v9.0).
- **HGNC collision groups verified (pre-2.0 uncommitted work).**
  The 4 collision groups are the same gene by Ensembl ID (SEM1 ←
  {C7orf76, SHFM1}, WASHC2A ← {FAM21A, FAM21B}, PTPN20 ←
  {PTPN20A, PTPN20B}, SCX ← {SCXA, SCXB}); the paralogue WASHC2C
  (FAM21C) maps to its own gene, never to WASHC2A. The builder now
  records `collision_ensembl` in the map and refuses to write it if
  a collision ever spans two different genes.
- **Calibration wording.** The "uncalibrated" claim for the Enrichr
  route is now scoped to 5- and 19-gene panels only, in `paper.md`
  and in `docs/validation/calibration-enrichr-route.md`: at 100
  genes (n=30) the result is indistinguishable from calibrated and
  nothing is claimed either way.

### Fixed
- **Base rate (step C).** Of 10,000 random predicted pairs, weak
  evidence is 12.90% in v10.0 (9.41% in v9.0) and CLIP-only is
  10.91% — under the 20% threshold, so the adoption passes the
  gate; but CLIP-only without PMID is 31.8% of all weak evidence
  in v10.0 (0% in v9.0), which is why the separate level exists.
  The rise in weak evidence comes entirely from the no-PMID mass;
  weak evidence with PMID actually falls (9.41% → 8.80%).

## [1.60] - 2026-11-21

### Fixed
- **HGNC map policy: previous symbol preferred over alias; Approved
  symbols never mapped.** The audit (`scripts/audit_gene_symbols.py`)
  and the builder (`scripts/build_hgnc_symbol_map.py`) merged
  `prev_symbol` and `alias_symbol` into one candidate set. Now the
  previous symbol is preferred: a unique prev candidate wins; an
  ambiguous previous symbol is discarded, never rescued by an alias;
  aliases are used only when no previous symbol reaches the
  libraries. And a gene that is itself an Approved HGNC symbol is
  never mapped — the string denotes that gene, not another gene's
  prev/alias entry. Without that guard, 45 current genes were
  wrongly renamed (e.g. BCAM -> BCAT2, while TargetScan's BCAM is
  ENSG00000187244 = HGNC BCAM itself, verified against the Ensembl
  IDs in the TargetScan source file). The map is 498 entries (487
  in GO/Reactome/WikiPathways + 11 KEGG-only), recovering AGPAT9 →
  GPAT3, LOR → LORICRIN (renamed 2019-10-15), MPP6 → PALS2,
  NOV → CCN3 and SEPT2 → SEPTIN2 — all five verified row by row
  against the HGNC file; 3 genes stay ambiguous and unmapped
  (B3GNT1, MUM1, QARS) and 1,504 current genes the libraries lack
  stay untouched. The universe is 13,016 distinct symbols
  (13,035 − 19: 15 merges into a symbol already present + 4
  collision groups; no gene is lost — every dropped count is two
  symbols of one gene). The audit JSON now records the HGNC
  download date (2026-11-21) and the policy. Reference panels
  unchanged; the cancer specificities are 0.84% / 0.36% / 0.51%
  (Reactome N 7,523). Local calibration re-run with the same
  seeds (7, 11, 3) and the corrected map: N-1 K5 GO 8.0%, Strict
  K2 GO 18.7%, thesis-like GO 0.0% (TargetScan column; the moves
  from the 543-entry map run are at most 3 panels, within sampling
  error). Enrichr calibration re-run with the corrected map (same
  seeds 7, 11, 13; population 13,016): K5 97/75/82/57, K19
  29/18/28/14 (53% any), K100 10/0/13/13 (23% any) — the verdict
  is unchanged (the route is uncalibrated at all sizes); at
  100 genes the any-library interval (0.10–0.42) includes the
  calibrated expectation of ≈0.185. The lists were redrawn from
  the new population, so the comparison with the previous run is
  not strictly paired; the largest move (Reactome K5, −10 pp) is
  under 2 standard errors for n=100.
- **Calibration wording.** The local-calibration doc no longer says
  the small rises "come from the symbol map": the re-runs use the
  same seeds, so the comparison is paired, and the differences are
  a handful of panels (at most 7 of 300) within sampling error — no
  rise is claimed.

### Changed
- `paper.md` trimmed 1,733 → 1,624 words (whole file; 1,577 body,
  1,662 splitting hyphenated words), so the count no longer depends
  on the counting method. Detail moved to `docs/`: the live
  convergence-test parameters (`docs/validation/convergence.md`),
  the 139,399 miRNA sets and DIANA's interaction types
  (`docs/comparison.md`), and the open question on Enrichr's
  adjusted-p mechanism (`docs/validation/calibration-enrichr-route.md`
  and the maintainers' email draft). Substance and AI disclosure
  untouched.

## [1.59] - 2026-11-21

### Added
- **HGNC gene-symbol map applied at load** (`local_db/hgnc_symbol_map.json`,
  538 entries, built by `scripts/build_hgnc_symbol_map.py` from the HGNC
  complete set, CC0, downloaded 2026-11-21; the builder refuses to write
  if the counts drift from the audit). TargetScan target sets and the
  miRTarBase index keys are mapped to current HGNC symbols when data is
  loaded (`load_local_data(mapped=True)` is the default; the pickle cache
  keeps the raw content, `mapped=False` returns it). Ambiguous genes (9)
  are never mapped; evidence merges when two symbols fall on the same
  current symbol, keeping the best level and the union of PMIDs and
  experiments. Test suite: `tests/test_gene_symbols.py` (4 tests).
- Gene-symbol audit (`scripts/audit_gene_symbols.py`,
  `docs/validation/gene-symbol-audit.json`, 2026-11-21): 2,496 of
  13,035 TargetScan genes (19.1%) appear in none of the three local
  enrichment libraries (GO, Reactome, WikiPathways; KEGG excluded —
  against all four the count is 2,424). C-digit-orf names: 318.
  KIAA/FAM/LOC: 232 (the audit's reported 246 could not be
  reproduced under any definition tried: KIAA 68, FAM 164, LOC 0).
  HGNC recovers 527 of the 2,496 with a unique current symbol in
  the libraries (21.1%; 11 more only via KEGG), 9 ambiguous, 311
  mapped to symbols the libraries do not hold, 1,649 unmapped.
  AARS → AARS1 test passes.

### Changed
- **Cancer panel: 18 → 19 genes (LCOR).** C10orf12 (a target of
  miR-21-5p and miR-155-5p) and LCOR (a target of let-7a-5p) are the
  same gene under two symbols; the Strict intersection missed it
  before the map. Lipid (5 genes) and control (6 genes) panels
  unchanged. Cancer n-1 unchanged (188 genes). Convergence (cancer,
  ±25% seed 7, 10,000 perms): p = 0.0816; at ±10% (seed 11, 5,000
  perms) p = 0.0422 — it crosses 0.05 at the stricter tolerance,
  so the panel holds at the reference tolerance, not as a general
  claim (`docs/validation/convergence.md` reworded). Enrichr live
  adjusted p for the 19-gene panel: 0.0262 / 0.0473 / 0.0473 (same
  three terms); exact specificities 0.83% / 0.36% / 0.51%. TargetScan
  universe 13,035 → 12,985 distinct symbols; Reactome ∩ universe
  7,337 → 7,528 (GO 10,095; WikiPathways 4,537; KEGG 5,117).
- Calibration re-run on the mapped universe (12,985 genes; see
  `docs/validation/calibration-enrichr-route.md` and
  `docs/validation/README.md`). Local route (TargetScan background):
  N-1 K5 GO 8.3% (was 7.3%), Reactome 11.7% (unchanged), WikiPathways
  3.7% (was 3.3%); Strict K2 GO 19.3% (was 17.0%), Reactome 13.3%
  (was 12.7%), WikiPathways 7.3% (unchanged); thesis-like GO 0.0%
  (was 1.5%), Reactome 3.0% (was 3.1%), WikiPathways 0.0% (unchanged).
  Live Enrichr route, engine criterion: 5-gene lists GO 0.96 (was
  0.98), Reactome 0.79 (was 0.75), WikiPathways 0.83 (was 0.86), KEGG
  0.63 (was 0.56), any 0.99 (was 1.00); the cancer-size run is now
  19 genes (`calibration-enrichr-K19-engine.json`, replacing K18):
  GO 0.34, Reactome 0.23, WikiPathways 0.21, KEGG 0.14, any 0.61
  (K18 was 0.30/0.18/0.19/0.10/0.51); 100-gene lists GO 0.13 (was
  0.10), Reactome 0.07 (was 0.03), WikiPathways 0.07 (was 0.17), KEGG
  0.07 (was 0.23), any 0.33 (was 0.37) — every library's 95%
  interval still includes 5%, so no size threshold can be claimed.
  The small rises come from the symbol map: more genes are
  recognised by the libraries, so a random list has slightly more
  chances to hit a term.

### Fixed
- **Historical effect of the arm fix (v1.58): results from versions
  up to 1.57 are not reliable for the 29 blocked names.** Until
  1.58, hsa-miR-519e-3p returned the miR-17 family's 1,385 genes,
  and 28 more names resolved to the wrong arm or the wrong family.
  Any result published or computed with ≤1.57 for those names must
  be re-run on 1.58. The thesis panel is unaffected. Recorded in
  paper.md (Validation and use) and `docs/validation/arm-resolution.md`.
- Arm audit, second check (`scripts/audit_arm_unverified_high.py`,
  `docs/validation/arm-unverified-high.json`, 2026-11-21): in the
  three reference panels (lipid, cancer, control), NO gene reaches
  "alta" only via validated pairs whose arm could not be verified —
  and more precisely, the bare-name fallback never even triggers
  there (every panel miRNA has an exact hit in the index: 185, 224,
  58, 1,141, 103, 666, 970, 673, 512, 283 pairs). The "alta" genes
  (ABCA1; AGO4, MEIS1; NFIA) rest on exact-name evidence. The
  fail-open policy for unverifiable keys (kept, with a warning) is
  locked in by a test, not just documented.

## [1.58] - 2026-11-21

### Added
- Lipid benchmark (`docs/validation/lipid-benchmark.md`): one PMID
  per miRNA–gene pair, each checked on PubMed 20 November 2026.
  ABCA1 is the only gene with abstract-level support (three 2010
  papers demonstrate miR-33 → ABCA1 repression); SNTB2 and TSC22D2
  rest on CLIP datasets only (the pair is in the dataset, not in
  the abstract); KPNA3 and SCN1A are predicted only, with no PMID
  invented for them. The paper now points to it.
- **ABCA1, all five panel miRNAs.** The benchmark now shows that
  every panel miRNA regulates ABCA1 in the published literature:
  miR-758 (21885853, Ramirez et al. 2011, luciferase), miR-144
  (23519695, Ramirez et al. 2013; 23519696, de Aguiar Vallim et
  al. 2013), miR-106b (22119192, Kim et al. 2012, functional with
  a rescue experiment) — all four PMIDs checked on PubMed. New
  miRNA × gene table. **Recovery gap, measured:** none of the four
  PMIDs appears anywhere in the official miRTarBase v9.0 index;
  2 of 5 functional pairs at strong level (40%, exact 95%
  Clopper–Pearson interval 5%–85%), 1 of 5 only as
  weak CLIP, 2 of 5 entirely missing; 4 of 7 functional PMIDs
  (57%) missing. The paper's Validation section now says so
  (1,667 words strict).
- **miRTarBase 2025 (v10.0) checked, not adopted.** The four
  ABCA1 PMIDs are absent from the whole v10.0 file (4,004,967
  rows; `scripts/check_mirtarbase_v10.py`), so the 57% gap is not
  a version problem. The builder now reads the v10.0 CSV as well
  as the v9.0 xlsx (UTF-8 BOM, PMID with ".0", human filter);
  the rebuilt v9.0 index is identical to the live one. The v10.0
  index (1,730,805 human pairs) leaves the five lipid genes'
  evidence identical (8/8 pairs) but changes cancer (AGO4 strong
  → weak; 2 IL6R pairs gone) and control (NFIA strong → weak,
  SOX6 none → strong, 2 new weak pairs). The live index stays
  v9.0: a version change is a milestone that needs the reference
  panels, the calibration and the benchmark re-run, and the
  user's approval.
- OpenAPI export (`docs/openapi.json`, regenerated by
  `scripts/export_openapi.py`, kept in sync by
  `tests/test_openapi.py`) and three verified usage examples
  (`docs/api-examples.md`; payloads in `docs/examples/`): lipid
  analysis (curl), cancer N-1 analysis (Python), compare two
  panels (curl).
- Comparison protocol (`docs/comparison-protocol.md`): 10 fixed
  panels (3 reference + 7 random, seeds 101–107,
  `docs/validation/comparison-panels.json`, generated by
  `scripts/generate_comparison_panels.py`), 4 measurable tasks
  (false positives, strong-target sensitivity, reproducibility,
  offline), run by hand under each tool's terms of use, report
  only what the run folder reproduces. KENRYU baseline:
  **0 of 7 random panels** with a significant local term —
  descriptive only (exact 95% Clopper–Pearson interval 0–41%),
  read with the positive controls (6 of 8) and the TargetScan
  sensitivity (54%)
  (`docs/validation/comparison-kenryu-baseline.json`).

### Fixed
- **Arm resolution (-5p/-3p): a miRNA name with an arm no longer
  silently receives the other arm's TargetScan targets.** Audit
  (2026-11-21, v1.56, `scripts/audit_arm_resolution.py`): 33 mature
  human miRBase 23.1 names resolve to a family whose key carries no
  arm; in 14 hairpins both arms return exactly the same gene set,
  which is impossible in a seed predictor. The fix verifies the arm
  by SEED: the TargetScan family seed (`miR_Family_Info.txt`) must
  occur in the query's mature sequence (miRBase 23.1) starting at
  offset 0-3, which tolerates the 1-2 nt of 5' annotation drift
  between miRBase versions (miR-144-3p: offset 3; miR-802-5p:
  offset 0) while rejecting the wrong arm, whose seed is absent.
  15 arms blocked (their predicted targets return empty, with an
  explicit warning, listed as `arm_unverified` in the analyze
  response, the CLI report and the UI). The seed check also catches
  14 letter-drop collapses (the `m_no_letter` fallback merged
  distinct paralogs, e.g. miR-519e-3p into the miR-17 family):
  29 blocked in total, 303 armed names keep their targets, 0 swaps
  (`docs/validation/arm_blocked.json`,
  `docs/validation/arm-resolution.md`). The same swap lived in the
  miRTarBase lookup (hsa-miR-107-5p received the -3p arm's 298
  validated pairs via the bare-name fallback): the engine now guards
  that fallback the same way, and the live Harmonizome fetch too.
  The thesis panel does not change (all five miRNAs resolve to armed
  family keys). Bare names and mouse are not checked (documented
  limits). When a seed file is missing the verification reports
  `arm_verification: unavailable` and the old behaviour returns,
  loudly. New files: `local_db/mir_family_info.txt` (tracked;
  downloaded by the updater with TargetScan 8.0) and
  `tests/test_arm_resolution.py` (11 tests, including the full
  332-name audit as a CI regression).
- v1.57 verified live on the Space (commit ef1edc1): /health 1.57,
  commit 1d720ee, mirtarbase_mode official, 285 families; lipid
  panel 5 genes, 19 Enrichr terms in 3 groups [10, 7, 2], local 0
  terms; cancer panel 18 genes, 3 terms, exact specificity
  [0.0061, 0.0032, 0.0048], local 0 terms. The live convergence
  matches the local engine bit for bit (seed 2547548954;
  p 0.126/0.123/0.097 at ±10/25/50%).
- CI failure (runs 44803067878, 44803478428):
  `test_cli_determinism_end_to_end` expected the prioritized gene
  order, but CI has no validation data (the official miRTarBase
  index is gitignored, the GMT is not in the repo) so the order
  falls back to alphabetical. The test now expects the order for
  the data actually available, the same pattern as
  `test_cli_offline`.
- **The manifest reported a stale commit.** `load_commit_hash()`
  read only `local_db/build_info.json`, which still said 1.54 /
  e7d6eb0 while HEAD had moved on. First fix: the file only won
  when its version matched APP_VERSION. That still had a hole: a
  file with the same version and an old commit would pass — and
  parsing the remote name fails on a copy hosted elsewhere (a
  GitLab fork). Final rule: in a Space — signalled by the
  `SPACE_ID` environment variable, which Hugging Face sets at
  runtime for every Space (docs: Spaces Overview, "Built-in
  environment variables") — `build_info.json` is the only record
  of the GitHub commit. Everywhere else the git HEAD is the code
  that runs and always wins. A checkout whose origin is
  huggingface.co is also treated as a Space (it tracks the Space
  repo, so its HEAD is the Space's commit). Tests: SPACE_ID set
  → the file wins; a GitLab remote → the HEAD wins.

## [1.57] - 2026-11-20

### Added
- Convergence test (`convergence_test.py`, repo root; written by
  Claude, reproduced bit-for-bit by Muse): the observed gene count
  against 10,000 random size-matched miRNA-family panels. Lipid
  panel: 4 effective seed families, 5 genes, empirical p = 0.14;
  cancer panel: 3 families, 18 genes, p = 0.10; neither is unusual
  at 0.05. The paper now says so and reports the effective size:
  the count is not evidence by itself; the value is in the
  per-gene evidence layer. Results in `docs/validation/convergence.md`.
- Convergence test in the report (`convergence` field of
  `/api/v1/analyze`, shown in the UI): 2,000 permutations per
  tolerance, p as a range across ±10–50%, resolution shown next to
  p, seed from the sorted family names, effective size with a
  same-family warning, capped at 10 effective families. Tests in
  `tests/test_convergence.py` (7 tests).
- Execution manifest (`execution_manifest` in every `/api/v1/analyze`
  response, in the CLI `analysis.json` and in the Markdown report):
  engine version, commit, data versions, parameters, convergence
  seed and a reproduce command. Same panel, same result — verified
  end to end by API and by CLI (`tests/test_manifest.py`, 5 tests;
  suite 249).
- Comparison table with only verified claims
  (`docs/comparison.md`, sources + consultation date 20 November
  2026), including a measurable four-point comparison (false
  positives on random panels, strong-target sensitivity,
  reproducibility, offline operation) and the exact citation of the
  miRPathDB 2027 evidence tiers. The paper's State of the field now
  points to it (paper.md 1,624 words).
- "How to use and cite" section in both READMEs, with the Zenodo
  concept DOI, and GitHub issue templates (bug report, feature
  request) that ask for the version/commit from the manifest.

### Fixed
- README (EN+ES): the lipid panel no longer says "Lipid pathways
  significant (q < 0.05)" and the cancer panel no longer says "3
  significant with p-adj < 0.05". Both now say the local route
  gives 0 terms and the Enrichr terms are exploratory, matching
  the paper.
- paper.md: test count 220 → 237.
- State of the field: DIANA-miRPath now cites v4.0 (one-sided
  Fisher exact test; Wallenius as the bias option) — permutation
  testing was the v3.0 alternative, not v4.0. The strong/weak split
  is now presented as miRTarBase's own classification, inherited,
  not a contribution. miRPathDB 2027 stays: the version exists
  (Flotho et al., NAR advance article 17 Nov 2026, doi
  10.1093/nar/gkag1178, PMID 43042465) — a correction claiming
  otherwise was checked against the Oxford Academic page and
  rejected; the episode is recorded in the paper review (§11).
  miRNet's "11 databases" (v1.0) verified in the primary source.

## [1.56] - 2026-11-19

### Added
- Exact specificity: `term_specificity` now returns the exact
  hypergeometric tail (scipy), equal to the unadjusted one-sided
  Fisher p-value with the TargetScan universe as population — no
  sampling, no seed, no 4-second cost. Legacy `n_panels`/`seed`
  arguments are accepted and ignored. The three cancer terms are now
  0.61%, 0.32%, 0.48% (was 0.0/0.0/0.01 from 100 draws); the reading
  does not change.
- Redundant-term grouping: terms resting on (almost) the same panel
  genes (overlap coefficient ≥ 0.75 on the panel genes in each term)
  collapse into the most significant one. The response keeps every
  term and adds `is_group_representative`, `group_size`,
  `group_members` and `group_of`; no p-value is modified. The main
  table shows only representatives, with a "+N similar terms" note.
  On the 19 lipid-panel terms this gives 3 groups (lipid terms on
  ABCA1, neuronal terms on SCN1A, nuclear-import terms on KPNA3).
- Space sync checklist (`RELEASING.md`): rollback step, /health
  commit must match the synced GitHub tag, reference panels run
  from a script, and a CI-driven sync as the medium-term goal.

### Fixed
- Calibration intervals are now exact Clopper–Pearson intervals
  (the K100 normal-approximation intervals understated the
  uncertainty; GO's interval 0.02–0.27 also includes 0.05, not only
  Reactome's). The any-library rate is not compared to 5%: with four
  calibrated libraries the expected value is ≈18.5% (1 − 0.95⁴).

## [1.55] - 2026-11-19

This release changes what users see as the enrichment result, so it
is cut as a release (see `RELEASING.md`): the local Fisher test is
now the method shown in the interface, and the live Enrichr terms
are a separate, labelled exploratory list.

### Fixed
- Specificity column: tooltip and docs now say "random gene lists of the same
  size drawn from the TargetScan universe" (not "random panels"); a value of
  0 is shown as "<1%" (100 draws).
- Attribution: the random-panel protocol (`validate_random_panels.py`) was
  written by Claude (Sonnet 5.5), not by the author — the [1.54] entry
  said "el protocolo del autor"; the published entry is left as-is and
  the correction is recorded here and in `docs/validation/README.md`.

### Added
- Enrichment: the LOCAL Fisher test (TargetScan-universe background)
  is now the method; the live Enrichr route is a separate
  `enrichment_enrichr` list in the analyze response, rendered in its
  own UI section as "p ajustado de Enrichr (exploratorio)". New helper
  `compute_enrichment_routes` (3 tests) + 2 frontend tests.
- docs/validation/: specificity of the 19 lipid-panel terms (version, date),
  re-estimated with 10,000 draws: 12 specific (below 0.03), 7 near the
  0.05 threshold (0.03–0.07, numbers only, no binary label).
- docs/validation/calibration-enrichr-route.md: K18-engine run (100
  18-gene lists): FPR 30% GO, 18% Reactome, 19% WikiPathways, 10%
  KEGG, 51% any library. Mechanism CONFIRMED: the Enrichr adjusted
  p is a BH correction over only the terms the API returns
  (`verify_enrichr_bh.py`: 1,077 of 1,077 rows reproduced across 10
  lists and 4 libraries, tolerance 1e-4). The local replica
  (`replica_bh_overlap.py`: 48 terms, BH-over-48 = 0.0342) was the
  same computation repeated, not an independent replica — it guards
  against a transcription error and adds no evidence. Same pattern
  KENRYU had before v1.45 (H-05). paper.md now carries the K18
  numbers and the confirmed mechanism (1,580 words).
- tests/fixtures/: synthetic miRTarBase official index (same schema, the
  7 thesis-panel pairs copied verbatim) so the official-format code path
  runs in CI; the 4 API tests that skipped in CI now run against it.
- CI: `smoke` job — clean-clone install (`pip install -e .`) plus the
  real CLI offline on the thesis panel (5 genes, ABCA1 strong).
- docs/validation/calibration-enrichr-route.md: live Enrichr route
  calibrated with random gene lists (100 x 5 genes: GO 98% false
  positives — the route is not calibrated at these sizes).
- paper.md: local enrichment is now the stated method (Enrichr terms
  reported as exploratory); State of the field as a comparison table
  (miRPathDB, miRNet, miRWalk, DIANA-miRPath, miEAA); new Software
  design and Research impact statement sections; trimmed to 1,580
  words (JOSS limit 750–1,750), detail moved to `docs/`.

### Fixed
- `term_specificity` returned a 2-decimal rounding that blurred the 5%
  generic threshold (a true 0.047 read as 0.05); it now returns the raw
  fraction and display layers round.
- The miRTarBase version string was hardcoded to the real file's counts;
  it is now built from the loaded index's own metadata, so a synthetic
  fixture reports itself as synthetic.
- Commit hash: `/health` and every analysis JSON now report `"commit"`,
  the short GitHub hash read from `local_db/build_info.json` (written
  by the sync step); `"unknown"` when the file is absent or malformed.

## [1.54] — 2026-11-19

### Legacy cap, fondo TargetScan, especificidad, declaración y calibración

Revisión del autor y de Claude sobre el diff de v1.53. (Nota: el diff
leído era el anterior a las enmiendas — los tres pendientes de la
declaración — modelos de ChatGPT y de Claude, autor de audit_v150.py —
ya estaban resueltos en v1.53; lo nuevo aquí es el papel completo de
Claude, abajo.)

**Modo legacy: tope en "media" + aviso visible.** Quien clone y arranque
sin el índice oficial de miRTarBase caía en el modo legacy (snapshot de
Harmonizome), donde "legacy" contaba como evidencia fuerte y
desbloqueaba "alta" sin nivel verificado: sus resultados no coincidirían
con el artículo. Ahora la evidencia legacy suma +30 puntos como antes,
pero NO desbloquea "alta" (el nivel queda en "media" con la explicación
"legacy"); la respuesta del análisis trae `"mirtarbase_mode": "legacy"`
y un aviso visible con el paso exacto para construir el índice, y la UI
lo muestra sobre la tabla de validación. `/health` también expone
`"mirtarbase_mode"`. El puntaje no cambia; solo el nivel. El panel de
la tesis no cambia (usa el índice oficial). Método test-first: el test
de tope y el test de modo/aviso fallaron primero contra v1.53.

**Construir el índice oficial, paso exacto (DATA_LICENSES.md).**
`python build_mirtarbase_official.py hsa_MTI.xlsx local_db/mirtarbase_official.json.gz`
con el `hsa_MTI.xlsx` (miRTarBase 9.0) descargado de la página oficial;
después, `GET /health` debe mostrar `"mirtarbase_mode": "official"`.

**Fondo de enriquecimiento: universo de TargetScan.** La ruta local
pasaba la unión del panel (v1.49): fondo dependiente del panel, casi
sin potencia (0 términos para el panel de cáncer N-1, 188 genes). Ahora
pasa el universo de TargetScan — todos los genes con algún sitio
predicho (13.035 humanos) — intersecado por biblioteca con los genes
de cada GMT (`local_enrichment`). El panel de cáncer N-1 pasa de 0 a
43 términos locales (14 GO + 29 Reactome). La ruta viva de Enrichr no
cambia (su API no admite fondo personalizado). Método test-first: el
test de integración (cáncer N-1 offline) falló primero contra v1.53.

**Columna de especificidad.** Cada término de enriquecimiento trae
`"specificity"`: la fracción de 100 paneles aleatorios del mismo
tamaño en los que el término aparece por azar (1,0 = genérico; 0,0 =
específico; `null` si el término no está en los GMT locales). La ruta
viva traduce el término al español, así que la búsqueda usa `RawTerm`
(el nombre en inglés, ahora preservado en la respuesta). Los 3
términos del panel de cáncer tienen especificidad 0,0 / 0,0 / 0,01.

**Calibración con paneles aleatorios (paper.md, docs/validation/).**
El protocolo del autor (`validate_random_panels.py`, 300 paneles N-1,
300 estrictos de 2 miRNAs, 65 tipo tesis; semillas 7, 11 y 3) se
reprodujo EXACTO con v1.53 (mismas cifras, mismos controles
positivos). Con el fondo de unión de GMT: 13,3%/16,3%/5,7% (N-1) y
28,7%/20,0%/12,7% (estricto) de falsos positivos; con el universo
TargetScan: 7,3%/11,7%/3,3% y 17,0%/12,7%/7,3%; con la unión del
panel: 2,3%/4,7%/1,0% y 2,0%/4,0%/1,3%. La unión del panel es la más
conservadora pero sin potencia; el universo TargetScan es el
compromiso con fundamento (población fija e independiente del panel).
Para paneles tipo tesis (4–6 genes) no hay inflación en ningún fondo.
Los controles positivos fueron débiles y NO se presentan como
validación (miR-17 recupera "ciclo celular"/E2F/TGF con términos
genéricos; let-7 solo "MAPK"; miR-34 no devuelve términos).

**Declaración de uso de IA.** Papel completo de Claude (Sonnet 5.5):
además de la auditoría independiente, la revisión de las correcciones
y la recomendación de prioridades, escribió el parche y las pruebas
del cambio "alta" y redactó los borradores de la entrada del
CHANGELOG y de la declaración; Muse los revisó y adaptó.

**H-12 verificado cerrado.** Las descripciones escritas a mano y el
primer resultado de PubMed van en bloques separados desde v1.44
("Notas editoriales sobre las rutas (contexto del equipo, sin cita)"
vs "Primeros resultados en PubMed (búsqueda por gen + término; no se
verificó que respalden la nota anterior)") y sin evidencia clínica el
texto dice "predicción computacional, pendiente de validación
experimental" (v1.47); la etiqueta es "significancia estadística",
nunca "funcional".

**Caché: sin cambio de versión (verificado).** La caché persistente
guarda traducciones, resultados de PubMed, resultados de Enrichr y
la investigación por gen (`v5|GENE|YEARS`); ningún nivel ni puntaje
de validación se cachea, y el enriquecimiento local se calcula en
fresco cada vez. Ni el cambio de "alta" (v1.53) ni el cambio de fondo
(v1.54) alteran ningún valor cacheado, así que `enr-v2` y `v5` no se
tocan.

**Casos de estudio recalculados con v1.54 (en vivo, 2026-11-19).**
Panel de la tesis: mismos 5 genes y puntajes (ABCA1 100 alta strong;
SNTB2/TSC22D2 65 media weak; KPNA3 50 media none; SCN1A 60 media
none; media 68,0; conteos 1/2/0/2), 19 términos. Panel de cáncer:
mismos 18 genes, niveles y 3 términos Reactome con los mismos p-adj
(0,0176/0,0424/0,0424); CHD7 sigue en "media" (75, débil).

**Suite: 206 tests** (198 + 8 nuevos: tope legacy, modo/aviso legacy,
fondo TargetScan — unidad e integración —, especificidad — unidad e
integración — y RawTerm).

## [1.53] — 2026-11-19

### Bloque 2 (resto) + "alta" exige evidencia fuerte + declaración de IA

Revisión del autor y de Claude sobre los diffs de v1.51/v1.52 (2026-11-19).
El puntaje no cambia en ningún gen; solo cambia la etiqueta de nivel
cuando la evidencia es débil. El parche y las pruebas del cambio "alta"
los escribió Claude (Sonnet 5.5); Muse los revisó, los adaptó a v1.52
y los verificó en vivo. El panel de la tesis no cambia.

**"Alta" exige evidencia fuerte.** La evidencia débil (CLIP-seq,
microarray, proteómica) es apoyo de alto rendimiento: suma +15 puntos
pero NO desbloquea "alta", porque no demuestra unión directa. "Alta"
exige evidencia fuerte (reporter, western blot o qPCR) de la
interacción, o el snapshot legacy (que no trae nivel de evidencia:
limitación documentada, no verificación). Un gen con evidencia débil
y literatura completa sigue con 85 puntos, pero queda en "media" con
la explicación "solo evidencia débil; 'alta' exige evidencia fuerte".
Panel de la tesis sin cambio (ABCA1 fuerte, 100 "alta"; SNTB2 y
TSC22D2 débiles pero con 65 y 55, por debajo de 70). En el panel de
cáncer, CHD7 pasa de "alta" (75, débil) a "media". Esto revierte la
decisión de diseño de v1.48 (la débil satisfacía el candado); el test
que la codificaba se actualizó en consecuencia. Método test-first: los
3 tests nuevos y el test actualizado fallaron primero contra v1.52.

**Interruptor "exigir evidencia experimental": sin cambio de criterio.**
Sigue conservando genes con CUALQUIER apoyo experimental (fuerte o
débil, o legacy) para al menos un miRNA del panel: su trabajo es
separar "con apoyo" de "solo predicho", no aplicar el criterio
estricto — ese ya lo codifica el nivel "alta". La etiqueta de la
interfaz ahora lo dice así ("fuerte o débil", con la aclaración de que
la débil no prueba unión directa); el comentario del motor ya lo decía.

**GMT de Harmonizome: descarga al arrancar.** La frase de
DATA_LICENSES.md que decía que la plataforma funcionaba "sin
redistribuir datos CC BY-NC-SA" no era cierta: el archivo sigue en el
historial de GitHub y el Space conserva una copia. Corregido: el
archivo NO está en el árbol actual del repo, pero permanece en el
historial (commits anteriores a v1.52) y el Space aún lleva una copia
que se elimina en el próximo sync. Además, el arranque ahora intenta
descargar el GMT de Harmonizome (best-effort, con timeout, nunca frena
el arranque) cuando no hay índice oficial ni GMT en disco, así que
ninguna copia nueva se distribuye ni en el repo ni en la imagen
desplegada. Desactivable con KENRYU_STARTUP_GMT_DOWNLOAD=0 (la suite
lo desactiva; tests nuevos en tests/test_startup_download.py).

**Declaración de uso de IA (paper.md).** Reescrita según la guía de
JOSS (título exacto "AI usage disclosure", verificado en
joss.readthedocs.io): nombra herramientas y modelos (Muse, Meta, Muse
Spark 1.3; Claude Sonnet 5.5, Anthropic; ChatGPT GPT-6, OpenAI, 2026,
nivel gratuito, según lo informado por la propia herramienta), cubre código,
documentación y redacción del artículo, describe el papel de Claude
(auditoría independiente, revisión de las correcciones, recomendación
de prioridades) y de ChatGPT (plan inicial y lectura de la auditoría),
y afirma la revisión y validación humana. El autor escribió su propio
script de verificación (audit_v150.py, no subido a GitHub) y lo ejecutó
en su máquina, con comprobaciones en vivo contra el servicio desplegado.
Eliminada la afirmación "con su aprobación en cada versión" (no
documentable: una versión se publicó antes de su confirmación) y la
traducción al español del paper (ahora solo en README.es.md).
También se corrigió "Decisiones de Claude" por "Recomendaciones de
Claude, aprobadas por el autor" en la sección [1.51], y las frases de
los READMEs que describían cuándo un gen llega a "alta".

## [1.52] — 2026-11-19

### Bloque 2: XSS residual, /health con lo cargado, GMT a descarga bajo demanda

Pendientes de la re-auditoría de v1.50 (2026-11-19, verificados contra el
código antes de corregirlos) + la declaración de uso de IA para JOSS,
aprobada por el autor. Ningún resultado de análisis cambia: los puntajes
y los genes del panel de la tesis son los mismos antes y después.

**XSS residual (auditoría v1.50).** La fila de la Sección II del informe
imprimible interpolaba `g`, `d.system`, `associated_routes` y `pathText`
sin escapar (línea 912 en v1.51); ahora los cuatro pasan por
`escapeHtml()`. Las tres imágenes del informe (`venn`, `volcano`, `ppi`)
interpolaban el base64 sin validar (líneas 937/938/957); ahora usan el
nuevo `b64Img()`, con el mismo guard `/^[A-Za-z0-9+/=]+$/` que `setImg`
— solo base64 válido llega al `src`. Guardas nuevas en
`tests/test_frontend_safety.py` (escritas primero, fallaron contra v1.51).

**/health refleja lo cargado.** `/health` usaba `load_local_versions()`
(el registro estático) y declaraba el índice oficial 9.0 y el snapshot
de KEGG en un clon fresco donde ninguno de los dos archivos existe.
Ahora usa `loaded_data_versions()` (v1.51): en un clon fresco dice
"not loaded". Contrato actualizado en
`tests/test_updater.py::test_health_reports_data_sources`.

**GMT de Harmonizome a descarga bajo demanda.** El snapshot
`local_db/mirtarbase_gene_set_library.gmt.gz` (CC BY-NC-SA 4.0) ya no
está en el repo (`git rm` + `.gitignore`); el actualizador ya sabía
descargarlo y ahora es la única vía para obtenerlo. El motor sigue
usándolo como respaldo cuando el índice oficial no existe (también en
el Space, que conserva su copia: modelo Galaxy — la plataforma mantiene
los datos, el repo no los redistribuye). El loader se prueba ahora
contra un GMT sintético (`_tiny_legacy_gmt`, mismo formato) en vez del
archivo real. Notas actualizadas en `DATA_LICENSES.md`, el actualizador
y los READMEs.

**Nota TargetScan.** `data_versions.json` decía `targetscan_full.json.zip`;
el archivo es `data/targetscan_full.txt.zip` desde v1.49. Corregido.

**`@app.on_event("startup")` → lifespan.** El decorador deprecado
generaba un `DeprecationWarning` en cada corrida; el arranque (caché +
bases locales) ahora vive en un handler `lifespan`. Todos los tests
usan `with TestClient(...)`, así que el lifespan se ejecuta en la
suite. De paso, el mensaje de arranque decía "v1.46"; ahora usa
`APP_VERSION`.

**Declaración de uso de IA (paper.md).** Se añade la sección
"AI use statement" con su traducción "Declaración de uso de IA",
aprobada por el autor: Muse ayudó con código y documentación; la
suite de pruebas la desarrolló el autor (aportada en estado avanzado)
y Muse ayudó con cobertura adicional para las correcciones de la
auditoría; Claude hizo la auditoría técnica independiente; el autor
verificó los resultados con sus propios scripts (audit_v150.py) y
comprobaciones en vivo.

## [1.51] — 2026-11-19

### Bloque 1: tabla miRNA × gen, versiones de datos y conteos de evidencia

Recomendaciones de Claude sobre la revisión de v1.50 (2026-11-19), aprobadas por el autor: las ideas 1
y 2 aprobadas como primer bloque, junto con el conteo de evidencia y
dos matices de etiquetado. Los pesos y el criterio de selección NO
cambian (el método de la tesis sigue intacto; el interruptor "exigir
evidencia experimental" sigue apagado por defecto).

**Idea 1 — tabla miRNA × gen.** Cada gen del dossier lleva ahora
`evidence_by_mirna`: una fila por cada miRNA del panel con su nivel,
sus PMIDs y sus experimentos ("none" cuando ese miRNA no tiene entrada
para el gen). Matiz 1: el nivel "fuerte" de un gen significa que AL
MENOS UN miRNA del panel tiene evidencia fuerte para ese gen — la
tabla muestra cuál. La ficha del gen en la interfaz muestra la tabla
("Evidencia por miRNA — qué miRNA valida este gen"), con los PMIDs
enlazados (solo enteros) y todo escapado (H-08).

**Idea 2 — versiones de datos en el JSON del análisis.** Nuevo
`loaded_data_versions()` en el motor: un diccionario plano
{fuente: versión} de lo que está cargado EN ESTE PROCESO (no de lo que
el JSON estático declara). Va en el nivel superior de la respuesta de
`/api/v1/analyze` (`"data_versions"`), en `validation.data_sources` y
en la exportación estructurada — los tres con el mismo valor. Sin el
índice oficial, miRTarBase dice "not loaded" en vez de la versión
declarada; igual para KEGG y TargetScan. (La corrección de /health
para que use lo cargado queda para el bloque 2.)

**Conteos de evidencia en el resumen.** `validation.summary` gana
`n_strong`, `n_weak`, `n_legacy` y `n_only_predicted` (por nivel de
evidencia de cada gen). Panel de la tesis: 1 fuerte, 2 débiles,
2 solo predichos — la cifra "N de 5" ahora sale del propio informe.
La insignia de la interfaz muestra los conteos.

**Matiz 2 — la evidencia débil no se llama "validado".** La evidencia
débil (Functional MTI Weak: CLIP-seq, microarray, proteómica) es apoyo
de alto rendimiento: respalda la interacción sin demostrar unión
directa. Cambios de etiqueta (peso +15 y el candado de "alta" NO
cambian, por la decisión de diseño de v1.48): el resumen del nivel en
el dossier dice "Apoyo de alto rendimiento (...)"; la capa en la
puntuación dice "miRTarBase (apoyo de alto rendimiento, evidencia
débil)"; la columna de la interfaz dice "apoyo de alto rendimiento";
el texto del método define ambas cosas.

**Números antes/después (corrida en vivo local, panel de la tesis):**
sin cambio en genes ni en puntajes — solo se añaden campos y se
reetiqueta. Suite: 180 → 189 tests (8 nuevos en
`tests/test_evidence_table.py` + 1 guarda en
`tests/test_frontend_safety.py`; los 3 tests que requieren el índice
oficial se omiten en CI, como antes — el archivo sintético para el CI
es el bloque 3).

**Inconsistencias señaladas por Claude, corregidas aquí:** SECURITY.md
listaba 1.49 (ahora 1.51) y paper.md decía 162 tests (ahora 189).

## [1.50] — 2026-11-19

### Lo que faltaba de la auditoría de Claude (H-10, H-20, H-19, H-22)

Al releer la auditoría completa contra v1.49, casi todo estaba ya
hecho (H-01 con su recarga de DB en memoria, H-03 con su test de
regresión, H-06 con la frase corregida en el caso de estudio). Quedaba
esto:

**H-10: límite de tasa en la API pública.** Ningún endpoint lo tenía:
un solo cliente podía inundar `/api/v1/analyze` (cada llamada quema
30–90 s de peticiones vivas a NCBI/Enrichr) y agotar la cuota de esas
APIs para todos. Middleware propio, sin dependencias nuevas: ventana
deslizante por IP (60 peticiones/minuto en `/api/v1/*`), 429 con
cabecera `Retry-After: 60`. `/health` queda fuera (los health checks
de Hugging Face no deben fallar) y también los preflight OPTIONS.
Detrás del proxy de HF todos los clientes llegan con la IP del proxy,
así que la clave es `X-Forwarded-For` (primer salto). Configurable
por entorno (`KENRYU_RATE_LIMIT=0` lo desactiva,
`KENRYU_RATE_LIMIT_PER_MINUTE` cambia el presupuesto); la suite de
tests lo desactiva por defecto en `tests/conftest.py`.

**H-20: cobertura de errores de red.** El código ya devolvía
vacíos en vez de lanzar, pero no había tests que lo probaran.
`tests/test_live_fetchers.py` (11 tests, con transportes mock de
httpx y `asyncio.sleep` falso): `safe_pubmed_request` reintenta los
429 y devuelve el JSON al tercer intento, y devuelve `None` tras
agotar reintentos, ante un 500 o ante un timeout; OMIM devuelve `[]`
ante 429, ClinVar `{"count": 0}` ante timeout, ClinicalTrials `[]`
ante 500 y timeout, y `enrich_genes_direct` devuelve `[]` cuando
`addList` da 429 o timeout (y sin red cuando la lista tiene menos
de 3 genes). `tests/test_frontend_safety.py` (3 tests): el `script.js`
servido define `escapeHtml`, no contiene la interpolación cruda
`openGenePanel('${...}')` y no tiene la versión vieja `v1.38`.

**H-19 (resto):** la fila de TargetScan en `DATA_LICENSES.md` aclara
ahora que el archivo es la lista de sitios conservados ya
prefiltrada por conservación y que el código no aplica ningún
umbral adicional de `context++` ni `PCT`.

**H-22 — verificado (con matices de esta máquina).** Docker instalado
en la máquina de desarrollo (docker.io 29.1.3; daemon con
`--iptables=false --bridge=none --storage-driver=vfs`): los pulls
funcionan, pero `docker run` y `docker build` NO pueden ejecutar
contenedores aquí — el kernel de la VM deniega `setns` (runc:
«operation not permitted») y el driver overlayfs no monta sobre el
overlay de la VM. La imagen se construyó y verificó con **buildah
`--isolation=chroot`** (sin namespaces): el `Dockerfile` del repo se
construye tal cual, usando un `Dockerfile.local-test` temporal que
solo añade el CA del proxy de la VM para que pip verifique pypi.org
(ese archivo temporal NO se commitea; el `Dockerfile` del repo no
cambia). Imagen `kenryu:1.50` verificada corriendo: /health 1.50, 285
familias, y el panel de la tesis offline con los mismos 5 genes y
puntajes que la corrida local (ABCA1 55, SNTB2/TSC22D2 40,
KPNA3/SCN1A 25; media 37,0 offline). Nota: buildah avisa que
HEALTHCHECK se ignora en formato OCI; el `Dockerfile` sigue siendo
válido para Docker/BuildKit, donde sí se respeta.

**Números antes/después (corrida en vivo local, 2026-11-19).**
Panel de la tesis (mismos 5 genes, mismos puntajes): ABCA1 100 alta,
SNTB2 65 media, TSC22D2 65 media, KPNA3 50 media, SCN1A 60 media;
media 68,0. El limitador no cambia ningún resultado: solo frena
ráfagas (verificado en vivo con `KENRYU_RATE_LIMIT_PER_MINUTE=3`:
3 peticiones 200, la 4.ª 429 con `Retry-After`).

**Tests: 162 → 180.**

## [1.49] — 2026-11-19

### Ronda 3 (resto): fondo del enriquecimiento, Open Targets + dossier de cuatro niveles, H-15 / H-19 / H-20 / H-22 / H-23

**Fondo del enriquecimiento: universo predicho.** El test de Fisher
local usaba como fondo la unión de cada biblioteca GMT (~20.000 genes
que el análisis nunca probó) y diluía todos los términos. Ahora el
fondo es el universo predicho: la unión de las dianas del panel
(TargetScan local, o el respaldo de miRTarBase cuando TargetScan no
tiene el miRNA). Es la población que el análisis probó. La ruta viva
de Enrichr usa el fondo propio de Enrichr (su API no admite un fondo
personalizado). Detalle: un fondo vacío ya no recae en silencio en la
unión de la GMT (era un `background or ...` que trataba el conjunto
vacío como «sin fondo»). Consecuencia honesta: el panel de control
pasaba de 1 término offline (GO, respuesta al TGF-beta, p-adj 0,0288)
a 0 — ese término era un artefacto del fondo viejo, porque las dianas
del panel ya están enriquecidas en esa ruta. Los paneles de la tesis y
de cáncer siguen en 0 términos offline: ningún término local reúne
≥2 de sus genes (el GMT de GO 2023, por ejemplo, no contiene ABCA1).

**Open Targets + dossier de cuatro niveles.** Cada gen investigado
trae ahora sus enfermedades asociadas desde Open Targets (GraphQL
gratuito, sin clave; el símbolo se resuelve a Ensembl vía mygene.info;
top 3 por puntuación de asociación). DisGeNET sigue siendo la fuente
cuando hay clave `DISGENET_API_KEY`; si no, Open Targets; el origen
queda en `diseases_source`. La enfermedad **no entra en el puntaje**:
dice que el gen importa, no que el microARN lo regule. El dossier es
ahora de cuatro niveles independientes por gen — predicción,
experimental, enfermedad, literatura/datos clínicos — con
`evidence_tiers` y `tiers_present` (n de 4) en la respuesta, una
columna «Dossier» en la tabla de validación y un bloque en la ficha
del gen. La caché de investigación sube de v4 a v5 (el resultado gana
`diseases` y `diseases_source`).

**H-15: caché atómica, acotada y con frescura.** `save_persistent_cache`
escribe a un temporal y renombra (`os.replace`) bajo un candado de
proceso; cada sección se recorta a un límite (traducciones 2.000,
PubMed 1.000, Enrichr 200, investigación 500; sobreviven las más
recientes); la investigación multi-fuente caduca a los 90 días
(`_saved_at` en el archivo; las traducciones y los metadatos de PubMed
son hechos estables y no caducan).

**H-20: descargas robustas.** `download_file` reintenta 429 (respeta
`Retry-After`, tope 60 s) y timeouts de conexión (3 intentos, backoff
2/4 s); un 404 no se reintenta; solo se admiten URLs https. Cobertura
nueva: `transform_html_fasta` escribe con `newline="\n"` (en Windows
el `\n` se traducía a `\r\n`: 161 vs 157 bytes) y `verify_fasta` ya no
exige dos cabeceras `>` — un FASTA válido de una sola entrada fallaba.

**H-22: Dockerfile.** Imagen base fijada por digest
(`python:3.11-slim@sha256:07d48240…`); `build-essential` se purga tras
`pip install` (no queda compilador en la imagen final); `local_db`
pasa de 777 a 755. No se pudo probar la construcción aquí (no hay
Docker en esta máquina).

**H-19: datos coherentes.** El zip de TargetScan se renombra a
`data/targetscan_full.txt.zip` (contiene un .txt, no JSON); el nombre
viejo queda como respaldo en el cargador para que el Space siga
funcionando hasta su próximo sync. La nota de `data_versions.json`
aclara el GMT de respaldo: 596 filas (576 hsa-miRNAs + 20 de otras
especies).

**H-23: documentación.** `SECURITY.md` con el correo real del autor
y sin la mención a NutriKen; los README corrigen «18k+ dianas» por
13.035 genes humanos verificados, STEP 02 ya no describe miRTarBase
como «Remote API» (es el índice local oficial 9.0), el ejemplo JSON
muestra base64 crudo (la API no añade el prefijo `data:image/png`) y
el roadmap marca hechos el miRTarBase oficial y el fondo; el
`index.html` duplicado y desactualizado de la raíz se elimina (la
interfaz vive en `static/`). `paper.md`: 31 → 162 tests, 18.000 →
13.035 genes; la declaración de uso de IA se elimina por indicación
del autor.

**Números antes/después (corrida en vivo local, 2026-11-19).**
Panel de la tesis (mismos 5 genes, mismos puntajes): ABCA1 100 alta,
SNTB2 65 media, TSC22D2 65 media, KPNA3 50 media, SCN1A 60 media;
media 68,0. Dossier: ABCA1/SNTB2/TSC22D2 4/4, KPNA3/SCN1A 3/4 (les
falta el nivel experimental). Enfermedades en vivo (Open Targets):
ABCA1 Tangier 0,838; SCN1A Dravet 0,877; KPNA3 paraplejia espástica
88 0,504; SNTB2 neurodegenerativa 0,312; TSC22D2 diabetes tipo 2
0,394. Enriquecimiento en vivo: 19 términos (sin cambio).
`require_experimental`: [ABCA1, SNTB2, TSC22D2] (sin cambio). Control
offline: 1 término → 0 (artefacto del fondo viejo, ver arriba).

**Tests: 139 → 162.**

## [1.48] — 2026-11-19

### Ronda 3: miRTarBase oficial 9.0 + dossier de tres niveles + KEGG por descarga

**miRTarBase oficial (el cambio central).** El snapshot de Harmonizome
queda como respaldo; el motor carga ahora el índice oficial de
miRTarBase 9.0, construido localmente desde `hsa_MTI.xlsx` con
`build_mirtarbase_official.py` (957.040 filas → 396.758 pares miRNA–gen,
2.947 miRNAs, 2,3 MB). Cada par lleva su nivel de evidencia
(fuerte = Functional MTI; débil = Functional MTI (Weak): CLIP-seq,
microarray, proteómica), sus PMIDs y sus experimentos. El archivo
derivado es un dato de ejecución y **no se sube a git**: para
reproducirlo, descarga `hsa_MTI.xlsx` desde la página oficial de
miRTarBase y ejecuta el constructor (ver `DATA_LICENSES.md`). Los
términos de redistribución de miRTarBase **no se pudieron verificar**:
el sitio oficial devolvió 404 en todas sus páginas al comprobarlo el
2026-11-19. Un detalle del archivo oficial: 1.863 símbolos génicos
(SEPT/MARCH/DEC1) venían convertidos en fechas por Excel; el
constructor los recupera a su símbolo original.

**Dossier de tres niveles.** Cada gen del dossier lleva
`evidence_level` (fuerte / débil / solo predicho) con sus
`mirtarbase_pmids` y `mirtarbase_experiments`; la interfaz muestra una
columna «Evidencia miRTarBase» con el nivel y los PMIDs enlazados a
PubMed. Puntaje: fuerte +30 (como antes), débil +15, sin evidencia +0.
La evidencia débil SÍ es validación experimental de la interacción
(CLIP muestra unión física), así que satisface el candado de «alta»;
su menor peso hace que «alta» sea más difícil de alcanzar con ella.
Los puntajes se calculan en vivo en cada análisis; no hay caché
estancada que subir de versión (la caché de investigación v4 guarda
solo PubMed/OMIM/ClinVar/ClinicalTrials, no la pertenencia a
miRTarBase).

**Interruptor opcional «exigir evidencia experimental».** Nuevo campo
`require_experimental` en `POST /api/v1/analyze` (apagado por defecto)
y casilla en la interfaz. Cuando está activo, el consenso solo
conserva genes con evidencia en miRTarBase (fuerte o débil) para al
menos un miRNA del panel. Es un filtro sobre el consenso, no un cambio
de método: por defecto todo sigue igual.

**KEGG vuelve, por descarga (decisión del usuario).** KEGG no se
redistribuye: el actualizador descarga la biblioteca KEGG_2021_Human
de Enrichr en ejecución (`POST /api/v1/updates/kegg/install`), el
archivo es gitignored y nunca se sube al repo. El enriquecimiento
offline lo carga si existe. El modelo es el de Galaxy: el repo trae el
código, no todos los datos; quien clone descarga los datos desde sus
páginas oficiales (ver `DATA_LICENSES.md`).

**Números antes/después (corrida en vivo local, 2026-11-19).**
Panel de la tesis (mismos 5 genes): ABCA1 100 alta (fuerte, PMIDs
20466882/20466885/20566875 + 22473208); SNTB2 50 → **65 media**
(débil: PAR-CLIP 22012620/21572407, HITS-CLIP 22473208); TSC22D2
50 → **65 media** (débil: PAR-CLIP 23592263, HITS-CLIP 24906430);
KPNA3 50 media (sin evidencia); SCN1A 60 media (sin evidencia).
Media del panel 62,0 → **68,0**. Dossier de tres niveles: 1 fuerte,
2 con apoyo CLIP/débil, 2 solo predichos. Con el interruptor activo,
el panel queda en [ABCA1, SNTB2, TSC22D2]. Enriquecimiento: 19
términos, máx. p-adj 0,0349, ahora 8 WikiPathways + 8 GO + **3 KEGG**
(«Transportadoras de ABC», «Metabolismo del colesterol», «Digestión y
absorción»). Panel de cáncer (18 genes, sin cambio): los 5
investigados cambian porque el índice oficial tiene más cobertura —
AAK1/CHD7/IL6R (validados en 2 miRNAs) entran por delante de
PURB/RPS6KA3/WNK3 (1 miRNA); AGO4 y MEIS1 siguen 80 alta (fuerte).
Panel de control (6 genes, sin cambio): NFIA sube al primer puesto
por prioridad (validación fuerte en 2 miRNAs), antes el orden era
alfabético por empate a 0.

**Versiones de datos (H-19, parcial).** `data_versions.json` registra
ahora el índice oficial de miRTarBase 9.0 y la instalación de KEGG
(snapshot con fecha), y la nota de Enrichr ya no dice «No local copy
yet».

**Tests: 139/139** (123 de v1.47 + 16 nuevos: niveles de puntaje y
dossier, carga del índice oficial y recuperación de símbolos,
interruptor, KEGG y su instalador).

## [1.47] — 2026-11-19

### Ronda 2, Commit A: afirmaciones y consistencia (H-12, H-18, H-11, H-17, H-16, H-21)

### H-04 (decisión del usuario, 2026-11-19): corregir el texto
- El manuscrito describía evidencia dual obligatoria; el código no lo
  hace. Decisión: el texto describe lo que el código hace de verdad.
  `paper.md`: el consenso se construye sobre predicciones de TargetScan
  y cada gen se verifica después contra miRTarBase y bases clínicas en
  el dossier (miRTarBase también es respaldo cuando TargetScan no tiene
  el miRNA). El caso de estudio y la introducción de la síntesis dicen
  lo mismo. La evidencia dual se mostrará por niveles en una ronda
  posterior (miRTarBase oficial + dossier de tres niveles), sin cambiar
  el método de la tesis.

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
- 123/123 (4 nuevos en Commit B: CORS rechaza otro Space, CORS acepta
  el origen exacto sin credenciales, cabeceras de seguridad presentes,
  escape en innerHTML; 3 nuevos en H-04: paper.md sin evidencia dual
  obligatoria, caso de estudio con el modo strict real, síntesis sin
  intersección TargetScan+miRTarBase).

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
### Fixed
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
