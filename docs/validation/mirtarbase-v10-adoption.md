# miRTarBase 10.0 adoption (v2.0 milestone)

Decision (2026-11-21): adopt miRTarBase 10.0 (Cui et al., *Nucleic
Acids Res* 2025;53(D1):D147–D156, doi 10.1093/nar/gkae1072, PMID
39578692) as a measured milestone with a major version (2.0). The
lipid panel is unchanged (8/8 pairs identical), so the risk to the
thesis is low; the gain is for other panels and for the future.

If step A or step C had failed, the adoption would have stayed on
v9.0 and been documented as a limitation instead. A documented
limitation is more defensible than an unmeasured update.

## Step A — Licence: terms UNVERIFIABLE (mitigation applies)

Checked 2026-11-21 via the Wayback CDX API
(`mirtarbase.cuhk.edu.cn`, the live site returns 404):

- No 200-status capture of the 2025 `download.php` (0 rows).
- No `LICENSE` capture under the `miRTarBase_2025` path (0 rows).
- The 2025 homepage has one 200 row (2025-03-03) that redirects to
  a 2025-05-02 capture returning 412 (blocked): not readable.
- The v9.0 LICENSE (2022 path) IS readable (capture 2024-05-30):
  data "LICNESED FOR ACADEMIC USE FREE OF CHARGE" (sic), no
  liability, "AS IS", no maintenance obligation. It does NOT mention
  redistribution or commercial use expressly.
- The v10.0 article is CC BY-NC 4.0 — that covers the article, not
  necessarily the data.

Mitigation (the criterion for continuing): the v10.0 derived index
is built on the user's machine from the user's own download of
`hsa_MTI.csv` (`validation-panels/mirtarbase_v10.json.gz`) and is
NEVER committed to the repository and NEVER included in the Space
sync package. The public Space keeps serving v9.0 by default. The
alternative — the Space downloading the official CSV itself at
runtime and building the index there — is a decision for the user,
not taken here.

## Step B — Level-change matrix: PASS

Matrix (from the derived indices, `mirtarbase_v10_change_matrix.json`):
8,190 of 320,428 shared pairs change level (2.6%). strong→strong
5,758, strong→weak 344, strong→none 1; weak→strong 337, weak→weak
306,128; none→strong 1,992, none→weak 5,516. 76,330 pairs disappear
(weak 65,635, none 8,641, **strong 2,054**); 1,410,377 appear
(weak 1,409,604, strong 769, none 4).

Verification (2026-11-21): 30 strong→weak and 30 disappeared strong
pairs, sampled (seed 11), each checked row by row against the raw
files. **30/30 downgrades and 30/30 disappearances are recuration in
the source data; 0 constructor errors.** Script:
`validation-panels/verify_mirtarbase_v10_stepB.py`; output
`mirtarbase_v10_stepB_verification.json`. (One script fix: the v9.0
xlsx column is "References (PMID)", not "Pubmed ID".)

## Step C — Base rate: PASS

10,000 random predicted (miRNA, gene) pairs, seed 7 (same sample as
the earlier investigation). Script:
`validation-panels/analyze_mirtarbase_v10_stepC.py`.

| Fraction of sampled pairs | v9.0 | v10.0 |
|---|---|---|
| strong | 1.61% (161) | 1.73% (173) |
| weak | 9.41% (941) | 12.90% (1,290) |
| weak with PMID | 9.41% (941) | 8.80% (880) |
| weak, CLIP-only | 7.48% (748) | 10.91% (1,091) |
| weak, CLIP-only, no PMID | 0.00% (0) | 4.10% (410) |
| weak, CLIP-only, with PMID | 7.48% (748) | 6.81% (681) |
| none | 0.46% (46) | 0.05% (5) |
| absent | 88.52% | 85.32% |

Criterion: if CLIP-only were very frequent (e.g. > 20% of sampled
pairs), it must be shown without points. CLIP-only is 10.91% of
sampled pairs in v10.0 — under the threshold, so step C passes.
But CLIP-only without PMID is 31.8% of ALL weak evidence in v10.0
(vs 0% in v9.0), so the separate level in step D is necessary: the
rise in weak evidence (9.41% → 12.90%) comes entirely from the
no-PMID mass (+4.10 pp); weak evidence with PMID actually falls
(9.41% → 8.80%). Strong evidence barely moves (1.61% → 1.73%).

## Step D — Levels: DONE (2026-11-21)

- Strong (with PMID): +30, unlocks "high".
- Weak with PMID: +15, never "high".
- **"clip_only"** (Functional MTI (Weak), NO PMID, all experiments
  CLIP): visible in the dossier ("Solo CLIP sin PMID … visible, sin
  puntos, no desbloquea «alta»"), +0 points, never unlocks "high".
- Bare miRNA names (no arm): marked "brazo no especificado" in the
  dossier table, whether the query fell back to the bare key or the
  query itself was bare (e.g. the control panel's hsa-miR-375).

Measured before implementing: v9.0 has 0 pairs without PMID at any
level; v10.0 has 0 strong pairs without PMID and 1,411,229 weak pairs
without PMID, ALL of them CLIP-only (non-CLIP 0, no experiments 0) —
so one new level covers the whole mass. The constructor is fail-first:
an unknown Support Type, a Functional MTI row without PMID, or a
weak row without PMID whose experiments are not all CLIP all raise
instead of being silently mapped. The per-support PMID sets are
tracked separately, so a PMID on a weak row never rescues a strong
claim.

Indices rebuilt with the new constructor: v9.0 is IDENTICAL
(396,758 pairs, mirnas dicts equal, 0 clip_only); v10.0
(`local_db/mirtarbase_v10.json.gz`, never committed, never shipped
to the Space) has 1,730,805 pairs — strong 8,856, weak 310,363, none
357, clip_only 1,411,229 — with 0 PMID/experiment diffs against the
old build and exactly 1,411,229 weak→clip_only transitions.

The engine loads the v10.0 index when its file is present and falls
back to v9.0 otherwise; `KENRYU_MIRTARBASE_VERSION=9.0` forces the
9.0 index (reproducing old results). The loaded version is fixed in
the manifest and in /health. Tests: `tests/test_mirtarbase_levels.py`
(11 tests, fail-first).

## Step E — Results: DONE (2026-11-21)

Script: `validation-panels/compare_mirtarbase_v10_stepE.py`; output
`mirtarbase_v10_stepE_panels.json`. Scores shown are the
miRTarBase-driven part only (TargetScan +25, level +30/+15/+0).

- **Lipids: 5 genes, IDENTICAL.** The thesis panel is untouched.
  The 8 benchmark pairs are identical, and the ABCA1 five-miRNA
  check is identical (all from armed keys).
- **Cancer: 19 genes.** AGO4 best stays strong (55/media) — the
  miR-155-5p pair downgraded strong→weak (same PMID 20584899) but
  let-7a-5p → AGO4 is still strong. **IL6R: weak (40/media) → none
  (25/baja)** — both its pairs are absent from the v10.0 raw file.
  LCOR best stays none (25/baja); its v9.0 "none" PMID (21572407)
  is gone in v10.0. Consequence: the 5 researched genes are now
  {AAK1, AGO4, CHD7, MEIS1, PURB} — PURB replaces IL6R on the
  alphabetical tiebreak.
- **Control: 6 genes.** NFIA best stays strong (55/media) — the
  miR-21-5p pair downgraded strong→weak (PMID 24980967) but
  miR-29a-3p → NFIA is still strong; new miR-375 → NFIA clip_only
  (bare name, marked). **SOX6: none (25/baja) → strong (55/media)**
  — miR-21-5p → SOX6 went from no Support Type to "Functional MTI"
  (PMID 32757994); new miR-375 → SOX6 clip_only (bare name,
  marked). Consequence: the prioritised order is now
  [NFIA, SOX6, GMFB, KLF12, PDE3A, SKI] (SOX6 second, GMFB third).
- **"High" gene set: no unexplained change.** No gene reaches
  "high" on clip_only evidence (scores are capped at "media" by the
  strong gate); the only level changes in the reference panels are
  the three recurations above (AGO4 pair, IL6R, SOX6, NFIA pair),
  each explained by the raw files.
- Side effect found and fixed in the tests: the AARS pair for
  hsa-miR-124-3p is absent from miRTarBase 10.0 entirely (the raw
  CSV has 0 AARS rows) — the gene is gone from the source, not from
  our index. The symbol-map test now asserts the 9.0 file
  explicitly for that pair.

## Step F — Performance: DONE (2026-11-21)

Measured on this VM (cold start, full app to first /health):

|  | v9.0 | v10.0 |
|---|---|---|
| index file (gz) | 2.3 MB | 8.3 MB |
| pairs | 396,758 | 1,730,805 |
| index parse alone | 2.2 s | 6.7 s |
| app startup to first /health | 6.2 s | 20.7 s |
| peak RSS (full app) | 477 MB | 1,397 MB |

The v10.0 index costs ~+14 s of cold start and ~+0.9 GB of RSS on
this small VM. That is stable for local/CLI use. The public Space is
unaffected by default: under the step A mitigation it keeps serving
v9.0 (no v10.0 file is shipped), so its startup and memory do not
change. If the user later chooses the runtime-download option for
the Space, these are the numbers that deployment would pay.

## Step G — Publication: PENDING

Version 2.0 with CHANGELOG and paper updated; option to use the
9.0 index and reproduce old results; data versions fixed in the
manifest.
