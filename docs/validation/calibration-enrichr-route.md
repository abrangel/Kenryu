# Calibration of the live Enrichr route

The paper's enrichment p-values (thesis and cancer panels) come from
Enrichr's own background, but the random-panel calibration measured
the LOCAL route, not Enrichr's. This document measures the live route
directly: random gene lists of the same size, drawn from the TargetScan
universe, are submitted to Enrichr, and the fraction with at least one
significant term (adjusted p < 0.05) is the route's false-positive
rate. A calibrated method should give ~5%.

- Script: `calibrate_enrichr.py` (repo root); results:
  `docs/validation/calibration-enrichr-K5-engine.json` (and
  `-K18-engine.json`)
- Method per list: one addList POST + four enrich GETs (GO Biological
  Process 2023, Reactome Pathways 2024, WikiPathways 2019 Human, KEGG
  2021 Human), 1 s delay between requests, backoff on 429.
- Criterion: EXACTLY the engine's (`enrich_genes_direct` in
  kenryu_engine.py): the first 8 rows of each library, row[6]
  (adjusted p) < 0.05. The false-positive rate is the fraction of
  lists with at least one such term, per library and overall. The
  raw-p counts (row[2]) and the all-rows counts are recorded as
  controls.
- Date: 2026-11-19. Library versions as above; note that the
  Enrichr website's default GO library has since moved to "GO
  Biological Process 2026" — the 2023 version used here is under
  the site's Legacy tab.

## One random list, verified three ways

List (first list of the K5 run, seed 7): ARL8B, COL4A2, IFITM3,
MAN2A2, SPATA2. GO Biological Process 2023, top 8 rows, as returned
by the API and as shown on the website (manual check, 2026-11-19,
no sign-in; the website's table header says the adjusted p-value is
computed by Benjamini–Hochberg):

| Rank | p-value (row[2]) | adjusted p (row[6]) | Overlap | Term |
|---|---|---|---|---|
| 1 | 0.001749 | 0.02515 | MAN2A2 | Protein Deglycosylation (GO:0006517) |
| 2 | 0.002498 | 0.02515 | SPATA2 | Regulation Of Programmed Necrotic Cell Death (GO:0062098) |
| 3 | 0.003246 | 0.02515 | IFITM3 | Negative Regulation Of Viral Transcription (GO:0032897) |
| 4 | 0.004243 | 0.02515 | IFITM3 | Regulation Of Viral Transcription (GO:0046782) |
| 5 | 0.004243 | 0.02515 | IFITM3 | Response To Interferon-Alpha (GO:0035455) |
| 6 | 0.004492 | 0.02515 | IFITM3 | Negative Regulation Of Viral Entry Into Host Cell (GO:0046597) |
| 7 | 0.004741 | 0.02515 | MAN2A2 | N-glycan Processing (GO:0006491) |
| 8 | 0.004990 | 0.02515 | SPATA2 | Regulation Of Necroptotic Process (GO:0060544) |

The API and the website agree exactly (same 8 terms, same p-values,
same adjusted p-values, same single-gene overlaps; the website shows
48 entries for this list, the same 48 rows the API returned). The
adjusted p is therefore really < 0.05 for one-gene overlaps — this
is not a column mix-up in the script. But the numbers are
incompatible with a Benjamini–Hochberg correction over the whole
library: with BH, the adjusted p at rank 8 is at most p₈ × m / 8, so
0.02515 ≤ 0.00499 × m / 8 gives m ≥ 40 — and with m = 5,400 terms an
adjusted p of 0.025 would require a term with p / rank ≈ 4.7 × 10⁻⁶,
which single-gene overlaps cannot produce (the smallest raw p here
is 0.0017). The values are compatible with a correction over far
fewer tests, without knowing how many. The mechanism is unexplained;
the paper reports the observation with the date and the library
versions and does not claim a mechanism.

## Consequence (holds whatever the exact rates turn out to be)

The Enrichr p-values in the two case studies are reported as
exploratory, not as significant findings. The calibrated method is
the local Fisher test with the TargetScan-universe background
(false-positive rate 0–3% at thesis-panel sizes; see
`docs/validation/`), which returns no terms for the 5-gene panel —
the panel is too small for that test to have power. The lipid
panel's support is at the gene level (ABCA1 with strong experimental
evidence). The Enrichr terms are reported with the per-term
specificity column as the guard: each term is rarely hit by a random
list (12 of 19 below 0.03, 7 near the threshold reported as numbers
only), even though some term or other is almost always hit.

## 5-gene lists (thesis-like): engine criterion — DONE 2026-11-19

100 random 5-gene lists, seed 7, re-run with the engine's criterion
(first 8 rows per library, adjusted p < 0.05). False-positive rate:

| Library | FPR (engine) | FPR (all rows) | Significant terms in the top 8 (total over 100 lists) |
|---|---|---|---|
| GO Biological Process | 0.98 | 0.98 | 765 (~7.7 per list) |
| Reactome | 0.75 | 0.75 | 450 |
| WikiPathways | 0.86 | 0.86 | 364 |
| KEGG | 0.56 | 0.56 | 173 |

**No list at all (0/100) came back without a significant term, under
either criterion** (any library: 1.0). The count correction was real
— 39 significant GO terms per list becomes ~7.7 in the top 8, the
most the engine can publish — but the false-positive rate is
identical, because the significant terms sit in the top rows. Control
with the raw p-value (row[2]) instead of the adjusted one: GO 1.00,
Reactome 0.98, WikiPathways 0.89, KEGG 0.95 — barely higher, which is
consistent with the adjusted p behaving like a correction over far
fewer tests than the whole library.

## 18-gene lists (cancer-like): PENDING

100 random 18-gene lists, seed 11, re-running 2026-11-19 with the
engine's criterion. To be filled in when the run finishes; it decides
the wording of the cancer panel's 3 terms.
