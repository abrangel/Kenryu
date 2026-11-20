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
is 0.0017).

### A consistent mechanism: BH over only the terms with an overlap

The author reproduced the list against the local GO library
(5,407 terms; universe 14,698 genes) and found a mechanism that
fits: 48 of the 5,407 terms have at least one overlapping gene; BH
over all 5,407 tests gives an adjusted p of 1.0 for the top 8, while
BH over only the 48 terms with an overlap gives 0.034, identical
across the top 8 — the same pattern Enrichr shows (0.02515) — and
the 8 terms come out in the same order as Enrichr's, except the
8th. The raw p differs slightly (0.0024 vs Enrichr's 0.0017)
because Enrichr's gene universe is larger, and the earlier bound
m ≥ 40 is consistent with m = 48. The same replica was run
independently against the repository's local GO library
(`docs/validation/replica_bh_overlap.py`, 2026-11-19): identical
result — 48 terms, raw p 0.002379, BH over the 48 = 0.0342 in the
top 8, same order for the first 7, 8th term differs (ours:
Negative Regulation Of Viral Life Cycle; Enrichr: Regulation Of
Necroptotic Process). Note also that the API returned 48 rows for
this list and the website shows 48 entries — the same count as the
terms with an overlap.

So the adjusted p behaves like a BH correction computed only over
the terms that have at least one gene in common with the list, not
over the whole library. This is a hypothesis consistent with one
list and two independent replicas; it is NOT confirmed, and the
definitive answer belongs to the maintainers (a draft question is
in `joss-resubmission/email-enrichr-maintainers-DRAFT.md`). It is
also the same pattern KENRYU itself had before v1.45: H-05 was
exactly this — BH over only the terms with an overlap — corrected
to BH over all non-empty terms. An easy mistake to make, which we
also made and corrected; the paper says so, without accusation.

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

## 18-gene lists (cancer-like): DONE 2026-11-19

100 random 18-gene lists, seed 11, engine criterion (first 8 rows,
adjusted p < 0.05):

| Library | FPR (engine) | Significant terms in the top 8 (total over 100 lists) |
|---|---|---|
| GO Biological Process | 0.30 | 102 (~1.0 per list) |
| Reactome | 0.18 | 50 |
| WikiPathways | 0.19 | 66 |
| KEGG | 0.10 | 21 |

**51 of 100 random lists came back with at least one "significant"
term in some library** — about ten times the 5% a calibrated method
should give. Control with the raw p-value instead of the adjusted
one: GO 1.00, Reactome 1.00, WikiPathways 0.99, KEGG 0.99. The
adjustment does cut that rate, but only to 0.10–0.30 per library
(0.51 in any library) — still 2–10× the 5% a calibrated method
should give. The cancer panel's three terms (Developmental Biology
0.0176, Recycling Pathway of L1 0.0424, Signaling by TGFBR3
0.0424) therefore cannot be presented as significant findings;
they are reported as exploratory, with their specificity values
(<1%, <1%, 0.01) as the guard.
