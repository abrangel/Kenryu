# Calibration of the live Enrichr route

The paper's enrichment p-values (thesis and cancer panels) come from
Enrichr's own background, but the random-panel calibration measured
the LOCAL route, not Enrichr's. This document measures the live route
directly: random gene lists of the same size, drawn from the TargetScan
universe, are submitted to Enrichr, and the fraction with at least one
significant term (adjusted p < 0.05) is the route's false-positive
rate. A calibrated method should give ~5%.

- Script: `calibrate_enrichr.py` (repo root); results:
  `docs/validation/calibration-enrichr-K5.json` (and `-K18.json`)
- Method per list: one addList POST + four enrich GETs (GO Biological
  Process 2023, Reactome Pathways 2024, WikiPathways 2019 Human, KEGG
  2021 Human), 1 s delay between requests, backoff on 429.
- Difference from the app's live path: the app keeps the top 8 terms
  per library; here ALL rows with adjusted p < 0.05 count. For the
  false-positive rate (at least one significant term) the difference is
  negligible — the top rows are already significant (checked on the
  first two lists).

## 5-gene lists (thesis-like): the route is not calibrated

100 random 5-gene lists, seed 7, 2026-11-19. False-positive rate:

| Library | FPR | Significant terms (total over 100 lists) |
|---|---|---|
| GO Biological Process | 0.98 | 3,891 (~39 per list) |
| Reactome | 0.75 | 1,681 |
| WikiPathways | 0.86 | 666 |
| KEGG | 0.56 | 192 |

**No list at all (0/100) came back without a significant term.**
Single-gene overlaps already reach adjusted p < 0.05: for example,
the first random list (ARL8B, COL4A2, IFITM3, MAN2A2, SPATA2) got 43
significant GO terms, the top ones each overlapping in one gene only
(e.g. "Protein Deglycosylation", one gene, raw p 0.0017, adjusted p
0.025). Enrichr's adjusted p-values do not behave like a
Benjamini–Hochberg correction over the whole library at these list
sizes.

Consequence: the 19 lipid-panel terms and their p-values cannot be
presented as significant findings. The calibrated method is the local
Fisher test with the TargetScan-universe background (false-positive
rate 0–3% at these sizes; see `docs/validation/`), which returns no
terms for the 5-gene panel — the panel is too small for that test to
have power. The Enrichr terms are reported as exploratory, with the
per-term specificity column as the guard: each term is rarely hit by a
random list (18 of 19 below 5%), even though some term or other is
almost always hit.

## 18-gene lists (cancer-like): PENDING

100 random 18-gene lists, seed 11, running 2026-11-19. To be filled
in when the run finishes; it decides the wording of the cancer
panel's 3 terms.
