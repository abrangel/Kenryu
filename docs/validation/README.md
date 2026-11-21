# Panel-level calibration

`validate_random_panels.py` (repo root) measures the false-positive rate
of the enrichment step: the fraction of random miRNA panels with at least
one significant term (adjusted p < 0.05). If the method were well
calibrated, the rate should be close to 5%.

The protocol and the script were written by Claude (Sonnet 5.5); Muse
reproduced the runs with v1.53 and verified the numbers matched.

Three backgrounds are compared in each run:

- **actual** — the union of each GMT library's genes (the pre-v1.49
  background; also what the live Enrichr path effectively uses).
- **targetscan** — all genes with any predicted TargetScan site,
  intersected per library with that library's genes (the v1.54
  background).
- **panel** — the union of the panel's own miRNAs' targets (the v1.49
  background).

## Reproduce

From the repo root, with the dependencies installed:

```bash
python validate_random_panels.py --panels 300 --size 5 --mode n-1 --seed 7 --out /tmp/r_n1.json
python validate_random_panels.py --panels 300 --size 2 --mode strict --seed 11 --out /tmp/r_s2.json
python validate_random_panels.py --panels 200 --size 5 --mode strict --min-genes 4 --max-genes 6 --seed 3 --out /tmp/r_thesis.json
```

## Results (2026-11-21, corrected HGNC symbol map: 498 entries; same TargetScan base and GMT libraries)

| Scenario | GO (actual / TargetScan) | Reactome (actual / TargetScan) | WikiPathways (actual / TargetScan) |
|---|---|---|---|
| N-1, 5 miRNAs, 300 panels, median 16 genes | 12.7% / 8.0% | 16.0% / 11.7% | 5.7% / 3.7% |
| Strict, 2 miRNAs, 300 panels, median 46 genes | 28.0% / 18.7% | 20.3% / 13.3% | 12.3% / 7.3% |
| Thesis-like: strict, 5 miRNAs, 4–6 genes, 67 panels | 3.0% / 0.0% | 3.0% / 3.0% | 0% / 0% |

The panel-union background (v1.49, third column in earlier versions of
this table) is no longer computed: the TargetScan universe has been the
background since v1.54. Previous values (2026-11-19, v1.53 data):
N-1 K5 13.3% / 7.3% (GO), 16.3% / 11.7% (Reactome), 5.7% / 3.3% (WP);
Strict K2 28.7% / 17.0%, 20.0% / 12.7%, 12.7% / 7.3%; thesis-like
(65 panels) 3.1% / 1.5%, 3.1% / 3.1%, 0% / 0%. The re-runs use the
same seeds (7, 11, 3), so the comparison is paired: the differences
are a handful of panels — at most 7 of 300 in the TargetScan column
(Strict K2 GO, 51 → 58), 3 of 300 in N-1 K5 GO (22 → 25), 1 panel in
thesis-like GO (1 of 65 → 0 of 66) — and lie within sampling error.
No rise is claimed. The runs were repeated 2026-11-21 with the
corrected map (498 entries, Approved-symbol guard): the numbers
moved by at most 3 panels (Strict K2 GO actual 87 → 84, TargetScan
58 → 56; N-1 K5 GO TargetScan 25 → 24; thesis-like 66 → 67 panels,
same rates), again within sampling error.
Positive controls were weak (miR-17 recovers "cell cycle", E2F and TGF
but its strongest terms are generic; let-7 recovers only "MAPK"; miR-34
returns no terms) and are not presented as validation: they use
predicted, not demonstrated, targets. See `paper.md` ("Panel-level
calibration") for the interpretation.
