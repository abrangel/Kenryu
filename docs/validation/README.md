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

## Results (2026-11-19, v1.53 data, same TargetScan base and GMT libraries)

| Scenario | GO (A / C / B) | Reactome (A / C / B) | WikiPathways (A / C / B) |
|---|---|---|---|
| N-1, 5 miRNAs, 300 panels, median 16 genes | 13.3% / 7.3% / 2.3% | 16.3% / 11.7% / 4.7% | 5.7% / 3.3% / 1.0% |
| Strict, 2 miRNAs, 300 panels, median 46 genes | 28.7% / 17.0% / 2.0% | 20.0% / 12.7% / 4.0% | 12.7% / 7.3% / 1.3% |
| Thesis-like: strict, 5 miRNAs, 4–6 genes, 65 panels | 3.1% / 1.5% / 0.0% | 3.1% / 3.1% / 3.1% | 0% / 0% / 0% |

The `validation_random_panels_*.json` files are the original runs; the
`bgB-*.json` files add the third (panel-union) column to the same runs.
Positive controls were weak (miR-17 recovers "cell cycle", E2F and TGF
but its strongest terms are generic; let-7 recovers only "MAPK"; miR-34
returns no terms) and are not presented as validation: they use
predicted, not demonstrated, targets. See `paper.md` ("Panel-level
calibration") for the interpretation.
