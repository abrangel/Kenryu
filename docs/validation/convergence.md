# Panel convergence test

Question: how many genes in common does this panel produce, and how
surprising is that number? The test compares the observed count with
10,000 random miRNA-family panels, matching each observed family by
target count (±25%), so a panel of "large" families is compared with
"large" families. It also reports the panel's effective size (distinct
seed families): miR-33a and miR-33b share a family and count once.

Script: `convergence_test.py` (repo root; written by Claude, Sonnet
5.5, reproduced bit-for-bit by Muse on 2026-11-20 — identical
empirical p-values to full precision). Run from the repo root:

```bash
python convergence_test.py --mirnas hsa-miR-33a-5p hsa-miR-144-3p hsa-miR-758-3p hsa-miR-106b-5p hsa-miR-33b-5p --mode strict --perms 10000 --seed 7
```

## Results (Strict mode, 10,000 permutations, seed 7; 2026-11-20, re-run 2026-11-21 after the v1.59 HGNC symbol map)

| Panel | Effective families | Genes observed | Null (median, p95) | Empirical p |
|---|---|---|---|---|
| Lipid (thesis) | 4 of 5 names | 5 | 2, 6 | 0.14 |
| Control (example) | 3 | 6 | 6, 13 | 0.57 |
| Cancer | 3 | 19 | 10, 21 | 0.08 |

No convergence is unusual at 0.05 at the reference tolerance. Sensitivity: varying the
size-matching tolerance between ±10% and ±50%, the cancer panel
stays between 0.04 and 0.09 and the lipid panel between 0.11 and
0.15 (cancer at ±10%, 5,000 permutations, seed 11: p = 0.042;
`docs/validation/convergence-cancer-tol010.json`). The cancer p
crosses 0.05 at the tightest tolerance, so the conclusion is stated
at the reference tolerance, not as a blanket claim.

## In the report (API)

The same test runs live in the report (`convergence` field of
`/api/v1/analyze`, `panel_convergence` in `kenryu_core.py`), with the
conditions approved 2026-11-20:

- Seed derived from the sorted family names (SHA-256 of the joined
  names), so two runs give the same number.
- 2,000 permutations per tolerance in the API; the resolution
  1/(perms+1) is shown next to the p.
- The p is reported as a RANGE across the matching tolerance
  (±10% to ±50%), not a single value.
- The panel's effective size is shown, with a warning when miRNAs
  share a family (lipid: 5 names, 4 seeds).
- Fixed label: "No mide relevancia biológica, solo si el número de
  genes en común es inusual." Panels with more than 10 effective
  families, a single family, or a miRNA without local TargetScan
  data are skipped with their reason.
- Same limit as above: the matching controls only the number of
  targets, not shared popular genes.

## Reading

The result does not invalidate the panels. It answers "is the number
of genes surprising?", not "is ABCA1 relevant?". The honest
conclusion is that the count is not evidence by itself, and that the
value is in the per-gene evidence layer. That is the argument the
paper sustains: the panel's support is at the gene level (ABCA1 with
strong experimental evidence), not in the convergence count and not
in pathway enrichment (the local route returns 0 terms).

## Limit

The size matching controls only the number of targets, not "popular"
genes with long 3′UTRs, which are shared by many miRNAs and inflate
the null for every panel. Treat the exact p with caution; the
conclusion at the reference tolerance (nothing unusual at 0.05) is
robust, but at ±10% the cancer panel is marginally below 0.05
(p = 0.042), so the claim is not tolerance-proof.
