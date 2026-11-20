# Lipid-panel terms: specificity

Rule fixed before applying it to these numbers: the specificity
number is reported for every term, and a binary label is used only
when the value is clearly away from the 0.05 threshold — below 0.03
"specific", above 0.07 "generic". Between 0.03 and 0.07 the term is
near the threshold and only the number is reported, with no label.
The rule is applied mechanically to all 19 terms below — no term was
picked or dropped after seeing the results.

Honesty note: the 100-draw values were already computed (live run on
the Space, 2026-11-19) before the rule was fixed. The 5% threshold
was proposed by Claude knowing only the range of the values
(0.00–0.08), not which terms exceeded it.

## Exact values (2026-11-19, Unreleased)

Specificity is now computed **exactly**, not estimated by sampling:
it is the hypergeometric tail P(overlap ≥ observed) for a random gene
list of the same size (5 genes) drawn from the TargetScan universe —
which equals the one-sided Fisher p-value **without** multiple-testing
adjustment, with the TargetScan universe as the population. It is a
guard against generic terms, not a significance test, and with 19
terms the 0.05 threshold must be read with that caution.

The table below recomputes the 19 terms with the exact function
(`term_specificity`, `kenryu_core.py`; script:
`validation-panels/regenerate_lipid_specificity_exact.py`, results:
`docs/validation/specificity-exact.json`). The 10,000-draw estimates
differ in the third digit; one verdict changed: Import Into Nucleus
(0.0311 → 0.0292 exact) moves from "near threshold" to "specific".

- Engine version: 1.55 (commit 3b68011; the exact function and the
  grouping are Unreleased work on top of it)
- Date: 2026-11-19
- Route: live Enrichr (the enrichment p-values below come from
  Enrichr's own background, not from the local Fisher test)
- Panel: the 5-miRNA thesis lipid panel (hsa-miR-33a-5p, hsa-miR-33b-5p,
  hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p), Strict mode
- Source file: `api-tests/panel-tesis-v154-space-2026-11-19.json`
- Specificity: exact hypergeometric tail (see above). A value below 1%
  is shown as "<1%".

## Result (exact): 13 specific, 6 near threshold (number only), 0 generic

| Source | Term | p-adj | Specificity (100) | Specificity (10,000) | Specificity (exact) | Verdict |
|---|---|---|---|---|---|---|
| WikiPathways | Metabolic pathway of LDL, HDL and TG, including diseases WP4522 | 0.0160 | 0.03 | 0.0129 | 0.0134 | specific  |
| WikiPathways | SREBF and miR33 in cholesterol and lipid homeostasis WP2011 | 0.0160 | 0.02 | 0.0186 | 0.0178 | specific  |
| WikiPathways | Statin Pathway WP430 | 0.0164 | 0.04 | 0.0216 | 0.0233 | specific  |
| WikiPathways | Nuclear Receptors in Lipid Metabolism and Toxicity WP299 | 0.0164 | 0.02 | 0.0252 | 0.0255 | specific  |
| WikiPathways | Rett syndrome causing genes WP4312 | 0.0166 | 0.06 | 0.0465 | 0.0451 | near threshold  |
| WikiPathways | Vitamin B12 Metabolism WP1533 | 0.0166 | 0.04 | 0.0339 | 0.0321 | near threshold  |
| WikiPathways | Folate Metabolism WP176 | 0.0187 | 0.06 | 0.0491 | 0.0483 | near threshold  |
| WikiPathways | Selenium Micronutrient Network WP15 | 0.0213 | 0.08 | 0.0576 | 0.0548 | near threshold  |
| GO | Detection Of Stimulus Involved In Sensory Perception Of Pain (GO:0062149) | 0.0213 | 0.01 | 0.0021 | 0.0026 | specific  |
| GO | Detection Of Mechanical Stimulus Involved In Sensory Perception (GO:0050974) | 0.0213 | 0.00 | 0.0033 | 0.0031 | specific  |
| GO | NLS-bearing Protein Import Into Nucleus (GO:0006607) | 0.0213 | 0.00 | 0.0079 | 0.0078 | specific  |
| GO | Cardiac Muscle Cell Contraction (GO:0086003) | 0.0213 | 0.00 | 0.0079 | 0.0088 | specific  |
| GO | Cardiac Muscle Cell Action Potential (GO:0086001) | 0.0237 | 0.00 | 0.0129 | 0.0139 | specific  |
| GO | Cardiac Muscle Cell Action Potential Involved In Contraction (GO:0086002) | 0.0237 | 0.00 | 0.0113 | 0.0139 | specific  |
| GO | Calcium Ion Import Across Plasma Membrane (GO:0098703) | 0.0250 | 0.02 | 0.0177 | 0.0160 | specific  |
| KEGG | Fat digestion and absorption | 0.0290 | 0.04 | 0.0195 | 0.0228 | specific  |
| KEGG | ABC transporters | 0.0290 | 0.03 | 0.0309 | 0.0316 | near threshold  |
| KEGG | Cholesterol metabolism | 0.0290 | 0.06 | 0.0350 | 0.0365 | near threshold  |
| GO | Import Into Nucleus (GO:0051170) | 0.0349 | 0.01 | 0.0311 | 0.0292 | specific  |

## Grouping: the 19 terms rest on 3 groups of panel genes

Terms that rest on (almost) the same panel genes collapse into the
most significant one (overlap coefficient ≥ 0.75 on the panel genes
in each term; no chaining; no p-value is modified). With a 5-gene
panel, terms that share a single gene collapse — which is correct:
it is the same evidence.

| Representative | Group size | Members (same panel genes) |
|---|---|---|
| Metabolic pathway of LDL, HDL and TG (WP4522) | 10 | SREBF and miR33; Nuclear Receptors; Statin Pathway; Vitamin B12 Metabolism; Folate Metabolism; Selenium Micronutrient Network; ABC transporters; Cholesterol metabolism; Fat digestion and absorption |
| Rett syndrome causing genes (WP4312) | 7 | Cardiac Muscle Cell Contraction; Detection Of Mechanical Stimulus; Detection Of Stimulus; Cardiac Muscle Cell Action Potential; Cardiac Muscle Cell Action Potential Involved In Contraction; Calcium Ion Import Across Plasma Membrane |
| NLS-bearing Protein Import Into Nucleus (GO:0006607) | 2 | Import Into Nucleus (GO:0051170) |

The 19 terms are therefore **not 19 independent findings**: they are
3 groups of evidence — the lipid terms on ABCA1, the neuronal terms
on SCN1A (with SNTB2), and the nuclear-import terms on KPNA3. The
paper presents them as groups, not as a list of 19.

## Reading notes

- Specificity **is** the unadjusted Fisher p-value with the TargetScan
  universe as population (the local Fisher test computes the same tail
  with the same population), so on the local path it adds little. On
  the live Enrichr path it adds value, because Enrichr uses its own
  background and specificity contrasts the result against the
  TargetScan universe.
- The three KEGG terms have a specificity here because the KEGG GMT
  was downloaded at runtime on the Space. In a fresh clone without
  that download, KEGG terms report specificity None.
- An Enrichr term may differ from the local GMT if the library
  versions do not match; in that case specificity is None.
- "Detection Of Stimulus", "Cardiac Muscle Cell Contraction" and the
  other GO terms driven by SCN1A/SNTB2 are statistically specific
  (they rarely appear in random lists) but their biological reading
  is the thesis's comorbidity hypothesis (dyslipidemia / neuronal
  excitability), not a lipid finding on its own. Specificity rules
  out "this term appears in any random list"; it does not establish
  biological relevance.
