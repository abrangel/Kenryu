# Lipid-panel terms: specificity

Rule fixed before applying it to these numbers: a term with
specificity >= 0.05 (5%) is marked "generic" and is NOT presented as
a finding in the paper. Terms below 5% are "specific". The rule is
applied mechanically to all 19 terms below, specific and generic
alike — no term was picked or dropped after seeing the results.

Honesty note: the values were already computed (live run on the Space,
2026-11-19) before the rule was fixed. The 5% threshold was proposed
by Claude knowing only the range of the values (0.00–0.08), not which
terms exceeded it. Four terms exceed it; they are reported as-is.

- Engine version: 1.54 (commit 704bc9b)
- Date: 2026-11-19
- Route: live Enrichr (the enrichment p-values below come from
  Enrichr's own background, not from the local Fisher test)
- Panel: the 5-miRNA thesis lipid panel (hsa-miR-33a-5p, hsa-miR-33b-5p,
  hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p), Strict mode
- Source file: `api-tests/panel-tesis-v154-space-2026-11-19.json`
- Specificity: fraction of 100 random gene lists of the same size
  (5 genes), drawn from the TargetScan universe, in which the term
  appears by chance. A value of 0 is shown as "<1%".

## Result: 15 specific, 4 generic

| Source | Term | p-adj | Specificity | Verdict |
|---|---|---|---|---|
| WikiPathways | Metabolic pathway of LDL, HDL and TG, including diseases WP4522 | 0.0160 | 0.03 | specific |
| WikiPathways | SREBF and miR33 in cholesterol and lipid homeostasis WP2011 | 0.0160 | 0.02 | specific |
| WikiPathways | Nuclear Receptors in Lipid Metabolism and Toxicity WP299 | 0.0164 | 0.02 | specific |
| WikiPathways | Statin Pathway WP430 | 0.0164 | 0.04 | specific |
| WikiPathways | Rett syndrome causing genes WP4312 | 0.0166 | 0.06 | generic |
| WikiPathways | Vitamin B12 Metabolism WP1533 | 0.0166 | 0.04 | specific |
| WikiPathways | Folate Metabolism WP176 | 0.0187 | 0.06 | generic |
| WikiPathways | Selenium Micronutrient Network WP15 | 0.0213 | 0.08 | generic |
| GO | Cardiac Muscle Cell Contraction (GO:0086003) | 0.0213 | 0.00 | specific |
| GO | Detection Of Mechanical Stimulus Involved In Sensory Perception (GO:0050974) | 0.0213 | 0.00 | specific |
| GO | Detection Of Stimulus Involved In Sensory Perception Of Pain (GO:0062149) | 0.0213 | 0.01 | specific |
| GO | NLS-bearing Protein Import Into Nucleus (GO:0006607) | 0.0213 | 0.00 | specific |
| GO | Cardiac Muscle Cell Action Potential (GO:0086001) | 0.0237 | 0.00 | specific |
| GO | Cardiac Muscle Cell Action Potential Involved In Contraction (GO:0086002) | 0.0237 | 0.00 | specific |
| GO | Calcium Ion Import Across Plasma Membrane (GO:0098703) | 0.0250 | 0.02 | specific |
| KEGG | ABC transporters | 0.0290 | 0.03 | specific |
| KEGG | Cholesterol metabolism | 0.0290 | 0.06 | generic |
| KEGG | Fat digestion and absorption | 0.0290 | 0.04 | specific |
| GO | Import Into Nucleus (GO:0051170) | 0.0349 | 0.01 | specific |

## Reading notes

- In the local path, specificity is almost exactly the unadjusted
  Fisher p-value (a Monte Carlo estimate of the same tail
  probability), so it adds little there. On the live Enrichr path it
  adds value, because Enrichr uses its own background and specificity
  contrasts the result against the TargetScan universe.
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
