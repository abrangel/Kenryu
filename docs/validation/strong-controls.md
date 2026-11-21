# Strong-evidence positive controls and TargetScan sensitivity

The earlier positive controls (miR-17, let-7, miR-34) used PREDICTED
targets and were weak. These use the STRONG-evidence target list of
each miRNA from the official miRTarBase 9.0 index (Functional MTI:
luciferase, qRT-PCR, western blot) and run the calibrated LOCAL
enrichment (TargetScan-universe background) on it.

- Engine version: 1.54 (commit 349ea47; the code touched after the
  sync did not change enrichment or validation)
- Date: 2026-11-19
- Script: `strong_controls.py` (repo root); results:
  `docs/validation/strong-controls.json`
- Pre-registration: the expected biology for each miRNA was fixed
  in the script's CANDIDATES table before the run (commit 9976b84
  contains the table with the expected terms), so the controls were
  not fitted to the output — the script prints the expected biology,
  then computes the enrichment and the verdict.

## TargetScan sensitivity (context, not a pass/fail)

Pooled over the 8 miRNAs: 404 of 754 strong targets are predicted by
TargetScan (0.536). Per miRNA it ranges from 0.35 (miR-21-5p) to 0.74
(miR-16-5p). Weak-evidence targets are recovered far less often
(0.13–0.53), as expected: weak evidence includes CLIP-seq and
microarray hits that are often non-conserved or indirect. TargetScan
only predicts conserved sites, so a sensitivity below 100% is
expected — it is a property of the predictor, not a failure of the
platform. The paper reports it as context for the consensus method:
the intersection with miRTarBase is what adds the validated targets
TargetScan misses.

## Positive controls: 6 of 8 recover their expected biology

| miRNA | strong | sens. strong | sens. weak | expected biology | verdict | top terms |
|---|---|---|---|---|---|---|
| hsa-miR-21-5p | 132 | 0.348 | 0.202 | apoptosis / PI3K-Akt / cancer (PTEN, PDCD4) | recovered | MicroRNAs in cancer (0.0000); Regulation Of Apoptotic Process (GO:0042981) (0.0000); Positive Regulation Of Apoptotic Process (GO: (0.0000) |
| hsa-miR-33a-5p | 31 | 0.581 | 0.257 | cholesterol / lipid metabolism (ABCA1) | recovered | PPARA Activates Gene Expression (0.0000); Regulation of Lipid Metabolism by PPARalpha (0.0000); Positive Regulation Of Cellular Catabolic Pro (0.0001) |
| hsa-let-7a-5p | 46 | 0.696 | 0.303 | cell cycle / development (HMGA2, LIN28) | recovered | MicroRNAs in cancer (0.0000); Cellular Senescence (0.0000); Human cytomegalovirus infection (0.0000) |
| hsa-miR-155-5p | 173 | 0.486 | 0.141 | immune response / inflammation (SOCS1) | partial | Regulation Of Transcription By RNA Polymerase (0.0000); Pathways in cancer (0.0000); Integrated Breast Cancer Pathway WP1984 (0.0000) |
| hsa-miR-34a-5p | 124 | 0.411 | 0.129 | cell cycle / p53 / apoptosis (BCL2, CDK4/6) | partial | Breast cancer (0.0000); MicroRNAs in cancer (0.0000); Pathways in cancer (0.0000) |
| hsa-miR-17-5p | 80 | 0.688 | 0.396 | cell cycle (E2F) | recovered | Negative Regulation Of Macromolecule Metaboli (0.0000); Human T-cell leukemia virus 1 infection (0.0000); Pancreatic cancer (0.0000) |
| hsa-miR-16-5p | 66 | 0.742 | 0.298 | cell cycle / apoptosis (BCL2) | recovered | MicroRNAs in cancer (0.0000); Pathways in cancer (0.0000); PI3K-Akt signaling pathway (0.0000) |
| hsa-miR-29a-3p | 102 | 0.676 | 0.528 | extracellular matrix / collagen | recovered | Extracellular Matrix Organization (0.0000); Small cell lung cancer (0.0000); PI3K-Akt signaling pathway (0.0000) |

## The two partial controls

- **miR-155-5p** (173 strong targets, 24 significant terms): the
  expected immune biology is essentially absent — no immune,
  inflammation, T-cell, B-cell or cytokine term reaches significance;
  only "Signaling by Interleukins" appears. The strong list is
  dominated by cancer genes, because that is what the validation
  literature studied.
- **miR-34a-5p** (124 strong targets, 24 significant terms): no p53,
  cell-cycle or apoptosis term reaches significance; only one G1/S
  transition term does. Same cause: the list's composition reflects
  the literature's cancer bias, not the miRNA's biology alone.

Both are reported as partial, not dropped. The lesson for the paper:
the controls validate the local enrichment route (it recovers real
biology from strong targets), but the KEGG/Reactome cancer terms
dominate because miRTarBase's validated literature is cancer-heavy.
A control built from a miRNA whose literature is not cancer-biased
would be a stronger test; none of the eight qualifies.

## Re-run on miRTarBase 10.0 (v2.0, 2026-11-21)

Same script, same pre-registered expected biology, against the
v10.0 index. Results: `docs/validation/strong-controls-v10.json`.
The verdicts use the same standard as the v9.0 run: the expected
biology must show up in the significant terms, not just one or two
stray terms among 24 cancer-dominated ones.

| miRNA | strong | sens. strong | sens. weak | verdict | top terms |
|---|---|---|---|---|---|
| hsa-miR-21-5p | 125 | 0.384 | 0.209 | recovered | MicroRNAs in cancer; Regulation Of Apoptotic Process (GO:0042981); Signaling Pathways in Glioblastoma WP2261 |
| hsa-miR-33a-5p | 25 | 0.600 | 0.261 | recovered | PPARA Activates Gene Expression; Regulation of Lipid Metabolism by PPARalpha; Metabolism of Lipids |
| hsa-let-7a-5p | 51 | 0.627 | 0.308 | recovered | MicroRNAs in cancer; Positive Regulation Of Cell Population Proliferation (GO:0008284); Regulation Of Cell Population Proliferation (GO:0042127) |
| hsa-miR-155-5p | 101 | 0.545 | 0.172 | partial | Regulation Of Transcription By RNA Polymerase II; Regulation Of DNA-templated Transcription; Positive Regulation Of DNA-templated Transcription |
| hsa-miR-34a-5p | 100 | 0.400 | 0.143 | partial | MicroRNAs in cancer; Regulation Of Cell Population Proliferation; Developmental Biology |
| hsa-miR-17-5p | 80 | 0.700 | 0.409 | recovered | Pancreatic cancer; MicroRNAs in cancer; Pathways in cancer |
| hsa-miR-16-5p | 64 | 0.750 | 0.310 | recovered | MicroRNAs in cancer; PI3K-Akt signaling pathway; Breast cancer |
| hsa-miR-29a-3p | 97 | 0.660 | 0.554 | recovered | Small cell lung cancer; Focal adhesion; AGE-RAGE signaling pathway in diabetic complications |

**Still 6 of 8 recovered; the same two partial (miR-155, miR-34a).**
miR-155 keeps exactly one immune term among its 24 significant
terms ("Signaling by Interleukins"); miR-34a keeps two
cell-cycle-related terms and no p53 or apoptosis term. Pooled
TargetScan sensitivity on strong targets: **358 of 643 (0.557)**,
vs 404 of 754 (0.536) on v9.0. The biggest single-miRNA move is
miR-155-5p: 173 → 101 strong targets (recuration in v10.0). The
paper's "6 of 8" claim stands under v10.0; its sensitivity figure
is updated to 358/643.
