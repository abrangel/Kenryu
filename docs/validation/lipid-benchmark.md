# Lipid benchmark: one PMID per miRNA–gene pair, checked on PubMed

Panel: hsa-miR-33a-5p, hsa-miR-33b-5p, hsa-miR-144-3p, hsa-miR-758-3p,
hsa-miR-106b-5p (Strict mode, 5 convergent genes).
Date checked: 20 November 2026.
Method: the pairs and their levels come from the official miRTarBase
v9.0 index as loaded by the engine (`get_validated_evidence`). Each
PMID below was opened on PubMed and its abstract read. A pair is
marked **abstract-level** when the abstract itself states the miRNA
regulates the gene, and **dataset-level** when the pair is only in the
paper's genome-wide CLIP dataset as curated by miRTarBase — the
abstract does not mention the pair. The method label (luciferase,
PAR-CLIP, HITS-CLIP) is miRTarBase's annotation; the abstract
confirms the pair, not the label.

## ABCA1 — strong (the control axis)

| miRNA | Level | PMID | What the paper shows |
| --- | --- | --- | --- |
| hsa-miR-33a-5p | strong | 20466882 | **Abstract-level.** Najafi-Shoushtari et al., *Science* 2010: "microRNAs (miR-33a/b) embedded within introns of the SREBP genes target ABCA1 ... for posttranscriptional repression." Antisense inhibition of miR-33 up-regulates ABCA1 and raises HDL. |
| hsa-miR-33a-5p | strong | 20466885 | **Abstract-level.** Rayner et al., *Science* 2010: "miR-33 inhibits the expression of ABCA1, thereby attenuating cholesterol efflux to apolipoprotein A1." Silencing miR-33 in vivo increases hepatic ABCA1 and plasma HDL. |
| hsa-miR-33a-5p | strong | 20566875 | **Abstract-level.** Marquart et al., *PNAS* 2010: "sequences in the 3' UTR of ABCA1 ... as targets for miR-33-mediated silencing." Hepatic miR-33 overexpression lowers ABCA1 mRNA, protein and plasma HDL; silencing raises them. |
| hsa-miR-33b-5p | strong | 20466882 | **Abstract-level.** Same paper as above; it names miR-33a/b together. |
| hsa-miR-106b-5p | weak | 22473208 | **Dataset-level.** Riley et al., *EMBO J* 2012: genome-wide HITS-CLIP in EBV-transformed B cells (310 human miRNAs). miR-106b is a miR-17-family miRNA; the ABCA1 pair is in the HITS-CLIP dataset, not in the abstract. |

## SNTB2 — weak, dataset-level only

| miRNA | Level | PMID | What the paper shows |
| --- | --- | --- | --- |
| hsa-miR-33a-5p | weak | 21572407 | **Dataset-level.** Kishore et al., *Nat Methods* 2011: a quantitative comparison of CLIP methods (HuR and Argonaute 2). The miR-33 → SNTB2 site is in the AGO2 PAR-CLIP dataset; the abstract does not mention SNTB2 or miR-33. |
| hsa-miR-33a-5p | weak | 22012620 | **Dataset-level.** Lipchina et al., *Genes Dev* 2011: AGO CLIP + miRNA perturbation in human ES cells, focused on the miR-302/367 cluster. The miR-33 → SNTB2 pair is in the target list; the abstract does not mention it. |
| hsa-miR-33b-5p | weak | 21572407, 22012620 | **Dataset-level.** Same two datasets as above. |
| hsa-miR-106b-5p | weak | 22473208 | **Dataset-level.** Same HITS-CLIP dataset as for ABCA1 above. |

## TSC22D2 — weak, dataset-level only

| miRNA | Level | PMID | What the paper shows |
| --- | --- | --- | --- |
| hsa-miR-144-3p | weak | 24906430 | **Dataset-level.** Pillai et al., *Breast Cancer Res Treat* 2014: HITS-CLIP in three breast cancer subtypes, focused on miR-9-5p and miR-193a/b-3p in steroid-receptor signalling. The miR-144 → TSC22D2 pair is in the HITS-CLIP dataset; the abstract does not mention it. |
| hsa-miR-758-3p | weak | 23592263 | **Dataset-level.** Whisnant et al., *mBio* 2013: PAR-CLIP of HIV-1 interaction with the cellular miRNA machinery. The miR-758 → TSC22D2 pair is in the PAR-CLIP dataset; the abstract is about HIV-1, not the pair. |

## KPNA3 and SCN1A — predicted only, no PMID

Neither gene has experimental evidence from any panel miRNA in the
official miRTarBase v9.0 index. Both are TargetScan predictions only.
The benchmark states this plainly and invents no PMID for them.

## What the benchmark establishes

- **ABCA1 is the only gene with abstract-level support**: three
  independent papers (2010) demonstrate miR-33 → ABCA1 repression
  in mouse and human cells and in vivo. The panel's control axis is
  real.
- **SNTB2 and TSC22D2 rest on CLIP datasets, not on papers about
  them.** A CLIP binding site in a supplementary table is evidence
  of binding, not of functional regulation. miRTarBase's "weak"
  level lumps both together; the benchmark separates them, and the
  report's per-gene dossier should keep that distinction visible.
- **KPNA3 and SCN1A have no experimental support at all.** They are
  convergent by prediction alone. Any claim about them is a
  hypothesis to test, not a finding.
- The honest summary of the panel is therefore: 1 gene with strong
  functional evidence, 2 genes with CLIP binding evidence only, 2
  genes predicted only. The convergence count (5 genes, p = 0.14
  against random panels) is not evidence by itself; the dossier is.
