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

All five panel miRNAs regulate ABCA1 in the published literature.
The papers say "miR-758", "miR-144" and "miR-106b"; the panel forms
are the guide arms miR-758-3p, miR-144-3p and miR-106b-5p.

| miRNA | Level | PMID | What the paper shows |
| --- | --- | --- | --- |
| hsa-miR-33a-5p | strong | 20466882 | **Abstract-level.** Najafi-Shoushtari et al., *Science* 2010: "microRNAs (miR-33a/b) embedded within introns of the SREBP genes target ABCA1 ... for posttranscriptional repression." Antisense inhibition of miR-33 up-regulates ABCA1 and raises HDL. |
| hsa-miR-33a-5p | strong | 20466885 | **Abstract-level.** Rayner et al., *Science* 2010: "miR-33 inhibits the expression of ABCA1, thereby attenuating cholesterol efflux to apolipoprotein A1." Silencing miR-33 in vivo increases hepatic ABCA1 and plasma HDL. |
| hsa-miR-33a-5p | strong | 20566875 | **Abstract-level.** Marquart et al., *PNAS* 2010: "sequences in the 3' UTR of ABCA1 ... as targets for miR-33-mediated silencing." Hepatic miR-33 overexpression lowers ABCA1 mRNA, protein and plasma HDL; silencing raises them. |
| hsa-miR-33b-5p | strong | 20466882 | **Abstract-level.** Same paper as above; it names miR-33a/b together. |
| hsa-miR-758-3p | **missing from the index** | 21885853 | **Abstract-level.** Ramirez et al., *Arterioscler Thromb Vasc Biol* 2011: "miR-758 directly targets the 3'-untranslated region of Abca1 as assessed by 3'-untranslated region luciferase reporter assays." In mouse and human cells miR-758 repressed ABCA1; anti-miR-758 increased it and raised cholesterol efflux to apoA1. The pair is **not in the miRTarBase index** (the index has 58 validated genes for miR-758-3p; ABCA1 is not among them). |
| hsa-miR-144-3p | **missing from the index** | 23519695 | **Abstract-level.** Ramirez et al., *Circ Res* 2013: "Overexpression of miR-144 reduces ABCA1 expression and attenuates cholesterol efflux to ApoA1 in macrophages. Delivery of miR-144 oligonucleotides to mice attenuates ABCA1 expression in the liver, reducing HDL levels. Conversely, silencing of miR-144 in mice increases the expression of ABCA1 and plasma HDL levels." The paper also tests the Abca1 3'UTR by luciferase (miR-33 + miR-144 together: 60% decrease). The pair is **not in the miRTarBase index** (the index has 224 validated genes for miR-144-3p; ABCA1 is not among them). |
| hsa-miR-144-3p | **missing from the index** | 23519696 | **Abstract-level.** de Aguiar Vallim et al., *Circ Res* 2013: "activation of the nuclear receptor FXR in vivo increases hepatic levels of miR-144, which in turn lowers hepatic ABCA1 and plasma HDL levels. We identified 2 complementary sequences to miR-144 in the 3' UTR of ABCA1 mRNA that are necessary for miR-144-dependent regulation." Same pair, second independent paper. |
| hsa-miR-106b-5p | weak (index); **strong in the literature** | 22473208 (index) / 22119192 (literature) | **Abstract-level (22119192).** Kim et al., *Exp Neurol* 2012: "miR-106b significantly decreased ABCA1 levels and impaired cellular cholesterol efflux in neuronal cells... rescued by expression of miR-106b-resistant ABCA1." The index carries this pair only as weak HITS-CLIP (22473208, Riley et al. 2012, genome-wide dataset); the functional paper 22119192 is **not in the index at all**. |

## miRNA × gene table (the panel's evidence at a glance)

`S` = strong, abstract-level; `W` = weak, CLIP dataset only;
`M` = published functional evidence missing from the index;
`—` = no evidence (predicted only).

| Gene | miR-33a-5p | miR-33b-5p | miR-144-3p | miR-758-3p | miR-106b-5p |
| --- | --- | --- | --- | --- | --- |
| ABCA1 | S | S | M (2 papers) | M (1 paper) | W in index, S in literature |
| SNTB2 | W | W | — | — | W |
| TSC22D2 | — | — | W | W | — |
| KPNA3 | — | — | — | — | — |
| SCN1A | — | — | — | — | — |

The strong message is not "ABCA1 is regulated by miR-33". It is:
**all five panel miRNAs converge on one gene, ABCA1, and each of
them has independent published functional evidence for it.** That is
the biological story of the case.

## Recovery gap: what the official miRTarBase index misses

Checked 20 November 2026 by loading the official miRTarBase v9.0
index and searching it for the four new PMIDs (21885853, 23519695,
23519696, 22119192): **none of the four appears anywhere in the
index, for any pair.** The index lists 20 miRNAs with ABCA1 evidence;
the panel's miR-758 and miR-144 are not among them.

- Functional ABCA1 pairs in the published literature, for the panel:
  5 of 5 miRNAs (miR-33a, miR-33b, miR-758, miR-144, miR-106b).
- In the index at strong level: 2 of 5 (miR-33a, miR-33b).
  **Recovery of functional pairs at strong level: 40%.**
- In the index in any form: 3 of 5 (miR-106b only as weak CLIP).
  **2 of 5 functional pairs are entirely missing.**
- Functional PMIDs for panel miRNAs → ABCA1: 7
  (20466882, 20466885, 20566875, 21885853, 23519695, 23519696,
  22119192). In the index: 3. **4 of 7 functional PMIDs (57%)
  are missing.**
- The gap is not that the miRNAs lack data: the index has 58
  validated genes for miR-758-3p and 224 for miR-144-3p. The
  specific ABCA1 pairs are missing. The earlier "miR-758-3p with 0
  targets" gap was the same shape: the miRNA is in the index, the
  pair is not.

This is an honest limitation for the paper, with numbers: the
per-gene dossier is only as complete as the index it reads. A user
who trusts the dossier alone would miss that miR-758 and miR-144
also regulate ABCA1, and would read miR-106b → ABCA1 as weak when
the literature has a functional paper with a rescue experiment.

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

- **ABCA1 is the convergence story of the panel**: all five panel
  miRNAs regulate it in independent published work — miR-33a/b
  (2010, three papers), miR-758 (2011, luciferase), miR-144 (2013,
  two papers, LXR and FXR pathways), miR-106b (2012, functional
  with a rescue experiment). The panel does not rest on one axis;
  five independent lines of literature meet on one gene.
- **The index recovers only part of that story**: 2 of 5 functional
  pairs at strong level, 1 of 5 only as weak CLIP, 2 of 5 missing
  entirely; 4 of 7 functional PMIDs missing. The dossier is only
  as complete as the index it reads, and the benchmark quantifies
  the gap instead of hiding it.
- **SNTB2 and TSC22D2 rest on CLIP datasets, not on papers about
  them.** A CLIP binding site in a supplementary table is evidence
  of binding, not of functional regulation. miRTarBase's "weak"
  level lumps both together; the benchmark separates them, and the
  report's per-gene dossier should keep that distinction visible.
- **KPNA3 and SCN1A have no experimental support at all.** They are
  convergent by prediction alone. Any claim about them is a
  hypothesis to test, not a finding.
- The honest summary of the panel is therefore: 1 gene with strong
  functional evidence from five miRNAs (two of them invisible to
  the index), 2 genes with CLIP binding evidence only, 2 genes
  predicted only. The convergence count (5 genes, p = 0.14
  against random panels) is not evidence by itself; the dossier is.
