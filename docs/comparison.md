# Comparison with existing tools

Full version of the State of the field table in `paper.md`, with only
verified claims. Each claim cites its primary source; all sources
were consulted on **20 November 2026**. No novelty adjectives are
used: the table states what each tool does, and the last row states
what KENRYU adds as a combination.

## The table

| Tool | Version cited | Input | What it does (verified) | Evidence per gene | Runs where |
|---|---|---|---|---|---|
| multiMiR | Ru et al. 2014, *Nucleic Acids Res* 42(17):e133, doi 10.1093/nar/gku631 | miRNA list (R) | 14 databases (8 predicted, 3 validated, 3 disease/drug); user-defined cutoffs for predicted binding strength | database of origin | R package / web server |
| DIANA-miRPath | v4.0, Tastsoglou et al. 2023, *Nucleic Acids Res* 51(W1):W154–W159, doi 10.1093/nar/gkad431 | miRNA list | target union / intersection (incl. miRNAs compensating each other on the same mRNAs); one-sided Fisher exact test; experimental (TarBase v8.0, miRTarBase, microCLIP) + predicted (microT-CDS, TargetScan); Wallenius non-central hypergeometric as the bias-accounting option; context-specific (cell/tissue) testing | pathway level | web |
| miRPathDB | 2027, Flotho et al., *Nucleic Acids Res*, advance article 17 Nov 2026, doi 10.1093/nar/gkag1178 | one miRNA or one pathway | precomputed miRNA–pathway associations; 3 cumulative evidence tiers (High: strong experimental; Medium: + broader experimental; Low: + predictions); miRTarBase 9.0 + TarBase v9; TargetScan v8 + miRDB v6.0; 4 species; 5,701 miRNA entries; 28,044 accessions | tier per association | web database |
| miRNet | Fan et al. 2016, *Nucleic Acids Res* 44(W1):W135–W141, doi 10.1093/nar/gkw288; v2.0 Chang et al. 2020, *Nucleic Acids Res* 48(W1):W244–W251 | miRNA/gene list | 11 miRNA databases (v1.0); network visualisation + functional analysis; v2.0 adds transcription factors, SNPs, >5-fold larger knowledgebases, in situ functional analysis | network edges | web |
| miRWalk | 3.0, Sticht et al. 2018, *PLOS ONE* 13(10):e0206239 | miRNA/gene | TarPmiR random-forest binding-site scan (5′UTR/CDS/3′UTR) + TargetScan/miRDB/miRTarBase | target lists | web form |
| miEAA | 2023, Aparicio-Puerta et al., *Nucleic Acids Res* 51(W1):W319–W325, doi 10.1093/nar/gkad392 | miRNA list (mature/precursor, 10 species) | enrichment over 139,399 miRNA sets (19 repositories); direct + indirect miRNA annotation; API + Python package | miRNA-set membership | web / API / package |
| **KENRYU** | v1.56 | **miRNA panel** | **TargetScan convergence (Strict/N-1/N-2) + miRTarBase evidence tiers (inherited) + per-gene literature/clinical dossier + calibrated local enrichment + convergence test + bilingual report** | **strong/weak/none per gene, with PMIDs** | **local/offline or hosted** |

## What is NOT novel (and is presented as inheritance, not contribution)

- **Intersection of several miRNAs** is covered by DIANA-miRPath
  v4.0 (union/intersection modules).
- **Strong/weak evidence levels** are miRTarBase's own
  classification: strong = reporter assay, Western blot, qPCR;
  weak = microarrays, next-generation sequencing, pSILAC
  (miRPathDB 2027 states it from miRTarBase's Support Type;
  miRPathDB 1.0, Backes et al. 2017, *Nucleic Acids Res*
  45(D1):D90–D96, states the same). KENRYU inherits it.
- **Cumulative evidence tiers**: miRPathDB 2027 has three, defined
  in its full text (not only the abstract): High = miRTarBase
  Functional MTI records + DIANA-TarBase confidence-1 records whose
  assays fall outside CLIP/AGO-based and nonpair-specific
  transcriptome- or proteome-wide classes; Medium = High plus the
  remaining experimental records (CLIP/AGO-based assays,
  nonpair-specific readouts after miRNA perturbation, Functional MTI
  (Weak) and Non-Functional MTI records, and miRTarBase entries
  without a Support Type); Low = Medium plus TargetScan and miRDB
  predictions without additional score filtering. The tiers describe
  the included target evidence, not the biological certainty of an
  enriched pathway (the paper says so).
- **Programmatic access** already exists in miRPathDB 2027 (read-only
  JSON API and Model Context Protocol endpoint, plus versioned bulk
  downloads). It is not claimed as a KENRYU differentiator.
- **Empirical testing against random sets** exists elsewhere:
  permutation testing was the DIANA-miRPath v3.0 alternative
  (Vlachos et al. 2015, *Nucleic Acids Res* 43(W1):W460–W466),
  replaced in v4.0 by Wallenius non-central hypergeometric testing;
  miRPathDB 2027 compares cross-species profiles against matched
  random target sets to account for shared target-gene biases.
- **Enrichment of miRNA lists** is covered by miEAA 2023 and
  DIANA-miRPath v4.0.

## What KENRYU adds, as a combination

Existing tools cover intersection, evidence by levels and target
enrichment; KENRYU contributes a reproducible offline flow, with
documented calibration against random lists and panels, and
documentation in Spanish. The pieces that no other tool in the
table combines are: documented calibration with random panels
(`docs/validation/`, false-positive rates per library and list
size), the panel convergence test (`docs/validation/convergence.md`),
the per-gene clinical dossier (PubMed, OMIM, ClinVar,
ClinicalTrials.gov with PMIDs), offline execution with version
records, exact term specificity, and documentation in Spanish.

miEAA 2023 is useful for positioning even though it answers another
question (it enriches miRNA lists, not gene lists): its paper states
that DIANA-miRPath and miRWalk "indirectly annotate miRNAs by
deriving Gene Ontology and pathway terms from their target genes,
although this approximation is not without caveats". That is the
same concern KENRYU's calibration addresses — the bias that target
lists carry.

## Not verified

- **miRWalk functions** are not directly verified by us. The table
  row rests on the cited paper (Sticht et al. 2018) and on the
  miEAA 2023 classification of miRWalk as indirect annotation via
  target genes.
- **miRNet 2.0 database count**: no count is cited for v2.0 (its
  paper says ">5-fold increased knowledgebases" without a number).
  The "11 databases" figure is v1.0 (Fan et al. 2016).
- **miRPathDB 2027**: this version exists — Flotho et al.,
  *Nucleic Acids Res*, advance article published 17 November 2026,
  doi 10.1093/nar/gkag1178. It is the version cited in `paper.md`;
  miRPathDB 2.0 (Backes et al. 2020, *Nucleic Acids Res*
  48(D1):D93–D98) used two subsets per miRNA (all validated +
  strong only); the 2027 release uses the three cumulative tiers.
