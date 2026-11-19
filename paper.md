---
title: 'KENRYU: A Platform for Convergent microRNA Target Analysis and Multi-Source Evidence Integration'
tags:
  - Python
  - bioinformatics
  - microRNA
  - TargetScan
  - miRTarBase
  - functional enrichment
authors:
  - name: Cesar Abrahan Manzo Carvajal
    orcid: 0009-0006-4438-7456
    affiliation: Independent Researcher
date: 18 November 2026
bibliography: paper.bib
---

# Summary
KENRYU is an open-source platform for the convergent analysis of microRNA (miRNA) target genes. Given a panel of miRNAs, it identifies the genes regulated by all — or nearly all — of them, combining thermodynamic target predictions from TargetScanHuman v8.0 [@lewis2015conserved], served locally, with experimentally validated interactions from miRTarBase [@huang2020mirtarbase]. Three consensus modes (Strict, N-1, N-2) let the user trade sensitivity against confidence. Each convergent gene is then annotated with functional enrichment (KEGG, Reactome, WikiPathways, and Gene Ontology Biological Process, via Enrichr [@kuleshov2016enrichr]) and with parallel evidence from PubMed, OMIM, ClinVar, and ClinicalTrials.gov [@sayers2022ncbi]. Results are shown as Venn diagrams, volcano plots, and a protein–protein interactome, and exported as a paginated academic report with a Vancouver-style bibliography. KENRYU is a FastAPI service with a web interface, containerized with Docker and deployed on Hugging Face Spaces, with a local database of 285 miRNA families and more than 18,000 targets. It is research-use-only software: it generates hypotheses and organizes evidence, and it does not provide clinical or diagnostic recommendations.

# Statement of need
miRNA target analysis is fragmented in practice: prediction and experimental validation live in separate tools, enrichment in another, and literature and clinical-trial evidence in yet others. A researcher studying a miRNA panel must therefore stitch together several resources by hand to answer one question — which genes does this panel converge on, and what is known about them? Existing tools typically export tables or isolated visualizations, or sum predictions from many predictors. KENRYU instead requires dual evidence: a target must be both predicted by TargetScan and experimentally validated in miRTarBase. It then integrates enrichment, literature, and clinical-trial evidence into a single report, so the convergent genes arrive with their evidence attached instead of as a bare list.

# State of the field
Several established tools cover parts of this workflow, and KENRYU differs from each of them in one specific respect. multiMiR [@ru2014multimir] is a Bioconductor R package that integrates validated (miRecords, miRTarBase, TarBase) and predicted (DIANA-microT, ElMMo, MicroCosm, miRanda, miRDB, PicTar, PITA, TargetScan) miRNA–target databases with miRNA–disease and miRNA–drug associations; it is a programming library, with no web interface, no consensus modes, and no report generation. miRWalk 3.0 [@sticht2018mirwalk] is a web resource that scans 5′UTR, CDS, and 3′UTR with a random-forest model (TarPmiR) and integrates TargetScan, miRDB, and miRTarBase; it returns target lists from a web form, with no dual-evidence consensus, no per-gene literature and clinical-trial dossier, and no local execution. DIANA-miRPath v3 [@vlachos2015diana] performs pathway analysis (KEGG and Gene Ontology) from microT-CDS predicted and TarBase validated targets; it works at pathway level, with no convergent-gene report and no PubMed, OMIM, ClinVar, or ClinicalTrials.gov evidence attached to each gene. miRNet [@fan2016mirnet] is a web tool that integrates 11 miRNA databases with network visualization and functional analysis; it is network-centred and web-only, with no Strict / N-1 / N-2 dual-evidence consensus and no bilingual report export. Earlier tools such as funcExplorer, SMiRK, CAP-miRSeq, and miRViz provide clustering, enrichment, or visualization, but each covers only part of the workflow and none requires the dual-evidence convergence KENRYU is built on. Broad multi-omics platforms (e.g. PandaOmics, BioJupies) address wider questions and do not focus on miRNA panel convergence. KENRYU fills the gap between a target predictor and a finished, evidence-backed report for a miRNA panel: dual-evidence convergence, multi-source evidence per gene, and an exportable bilingual report, runnable locally or as a hosted service.

# Validation and use
KENRYU was validated on a five-miRNA lipid panel (hsa-miR-33a-5p, hsa-miR-33b-5p, hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p). In Strict mode it recovered the known miR-33 → ABCA1 axis and returned five convergent genes — ABCA1, KPNA3, SCN1A, SNTB2, and TSC22D2 — with significant enrichment in lipid and cholesterol pathways. This analysis is part of the author's master's thesis in bioinformatics (Universidad Internacional de La Rioja, 2026). A second, independent case study used a cancer panel never tried before (hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p): run live against the deployed service on 2026-11-11, Strict mode returned 18 convergent genes — AAK1, ACVR2A, AGO4, BNC2, BRWD3, CHD7, DCUN1D3, FOXP1, IL6R, KPNA4, MEIS1, PURB, RDX, RPS6KA3, SESN3, UBN2, WNK3, XKR4 — with three significant enrichment terms (Developmental Biology, Reactome, p-adj 0.0176; Recycling Pathway of L1, 0.0424; Signaling by TGFBR3, 0.0424), a biologically coherent result (AGO4 is a RISC component; ACVR2A sits in the TGF-beta pathway). The full case study is documented in `docs/case-study-cancer-panel.md`. Both cases are computational validations: KENRYU generates hypotheses and organizes evidence, and neither substitutes for experimental validation. The software is covered by an offline test suite (22 tests, including the exact reproduction of the five-gene panel result) and continuous integration, and it is archived on Zenodo (DOI 10.5281/zenodo.20564689).

# Acknowledgements
The author thanks Sergio Buenestado Serrano, thesis director at Universidad Internacional de La Rioja, for his guidance on the thesis work behind this software.

# AI use statement
The author wrote the software himself, with assistance from the AI tools Google Gemini and Anthropic Claude (code suggestions, debugging, and documentation drafts). All AI-assisted contributions were reviewed, tested, and verified by the author, who takes full responsibility for the software.

# References
