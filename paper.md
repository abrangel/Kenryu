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
Tools such as funcExplorer, SMiRK, CAP-miRSeq, and miRViz provide clustering, enrichment, or visualization, but each covers only part of the workflow and none requires the dual-evidence convergence KENRYU is built on. Broad multi-omics platforms (e.g. PandaOmics, BioJupies) address wider questions and do not focus on miRNA panel convergence. KENRYU fills the gap between a target predictor and a finished, evidence-backed report for a miRNA panel.

# Validation and use
KENRYU was validated on a five-miRNA lipid panel (hsa-miR-33a-5p, hsa-miR-33b-5p, hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p). In Strict mode it recovered the known miR-33 → ABCA1 axis and returned five convergent genes — ABCA1, KPNA3, SCN1A, SNTB2, and TSC22D2 — with significant enrichment in lipid and cholesterol pathways. This analysis is part of the author's master's thesis in bioinformatics (Universidad Internacional de La Rioja, 2026). The software is covered by an offline test suite (21 tests, including the exact reproduction of the five-gene panel result) and continuous integration, and it is archived on Zenodo (DOI 10.5281/zenodo.20564689).

# Acknowledgements
The author thanks Sergio Buenestado Serrano, thesis director at Universidad Internacional de La Rioja, for his guidance on the thesis work behind this software.

# AI use statement
The author wrote the software himself, with assistance from the AI tools Google Gemini and Anthropic Claude (code suggestions, debugging, and documentation drafts). All AI-assisted contributions were reviewed, tested, and verified by the author, who takes full responsibility for the software.

# References
