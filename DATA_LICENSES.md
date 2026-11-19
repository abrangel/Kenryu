# Data licenses — what the MIT license does and does not cover

The MIT license in `LICENSE` applies **only to the source code** in this
repository (the Python modules, the frontend, the tests and the
documentation). It does **not** cover the datasets in `data/` and
`local_db/`, which keep their own licenses and terms:

| Dataset | File(s) | License / terms |
|---|---|---|
| **TargetScan v8.0** (predictions) | `data/targetscan_full.json.zip` | Academic use with citation (Agarwal et al. 2015; McGeary et al. 2019). No OSI/CC redistribution license is published; the file is bundled for academic use with the required citation. |
| **miRBase v22/v23.1** (miRNA names) | `local_db/mirbase_*.fa` | Public domain / CC0. |
| **miRTarBase snapshot (Harmonizome)** | `local_db/mirtarbase_gene_set_library.gmt.gz` | Harmonizome (Ma'ayan Lab) is distributed under **CC BY-NC-SA 4.0** (non-commercial, share-alike). This snapshot is bundled as downloaded; replacing it with the official miRTarBase release is planned. |
| **Reactome** | `local_db/gmt_Reactome_Pathways_2024.gmt` | **CC0 1.0** (via Enrichr gene-set library). |
| **WikiPathways** | `local_db/gmt_WikiPathways_2019_Human.gmt` | **CC0 1.0** (via Enrichr gene-set library). |
| **Gene Ontology** | `local_db/gmt_GO_Biological_Process_2023.gmt` | **CC BY 4.0** (via Enrichr gene-set library), with attribution. |

## What is NOT in this repository (and why)

**KEGG** (`gmt_KEGG_2021_Human.gmt`) was removed from the repository on
2026-11-19: Kanehisa Laboratories prohibits redistributing or packaging
local copies of KEGG in third-party repositories or services without an
agreement with Pathway Solutions. KEGG enrichment is still available
through the **live Enrichr path** (`POST /api/v1/analyze` with network
access); the offline path uses Reactome, WikiPathways and GO only.

## Live data

When the service runs with network access it also queries, but does not
redistribute: PubMed/NCBI E-utilities, Enrichr (Ma'ayan Lab), STRING-DB,
UniProt, OMIM, ClinVar and ClinicalTrials.gov. Each is used under its
own terms of service; see the source table in the README.
