# Data licenses — what the MIT license does and does not cover

The MIT license in `LICENSE` applies **only to the source code** in this
repository (the Python modules, the frontend, the tests and the
documentation). It does **not** cover the datasets in `data/` and
`local_db/`, which keep their own licenses and terms.

## How data is obtained (the Galaxy model)

KENRYU does **not** ship every dataset in the repository. Like Galaxy and
other production bioinformatics platforms, it ships the code plus an
updater (`kenryu_updater.py`) that downloads each dataset from its
official source on demand, for local use. Cloning the repository gives
you the engine; running the updater (or the first analysis, which warns
about missing data) fetches the reference data. To reproduce a published
result, download each dataset below from its official page, or run
`POST /api/v1/updates/{source}/install` (admin token required).

Datasets that ARE bundled in the repository are listed as such below.
Datasets that are NOT bundled are marked "download on demand" and are
never committed: their files are gitignored.

| Dataset | File(s) | Bundled? | License / terms |
|---|---|---|---|
| **TargetScan v8.0** (predictions) | `data/targetscan_full.txt.zip` | Yes | Academic use with citation (Agarwal et al. 2015; McGeary et al. 2019). No OSI/CC redistribution license is published; the file is bundled for academic use with the required citation. The zip holds `Predicted_Targets_Info.default_predictions.txt` (a TSV, not JSON — renamed from `targetscan_full.json.zip` in v1.49): the conserved sites of conserved families, already pre-filtered by conservation — no additional `context++` or `PCT` threshold is applied by the code. |
| **miRBase v22/v23.1** (miRNA names) | `local_db/mirbase_*.fa` | Yes | Public domain / CC0. |
| **miRTarBase 9.0 official** (validated targets, with evidence levels) | `local_db/mirtarbase_official.json.gz` (derived index) | No — download on demand | Built locally by `build_mirtarbase_official.py` from the user's own download of `hsa_MTI.xlsx` (miRTarBase 9.0, 957,040 human interactions). The derived index (396,758 pairs with strong/weak level + PMID + experiment) is runtime data, never committed. **Redistribution terms could not be verified: the official site (miRTarBase.cuhk.edu.cn) returned 404 on every page when checked on 2026-11-19.** The raw file is never committed; to reproduce, download `hsa_MTI.xlsx` from the official miRTarBase page and run the builder. Cite Huang et al. 2022 (Nucleic Acids Res. 50:D222–D230, doi:10.1093/nar/gkab1079). |
| **miRTarBase snapshot (Harmonizome)** | `local_db/mirtarbase_gene_set_library.gmt.gz` | No — download on demand | Harmonizome (Ma'ayan Lab) is distributed under **CC BY-NC-SA 4.0** (non-commercial, share-alike). Not bundled since v1.52: the updater downloads the GMT snapshot on demand (see the `mirtarbase` entry in `kenryu_updater.py`); used only when the official index is absent; it carries no evidence levels. |
| **KEGG** | `local_db/gmt_KEGG_2021_Human.gmt` | No — download on demand | **Not bundled.** Kanehisa Laboratories prohibits redistributing or packaging local copies of KEGG in third-party repositories or services without an agreement with Pathway Solutions. The updater downloads the Enrichr KEGG_2021_Human gene-set library on demand, for local use only; the file is gitignored and never committed. KEGG enrichment is also available through the **live Enrichr path** (`POST /api/v1/analyze` with network access). |
| **Reactome** | `local_db/gmt_Reactome_Pathways_2024.gmt` | Yes | **CC0 1.0** (via Enrichr gene-set library). |
| **WikiPathways** | `local_db/gmt_WikiPathways_2019_Human.gmt` | Yes | **CC0 1.0** (via Enrichr gene-set library). |
| **Gene Ontology** | `local_db/gmt_GO_Biological_Process_2023.gmt` | Yes | **CC BY 4.0** (via Enrichr gene-set library), with attribution. |

## What is NOT in this repository (and why)

**KEGG** was removed from the repository on 2026-11-19 (see above) and
restored on 2026-11-19 as a **download-on-demand** dataset: the updater
fetches it from Enrichr at runtime, so the platform keeps KEGG
enrichment without redistributing KEGG data.

**miRTarBase 9.0 (official)** is not in the repository either: the raw
`hsa_MTI.xlsx` and the derived index are runtime data. The repository
ships the builder script (`build_mirtarbase_official.py`) and the loader,
so anyone can rebuild the index from their own official download.

**miRTarBase snapshot (Harmonizome)** was removed from the repository
on 2026-11-19 (v1.52) and restored as a **download-on-demand** dataset:
the updater fetches the GMT snapshot from Harmonizome at runtime when
the official index is absent, so the platform keeps the legacy fallback
without redistributing CC BY-NC-SA data.

## Live data

When the service runs with network access it also queries, but does not
redistribute: PubMed/NCBI E-utilities, Enrichr (Ma'ayan Lab), STRING-DB,
UniProt, OMIM, ClinVar and ClinicalTrials.gov. Each is used under its
own terms of service; see the source table in the README.
