# Case study: cancer miRNA panel (second computational validation)

Second independent use of KENRYU beyond the lipid panel of the master's
thesis. Run live against the deployed Hugging Face Space on 2026-11-11
(`POST /api/v1/analyze`, `years: 10`, `mode: strict`). Full JSON response
is kept with the project's API test records
(`api-tests/panel-cancer-21-155-let7a.json` in the working notes).

## Panel

- **miRNAs:** hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p
- **Mode:** Strict (a gene must be predicted by TargetScan for every miRNA in the panel; miRTarBase acts as a fallback when TargetScan has no data for a miRNA, and as a post-hoc validation layer in the dossier)
- **Runtime:** 51 s (partial cache)

## Result

**18 convergent genes:** AAK1, ACVR2A, AGO4, BNC2, BRWD3, CHD7, DCUN1D3,
FOXP1, IL6R, KPNA4, MEIS1, PURB, RDX, RPS6KA3, SESN3, UBN2, WNK3, XKR4.

**Enrichment:** 24 terms (8 Reactome, 8 Gene Ontology Biological Process,
8 WikiPathways — balanced across the four sources). Three significant with
p-adj < 0.05 (Benjamini-Hochberg):

| Term | Source | p-adj |
|------|--------|-------|
| Developmental Biology | Reactome | 0.0176 |
| Recycling Pathway of L1 | Reactome | 0.0424 |
| Signaling by TGFBR3 | Reactome | 0.0424 |

## Biological coherence

The result is internally coherent, which is what a second case study must
show beyond the count of genes:

- **AGO4** is a component of the RISC complex — the very machinery
  microRNAs use — so its presence in a miRNA convergence core is expected,
  not surprising.
- **ACVR2A** is an activin receptor in the TGF-beta pathway, coherent with
  the significant "Signaling by TGFBR3" enrichment term.
- **IL6R** and **FOXP1** are established cancer and inflammation genes,
  coherent with the choice of panel (miR-21-5p and miR-155-5p are
  well-known oncomiRs; let-7a-5p is a tumour suppressor).

Only the genes investigated in the run carry PubMed (PMID) and OMIM evidence in the JSON; the remaining genes are reported without a literature dossier.

## Reproduce it

```bash
curl -X POST "https://kenryu007-bioinformatica.hf.space/api/v1/analyze" \
  -H "Content-Type: application/json" \
  -d '{"mirnas": ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"], "years": 10, "mode": "strict"}'
```

## Scope and limits (read this before citing it)

- This is a **second computational case**, not external adoption and not
  experimental validation. KENRYU generates hypotheses and organizes
  evidence; it does not substitute for in vitro or in vivo validation.
- The panel was chosen by the author (cancer biology), run once, and
  reported as run. No cherry-picking: the 18 genes and the three
  significant terms are the complete Strict-mode output.
- The original run predates the let-7 naming fix in `clean_mirna_name_cloud`
  (see CHANGELOG): the remote miRTarBase lookup for the let-7 family now
  queries the correct name (`hsa-let-7a-5p` instead of `hsa-miR-7a-5p`).
  Re-run live against the deployed Space on 2026-11-18 (52 s, after the
  v1.38 sync): the Strict core is **identical — the same 18 genes**, 24
  enrichment terms, and the response contains `hsa-let-7a-5p` with no
  occurrence of the wrong `hsa-miR-7a-5p`. The fix changed the lookup name
  without changing the validated result.
