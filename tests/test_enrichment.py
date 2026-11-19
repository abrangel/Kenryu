"""Enrichment background: the predicted universe.

The Fisher test's background must be the set of genes that could have
been found: the union of the panel's predicted targets (TargetScan
local, or the miRTarBase fallback when TargetScan has no data for a
miRNA). Using the union of each GMT library as the background is the
wrong population: it is ~20k genes the analysis never tested, so it
dilutes every term and makes the offline enrichment return nothing
for small panels (the thesis panel returned 0 terms offline).
"""
from kenryu_core import fisher_enrichment, local_enrichment


def test_fisher_empty_background_returns_nothing():
    # Un fondo vacío es una población vacía, no "sin fondo": nunca debe
    # recaer en silencio en la unión de la GMT.
    lib = {"Term A": {"A", "B", "C", "D", "E"},
           "Filler": {f"F{i}" for i in range(95)}}
    assert fisher_enrichment({"A", "B"}, lib, background=set()) == []
    assert fisher_enrichment({"A", "B"}, lib, background=None) != []


def test_fisher_background_excludes_genes_outside_it():
    # Z is in the gene list but NOT in the background: it cannot have
    # been found, so it must be dropped before the 2x2 table.
    lib = {"Term A": {"A", "B", "C", "D", "E"}}
    bg = {"A", "B", "C", "D", "E"} | {f"F{i}" for i in range(95)}
    with_bg = fisher_enrichment({"A", "B", "Z"}, lib, background=bg)
    without_bg = fisher_enrichment({"A", "B"}, lib, background=bg)
    assert with_bg == without_bg
    assert with_bg and with_bg[0]["Overlap"] == 2


def test_local_enrichment_passes_background_through():
    libs = {
        "KEGG": {"Term A": {"A", "B", "C", "D", "E"}},
        "GO": {"Term B": {"A", "F", "G", "H"}},
    }
    bg = {"A", "B", "C", "D", "E", "F", "G", "H"}
    got = local_enrichment({"A", "B"}, libs, background=bg)
    expected = []
    for label, lib in libs.items():
        for row in fisher_enrichment({"A", "B"}, lib, background=bg):
            expected.append({"Term": row["Term"], "Source": label,
                             "Pval": row["Pval"], "RawTerm": row["Term"]})
    expected.sort(key=lambda r: r["Pval"])
    assert got == expected[:24]


def test_predicted_background_makes_real_overlap_significant():
    # Term A is 5,000 of the 20,000 genes in the GMT union (25%): 3 of
    # 5 common genes in it is unremarkable there. But the panel's
    # predicted targets are enriched for Term A: only 10 of the 100
    # predicted genes are in it, and the same 3-of-5 overlap is a
    # large share of that small term. This is the exact failure the
    # thesis panel hit offline (0 terms).
    term_a = {f"T{i}" for i in range(5000)}
    filler = {f"F{i}" for i in range(15000)}
    libs = {"KEGG": {"Term A": term_a, "Filler": filler}}
    genes = {"T1", "T2", "T3", "F1", "F2"}
    predicted = {f"T{i}" for i in range(1, 11)} | {f"F{i}" for i in range(1, 91)}
    union_rows = local_enrichment(genes, libs)
    assert union_rows == [], "GMT-union background must NOT call this significant"
    predicted_rows = local_enrichment(genes, libs, background=predicted)
    assert predicted_rows, "predicted background must call Term A significant"
    assert predicted_rows[0]["Term"] == "Term A"
    # The overlap behind that row, checked directly (local_enrichment
    # rows carry Term/Source/Pval only).
    direct = fisher_enrichment(genes, libs["KEGG"], background=predicted)
    assert direct[0]["Term"] == "Term A" and direct[0]["Overlap"] == 3
