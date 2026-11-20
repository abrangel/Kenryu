"""Enrichment background: the TargetScan universe.

The Fisher test's background must be the set of genes that could have
been found: all genes with any predicted TargetScan site (13,035 human
genes), intersected per library with that library's genes. Using the
union of each GMT library as the background is the wrong population:
it is ~20k genes the analysis never tested, so it inflates
significance (the random-panel protocol: 13-29% false positives with
the library union vs 7-17% with the TargetScan universe). The panel's
own union (v1.49) was the opposite failure: a panel-dependent
background with almost no power (0 terms for the 188-gene N-1 cancer
panel, vs 43 with the TargetScan universe).
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
    # v1.54: the background reaches each library intersected with that
    # library's genes (the effective population). KEGG's union is
    # {A..E}, GO's is {A,F,G,H}; the 8-gene background becomes 5 and 4.
    libs = {
        "KEGG": {"Term A": {"A", "B", "C", "D", "E"}},
        "GO": {"Term B": {"A", "F", "G", "H"}},
    }
    bg = {"A", "B", "C", "D", "E", "F", "G", "H"}
    got = local_enrichment({"A", "B"}, libs, background=bg)
    expected = []
    for label, lib in libs.items():
        lib_bg = bg & set().union(*lib.values())
        for row in fisher_enrichment({"A", "B"}, lib, background=lib_bg):
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


def test_local_enrichment_uses_targetscan_universe_per_library():
    """v1.54: the enrichment background is the TargetScan universe (all
    genes with any predicted TargetScan site), intersected per library
    with that library's genes — NOT the panel's own union and NOT the
    library union. Synthetic case: 5 genes all inside an 8-gene term.
    With the 16-gene library union as background the term is significant
    (p_adj 0.026); with the 10-gene TargetScan universe it is not
    (p_adj 0.44), because the term is a large share of the small
    universe. local_enrichment must match the direct per-library
    computation with the intersected background."""
    lib = {"T1": {f"{c}" for c in "abcdefgh"},
           "T2": {f"{c}" for c in "ijklmnop"}}
    genes = {"a", "b", "c", "d", "e"}
    universe = {f"{c}" for c in "abcdefghij"}  # 10 genes
    with_union = local_enrichment(genes, {"L": lib}, background=None)
    assert [r["Term"] for r in with_union] == ["T1"]
    with_ts = local_enrichment(genes, {"L": lib}, background=universe)
    assert with_ts == []
    # And it equals the direct per-library computation with the
    # intersected background (universe ∩ library genes).
    expected = []
    for label, l in {"L": lib}.items():
        bg = universe & set().union(*l.values())
        for row in fisher_enrichment(genes, l, background=bg):
            expected.append({"Term": row["Term"], "Source": label,
                             "Pval": row["Pval"], "RawTerm": row["Term"]})
    expected.sort(key=lambda r: r["Pval"])
    assert with_ts == expected[:24]


def test_term_specificity_generic_vs_specific():
    """v1.54: a term covering the whole background appears in every
    random panel (specificity 1.0); a term with one gene in 10,000
    never does (0.0); a 100-in-1000 term with 3-of-10 observed is
    rare by chance (0.04 with the fixed seed)."""
    from kenryu_core import term_specificity
    assert term_specificity(set(range(100)), set(range(100)), 5, 5) == 1.0
    assert term_specificity({"g0"}, {f"g{i}" for i in range(10000)}, 5, 1) == 0.0
    assert term_specificity({f"t{i}" for i in range(100)},
                            {f"t{i}" for i in range(100)} | {f"f{i}" for i in range(900)},
                            10, 3) == 0.04
    assert term_specificity({"A"}, set(), 5, 1) is None
    assert term_specificity(set(), {"A", "B"}, 5, 1) is None


def test_term_specificity_not_rounded_at_10000_draws():
    """Unreleased: with 10,000 draws the estimate keeps its full
    precision — rounding to 2 decimals would blur the 5% generic
    threshold (a true 0.068 must not read as 0.07, nor 0.047 as 0.05)."""
    import random
    from kenryu_core import term_specificity
    term = {f"t{i}" for i in range(100)}
    bg = term | {f"f{i}" for i in range(900)}
    rng = random.Random(7)
    hits = sum(1 for _ in range(10000)
               if len(set(rng.sample(sorted(bg), 10)) & term) >= 3)
    expected = hits / 10000
    assert expected != round(expected, 2)  # the case discriminates
    assert term_specificity(term, bg, 10, 3, n_panels=10000, seed=7) == expected
