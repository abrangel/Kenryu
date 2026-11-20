"""Unreleased (2026-11-19): exact specificity and redundant-term grouping.

Specificity = P(overlap >= observed) for a random gene list of the same size
drawn from the population: the hypergeometric survival function. It replaces
the Monte Carlo estimate (100 draws, 1% resolution) with the exact value.
Grouping collapses enrichment terms that rest on (almost) the same panel
genes, without touching any p-value.
"""
import pytest
from scipy.stats import hypergeom

from kenryu_core import term_specificity, group_redundant_terms


def _exact(N, K, n, k):
    return float(hypergeom.sf(k - 1, N, K, n))


# ── exact specificity ──────────────────────────────────────────────────────
@pytest.mark.parametrize("N,K,n,k", [(1000, 100, 10, 3), (7337, 977, 18, 7),
                                      (7337, 35, 18, 2), (10000, 1, 5, 1)])
def test_specificity_is_the_exact_hypergeometric_tail(N, K, n, k):
    term = {f"t{i}" for i in range(K)}
    bg = term | {f"f{i}" for i in range(N - K)}
    assert term_specificity(term, bg, n, k) == pytest.approx(_exact(N, K, n, k), abs=1e-12)


def test_specificity_edges():
    bg = {f"g{i}" for i in range(50)}
    term = {f"g{i}" for i in range(10)}
    assert term_specificity(term, bg, 5, 0) == 1.0          # overlap 0: always reached
    assert term_specificity(term, bg, 5, 11) == 0.0         # impossible overlap
    assert term_specificity(term, bg, 500, 10) == 1.0       # list > population: clipped
    assert term_specificity(set(), bg, 5, 1) is None        # no term genes in population
    assert term_specificity({"zzz"}, bg, 5, 1) is None      # term outside population
    assert term_specificity(term, set(), 5, 1) is None      # no population
    assert term_specificity(term, bg, 0, 1) is None         # empty list


def test_specificity_is_deterministic_and_ignores_legacy_draw_arguments():
    bg = {f"g{i}" for i in range(200)}
    term = {f"g{i}" for i in range(30)}
    base = term_specificity(term, bg, 10, 3)
    assert term_specificity(term, bg, 10, 3, n_panels=10, seed=1) == base
    assert term_specificity(term, bg, 10, 3, n_panels=10000, seed=99) == base


def test_specificity_is_not_zero_for_rare_terms():
    # The exact value is tiny but positive: Monte Carlo would have said 0.
    bg = {f"g{i}" for i in range(10000)}
    v = term_specificity({"g0"}, bg, 5, 1)
    assert v == pytest.approx(0.0005, abs=1e-9) and v > 0


# ── grouping ───────────────────────────────────────────────────────────────
def _item(term, p, source="GO"):
    return {"Term": term, "RawTerm": term, "Source": source, "Pval": p}


def _hits(**kw):
    return {("GO", k): set(v) for k, v in kw.items()}


def test_identical_hit_sets_collapse_into_the_most_significant_term():
    items = [_item("B", 0.02), _item("A", 0.01), _item("C", 0.03)]
    hits = _hits(A="xy", B="xy", C="xy")
    groups = group_redundant_terms(items, hits)
    assert len(groups) == 1
    assert groups[0]["representative"]["Term"] == "A"
    assert [m["Term"] for m in groups[0]["members"]] == ["B", "C"]


def test_overlap_coefficient_threshold():
    items = [_item("A", 0.01), _item("B", 0.02), _item("C", 0.03)]
    # B subset of A -> coefficient 1.0 (groups); C shares 1 of 2 -> 0.5 (stays)
    hits = _hits(A="abcd", B="abc", C="ax")
    groups = group_redundant_terms(items, hits, threshold=0.75)
    assert [g["representative"]["Term"] for g in groups] == ["A", "C"]
    assert [m["Term"] for m in groups[0]["members"]] == ["B"]


def test_no_chaining_members_must_match_the_representative():
    # B~A and C~B but C is not close to A: C must not be dragged in via B.
    items = [_item("A", 0.01), _item("B", 0.02), _item("C", 0.03)]
    hits = _hits(A="abcd", B="abcdef"[:4] + "e", C="de")
    groups = group_redundant_terms(items, hits, threshold=0.75)
    assert [g["representative"]["Term"] for g in groups] == ["A", "C"]


def test_unknown_or_empty_hits_stay_alone():
    items = [_item("A", 0.01), _item("B", 0.02), _item("C", 0.03)]
    hits = {("GO", "A"): set(), ("GO", "B"): set()}      # C has no entry
    groups = group_redundant_terms(items, hits)
    assert [g["representative"]["Term"] for g in groups] == ["A", "B", "C"]
    assert all(g["members"] == [] for g in groups)


def test_ties_break_by_term_name_and_input_is_not_mutated():
    items = [_item("Zeta", 0.02), _item("Alfa", 0.02)]
    snapshot = [dict(i) for i in items]
    groups = group_redundant_terms(items, _hits(Zeta="ab", Alfa="ab"))
    assert groups[0]["representative"]["Term"] == "Alfa"
    assert items == snapshot
    assert group_redundant_terms([], {}) == []


def test_pvalues_are_never_modified():
    items = [_item("A", 0.01), _item("B", 0.0424)]
    groups = group_redundant_terms(items, _hits(A="ab", B="ab"))
    assert groups[0]["members"][0]["Pval"] == 0.0424
