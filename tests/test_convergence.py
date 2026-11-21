"""Panel convergence test in the report (Unreleased).

Approved 2026-11-20 with conditions: seed derived from the sorted
family names (same input -> same number); ~2,000 permutations in the
API with the resolution shown next to p; p reported as a range across
the matching tolerance (+-10% to +-50%), not a single value; fixed
label; effective panel size shown with a warning when miRNAs share a
family; the matching controls only target count, not shared popular
genes.
"""
import pytest

import kenryu_engine as ke
from kenryu_core import (
    CONVERGENCE_LABEL,
    CONVERGENCE_LIMIT,
    panel_convergence,
)

LIPID = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
         "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
CANCER = ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"]


def _panel(mirnas):
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.TARGETSCAN_DB
    sets = [ke.get_targets_local(m, db) for m in mirnas]
    labels = ke.family_labels_for_targets(mirnas, sets, db)
    return sets, labels, db


def test_seed_is_deterministic_and_order_independent():
    sets, labels, db = _panel(LIPID)
    a = panel_convergence(sets, labels, db)
    b = panel_convergence(sets, list(reversed(labels)), db)
    assert a["p_by_tolerance"] == b["p_by_tolerance"]
    assert a["p_min"] == b["p_min"] and a["p_max"] == b["p_max"]


def test_lipid_panel_reproduces_reference():
    sets, labels, db = _panel(LIPID)
    r = panel_convergence(sets, labels, db, perms=10000, tolerances=(0.25,))
    assert r["k_names"] == 5
    assert r["k_effective"] == 4  # miR-33a and miR-33b share a family
    assert r["same_family_warning"] is True
    assert r["genes_observed"] == 5
    # Reference (convergence_test.py, seed 7, 10,000 perms): p = 0.137.
    # The seed here derives from the family names, so the draw differs;
    # the conclusion (nothing unusual at 0.05) must hold.
    assert 0.05 < r["p_by_tolerance"]["0.25"] < 0.25
    assert r["resolution"] == pytest.approx(1 / 10001)


def test_cancer_panel_reproduces_reference():
    sets, labels, db = _panel(CANCER)
    r = panel_convergence(sets, labels, db, perms=10000, tolerances=(0.25,))
    assert r["k_names"] == 3
    assert r["k_effective"] == 3
    assert r["same_family_warning"] is False
    # v1.59: 19, not 18 — the HGNC map unifies C10orf12 (miR-21-5p,
    # miR-155-5p) with LCOR (let-7-5p): the same gene under two symbols.
    assert r["genes_observed"] == 19
    assert 0.05 < r["p_by_tolerance"]["0.25"] < 0.20  # reference: 0.082


def test_p_is_a_range_across_tolerances_with_resolution():
    sets, labels, db = _panel(CANCER)
    r = panel_convergence(sets, labels, db)  # API defaults: 2,000 perms
    assert set(r["p_by_tolerance"]) == {"0.10", "0.25", "0.50"}
    assert r["perms"] == 2000
    assert r["resolution"] == pytest.approx(1 / 2001)
    assert r["p_min"] <= r["p_max"]
    assert all(0 < p < 1 for p in r["p_by_tolerance"].values())
    assert r["label"] == CONVERGENCE_LABEL
    assert r["limit"] == CONVERGENCE_LIMIT
    assert r["skipped"] is None


def test_too_many_families_is_skipped():
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.TARGETSCAN_DB
    # First 12 DISTINCT target sets (composite families share sets).
    names, seen = [], set()
    for n in sorted(db):
        if frozenset(db[n]) not in seen:
            seen.add(frozenset(db[n]))
            names.append(n)
        if len(names) == 12:
            break
    sets = [set(db[n]) for n in names]
    r = panel_convergence(sets, names, db)
    assert r["k_effective"] == 12
    assert r["skipped"] == "too_many_families"
    assert r["p_by_tolerance"] == {}
    assert r["genes_observed"] == len(set.intersection(*[set(s) for s in sets]))


def test_single_family_is_skipped():
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.TARGETSCAN_DB
    name = sorted(db)[0]
    r = panel_convergence([set(db[name])], [name], db)
    assert r["k_effective"] == 1
    assert r["skipped"] == "needs_at_least_two_families"
    assert r["p_by_tolerance"] == {}


def test_analyze_returns_convergence_for_lipid_panel():
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze",
                        json={"mirnas": LIPID, "mode": "strict",
                              "offline": True})
    assert r.status_code == 200
    conv = r.json()["convergence"]
    assert conv["k_names"] == 5
    assert conv["k_effective"] == 4
    assert conv["same_family_warning"] is True
    assert conv["genes_observed"] == 5
    assert set(conv["p_by_tolerance"]) == {"0.10", "0.25", "0.50"}
    assert conv["resolution"] == pytest.approx(1 / 2001)
    assert conv["label"] == CONVERGENCE_LABEL
    assert conv["limit"] == CONVERGENCE_LIMIT
