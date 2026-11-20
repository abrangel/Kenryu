"""Execution manifest (Unreleased).

One command must reproduce a run: engine version, commit, data
versions, parameters, seeds (including the convergence seed) and a
reproduce command. Same panel, same result — by API and by CLI.
"""
import json

import kenryu_engine as ke
from kenryu_core import panel_convergence

LIPID = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
         "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
EMPTY_PAIR = ["hsa-miR-412-5p", "hsa-miR-551a-3p"]  # strict intersection: 0 genes


def _analyze(mirnas, mode="strict"):
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze",
                        json={"mirnas": mirnas, "mode": mode, "offline": True})
    assert r.status_code == 200
    return r.json()


def test_convergence_seed_is_exposed_and_used():
    data = _analyze(LIPID)
    conv = data["convergence"]
    assert isinstance(conv["seed"], int)
    # Recomputing with the same labels must give the same p: the exposed
    # seed is the one the report used.
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.TARGETSCAN_DB
    sets = [ke.get_targets_local(m, db) for m in LIPID]
    labels = ke.family_labels_for_targets(LIPID, sets, db)
    again = panel_convergence(sets, labels, db)
    assert again["seed"] == conv["seed"]
    assert again["p_by_tolerance"] == conv["p_by_tolerance"]


def test_api_determinism_end_to_end():
    a = _analyze(LIPID)
    b = _analyze(LIPID)
    assert a["convergence"] == b["convergence"]
    assert a["execution_manifest"] == b["execution_manifest"]
    assert a["common_genes"] == b["common_genes"]


def test_manifest_has_everything_to_reproduce():
    data = _analyze(LIPID)
    m = data["execution_manifest"]
    assert m["engine_version"] == ke.APP_VERSION
    assert isinstance(m["commit"], str) and len(m["commit"]) >= 7
    assert "targetscan" in str(m["data_versions"]).lower() or m["data_versions"]
    p = m["parameters"]
    assert p["mirnas"] == LIPID
    assert p["mode"] == "strict"
    assert p["offline"] is True
    assert p["species"] == "human"
    assert p["require_experimental"] is False
    assert isinstance(m["seeds"]["convergence"], int)
    assert isinstance(m["reproduce"], str) and "hsa-miR-33a-5p" in m["reproduce"]


def test_cli_determinism_end_to_end(tmp_path):
    import kenryu_cli as cli
    panel = tmp_path / "panel.txt"
    panel.write_text("\n".join(LIPID) + "\n")
    out1, out2 = tmp_path / "out1", tmp_path / "out2"
    assert cli.main([str(panel), "--offline", "--out", str(out1)]) == 0
    assert cli.main([str(panel), "--offline", "--out", str(out2)]) == 0
    d1 = json.loads((out1 / "analysis.json").read_text())
    d2 = json.loads((out2 / "analysis.json").read_text())
    assert d1["convergence"] == d2["convergence"]
    assert d1["execution_manifest"] == d2["execution_manifest"]
    assert d1["common_genes"] == d2["common_genes"]
    # Prioritized order, not alphabetical: ABCA1 (validated by 2 panel
    # miRNAs) first, SNTB2/TSC22D2 (weak support) next, KPNA3/SCN1A
    # (predicted only) last.
    assert d1["common_genes"] == ["ABCA1", "SNTB2", "TSC22D2", "KPNA3", "SCN1A"]
    # The kept report carries the manifest, not just the JSON.
    assert "Manifiesto de ejecución" in (out1 / "report.md").read_text()


def test_manifest_in_empty_intersection_response():
    data = _analyze(EMPTY_PAIR)
    assert data["common_genes"] == []
    assert data["convergence"]["genes_observed"] == 0
    assert data["execution_manifest"]["parameters"]["mirnas"] == EMPTY_PAIR
