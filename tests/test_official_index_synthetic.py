"""Official-format miRTarBase tests that run in CI against the synthetic
index (tests/fixtures/mirtarbase_official_synthetic.json.gz).

The real 2.3 MB index is gitignored, so the equivalent tests against the
real file (test_engine.py, test_evidence_table.py, test_local_data.py)
skip in CI and only run where the index was built. These run everywhere:
the fixture has the same schema and the thesis panel's 7 pairs copied
verbatim from the real index, so the behaviour under test is identical
for the panel — and the version string says "synthetic" instead of
claiming the real file's counts.
"""
from pathlib import Path

import kenryu_engine as ke

_THESIS = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
           "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
_FIXTURE = Path(__file__).resolve().parent / "fixtures" / "mirtarbase_official_synthetic.json.gz"
_SYNTHETIC_VERSION = "miRTarBase official index (synthetic CI fixture, 5 miRNAs, 7 pairs)"


def test_synthetic_index_loads_with_thesis_values():
    """The fixture parses with the real loader and carries the thesis
    panel's values verbatim (ABCA1 strong with its 3 PMIDs, SNTB2 and
    TSC22D2 weak, KPNA3/SCN1A absent)."""
    import kenryu_core as kc
    db = kc.load_mirtarbase_official(_FIXTURE)
    assert len(db) == 5
    assert db["hsa-miR-33a-5p"]["ABCA1"]["level"] == "strong"
    assert db["hsa-miR-33a-5p"]["ABCA1"]["pmids"] == [20466882, 20466885, 20566875]
    assert db["hsa-miR-33a-5p"]["SNTB2"]["level"] == "weak"
    assert db["hsa-miR-758-3p"]["TSC22D2"]["level"] == "weak"
    assert "KPNA3" not in db["hsa-miR-33a-5p"]
    assert "SCN1A" not in db["hsa-miR-33a-5p"]


def test_synthetic_index_version_string_is_honest(synthetic_official_index):
    """The version string for the fixture must say it is synthetic —
    it must never claim the real file's 957,040 rows / 2,947 miRNAs."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["data_sources"]["mirtarbase"] == _SYNTHETIC_VERSION


def test_analyze_require_experimental_filters_thesis_panel(synthetic_official_index):
    """v1.48: con require_experimental=true, el panel de la tesis conserva
    solo los genes con evidencia en miRTarBase (ABCA1 fuerte, SNTB2 y
    TSC22D2 débil); KPNA3 y SCN1A (solo predichos) salen del consenso."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True, "require_experimental": True})
    assert r.status_code == 200
    assert sorted(r.json()["common_genes"]) == ["ABCA1", "SNTB2", "TSC22D2"]


def test_analyze_validation_report_has_evidence_levels(synthetic_official_index):
    """v1.48: el dossier trae el nivel de evidencia por gen (fuerte /
    débil / solo predicho) con sus PMIDs."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    genes = r.json()["validation"]["genes"]
    assert genes["ABCA1"]["evidence_level"] == "strong"
    assert 20466885 in genes["ABCA1"]["mirtarbase_pmids"]
    assert genes["SNTB2"]["evidence_level"] == "weak"
    assert genes["TSC22D2"]["evidence_level"] == "weak"


def test_analyze_response_has_data_versions_top_level(synthetic_official_index):
    """Idea 2: las versiones van en el nivel superior del JSON del análisis,
    y el dossier y la exportación llevan las mismas. Con el índice
    sintético, la versión dice que es sintético."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    body = r.json()
    assert body["data_versions"]["mirtarbase"] == _SYNTHETIC_VERSION
    assert body["data_versions"]["targetscan"].startswith("8.0")
    assert body["validation"]["data_sources"] == body["data_versions"]
    assert body["export"]["data_sources"] == body["data_versions"]


def test_thesis_panel_evidence_counts_and_table_live(synthetic_official_index):
    """El panel de la tesis offline: 1 fuerte, 2 débiles, 2 solo predichos,
    y la tabla por pareja con los valores del índice sintético."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": _THESIS, "years": 10, "mode": "strict",
            "offline": True})
    assert r.status_code == 200
    v = r.json()["validation"]
    s = v["summary"]
    assert (s["n_strong"], s["n_weak"], s["n_only_predicted"]) == (1, 2, 2)
    abca1 = v["genes"]["ABCA1"]["evidence_by_mirna"]
    assert abca1["hsa-miR-33a-5p"]["level"] == "strong"
    assert 20466885 in abca1["hsa-miR-33a-5p"]["pmids"]
    assert abca1["hsa-miR-33b-5p"]["level"] == "strong"
    assert abca1["hsa-miR-106b-5p"]["level"] == "weak"
    assert abca1["hsa-miR-144-3p"]["level"] == "none"
    assert abca1["hsa-miR-758-3p"]["level"] == "none"
    kpna3 = v["genes"]["KPNA3"]["evidence_by_mirna"]
    assert all(e["level"] == "none" for e in kpna3.values())
