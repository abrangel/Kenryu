"""OpenAPI export and usage examples (Unreleased).

docs/openapi.json must match the app's live schema, and the three
example payloads in docs/examples/ must run against the API.
"""
import json
from pathlib import Path

import kenryu_engine as ke

ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = ROOT / "docs" / "examples"


def _client():
    from fastapi.testclient import TestClient
    return TestClient(ke.app)


def test_openapi_export_is_in_sync():
    exported = json.loads((ROOT / "docs" / "openapi.json").read_text(encoding="utf-8"))
    with _client() as client:
        live = client.get("/openapi.json").json()
    assert exported["openapi"] == live["openapi"]
    assert exported["info"]["version"] == live["info"]["version"] == ke.APP_VERSION
    for path in ("/api/v1/analyze", "/api/v1/compare", "/api/v1/batch",
                 "/api/v1/reverse", "/health"):
        assert path in exported["paths"]
    assert exported == live


def test_example_payloads_run():
    cases = [
        ("analyze-lipid.json", "/api/v1/analyze",
         ["common_genes", "convergence", "execution_manifest"]),
        ("analyze-cancer-n1.json", "/api/v1/analyze",
         ["common_genes", "convergence", "execution_manifest"]),
        ("compare-lipid-cancer.json", "/api/v1/compare",
         ["panels", "n_union"]),
    ]
    with _client() as client:
        for name, path, fields in cases:
            payload = json.loads((EXAMPLES / name).read_text(encoding="utf-8"))
            r = client.post(path, json=payload)
            assert r.status_code == 200, name
            for f in fields:
                assert f in r.json(), (name, f)
