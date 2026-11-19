"""Frontend safety guards (auditoría externa, H-08 / H-21).

The XSS mitigation (escapeHtml) and the dynamic version live in
static/script.js, served by the app itself. These tests fetch the
served file and assert the guards are present, so a future edit of
the frontend cannot silently reintroduce raw innerHTML interpolation
or a hardcoded version.
"""
from fastapi.testclient import TestClient

import kenryu_engine as ke


def _script_text():
    with TestClient(ke.app) as client:
        r = client.get("/script.js")
    assert r.status_code == 200
    return r.text


def test_script_js_defines_escape_html():
    """H-08: every external/user value goes through escapeHtml."""
    assert "function escapeHtml" in _script_text()


def test_script_js_has_no_raw_gene_interpolation():
    """H-08: the old pattern openGenePanel('${g}') injected the gene
    name raw into an inline handler. The safe pattern uses data-gene
    + this.dataset.gene."""
    text = _script_text()
    assert "openGenePanel('${" not in text
    assert "data-gene" in text


def test_script_js_has_no_hardcoded_old_version():
    """H-21: the UI used to print v1.38 while the backend was v1.44.
    The version now comes from the API response (engine_version)."""
    text = _script_text()
    assert "v1.38" not in text
    assert "engine_version" in text
