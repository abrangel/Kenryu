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


def test_script_js_has_evidence_by_mirna_table():
    """v1.51 (idea 1): the gene panel shows the per-miRNA evidence table,
    with the weak level labelled high-throughput support, not validation."""
    text = _script_text()
    assert "evidence_by_mirna" in text
    assert "apoyo de alto rendimiento" in text
    assert "Evidencia por miRNA" in text


def test_report_table_row_escapes_values():
    """v1.52 (bloque 2, XSS residual): the printable report's Section II
    row interpolated g, d.system, associated_routes and pathText raw.
    Every value must go through escapeHtml before innerHTML."""
    text = _script_text()
    row = next(l for l in text.splitlines() if "tableRows +=" in l)
    assert "escapeHtml(g)" in row
    assert "escapeHtml(d.system" in row
    assert "escapeHtml((d.associated_routes" in row
    assert "escapeHtml(pathText)" in row


def test_report_images_validate_base64():
    """v1.52 (bloque 2): the report's venn/volcano/ppi images interpolated
    data.*_plot raw into data:image/png;base64,.... Only validated
    base64 may reach an <img> src (same guard as setImg)."""
    text = _script_text()
    assert "data:image/png;base64,${data.venn_plot}" not in text
    assert "data:image/png;base64,${data.volcano_plot}" not in text
    assert "data:image/png;base64,${data.ppi_plot}" not in text
    assert "b64Img(" in text
