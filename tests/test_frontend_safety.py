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


def test_legacy_warnings_render_escaped():
    """v1.54: the visible legacy-mode warning banner renders data.warnings
    through escapeHtml, never raw (the warning text comes from the server,
    but the same rule applies to every interpolated value)."""
    text = _script_text()
    line = next(l for l in text.splitlines() if "data.warnings" in l)
    assert "escapeHtml(w)" in line


def test_specificity_tooltip_says_random_gene_lists():
    # La especificidad muestrea LISTAS aleatorias de genes del universo
    # de TargetScan, no paneles de miRNAs: el texto debe decirlo así.
    js = _script_text()
    assert "listas aleatorias de genes" in js
    assert "paneles aleatorios" not in js


def test_specificity_zero_renders_less_than_one_percent():
    # Con 100 muestras la resolución es del 1%: un 0 observado es "<1%",
    # nunca "0%".
    js = _script_text()
    assert "'<1%'" in js
    assert "0% = específico" not in js


def test_enrichr_section_labels_adjusted_p_exploratory():
    """Unreleased (2026-11-19): the live Enrichr adjusted p is not a
    standard FDR at small list sizes, so the UI must label it as
    exploratory and render the Enrichr terms in their own section,
    separate from the local (method) enrichment."""
    text = _script_text()
    assert "p ajustado de Enrichr (exploratorio)" in text
    assert "enrichment_enrichr" in text
    assert "Enrichr (exploratorio)" in text


def test_index_html_has_enrichr_exploratory_container():
    """Unreleased (2026-11-19): the page has a separate container for
    the Enrichr cross-check, so the two routes never share a table."""
    from fastapi.testclient import TestClient
    with TestClient(ke.app) as client:
        r = client.get("/")
    assert r.status_code == 200
    assert "enrich-container-enrichr" in r.text
    assert "Enrichr (exploratorio)" in r.text


def test_specificity_formatter_keeps_the_five_percent_threshold_readable():
    # Unreleased: valor exacto. <1% por debajo del 1%, un decimal entre
    # 1% y 10% (un 4,7% no debe leerse como 5%), entero a partir del 10%.
    js = _script_text()
    assert "function fmtSpec(" in js
    assert "toFixed(1)" in js and "'<1%'" in js
    assert "Math.round(item.specificity * 100)" not in js


def test_specificity_header_and_tooltip_use_the_exact_wording():
    js = _script_text()
    assert js.count('title="${SPEC_TITLE}"') == 2          # las dos tablas
    assert "cola hipergeométrica" in js and "no es una prueba de significancia" in js


def test_group_badge_is_escaped_and_members_are_hidden_from_the_main_table():
    js = _script_text()
    assert js.count("is_group_representative !== false") == 2
    assert "escapeHtml(gm.join('; '))" in js          # el título del badge se escapa
    assert "Array.isArray(item.group_members)" in js   # y solo si es una lista
