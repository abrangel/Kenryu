"""v1.51 (bloque 1, decisiones de Claude 2026-11-19): la tabla miRNA x gen
en el dossier, los conteos de evidencia en el resumen, las versiones de
datos en el JSON del análisis y el reetiquetado de la evidencia débil
como apoyo de alto rendimiento (no "validado").

Matiz 1 de Claude: "fuerte" significa que AL MENOS UN miRNA del panel
tiene evidencia fuerte para ese gen; la tabla por pareja cierra esa
ambigüedad. Matiz 2: la evidencia débil incluye CLIP-seq y otros datos
de alto rendimiento que no demuestran unión directa.
"""

import pytest

import kenryu_core as kc
import kenryu_engine as ke

_THESIS = ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
           "hsa-miR-106b-5p", "hsa-miR-33b-5p"]
_OFFICIAL = ke.BASE_DIR / "local_db" / "mirtarbase_official.json.gz"

_RESEARCH = {
    "pubmed": [{"pmid": "111", "term": "x", "source": "PubMed"}],
    "omim": [{"id": "600046"}],
    "clinvar": {"count": 3, "url": "https://example.com"},
    "trials": [],
}
_EMPTY = {"pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""},
          "trials": []}


def _synthetic_report(evidence_by_mirna, gene_research=None):
    return kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p", "hsa-miR-144-3p"],
        common_genes=["ABCA1", "GENE2", "GENE3"],
        gene_research=({"ABCA1": _RESEARCH, "GENE2": _EMPTY,
                        "GENE3": _RESEARCH} if gene_research is None
                       else gene_research),
        local_by_mirna={"hsa-miR-33a-5p": {"ABCA1", "GENE2", "GENE3"},
                        "hsa-miR-144-3p": {"ABCA1", "GENE2"}},
        remote_by_mirna={"hsa-miR-33a-5p": {"ABCA1"},
                         "hsa-miR-144-3p": {"GENE2"}},
        data_versions={"mirbase": "22", "targetscan": "8.0"},
        evidence_by_mirna=evidence_by_mirna,
    )


_SYNTHETIC_EVIDENCE = {
    "hsa-miR-33a-5p": {"ABCA1": {"level": "strong", "pmids": [20466885],
                                 "experiments": ["Luciferase reporter assay"]}},
    "hsa-miR-144-3p": {"GENE2": {"level": "weak", "pmids": [22012620],
                                 "experiments": ["PAR-CLIP"]}},
}


def test_weak_tier_is_high_throughput_support_not_validation():
    """Matiz 2 de Claude: la evidencia débil no se llama 'validado'."""
    weak = kc.build_evidence_tiers(True, "weak", {})["experimental"]["summary"]
    assert "Apoyo de alto rendimiento" in weak
    assert "Validado experimentalmente" not in weak
    assert "no demuestra unión directa" in weak
    strong = kc.build_evidence_tiers(True, "strong", {})["experimental"]["summary"]
    assert "Validado experimentalmente (miRTarBase, evidencia fuerte)" in strong


def test_summary_counts_evidence_levels():
    report = _synthetic_report(_SYNTHETIC_EVIDENCE)
    s = report["summary"]
    assert s["n_strong"] == 1  # ABCA1
    assert s["n_weak"] == 1  # GENE2
    assert s["n_legacy"] == 0
    assert s["n_only_predicted"] == 1  # GENE3


def test_evidence_by_mirna_table_covers_every_panel_mirna():
    """Matiz 1 de Claude: cada gen lleva una fila por cada miRNA del panel,
    con su nivel, sus PMIDs y sus experimentos (vacío cuando no hay)."""
    report = _synthetic_report(_SYNTHETIC_EVIDENCE)
    abca1 = report["genes"]["ABCA1"]["evidence_by_mirna"]
    assert set(abca1) == {"hsa-miR-33a-5p", "hsa-miR-144-3p"}
    assert abca1["hsa-miR-33a-5p"]["level"] == "strong"
    assert abca1["hsa-miR-33a-5p"]["pmids"] == [20466885]
    assert abca1["hsa-miR-33a-5p"]["experiments"] == [
        "Luciferase reporter assay"]
    assert abca1["hsa-miR-144-3p"] == {"level": "none", "pmids": [],
                                       "experiments": [], "bare_name": False}
    gene3 = report["genes"]["GENE3"]["evidence_by_mirna"]
    assert all(v["level"] == "none" for v in gene3.values())


def test_method_text_defines_strong_as_at_least_one_mirna():
    report = _synthetic_report(_SYNTHETIC_EVIDENCE)
    assert "al menos un miRNA del panel" in report["method"]
    assert "no demuestra unión directa" in report["method"]


def test_not_evaluated_layer_uses_best_level():
    """Un gen no evaluado con evidencia débil no se etiqueta 'validado'."""
    report = _synthetic_report(_SYNTHETIC_EVIDENCE,
                               gene_research={"ABCA1": _RESEARCH,
                                              "GENE3": _RESEARCH})
    layers = report["genes"]["GENE2"]["layers"]
    assert "miRTarBase (apoyo de alto rendimiento)" in layers
    assert not any("validado experimentalmente" in layer for layer in layers)


def test_loaded_data_versions_reports_not_loaded():
    """Las versiones del análisis dicen lo que está cargado, no lo que el
    JSON estático declara: sin índice oficial, sin KEGG y sin TargetScan,
    los tres dicen 'not loaded'."""
    import kenryu_engine as ke2
    orig_format = ke2.MIRTARBASE_FORMAT
    orig_db = ke2.TARGETSCAN_DB
    orig_gmt = dict(ke2.GMT_LIBRARIES)
    try:
        ke2.MIRTARBASE_FORMAT = None
        ke2.TARGETSCAN_DB = {}
        ke2.GMT_LIBRARIES.pop("KEGG", None)
        v = ke2.loaded_data_versions()
    finally:
        ke2.MIRTARBASE_FORMAT = orig_format
        ke2.TARGETSCAN_DB = orig_db
        ke2.GMT_LIBRARIES.clear()
        ke2.GMT_LIBRARIES.update(orig_gmt)
    assert v["mirtarbase"] == "not loaded"
    assert v["kegg"] == "not loaded"
    assert v["targetscan"] == "not loaded"
