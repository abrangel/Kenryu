"""Tests for the validation core (kenryu_core.validation_score and
build_validation_report): every gene gets a score from independent,
verifiable evidence layers. Pure functions, no network.
"""

import kenryu_core as kc

FULL_RESEARCH = {
    "pubmed": [
        {"pmid": "111", "term": "x", "source": "PubMed"},
        {"pmid": "222", "term": "y", "source": "PubMed"},
        {"pmid": "333", "term": "z", "source": "PubMed"},
    ],
    "omim": [{"id": "600046"}],
    "clinvar": {"count": 214, "url": "https://example.com"},
    "trials": [{"nct": "NCT01456650"}],
}
EMPTY_RESEARCH = {
    "pubmed": [], "omim": [], "clinvar": {"count": 0, "url": ""}, "trials": [],
}


def test_full_evidence_scores_100():
    s = kc.validation_score("ABCA1", FULL_RESEARCH, predicted=True, validated=True)
    assert s["score"] == 100  # 25 + 30 + 15 + 10 + 10 + 10
    assert s["level"] == "alta"
    assert s["missing"] == []
    assert "miRTarBase" in ", ".join(s["layers"])
    assert "ABCA1" in s["explanation"]


def test_partial_evidence_scores_low():
    s = kc.validation_score("GENE1", EMPTY_RESEARCH, predicted=True, validated=False)
    assert s["score"] == 25
    assert s["level"] == "baja"
    assert "miRTarBase (validado experimentalmente)" in s["missing"]
    assert "PubMed" in s["missing"]


def test_no_evidence_scores_zero():
    s = kc.validation_score("GENE2", None, predicted=False, validated=False)
    assert s["score"] == 0
    assert s["level"] == "baja"
    assert s["layers"] == []


def test_pubmed_capped_at_15():
    research = dict(FULL_RESEARCH, pubmed=[{"pmid": str(i)} for i in range(10)])
    s = kc.validation_score("G", research, predicted=False, validated=False)
    assert s["score"] == 15 + 10 + 10 + 10  # pubmed capped, no pred/val
    assert s["level"] == "media"


def test_build_validation_report():
    report = kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p", "hsa-miR-144-3p"],
        common_genes=["ABCA1", "GENE2"],
        gene_research={"ABCA1": FULL_RESEARCH, "GENE2": EMPTY_RESEARCH},
        local_by_mirna={"hsa-miR-33a-5p": {"ABCA1", "GENE2"},
                        "hsa-miR-144-3p": {"ABCA1"}},
        remote_by_mirna={"hsa-miR-33a-5p": {"ABCA1"},
                         "hsa-miR-144-3p": {"ABCA1"}},
        data_versions={"mirbase": "22", "targetscan": "8.0"},
    )
    assert report["genes"]["ABCA1"]["score"] == 100
    assert report["genes"]["ABCA1"]["level"] == "alta"
    # GENE2 is predicted (in a local set) but not validated, no research.
    assert report["genes"]["GENE2"]["score"] == 25
    assert report["genes"]["GENE2"]["level"] == "baja"
    assert report["summary"]["n_genes"] == 2
    assert report["summary"]["high"] == 1
    assert report["summary"]["low"] == 1
    assert report["summary"]["mean_score"] == 62.5
    assert report["summary"]["not_validated"] == ["GENE2"]
    assert report["data_sources"]["mirbase"] == "22"
    assert "verificables" in report["method"]


def test_high_requires_experimental_validation():
    # 25 (predicho) + 15 (PubMed) + 10 + 10 + 10 = 70 puntos SIN miRTarBase:
    # la literatura y la clínica no pueden compensar la falta de evidencia
    # sobre la interacción. El puntaje se mantiene (70), el nivel queda en
    # "media" y la explicación dice por qué.
    s = kc.validation_score("GENE3", FULL_RESEARCH, predicted=True, validated=False)
    assert s["score"] == 70
    assert s["level"] == "media"
    assert "sin validación experimental" in s["explanation"]
    # Con la misma evidencia MÁS miRTarBase, sí es "alta".
    s2 = kc.validation_score("GENE3", FULL_RESEARCH, predicted=True, validated=True)
    assert s2["score"] == 100
    assert s2["level"] == "alta"


def test_build_validation_report_empty():
    report = kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p"], common_genes=[], gene_research={},
        local_by_mirna={}, remote_by_mirna={}, data_versions={})
    assert report["summary"]["n_genes"] == 0
    assert report["summary"]["mean_score"] is None
    assert report["genes"] == {}


# ── Ronda 3 (v1.48) — miRTarBase oficial + dossier de tres niveles ──────────
# Niveles de evidencia por (miRNA, gen): strong (Functional MTI: reporter,
# western, qPCR) +30; weak (Functional MTI (Weak): CLIP-seq, microarray,
# proteómica) +15; none (sin evidencia funcional) +0. La evidencia débil
# SÍ es validación experimental de la interacción (CLIP muestra unión
# física), así que satisface el candado de "alta" igual que la fuerte;
# el puntaje menor hace que "alta" sea más difícil de alcanzar con ella.


def test_validation_score_weak_is_half_of_strong():
    # Misma evidencia que el caso 100, pero validación débil: 25 + 15 +
    # 15 + 10 + 10 + 10 = 85 (no 100).
    s = kc.validation_score("ABCA1", FULL_RESEARCH, predicted=True, validated="weak")
    assert s["score"] == 85
    assert s["level"] == "alta"  # débil SÍ es evidencia experimental
    assert "débil" in ", ".join(s["layers"])


def test_validation_score_strong_equals_legacy_true():
    s_strong = kc.validation_score("ABCA1", FULL_RESEARCH, predicted=True, validated="strong")
    s_true = kc.validation_score("ABCA1", FULL_RESEARCH, predicted=True, validated=True)
    s_legacy = kc.validation_score("ABCA1", FULL_RESEARCH, predicted=True, validated="legacy")
    assert s_strong["score"] == s_true["score"] == s_legacy["score"] == 100


def test_validation_score_none_gives_no_points_and_caps_level():
    s = kc.validation_score("GENE3", FULL_RESEARCH, predicted=True, validated="none")
    assert s["score"] == 70
    assert s["level"] == "media"
    assert "sin validación experimental" in s["explanation"]


def test_build_validation_report_carries_evidence_levels():
    report = kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p", "hsa-miR-144-3p"],
        common_genes=["ABCA1", "SNTB2", "KPNA3"],
        gene_research={"ABCA1": FULL_RESEARCH, "SNTB2": FULL_RESEARCH,
                       "KPNA3": FULL_RESEARCH},
        local_by_mirna={"hsa-miR-33a-5p": {"ABCA1", "SNTB2", "KPNA3"},
                        "hsa-miR-144-3p": {"ABCA1", "SNTB2", "KPNA3"}},
        remote_by_mirna={"hsa-miR-33a-5p": {"ABCA1", "SNTB2"},
                         "hsa-miR-144-3p": {"SNTB2"}},
        evidence_by_mirna={
            "hsa-miR-33a-5p": {
                "ABCA1": {"level": "strong", "pmids": [20466885], "experiments": ["Luciferase reporter assay"]},
                "SNTB2": {"level": "weak", "pmids": [22012620], "experiments": ["PAR-CLIP"]}},
            "hsa-miR-144-3p": {
                "SNTB2": {"level": "weak", "pmids": [24906430], "experiments": ["HITS-CLIP"]}}},
        data_versions={"mirbase": "22", "targetscan": "8.0"},
    )
    # Fuerte: 100, alta. Débil: 85, alta (débil + resto de capas llega a 70+).
    assert report["genes"]["ABCA1"]["evidence_level"] == "strong"
    assert report["genes"]["ABCA1"]["score"] == 100
    assert report["genes"]["SNTB2"]["evidence_level"] == "weak"
    assert report["genes"]["SNTB2"]["score"] == 85
    assert report["genes"]["SNTB2"]["mirtarbase_pmids"] == [22012620, 24906430]
    assert report["genes"]["KPNA3"]["evidence_level"] == "none"
    assert report["genes"]["KPNA3"]["score"] == 70
    assert report["genes"]["KPNA3"]["level"] == "media"


def test_build_validation_report_not_evaluated_carries_evidence_level():
    # Los genes fuera del tope muestran su nivel local sin puntaje.
    report = kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p"],
        common_genes=["ABCA1", "KPNA3"],
        gene_research={"ABCA1": FULL_RESEARCH},
        local_by_mirna={"hsa-miR-33a-5p": {"ABCA1", "KPNA3"}},
        remote_by_mirna={"hsa-miR-33a-5p": {"ABCA1"}},
        evidence_by_mirna={"hsa-miR-33a-5p": {
            "ABCA1": {"level": "strong", "pmids": [20466885], "experiments": []}}},
        data_versions={},
    )
    assert report["genes"]["KPNA3"]["level"] == "not_evaluated"
    assert report["genes"]["KPNA3"]["evidence_level"] == "none"
    assert report["genes"]["ABCA1"]["evidence_level"] == "strong"
