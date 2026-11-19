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


def test_build_validation_report_empty():
    report = kc.build_validation_report(
        mirnas=["hsa-miR-33a-5p"], common_genes=[], gene_research={},
        local_by_mirna={}, remote_by_mirna={}, data_versions={})
    assert report["summary"]["n_genes"] == 0
    assert report["summary"]["mean_score"] is None
    assert report["genes"] == {}
