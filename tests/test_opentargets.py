"""Open Targets disease associations + the four-level dossier.

Open Targets (free GraphQL, no key) is the disease layer of the
dossier: symbol -> Ensembl ID via mygene.info, then associated
diseases with their association scores. It fills the "diseases" key
of the gene research ONLY when DisGeNET has no key configured
(DisGeNET stays the source when it has one); the origin is recorded
in "diseases_source". It never enters the validation score: a
disease association says the gene matters, not that the miRNA
regulates it.

The four dossier tiers, per gene:
  1. prediction — TargetScan predicted the gene for the panel.
  2. experimental — miRTarBase level (strong / weak / none).
  3. disease — Open Targets / DisGeNET associations.
  4. clinical_literature — PubMed / OMIM / ClinVar / ClinicalTrials.
"""
import asyncio

import httpx
import pytest

import kenryu_engine as ke
from kenryu_core import build_evidence_tiers, build_validation_report

MYGENE_OK = {"hits": [{"ensembl": {"gene": "ENSG00000165029"}}]}
MYGENE_LIST = {"hits": [{"ensembl": [{"gene": "ENSG00000165029"}]}]}
MYGENE_NONE = {"hits": []}
OT_ROWS = {"data": {"target": {"associatedDiseases": {"rows": [
    {"disease": {"name": "Tangier disease"}, "score": 0.838},
    {"disease": {"name": "Hypoalphalipoproteinemia"}, "score": 0.748},
    {"disease": {"name": "Metabolic syndrome"}, "score": 0.584},
    {"disease": {"name": "Hypercholesterolemia"}, "score": 0.571},
]}}}}


def _client(handler):
    return httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False)


def _mock_ot(request):
    url = str(request.url)
    if "mygene.info" in url:
        return httpx.Response(200, json=MYGENE_OK)
    if "opentargets.org" in url:
        return httpx.Response(200, json=OT_ROWS)
    return httpx.Response(404, json={})


async def _run(coro):
    return await coro


def test_fetch_opentargets_diseases_mocked():
    async def _go():
        async with _client(_mock_ot) as client:
            return await ke.fetch_opentargets_diseases("ABCA1", client)

    res = asyncio.run(_go())
    assert res["ensembl_id"] == "ENSG00000165029"
    assert res["diseases"] == [
        {"name": "Tangier disease", "score": 0.838},
        {"name": "Hypoalphalipoproteinemia", "score": 0.748},
        {"name": "Metabolic syndrome", "score": 0.584},
    ]  # top 3 by score; the 4th row is dropped


def test_fetch_opentargets_handles_ensembl_list_shape():
    def _mock(request):
        url = str(request.url)
        if "mygene.info" in url:
            return httpx.Response(200, json=MYGENE_LIST)
        if "opentargets.org" in url:
            return httpx.Response(200, json=OT_ROWS)
        return httpx.Response(404, json={})

    async def _go():
        async with _client(_mock) as client:
            return await ke.fetch_opentargets_diseases("ABCA1", client)

    res = asyncio.run(_go())
    assert res["ensembl_id"] == "ENSG00000165029"
    assert len(res["diseases"]) == 3


def test_fetch_opentargets_handles_missing_ensembl():
    def _mock(request):
        url = str(request.url)
        if "mygene.info" in url:
            return httpx.Response(200, json=MYGENE_NONE)
        return httpx.Response(404, json={})

    async def _go():
        async with _client(_mock) as client:
            return await ke.fetch_opentargets_diseases("NOPE", client)

    res = asyncio.run(_go())
    assert res == {"diseases": []}


def test_fetch_opentargets_never_raises():
    def _mock(request):
        raise httpx.ConnectError("boom")

    async def _go():
        async with _client(_mock) as client:
            return await ke.fetch_opentargets_diseases("ABCA1", client)

    res = asyncio.run(_go())
    assert res == {"diseases": []}


def test_research_gene_complete_uses_opentargets_when_disgenet_empty(monkeypatch):
    """v1.49: sin clave DISGENET_API_KEY, las enfermedades del dossier
    vienen de Open Targets y la caché sube a v5 (el resultado gana un
    campo nuevo)."""
    monkeypatch.delenv("DISGENET_API_KEY", raising=False)
    ke.PERSISTENT_CACHE["gene_research"].pop("v5|ABCA1|10", None)
    ke.PERSISTENT_CACHE["gene_research"].pop("v4|ABCA1|10", None)

    def _mock(request):
        url = str(request.url)
        params = request.url.params
        if "clinicaltrials.gov" in url:
            return httpx.Response(200, json={"studies": []})
        if "esearch.fcgi" in url:
            return httpx.Response(200, json={"esearchresult": {"idlist": [], "count": "0"}})
        if "mygene.info" in url:
            return httpx.Response(200, json=MYGENE_OK)
        if "opentargets.org" in url:
            return httpx.Response(200, json=OT_ROWS)
        return httpx.Response(200, json={"result": {}})

    async def _go():
        async with _client(_mock) as client:
            return await ke.research_gene_complete("ABCA1", {}, 10, client)

    res = asyncio.run(_go())
    assert res["diseases_source"] == "opentargets"
    assert res["diseases"][0] == {"name": "Tangier disease", "score": 0.838}
    assert "v5|ABCA1|10" in ke.PERSISTENT_CACHE["gene_research"]
    assert "v4|ABCA1|10" not in ke.PERSISTENT_CACHE["gene_research"]


def test_evidence_tiers_four_levels():
    research = {
        "pubmed": [{"pmid": "123"}],
        "omim": [{"omim_id": "600046"}],
        "clinvar": {"count": 3},
        "trials": [{"nct_id": "NCT00000000"}],
        "diseases": [{"name": "Tangier disease", "score": 0.838}],
    }
    tiers = build_evidence_tiers(True, "strong", research)
    assert tiers["prediction"]["present"] is True
    assert tiers["experimental"]["present"] is True
    assert tiers["disease"]["present"] is True
    assert tiers["clinical_literature"]["present"] is True
    assert sum(1 for t in tiers.values() if t["present"]) == 4


def test_evidence_tiers_prediction_only():
    tiers = build_evidence_tiers(True, "none", {})
    assert tiers["prediction"]["present"] is True
    assert tiers["experimental"]["present"] is False
    assert tiers["disease"]["present"] is False
    assert tiers["clinical_literature"]["present"] is False
    assert sum(1 for t in tiers.values() if t["present"]) == 1


def test_validation_report_carries_tiers():
    research = {"ABCA1": {
        "pubmed": [{"pmid": "123"}], "omim": [], "clinvar": {"count": 0},
        "trials": [], "diseases": [{"name": "Tangier disease", "score": 0.838}],
    }}
    report = build_validation_report(
        ["hsa-miR-33a-5p"], ["ABCA1"], research,
        {"hsa-miR-33a-5p": {"ABCA1"}}, {"hsa-miR-33a-5p": {"ABCA1"}}, {},
        evidence_by_mirna={"hsa-miR-33a-5p": {"ABCA1": {
            "level": "strong", "pmids": [20466885], "experiments": ["Luciferase reporter assay"]}}})
    tiers = report["genes"]["ABCA1"]["evidence_tiers"]
    assert tiers["experimental"]["present"] is True
    assert tiers["disease"]["present"] is True
    assert tiers["disease"]["summary"].startswith("1 enfermedad")
    assert report["genes"]["ABCA1"]["tiers_present"] == 4
