"""Enrichment routes (Unreleased, 2026-11-19).

The LOCAL Fisher test with the TargetScan-universe background is the
method; the live Enrichr route is an exploratory cross-check whose
adjusted p-values are not a standard FDR at small list sizes (see
docs/validation/calibration-enrichr-route.md). The API must return
both lists separately, and the UI must label the Enrichr adjusted
p-values as exploratory — presenting them as significance in the
interface would mislead.
"""
import asyncio

import kenryu_engine as ke


class _FakeEnrichrResp:
    def __init__(self, status, payload):
        self.status_code = status
        self._payload = payload

    def json(self):
        return self._payload


class _FakeEnrichrClient:
    """Enrichr falso: un término significativo (p-adj 0,01) por biblioteca."""

    def __init__(self, rows=None):
        self._rows = rows if rows is not None else [
            [1, "Term A", 0.001, 5.0, 10.0, "GENE", 0.01],
        ]

    async def post(self, url, files=None, timeout=None):
        return _FakeEnrichrResp(200, {"userListId": 1})

    async def get(self, url, params=None, timeout=None):
        return _FakeEnrichrResp(200, {params["backgroundType"]: list(self._rows)})


def _term_genes():
    """A gene list the local route enriches: the genes of one Reactome
    term (20-50 genes in the TargetScan universe), so the term itself
    comes back with a tiny adjusted p."""
    ke.TARGETSCAN_DB = ke.load_local_data()
    ke.load_gmt_libraries()
    universe = set().union(*ke.TARGETSCAN_DB.values())
    for name in sorted(ke.GMT_LIBRARIES["Reactome"]):
        genes = sorted(set(ke.GMT_LIBRARIES["Reactome"][name]) & universe)
        if 20 <= len(genes) <= 50:
            return genes, name, universe
    raise AssertionError("no Reactome term with 20-50 genes in the universe")


def test_compute_enrichment_routes_offline_returns_local_only():
    genes, term_name, _ = _term_genes()
    local_terms, enrichr_terms = asyncio.run(
        ke.compute_enrichment_routes(genes, offline=True, client=None))
    assert enrichr_terms == []
    assert any(t["Term"] == term_name for t in local_terms)
    assert all(t["Pval"] < 0.05 for t in local_terms)


def test_compute_enrichment_routes_online_splits_routes():
    genes, term_name, _ = _term_genes()
    # La caché de Enrichr (enr-v2) es por lista de genes: limpiarla para
    # que el cliente falso de este test sea el que responda, no la
    # entrada que dejó el test anterior con la misma lista.
    ke.PERSISTENT_CACHE["enrichr"].clear()
    local_off, _ = asyncio.run(
        ke.compute_enrichment_routes(genes, offline=True, client=None))
    local_on, enrichr_on = asyncio.run(
        ke.compute_enrichment_routes(genes, offline=False,
                                     client=_FakeEnrichrClient()))
    # The local route is the same online and offline: it is the method.
    assert [t["Term"] for t in local_on] == [t["Term"] for t in local_off]
    assert any(t["Term"] == term_name for t in local_on)
    # The Enrichr terms are a separate list, one per library.
    assert len(enrichr_on) == 4
    assert all(t["Term"] == "Term A" for t in enrichr_on)
    assert all(t["Pval"] == 0.01 for t in enrichr_on)


def test_compute_enrichment_routes_online_empty_enrichr():
    genes, term_name, _ = _term_genes()
    ke.PERSISTENT_CACHE["enrichr"].clear()
    local_terms, enrichr_terms = asyncio.run(
        ke.compute_enrichment_routes(genes, offline=False,
                                     client=_FakeEnrichrClient(rows=[])))
    assert any(t["Term"] == term_name for t in local_terms)
    assert enrichr_terms == []
