"""Error handling of the live research fetchers (auditoría externa, H-20).

H-20 listed as untested: HTTP 429 / timeouts / network errors in
`safe_pubmed_request`, `enrich_genes_direct`, `fetch_omim_for_gene`,
`fetch_clinvar_for_gene` and `fetch_clinicaltrials_for_gene`. The
code already returns empty fallbacks instead of raising; these tests
prove it with mocked transports (no network), so a future refactor
cannot turn a dead API into a 500.

Pattern follows tests/test_opentargets.py: httpx.MockTransport +
asyncio.run() inside sync test functions. asyncio.sleep is faked so
the 429 backoffs (2-4-6 s) cost nothing.
"""
import asyncio

import httpx
import pytest

import kenryu_engine as ke


def _client(handler):
    return httpx.AsyncClient(transport=httpx.MockTransport(handler), trust_env=False)


@pytest.fixture
def no_sleep(monkeypatch):
    async def _fake(*args, **kwargs):
        return None

    monkeypatch.setattr(asyncio, "sleep", _fake)


# ── safe_pubmed_request ──────────────────────────────────────────────────────

def test_safe_pubmed_request_retries_429_then_succeeds(no_sleep):
    """429 twice, then 200: the retry returns the payload."""
    calls = {"n": 0}

    def handler(request):
        calls["n"] += 1
        if calls["n"] < 3:
            return httpx.Response(429, json={})
        return httpx.Response(200, json={"esearchresult": {"count": "7"}})

    async def _go():
        async with _client(handler) as client:
            return await ke.safe_pubmed_request(client, "https://x", {}, retries=3)

    assert asyncio.run(_go()) == {"esearchresult": {"count": "7"}}
    assert calls["n"] == 3


def test_safe_pubmed_request_gives_up_after_retries(no_sleep):
    """429 on every attempt: None after exactly `retries` calls."""
    calls = {"n": 0}

    def handler(request):
        calls["n"] += 1
        return httpx.Response(429, json={})

    async def _go():
        async with _client(handler) as client:
            return await ke.safe_pubmed_request(client, "https://x", {}, retries=3)

    assert asyncio.run(_go()) is None
    assert calls["n"] == 3


def test_safe_pubmed_request_returns_none_on_500(no_sleep):
    def handler(request):
        return httpx.Response(500, text="boom")

    async def _go():
        async with _client(handler) as client:
            return await ke.safe_pubmed_request(client, "https://x", {}, retries=3)

    assert asyncio.run(_go()) is None


def test_safe_pubmed_request_returns_none_on_timeout(no_sleep):
    def handler(request):
        raise httpx.ConnectTimeout("no route", request=request)

    async def _go():
        async with _client(handler) as client:
            return await ke.safe_pubmed_request(client, "https://x", {}, retries=2)

    assert asyncio.run(_go()) is None


# ── OMIM / ClinVar / ClinicalTrials ──────────────────────────────────────────

def test_fetch_omim_returns_empty_on_429(no_sleep):
    """OMIM esearch 429s: [] instead of an exception."""
    def handler(request):
        return httpx.Response(429, json={})

    async def _go():
        async with _client(handler) as client:
            return await ke.fetch_omim_for_gene("ABCA1", client)

    assert asyncio.run(_go()) == []


def test_fetch_clinvar_returns_zero_on_timeout(no_sleep):
    """ClinVar esearch times out: {"count": 0, "url": ""}."""
    def handler(request):
        raise httpx.ReadTimeout("slow", request=request)

    async def _go():
        async with _client(handler) as client:
            return await ke.fetch_clinvar_for_gene("ABCA1", client)

    assert asyncio.run(_go()) == {"count": 0, "url": ""}


def test_fetch_clinicaltrials_returns_empty_on_500(no_sleep):
    def handler(request):
        return httpx.Response(500, json={})

    async def _go():
        async with _client(handler) as client:
            return await ke.fetch_clinicaltrials_for_gene("ABCA1", client)

    assert asyncio.run(_go()) == []


def test_fetch_clinicaltrials_returns_empty_on_timeout(no_sleep):
    def handler(request):
        raise httpx.ConnectTimeout("no route", request=request)

    async def _go():
        async with _client(handler) as client:
            return await ke.fetch_clinicaltrials_for_gene("ABCA1", client)

    assert asyncio.run(_go()) == []


# ── enrich_genes_direct (live Enrichr path) ───────────────────────────────────

def _clear_enrichr_cache():
    ke.PERSISTENT_CACHE["enrichr"].clear()


def test_enrich_genes_direct_short_list_returns_empty_without_network():
    """Fewer than 3 genes: [] without touching the network."""
    _clear_enrichr_cache()
    calls = {"n": 0}

    def handler(request):
        calls["n"] += 1
        return httpx.Response(200, json={})

    async def _go():
        async with _client(handler) as client:
            return await ke.enrich_genes_direct(["ABCA1", "KPNA3"], client)

    assert asyncio.run(_go()) == []
    assert calls["n"] == 0


def test_enrich_genes_direct_returns_empty_on_addlist_429(no_sleep):
    """Enrichr addList 429: [] instead of an exception."""
    _clear_enrichr_cache()

    def handler(request):
        return httpx.Response(429, json={})

    async def _go():
        async with _client(handler) as client:
            return await ke.enrich_genes_direct(
                ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"], client)

    assert asyncio.run(_go()) == []


def test_enrich_genes_direct_returns_empty_on_timeout(no_sleep):
    """Enrichr addList times out twice: [] instead of an exception."""
    _clear_enrichr_cache()

    def handler(request):
        raise httpx.ConnectTimeout("no route", request=request)

    async def _go():
        async with _client(handler) as client:
            return await ke.enrich_genes_direct(
                ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"], client)

    assert asyncio.run(_go()) == []
