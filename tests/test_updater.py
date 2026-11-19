"""Tests for the self-update checker (kenryu_updater.py).

All network is mocked (httpx.MockTransport) except the endpoint smoke
test, which runs against the real code path with a monkeypatched checker.
"""

import json

import httpx
import pytest

import kenryu_updater as ku

MIRBASE_README = (
    "The miRBase Sequence Database -- Release 23.1 (Sep 2026)\n"
    "--------------------------------------------------------\n"
)
TARGETSCAN_SNIPPET = "<title>TargetScanHuman 8.0</title>"
MIRTARBASE_SNIPPET = "<h1>miRTarBase 10.0</h1>"
ENRICHR_STATS = json.dumps({"statistics": [{"libraryName": "KEGG_2021_Human"}]})


def _mock_client(handler):
    return httpx.Client(transport=httpx.MockTransport(handler), trust_env=False)


# ── Version comparison ───────────────────────────────────────────────────

def test_compare_versions():
    assert ku.compare_versions("22", "23.1") == -1
    assert ku.compare_versions("8.0", "8.0") == 0
    assert ku.compare_versions("10.0", "9.0") == 1
    assert ku.compare_versions("22", "22.1") == -1
    # Unpinned / live / unparseable can never be compared.
    assert ku.compare_versions("live via Harmonizome (unpinned)", "10.0") is None
    assert ku.compare_versions("live API", "live") is None


def test_version_tuple():
    assert ku.version_tuple("23.1") == (23, 1)
    assert ku.version_tuple("22") == (22,)
    assert ku.version_tuple("nonsense") is None


# ── Parsers (real page shapes, verified 2026-11-18) ─────────────────────

def test_parsers():
    assert ku.parse_mirbase_readme(MIRBASE_README) == "23.1"
    assert ku.parse_mirbase_readme("nothing here") is None
    assert ku.parse_targetscan_page(TARGETSCAN_SNIPPET) == "8.0"
    assert ku.parse_mirtarbase_page(MIRTARBASE_SNIPPET) == "10.0"
    assert ku.parse_enrichr_stats(ENRICHR_STATS) == "live"
    assert ku.parse_enrichr_stats("not json") is None
    assert ku.parse_enrichr_stats("{}") is None


# ── Full check with mocked network ───────────────────────────────────────

def _handler(request):
    url = str(request.url)
    if "mirbase.org" in url:
        return httpx.Response(200, text=MIRBASE_README)
    if "targetscan.org" in url:
        return httpx.Response(200, text=TARGETSCAN_SNIPPET)
    if "mirtarbase.cuhk.edu.cn" in url:
        return httpx.Response(200, text=MIRTARBASE_SNIPPET)
    if "Enrichr" in url:
        return httpx.Response(200, text=ENRICHR_STATS)
    return httpx.Response(404, text="not found")


def test_check_updates_mocked():
    report = ku.check_updates(client=_mock_client(_handler))
    by_key = {s["key"]: s for s in report["sources"]}
    # miRBase local 22 vs live 23.1 -> outdated (the real finding).
    assert by_key["mirbase"]["status"] == "outdated"
    assert by_key["mirbase"]["latest_checked"] == "23.1"
    # TargetScan 8.0 == 8.0 -> up to date.
    assert by_key["targetscan"]["status"] == "up-to-date"
    # miRTarBase local version is not pinned -> unpinned, never "current".
    assert by_key["mirtarbase"]["status"] == "unpinned"
    assert by_key["mirtarbase"]["latest_checked"] == "10.0"
    # Enrichr answers -> live.
    assert by_key["enrichr"]["status"] == "live"
    assert report["updates_available"] == 1
    assert report["unpinned"] == 1


def test_check_updates_offline_never_raises():
    def _down(request):
        raise httpx.ConnectError("no network", request=request)

    report = ku.check_updates(client=_mock_client(_down))
    assert all(s["status"] == "unknown" for s in report["sources"])
    assert report["updates_available"] == 0
    assert report["unknown"] == len(report["sources"])
    # And the human summary still renders.
    assert len(ku.summarize(report)) == len(report["sources"])


def test_make_client_constructs_with_broken_no_proxy():
    # This VM's NO_PROXY holds bracketed IPv6 ([::1]); httpx 0.26 with
    # trust_env=True raises InvalidURL on construction. The factory must
    # never raise, whatever the environment holds.
    client = ku._make_client(timeout=5.0)
    client.close()


def test_load_local_versions():
    versions = ku.load_local_versions()
    assert versions["mirbase"] == "22"
    assert versions["targetscan"] == "8.0"
    assert versions["mirtarbase"] == "live via Harmonizome (unpinned)"
    assert versions["enrichr"] == "live API"
    # Missing file falls back to the registry defaults.
    assert ku.load_local_versions("/nonexistent/path.json")["mirbase"] == "22"


# ── API endpoint ─────────────────────────────────────────────────────────

def test_updates_endpoint(monkeypatch):
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    canned = ku.check_updates(client=_mock_client(_handler))
    monkeypatch.setattr(ke, "check_updates", lambda timeout=8.0: canned)
    with TestClient(ke.app) as client:
        r = client.get("/api/v1/updates")
    assert r.status_code == 200
    body = r.json()
    assert body["updates_available"] == 1
    assert {s["key"]: s["status"] for s in body["sources"]}["mirbase"] == "outdated"


def test_health_reports_data_sources():
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    with TestClient(ke.app) as client:
        r = client.get("/health")
    assert r.status_code == 200
    assert r.json()["data_sources"]["mirbase"] == "22"
    assert r.json()["data_sources"]["targetscan"] == "8.0"
