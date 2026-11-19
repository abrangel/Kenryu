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
    # Pin the old local versions explicitly: the real data_versions.json
    # may already record a newer install, and the check must compare
    # against what the caller pins, not the registry default.
    pinned = {"mirbase": "22", "targetscan": "8.0",
              "mirtarbase": "live via Harmonizome (unpinned)", "enrichr": "live API",
              "kegg": "not installed"}
    report = ku.check_updates(client=_mock_client(_handler), local_versions=pinned)
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
    # KEGG is not installed -> unpinned (available to download on demand).
    assert by_key["kegg"]["status"] == "unpinned"
    assert report["updates_available"] == 1
    assert report["unpinned"] == 2


def test_check_updates_uses_pinned_versions():
    # After an install records 23.1, the same live 23.1 must read as
    # up-to-date, not outdated (regression: the checker used to read the
    # registry default "22" forever).
    pinned = {"mirbase": "23.1", "targetscan": "8.0",
              "mirtarbase": "live via Harmonizome (unpinned)", "enrichr": "live API"}
    report = ku.check_updates(client=_mock_client(_handler), local_versions=pinned)
    by_key = {s["key"]: s for s in report["sources"]}
    assert by_key["mirbase"]["status"] == "up-to-date"
    assert by_key["mirbase"]["local"] == "23.1"
    assert report["updates_available"] == 0


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


def test_make_async_client_constructs_with_broken_no_proxy():
    # Same trap, async twin (the analysis path uses AsyncClient).
    import asyncio

    async def _build():
        async with ku.make_async_client(timeout=5.0):
            pass

    asyncio.run(_build())


def test_load_local_versions(tmp_path):
    # Against an explicit file: the pinned record is whatever the file
    # says (the real file now records 23.1 after the live install).
    versions_file = tmp_path / "data_versions.json"
    versions_file.write_text(json.dumps({
        "mirbase": {"local": "22"},
        "targetscan": "8.0",
    }))
    versions = ku.load_local_versions(versions_file)
    assert versions["mirbase"] == "22"
    assert versions["targetscan"] == "8.0"
    # Missing file falls back to the registry defaults.
    assert ku.load_local_versions("/nonexistent/path.json")["mirbase"] == "22"


# ── API endpoint ─────────────────────────────────────────────────────────

def test_updates_endpoint(monkeypatch):
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    canned = ku.check_updates(
        client=_mock_client(_handler),
        local_versions={"mirbase": "22", "targetscan": "8.0",
                        "mirtarbase": "live via Harmonizome (unpinned)",
                        "enrichr": "live API"})
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
    # /health reports whatever data_versions.json pins (23.1 after the
    # live install); the contract is that it matches the pinned record.
    assert r.json()["data_sources"] == ku.load_local_versions()
    assert set(r.json()["data_sources"]) == {"mirbase", "targetscan",
                                             "mirtarbase", "enrichr", "kegg"}


# ── Download + install ─────────────────────────────────────────────────────

FASTA = ">hsa-miR-33a-5p MIMAT0000091\nGUGCAUUGUAGUUGCAUUGCA\n>hsa-miR-144-3p MIMAT0000436\nUACAGUAUAGAUGAUGUACU\n"


def _install_handler(request):
    url = str(request.url)
    if "mirbase.org/download/CURRENT/README" in url:
        return httpx.Response(200, text=MIRBASE_README)
    if url.endswith("mature.fa"):
        return httpx.Response(200, text=FASTA)
    if url.endswith("hairpin.fa"):
        return httpx.Response(200, text=FASTA)
    return httpx.Response(404, text="not found")


def _seed_versions(tmp_path):
    versions = {
        "mirbase": {"local": "22", "latest_known": "23.1"},
        "targetscan": {"local": "8.0", "latest_known": "8.0"},
    }
    d = tmp_path / "local_db"
    d.mkdir()
    (d / "data_versions.json").write_text(json.dumps(versions))
    return d


def test_download_file_verifies_and_is_atomic(tmp_path):
    dest = tmp_path / "mature.fa"
    client = _mock_client(_install_handler)
    path, sha, size = ku.download_file(
        "https://www.mirbase.org/download/CURRENT/mature.fa",
        dest, client, verify="fasta")
    assert path.read_text() == FASTA
    assert size == len(FASTA.encode())
    assert len(sha) == 64
    # A bad payload raises and leaves nothing at dest.
    bad = tmp_path / "bad.fa"
    with pytest.raises(ValueError):
        ku.download_file("https://example.com/bad.fa", bad,
                         _mock_client(lambda r: httpx.Response(200, text="nope")),
                         verify="fasta")
    assert not bad.exists()
    assert not (tmp_path / "bad.fa.tmp").exists()


def test_transform_html_fasta(tmp_path):
    # miRBase serves .fa files wrapped in HTML (<p>, <br>, &gt;). The real
    # 23.1 mature.fa unwraps to exactly 74,463 entries; here a small sample.
    wrapped = ("<p>&gt;cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p<br>"
               "UGAGGUAGUAGGUUGUAUAGUU<br>"
               "&gt;hsa-miR-33a-5p MIMAT0000091 Homo sapiens miR-33a-5p<br>"
               "GUGCAUUGUAGUUGCAUUGCA<br></p>\n")
    dest = tmp_path / "mature.fa"
    client = _mock_client(lambda r: httpx.Response(200, text=wrapped))
    path, sha, size = ku.download_file(
        "https://www.mirbase.org/download/CURRENT/mature.fa",
        dest, client, verify="fasta", transform="html_fasta")
    text = path.read_text()
    assert text.startswith(">cel-let-7-5p MIMAT0000001")
    assert ">hsa-miR-33a-5p MIMAT0000091" in text
    assert "<" not in text and "&gt;" not in text
    assert size == len(text.encode())
    # And the raw HTML-wrapped payload fails verification without the transform.
    with pytest.raises(ValueError):
        ku.download_file("https://www.mirbase.org/download/CURRENT/mature.fa",
                         tmp_path / "raw.fa", client, verify="fasta")


def test_install_source_update_end_to_end(tmp_path):
    d = _seed_versions(tmp_path)
    old = d / "mirbase_mature.fa"
    old.write_text(">old entry\nACGT\n")
    report = ku.install_source_update(
        "mirbase", client=_mock_client(_install_handler), base_dir=tmp_path)
    assert report["old_version"] == "22"
    assert report["new_version"] == "23.1"
    assert len(report["installed"]) == 2
    assert (d / "mirbase_mature.fa").read_text() == FASTA
    assert (d / "mirbase_hairpin.fa").read_text() == FASTA
    backups = [tmp_path / b for b in report["backups"]]
    assert len(backups) == 1  # only mature.fa existed before
    assert backups[0].read_text().startswith(">old entry")
    assert "bak-22" in backups[0].name
    recorded = json.loads((d / "data_versions.json").read_text())
    assert recorded["mirbase"]["local"] == "23.1"
    assert "updated_at" in recorded["mirbase"]
    assert "23.1" in report["notice"]


def test_install_refuses_without_files_or_unknown_source(tmp_path):
    # enrichr has no downloadable files (live API only).
    with pytest.raises(ValueError, match="no downloadable files"):
        ku.install_source_update("enrichr", client=_mock_client(_install_handler),
                                 base_dir=tmp_path)
    with pytest.raises(ValueError, match="unknown source"):
        ku.install_source_update("nope", client=_mock_client(_install_handler),
                                 base_dir=tmp_path)


def test_install_refuses_when_latest_cannot_be_verified(tmp_path):
    d = _seed_versions(tmp_path)
    live = d / "mirbase_mature.fa"
    live.write_text(">old entry\nACGT\n")

    def _down(request):
        raise httpx.ConnectError("no network", request=request)

    with pytest.raises(RuntimeError, match="could not verify"):
        ku.install_source_update("mirbase", client=_mock_client(_down),
                                 base_dir=tmp_path)
    # The live file is untouched.
    assert live.read_text().startswith(">old entry")


def test_install_endpoint(monkeypatch):
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    canned = {"source": "mirbase", "old_version": "22", "new_version": "23.1",
              "installed": [], "backups": [], "notice": "ok"}
    monkeypatch.setattr(ke, "install_source_update", lambda source: canned)
    # v1.45 (H-01): el endpoint exige KENRYU_ADMIN_TOKEN en X-Admin-Token.
    monkeypatch.setenv("KENRYU_ADMIN_TOKEN", "test-token")
    headers = {"X-Admin-Token": "test-token"}
    with TestClient(ke.app) as client:
        assert client.post("/api/v1/updates/mirbase/install").status_code == 403
        r = client.post("/api/v1/updates/mirbase/install", headers=headers)
    assert r.status_code == 200
    assert r.json()["new_version"] == "23.1"

    monkeypatch.setattr(ke, "install_source_update",
                        lambda source: (_ for _ in ()).throw(ValueError("unknown source: x")))
    with TestClient(ke.app) as client:
        r = client.post("/api/v1/updates/x/install", headers=headers)
    assert r.status_code == 400


def test_reverse_endpoint():
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    with TestClient(ke.app) as client:
        r = client.get("/api/v1/reverse", params={"gene": "abca1"})
    assert r.status_code == 200
    body = r.json()
    assert body["gene"] == "ABCA1"
    # DB keys are TargetScan family names; the miR-33 family is mir-33-5p.
    assert "mir-33-5p" in body["mirnas"]
    assert body["n"] == len(body["mirnas"]) > 0
    assert body["n"] == len(body["mirnas"])
    with TestClient(ke.app) as client:
        r = client.get("/api/v1/reverse", params={"gene": "ZZZZNOTAGENE"})
    assert r.json()["mirnas"] == []


def test_analyze_offline_end_to_end():
    """The full pipeline with offline=true makes NO network calls and
    still returns genes, validation, export, heatmap (local enrichment
    is empty for this panel since v1.49: the predicted universe correctly
    discounts its only GMT-union term)."""
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
            "mode": "strict", "years": 10, "offline": True})
    assert r.status_code == 200
    body = r.json()
    assert body["offline"] is True
    assert len(body["common_genes"]) == 6  # same core as the live run
    # v1.48 (Ronda 3): con el índice oficial de miRTarBase, NFIA tiene
    # validación experimental (n_val=2) y sube al primer puesto por
    # prioridad; con el GMT anterior el orden era alfabético.
    official = ke.BASE_DIR / "local_db" / "mirtarbase_official.json.gz"
    expected = (["NFIA", "GMFB", "KLF12", "PDE3A", "SKI", "SOX6"] if official.exists()
                else ["GMFB", "KLF12", "NFIA", "PDE3A", "SKI", "SOX6"])
    assert body["common_genes"] == expected
    # v1.49 (Ronda 3): con el universo predicho como fondo, el panel de
    # control ya no devuelve términos offline. Su único término (GO,
    # respuesta al TGF-beta, p-adj 0,0288) era un artefacto del fondo de
    # la unión de GMTs: las dianas del panel ya están enriquecidas en esa
    # ruta, así que el fondo correcto la descuenta. El enriquecimiento
    # local sigue funcionando — lo prueban los tests de test_enrichment.py
    # y el cálculo directo con el universo predicho.
    assert body["enrichment"] == []
    assert body["validation"]["summary"]["n_genes"] == 6
    assert body["validation"]["data_sources"]["mirbase"] == "23.1"
    assert body["export"]["schema_version"] == "1.0"
    assert len(body["export"]["edges"]) > 0
    assert body["heatmap_plot"]
    assert body["venn_plot"]
    # Volcano and PPI need live data (PubMed counts, STRING): offline they
    # are honestly absent instead of simulated.
    assert body["volcano_plot"] is None
    assert body["ppi_plot"] is None
    # No live research in offline mode.
    for res in body["gene_research"].values():
        assert res["year_window_used"] == "offline"
        assert res["pubmed"] == []


# ── Ronda 3 (v1.48) — KEGG como descarga en ejecución, nunca en el repo ─────
KEGG_GMT = ("ABC transporters\t\tABCA1\tABCB1\n"
            "Glycolysis / Gluconeogenesis\t\tGAPDH\tHK1\n")


def _kegg_handler(request):
    url = str(request.url)
    if "datasetStatistics" in url:
        return httpx.Response(200, text=ENRICHR_STATS)
    if "geneSetLibrary" in url and "KEGG_2021_Human" in url:
        return httpx.Response(200, text=KEGG_GMT)
    return httpx.Response(404, text="not found")


def test_verify_gmt(tmp_path):
    good = tmp_path / "good.gmt"
    good.write_text(KEGG_GMT)
    assert ku.verify_gmt(good) is True
    bad = tmp_path / "bad.gmt"
    bad.write_text("not a gmt file\n")
    with pytest.raises(ValueError):
        ku.verify_gmt(bad)


def test_kegg_source_registered():
    entry = next(s for s in ku.DATA_SOURCES if s["key"] == "kegg")
    assert entry["files"], "KEGG must have a downloadable file"
    spec = entry["files"][0]
    assert "KEGG_2021_Human" in spec["url"]
    assert spec["local"] == "local_db/gmt_KEGG_2021_Human.gmt"
    assert spec["verify"] == "gmt"
    assert "never bundled" in entry["note"].lower()


def test_install_kegg_end_to_end(tmp_path):
    d = _seed_versions(tmp_path)
    report = ku.install_source_update(
        "kegg", client=_mock_client(_kegg_handler), base_dir=tmp_path)
    installed = d / "gmt_KEGG_2021_Human.gmt"
    assert installed.read_text() == KEGG_GMT
    assert report["new_version"].startswith("snapshot-")
    recorded = json.loads((d / "data_versions.json").read_text())
    assert recorded["kegg"]["local"].startswith("snapshot-")
    assert len(report["installed"]) == 1


# ── Ronda 3 (v1.49) — H-20: descargas robustas ─────────────────────────────
# Método: pruebas escritas primero; download_file no reintentaba 429 ni
# timeouts, y transform_html_fasta escribía \r\n en Windows.


def test_transform_html_fasta_writes_lf_bytes(tmp_path, monkeypatch):
    """H-20: simula la traducción de saltos de línea de Windows
    (write_text sin newline convierte \n en \r\n: 161 vs 157 bytes).
    El transform debe pasar newline="\n" para que el archivo quede con
    bytes exactos y size == len(text.encode())."""
    import pathlib
    real_write_text = pathlib.Path.write_text

    def _windows_write_text(self, data, encoding=None, errors=None, newline=None):
        if newline is None:
            data = data.replace("\n", "\r\n")
        return real_write_text(self, data, encoding=encoding, errors=errors, newline=newline)

    monkeypatch.setattr(pathlib.Path, "write_text", _windows_write_text)
    dest = tmp_path / "mature.fa"
    dest.write_text("<p>&gt;hsa-miR-33a-5p MIMAT0000091<br>GUGCAUUGUAGUUGCAUUGCA<br></p>\n",
                    encoding="utf-8")
    assert ku.transform_html_fasta(dest) is True
    data = dest.read_bytes()
    assert b"\r\n" not in data
    assert data.startswith(b">hsa-miR-33a-5p MIMAT0000091")


def test_download_file_retries_on_429(tmp_path, monkeypatch):
    monkeypatch.setattr(ku, "_DOWNLOAD_RETRY_DELAYS", (0.0, 0.0))
    calls = {"n": 0}

    def _handler(request):
        calls["n"] += 1
        if calls["n"] <= 2:
            return httpx.Response(429, text="slow down")
        return httpx.Response(200, content=b">x\nACGT\n")

    client = _mock_client(_handler)
    path, sha, size = ku.download_file(
        "https://example.org/x.fa", tmp_path / "x.fa", client, verify="fasta")
    assert calls["n"] == 3
    assert size == 8
    assert path.read_bytes() == b">x\nACGT\n"


def test_download_file_retries_on_timeout_then_raises(tmp_path, monkeypatch):
    monkeypatch.setattr(ku, "_DOWNLOAD_RETRY_DELAYS", (0.0, 0.0))
    calls = {"n": 0}

    def _handler(request):
        calls["n"] += 1
        raise httpx.ConnectTimeout("timed out")

    client = _mock_client(_handler)
    dest = tmp_path / "x.fa"
    with pytest.raises(httpx.ConnectTimeout):
        ku.download_file("https://example.org/x.fa", dest, client)
    assert calls["n"] == 3
    assert not dest.exists()
    assert not dest.with_suffix(".fa.tmp").exists()


def test_download_file_no_retry_on_404(tmp_path, monkeypatch):
    monkeypatch.setattr(ku, "_DOWNLOAD_RETRY_DELAYS", (0.0, 0.0))
    calls = {"n": 0}

    def _handler(request):
        calls["n"] += 1
        return httpx.Response(404, text="not found")

    client = _mock_client(_handler)
    with pytest.raises(httpx.HTTPStatusError):
        ku.download_file("https://example.org/x.fa", tmp_path / "x.fa", client)
    assert calls["n"] == 1


def test_download_file_rejects_http(tmp_path):
    calls = {"n": 0}

    def _handler(request):
        calls["n"] += 1
        return httpx.Response(200, content=b">x\nACGT\n")

    client = _mock_client(_handler)
    with pytest.raises(ValueError, match="https"):
        ku.download_file("http://example.org/x.fa", tmp_path / "x.fa", client)
    assert calls["n"] == 0


def test_verify_fasta_single_entry(tmp_path):
    # H-20 (cobertura): el verificador exigía DOS cabeceras '>' — un
    # FASTA válido de una sola entrada fallaba con "no entries".
    single = tmp_path / "one.fa"
    single.write_text(">hsa-miR-33a-5p MIMAT0000091\nGUGCAUUGUAGUUGCAUUGCA\n")
    assert ku.verify_fasta(single) is True
    empty = tmp_path / "empty.fa"
    empty.write_text("")
    with pytest.raises(ValueError):
        ku.verify_fasta(empty)
    not_fasta = tmp_path / "nope.fa"
    not_fasta.write_text("ACGT\n")
    with pytest.raises(ValueError):
        ku.verify_fasta(not_fasta)
