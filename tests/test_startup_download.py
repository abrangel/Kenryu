"""Startup download of the Harmonizome GMT (v1.53, bloque 2 / licencia).

The CC BY-NC-SA snapshot is no longer distributed in the repo (v1.52) nor
in the deployed image: at startup, when neither the official miRTarBase
index nor the GMT is on disk, the service tries to fetch the GMT from
Harmonizome (best-effort, never blocks startup). Tests use a fake
downloader and tmp paths; the real network is never touched here.
"""
import asyncio
import gzip

import kenryu_engine as ke


def _fake_gmt_bytes():
    return gzip.compress(b"hsa-miR-21-5p\tMIRT000001\tPTEN\tPDCD4\n")


def test_startup_downloads_legacy_gmt_when_missing(tmp_path, monkeypatch):
    monkeypatch.setenv("KENRYU_STARTUP_GMT_DOWNLOAD", "1")
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(tmp_path / "official.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(tmp_path / "mirtarbase_gene_set_library.gmt.gz"))

    async def fake_download(url):
        assert url.startswith("https://maayanlab.cloud/")
        return _fake_gmt_bytes()

    monkeypatch.setattr(ke, "_download_gmt_bytes", fake_download)
    assert asyncio.run(ke._ensure_legacy_gmt()) is True
    assert (tmp_path / "mirtarbase_gene_set_library.gmt.gz").exists()


def test_startup_download_failure_never_blocks_startup(tmp_path, monkeypatch):
    monkeypatch.setenv("KENRYU_STARTUP_GMT_DOWNLOAD", "1")
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(tmp_path / "official.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(tmp_path / "mirtarbase_gene_set_library.gmt.gz"))

    async def failing_download(url):
        raise ConnectionError("Harmonizome down")

    monkeypatch.setattr(ke, "_download_gmt_bytes", failing_download)
    assert asyncio.run(ke._ensure_legacy_gmt()) is False
    assert not (tmp_path / "mirtarbase_gene_set_library.gmt.gz").exists()


def test_startup_download_skipped_when_official_present(tmp_path, monkeypatch):
    monkeypatch.setenv("KENRYU_STARTUP_GMT_DOWNLOAD", "1")
    (tmp_path / "official.json.gz").write_bytes(b"fake-official")
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(tmp_path / "official.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(tmp_path / "mirtarbase_gene_set_library.gmt.gz"))

    async def must_not_be_called(url):
        raise AssertionError("download must not run when the official index exists")

    monkeypatch.setattr(ke, "_download_gmt_bytes", must_not_be_called)
    assert asyncio.run(ke._ensure_legacy_gmt()) is True
    assert not (tmp_path / "mirtarbase_gene_set_library.gmt.gz").exists()


def test_startup_download_disabled_by_env_flag(tmp_path, monkeypatch):
    monkeypatch.setenv("KENRYU_STARTUP_GMT_DOWNLOAD", "0")
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(tmp_path / "official.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(tmp_path / "mirtarbase_gene_set_library.gmt.gz"))

    async def must_not_be_called(url):
        raise AssertionError("download must not run when disabled by env flag")

    monkeypatch.setattr(ke, "_download_gmt_bytes", must_not_be_called)
    assert asyncio.run(ke._ensure_legacy_gmt()) is False
