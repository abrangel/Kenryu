"""H-15 (Ronda 3): the persistent cache must be atomic, bounded and
fresh. A crash mid-write must never leave a partial analysis_cache.json;
each section is capped (oldest entries dropped first); gene-research
entries expire after CACHE_TTL_DAYS (translations and PubMed metadata
are stable facts and do not expire).

conftest redirects ke.CACHE_FILE to a temp dir and starts the in-memory
cache empty, so these tests never touch the real cache file.
"""
import copy
import datetime
import json

import pytest

import kenryu_engine as ke


@pytest.fixture
def _snapshot():
    snap = copy.deepcopy(ke.PERSISTENT_CACHE)
    yield
    ke.PERSISTENT_CACHE.clear()
    ke.PERSISTENT_CACHE.update(snap)


def test_save_cache_is_atomic_on_failure(_snapshot, monkeypatch):
    ke.PERSISTENT_CACHE["trans"] = {"hola": "hello"}
    ke.save_persistent_cache()
    assert json.loads(ke.CACHE_FILE.read_text())["trans"] == {"hola": "hello"}

    def _boom(*args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(json, "dump", _boom)
    ke.save_persistent_cache()  # must not raise, must not corrupt
    assert json.loads(ke.CACHE_FILE.read_text())["trans"] == {"hola": "hello"}
    assert not ke.CACHE_FILE.with_suffix(".tmp").exists()


def test_save_cache_enforces_limits(_snapshot):
    ke.PERSISTENT_CACHE["trans"] = {f"k{i}": f"v{i}" for i in range(3000)}
    ke.save_persistent_cache()
    data = json.loads(ke.CACHE_FILE.read_text())
    assert len(data["trans"]) <= ke._CACHE_LIMITS["trans"]
    # The newest entries survive; the oldest are dropped first.
    assert f"k{2999}" in data["trans"]
    assert "k0" not in data["trans"]


def test_load_cache_drops_expired_gene_research(_snapshot):
    old = (datetime.datetime.now(datetime.timezone.utc)
           - datetime.timedelta(days=ke.CACHE_TTL_DAYS + 10))
    payload = {
        "_saved_at": old.isoformat(),
        "trans": {"hola": "hello"},
        "pubmed": {},
        "enrichr": {},
        "gene_research": {"v5|ABCA1|10": {"pubmed": []}},
    }
    ke.CACHE_FILE.write_text(json.dumps(payload))
    ke.PERSISTENT_CACHE["gene_research"] = {}
    ke.load_persistent_cache()
    assert ke.PERSISTENT_CACHE["trans"].get("hola") == "hello"
    assert ke.PERSISTENT_CACHE["gene_research"] == {}


def test_load_cache_keeps_fresh_gene_research(_snapshot):
    payload = {
        "_saved_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "trans": {},
        "pubmed": {},
        "enrichr": {},
        "gene_research": {"v5|ABCA1|10": {"pubmed": []}},
    }
    ke.CACHE_FILE.write_text(json.dumps(payload))
    ke.PERSISTENT_CACHE["gene_research"] = {}
    ke.load_persistent_cache()
    assert "v5|ABCA1|10" in ke.PERSISTENT_CACHE["gene_research"]
