"""Test setup: run from the repository root so relative paths
(static/, local_db/) resolve exactly as in production."""
import os
import sys
import tempfile
import types
from pathlib import Path

try:
    import gseapy  # noqa: F401
except ImportError:  # pinned gseapy==1.1.1 needs Rust to build; engine never uses it
    sys.modules["gseapy"] = types.ModuleType("gseapy")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
os.chdir(REPO_ROOT)

# v1.50 (auditoría externa, H-10): el limitador de tasa está desactivado
# dentro de la suite (KENRYU_RATE_LIMIT=0) para que los cientos de
# peticiones de los tests nunca den 429. Los tests de tests/test_rate_limit.py
# lo activan explícitamente con monkeypatch y resetean su estado.
os.environ.setdefault("KENRYU_RATE_LIMIT", "0")

# v1.53: la descarga del GMT al arrancar está desactivada dentro de la
# suite (KENRYU_STARTUP_GMT_DOWNLOAD=0) para que ningún test toque la
# red; tests/test_startup_download.py la activa explícitamente con
# monkeypatch y un descargador falso.
os.environ.setdefault("KENRYU_STARTUP_GMT_DOWNLOAD", "0")

# v1.45 (auditoría externa, H-15): a regression test once overwrote the
# committed `local_db/analysis_cache.json` with a single TESTGENE entry.
# Tests must never read from or write to the real cache file: the cache
# path is redirected to a temp dir and the in-memory cache starts empty.
import kenryu_engine as ke

_TMP = Path(tempfile.mkdtemp(prefix="kenryu-test-cache-"))


def pytest_configure():
    ke.CACHE_FILE = _TMP / "analysis_cache.json"
    ke.PERSISTENT_CACHE = {"trans": {}, "pubmed": {}, "enrichr": {}, "gene_research": {}}


# Unreleased: the real official miRTarBase index (2.3 MB) is gitignored
# and never exists in a fresh clone, so the tests that need the official
# format skipped in CI. This fixture points the engine at the committed
# synthetic index (tests/fixtures/mirtarbase_official_synthetic.json.gz:
# same schema, the 7 thesis-panel pairs copied verbatim from the real
# index) so the official-format code path runs in CI too.
import pytest as _pytest


@_pytest.fixture
def synthetic_official_index(monkeypatch, tmp_path):
    fixture = REPO_ROOT / "tests" / "fixtures" / "mirtarbase_official_synthetic.json.gz"
    monkeypatch.setattr(ke, "MIRTARBASE_OFFICIAL", str(fixture))
    monkeypatch.setattr(ke, "MIRTARBASE_V10", str(tmp_path / "nope-v10.json.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_GMT", str(tmp_path / "nope.gmt.gz"))
    monkeypatch.setattr(ke, "MIRTARBASE_DB", {})
    monkeypatch.setattr(ke, "MIRTARBASE_FORMAT", None)
    monkeypatch.setattr(ke, "MIRTARBASE_META", None)
    yield fixture
