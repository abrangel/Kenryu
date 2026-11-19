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

# v1.45 (auditoría externa, H-15): a regression test once overwrote the
# committed `local_db/analysis_cache.json` with a single TESTGENE entry.
# Tests must never read from or write to the real cache file: the cache
# path is redirected to a temp dir and the in-memory cache starts empty.
import kenryu_engine as ke

_TMP = Path(tempfile.mkdtemp(prefix="kenryu-test-cache-"))


def pytest_configure():
    ke.CACHE_FILE = _TMP / "analysis_cache.json"
    ke.PERSISTENT_CACHE = {"trans": {}, "pubmed": {}, "enrichr": {}, "gene_research": {}}
