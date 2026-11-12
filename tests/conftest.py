"""Test setup: run from the repository root so relative paths
(static/, local_db/) resolve exactly as in production."""
import os
import sys
import types
from pathlib import Path

try:
    import gseapy  # noqa: F401
except ImportError:  # pinned gseapy==1.1.1 needs Rust to build; engine never uses it
    sys.modules["gseapy"] = types.ModuleType("gseapy")

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
os.chdir(REPO_ROOT)
