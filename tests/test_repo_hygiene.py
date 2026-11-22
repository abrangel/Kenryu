"""Repo-level hygiene guards (H-02, H-09, H-18) — the three regression
gaps found by the 2026-11-22 audit review.

These tests guard things a unit test cannot: that a licensed data file
is never committed again, that the dependency scan runs, and that the
fixed "confidence" literals never return. They run in the normal
suite on all three OSes; pip-audit itself is a separate CI job
(see .github/workflows/ci.yml) because it needs the network.
"""
import re
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# H-02: files that must never be tracked in git. KEGG and the
# Harmonizome snapshot are third-party data under licences that forbid
# redistribution in an MIT repo; the derived miRTarBase indexes are
# runtime data, built locally, never committed (step A of the v10.0
# adoption). The GMT allowlist is the three redistributable libraries.
FORBIDDEN_TRACKED = [
    re.compile(r"local_db/gmt_KEGG", re.IGNORECASE),
    re.compile(r"mirtarbase_gene_set_library", re.IGNORECASE),
    re.compile(r"local_db/mirtarbase_(official|v10)\.json\.gz$"),
]
ALLOWED_GMTS = {
    "local_db/gmt_GO_Biological_Process_2023.gmt",
    "local_db/gmt_Reactome_Pathways_2024.gmt",
    "local_db/gmt_WikiPathways_2019_Human.gmt",
}

# H-18: the fixed "confidence" values (85/88/92/95) that used to sit
# in gene_details and in the UI's "Rigor Científico" bar fallback.
# Any line that mentions confidence next to one of those literals is
# a regression; the numbers may appear elsewhere (ports, counts),
# but never next to "confidence".
CONFIDENCE_RE = re.compile(
    r"(?i)confid[^\n]*\b(85|88|92|95)\b|\b(85|88|92|95)\b[^\n]*confid"
)
CONFIDENCE_FILES = [
    REPO_ROOT / "kenryu_engine.py",
    REPO_ROOT / "kenryu_core.py",
    REPO_ROOT / "static" / "script.js",
]


def _tracked_files():
    out = subprocess.run(
        ["git", "ls-files"], cwd=REPO_ROOT, capture_output=True,
        text=True, check=True,
    )
    return out.stdout.splitlines()


def test_no_forbidden_data_file_is_tracked():
    tracked = _tracked_files()
    bad = [f for f in tracked
           if any(p.search(f) for p in FORBIDDEN_TRACKED)]
    assert not bad, (
        "forbidden data file(s) tracked in git (H-02): " + ", ".join(bad)
    )


def test_only_allowlisted_gmts_are_tracked():
    tracked = _tracked_files()
    gmts = [f for f in tracked if f.endswith(".gmt")]
    bad = [f for f in gmts if f not in ALLOWED_GMTS]
    assert not bad, (
        "tracked GMT outside the allowlist (H-02): " + ", ".join(bad)
        + " — add it to ALLOWED_GMTS only with a licence check"
    )


def test_no_fixed_confidence_literals():
    bad = []
    for path in CONFIDENCE_FILES:
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if CONFIDENCE_RE.search(line):
                bad.append(f"{path.name}:{i}: {line.strip()[:100]}")
    assert not bad, (
        "fixed confidence literal(s) back in the code (H-18):\n"
        + "\n".join(bad)
    )
