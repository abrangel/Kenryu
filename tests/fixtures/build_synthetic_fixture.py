"""Build the synthetic miRTarBase official index used by the CI suite.

The real index (local_db/mirtarbase_official.json.gz, 2.3 MB, built from
hsa_MTI.xlsx) is gitignored and never exists in a fresh clone, so the
tests that need the official format skipped in CI. This fixture has the
SAME schema (payload with version/source/built_from/built_at/rows_read/
pairs/mirnas; each gene entry has level/pmids/experiments) but only the
7 pairs of the thesis panel, with values COPIED VERBATIM from the real
index (2026-11-19) — no invented PMIDs. It is labelled synthetic in its
own metadata, and the engine reports it as such in /health.

Regenerate: python tests/fixtures/build_synthetic_fixture.py
"""
import gzip
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent

# Values copied verbatim from the real miRTarBase 9.0 index.
MIRNAS = {
    "hsa-miR-33a-5p": {
        "ABCA1": {"level": "strong", "pmids": [20466882, 20466885, 20566875],
                  "experiments": ["Luciferase reporter assay"]},
        "SNTB2": {"level": "weak", "pmids": [21572407, 22012620],
                  "experiments": ["PAR-CLIP"]},
    },
    "hsa-miR-33b-5p": {
        "ABCA1": {"level": "strong", "pmids": [20466882],
                  "experiments": ["Luciferase reporter assay"]},
        "SNTB2": {"level": "weak", "pmids": [21572407, 22012620],
                  "experiments": ["PAR-CLIP"]},
    },
    "hsa-miR-144-3p": {
        "TSC22D2": {"level": "weak", "pmids": [24906430],
                    "experiments": ["HITS-CLIP"]},
    },
    "hsa-miR-758-3p": {
        "TSC22D2": {"level": "weak", "pmids": [23592263],
                    "experiments": ["PAR-CLIP"]},
    },
    "hsa-miR-106b-5p": {
        "ABCA1": {"level": "weak", "pmids": [22473208],
                  "experiments": ["HITS-CLIP"]},
    },
}

PAIRS = sum(len(v) for v in MIRNAS.values())


def main():
    payload = {
        "version": "miRTarBase 9.0",
        "source": ("synthetic fixture for CI: same schema as "
                   "build_mirtarbase_official.py output; the 7 thesis-panel "
                   "pairs copied verbatim from the real index"),
        "built_from": "synthetic",
        "built_at": "2026-11-19",
        "rows_read": PAIRS,
        "pairs": PAIRS,
        "synthetic": True,
        "mirnas": MIRNAS,
    }
    out = HERE / "mirtarbase_official_synthetic.json.gz"
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, ensure_ascii=False)
    print(f"wrote {out} ({PAIRS} pairs, {len(MIRNAS)} miRNAs)")


if __name__ == "__main__":
    main()
