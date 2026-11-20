#!/usr/bin/env python3
"""Generate the fixed 10 panels for the tool-comparison protocol.

Seven random panels (R1–R7, fixed seeds 101–107, sizes 2, 3, 5, 5, 5,
8, 10) plus the three reference panels (lipid, cancer, control).
One miRNA per TargetScan family, so each panel has as many effective
seed families as names. The panel list is committed before any tool
is run — nobody picks panels after seeing results.

Run from the repo root:  python scripts/generate_comparison_panels.py
Output: docs/validation/comparison-panels.json
"""
import json
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import kenryu_engine as ke  # noqa: E402

SIZES = [2, 3, 5, 5, 5, 8, 10]
SEEDS = [101, 102, 103, 104, 105, 106, 107]

REFERENCE = {
    "lipid": ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
              "hsa-miR-106b-5p", "hsa-miR-33b-5p"],
    "cancer": ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"],
    "control": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
}


def mature_names():
    """Human mature miRNA names from the miRBase fasta in local_db."""
    names = []
    for line in (ROOT / "local_db" / "mirbase_mature.fa").read_text().splitlines():
        if line.startswith(">hsa-"):
            names.append(line[1:].split()[0])
    return sorted(set(names))


def main():
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.TARGETSCAN_DB
    # Group miRNA names by their exact target set (one group = one family).
    by_set = {}
    for name in mature_names():
        genes = ke.get_targets_local(name, db)
        if genes:
            by_set.setdefault(frozenset(genes), []).append(name)
    sets = sorted(by_set, key=lambda s: sorted(by_set[s])[0])
    print(f"{len(by_set)} distinct target sets from {len(db)} families")
    out = {"reference": REFERENCE, "random": []}
    for i, (seed, size) in enumerate(zip(SEEDS, SIZES), start=1):
        rng = random.Random(seed)
        chosen = rng.sample(sets, size)
        mirnas = [sorted(by_set[s])[0] for s in chosen]
        assert len(set(mirnas)) == size
        out["random"].append({"id": f"R{i}", "seed": seed, "size": size,
                              "mirnas": sorted(mirnas)})
        print(f"R{i} (seed {seed}, n={size}): {', '.join(sorted(mirnas))}")
    path = ROOT / "docs" / "validation" / "comparison-panels.json"
    path.write_text(json.dumps(out, indent=2) + "\n", encoding="utf-8")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
