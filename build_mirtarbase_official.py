#!/usr/bin/env python3
"""Build the local miRTarBase official index from the downloaded hsa_MTI.xlsx.

miRTarBase 9.0 (hsa_MTI.xlsx, 957.040 rows) -> local_db/mirtarbase_official.json.gz

The output is a COMPACT DERIVED index, not the raw file:
    {"version": "miRTarBase 9.0", "source": ..., "mirnas": {
        "hsa-miR-33a-5p": {
            "ABCA1": {"level": "strong",
                      "pmids": [20466882, 20466885, 20566875],
                      "experiments": ["Luciferase reporter assay"]},
            "SNTB2": {"level": "weak",
                      "pmids": [21572407, 22012620],
                      "experiments": ["PAR-CLIP"]}}}}

Level per (miRNA, gene): "strong" if any row has Support Type
"Functional MTI"; "weak" if any has "Functional MTI (Weak)"; else "none"
(only "Non-Functional MTI" rows, or rows with no Support Type).
The xlsx repeats each interaction twice (once with a Support Type, once
with None and "//Other" appended to Experiments); this script deduplicates
by (miRNA, gene, PMID, experiment) and only takes experiment names from
rows that carry a Support Type.

The derived file is RUNTIME DATA, not source: it is NOT committed to git
(see .gitignore). To reproduce it, download hsa_MTI.xlsx from the official
miRTarBase page (https://mirtarBase.cuhk.edu.cn/) and run:
    python3 build_mirtarbase_official.py /path/to/hsa_MTI.xlsx
The redistribution terms of miRTarBase could not be verified: the official
site returned 404 on every page when checked on 2026-11-19 (see
DATA_LICENSES.md). The raw xlsx is never committed; only this derived
index, built locally from the user's own download, is used by the engine.

Usage: python3 build_mirtarbase_official.py [xlsx_path] [output_path]
Defaults: the browser-downloaded copy in ~/workspace/browser_downloads/
and local_db/mirtarbase_official.json.gz next to this script.
"""

import gzip
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

STRONG = "Functional MTI"
WEAK = "Functional MTI (Weak)"

DEFAULT_XLSX = sorted(
    Path.home().glob("workspace/browser_downloads/*/browser-download-*-hsa_MTI.xlsx")
)
DEFAULT_XLSX = str(DEFAULT_XLSX[-1]) if DEFAULT_XLSX else "hsa_MTI.xlsx"


# Gene symbols Excel mangled into dates when the xlsx was generated
# (SEPT1 -> 1-Sep, MARCH1 -> 1-Mar, DEC1 -> 1-Dec). The old Harmonizome
# snapshot uses the same short forms, so recovering them keeps the
# matching identical. Any other non-string value is skipped and counted.
_MANGLED_MONTH = {3: "MARCH", 9: "SEPT", 12: "DEC"}


def _clean_gene(value):
    if isinstance(value, str):
        return value.strip() or None
    if hasattr(value, "month") and hasattr(value, "day"):
        prefix = _MANGLED_MONTH.get(value.month)
        if prefix:
            return f"{prefix}{value.day}"
    return None


def build(xlsx_path, out_path):
    import openpyxl

    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = ws.iter_rows(values_only=True)
    header = [str(c) for c in next(rows)]
    i_mir = header.index("miRNA")
    i_gene = header.index("Target Gene")
    i_exp = header.index("Experiments")
    i_sup = header.index("Support Type")
    i_pmid = header.index("References (PMID)")

    # (mirna, gene) -> {"level": str, "pmids": set, "experiments": set}
    index = {}
    n_rows = 0
    n_recovered = 0
    n_skipped = 0
    level_counts = Counter()
    for r in rows:
        n_rows += 1
        mirna = r[i_mir].strip() if isinstance(r[i_mir], str) else None
        gene = _clean_gene(r[i_gene])
        if gene and hasattr(r[i_gene], "month"):
            n_recovered += 1
        if not mirna or not gene:
            n_skipped += 1
            continue
        support = (r[i_sup] or "").strip()
        pmid = r[i_pmid]
        experiments = (r[i_exp] or "").strip()
        entry = index.setdefault(mirna, {}).setdefault(
            gene, {"level": "none", "pmids": set(), "experiments": set()})
        if support == STRONG:
            entry["level"] = "strong"
        elif support == WEAK and entry["level"] != "strong":
            entry["level"] = "weak"
        if isinstance(pmid, int):
            entry["pmids"].add(pmid)
        elif isinstance(pmid, str) and pmid.strip().isdigit():
            entry["pmids"].add(int(pmid.strip()))
        # Only rows with a Support Type carry clean experiment names;
        # the duplicated None rows append "//Other". Experiments are
        # separated by "//" or ";" in the source file.
        if support:
            for exp in experiments.replace(";", "//").split("//"):
                exp = exp.strip()
                if exp and exp != "Other":
                    entry["experiments"].add(exp)

    mirnas = {}
    n_pairs = 0
    for mirna, genes in index.items():
        mirnas[mirna] = {}
        for gene, ev in genes.items():
            n_pairs += 1
            level_counts[ev["level"]] += 1
            mirnas[mirna][gene] = {
                "level": ev["level"],
                "pmids": sorted(ev["pmids"]),
                "experiments": sorted(ev["experiments"]),
            }

    payload = {
        "version": "miRTarBase 9.0",
        "source": "https://mirtarbase.cuhk.edu.cn/",
        "built_from": "hsa_MTI.xlsx",
        "built_at": datetime.now(timezone.utc).isoformat(),
        "rows_read": n_rows,
        "pairs": n_pairs,
        "mirnas": mirnas,
    }
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"rows: {n_rows}  pairs: {n_pairs}  miRNAs: {len(mirnas)}")
    print(f"levels: {dict(level_counts)}")
    print(f"recovered Excel-mangled gene symbols: {n_recovered}  skipped rows: {n_skipped}")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
    return payload


if __name__ == "__main__":
    xlsx = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX
    out = (sys.argv[2] if len(sys.argv) > 2 else
           str(Path(__file__).parent / "local_db" / "mirtarbase_official.json.gz"))
    build(xlsx, out)
