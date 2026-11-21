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
"Functional MTI" (and at least one of those rows carries a PMID);
"weak" if any has "Functional MTI (Weak)" (and at least one of those
rows carries a PMID); "clip_only" (v2.0) if the weak rows carry NO PMID
and every experiment is a CLIP assay — visible in the dossier, +0
points, never unlocks "alta"; else "none" (only "Non-Functional MTI"
rows, or rows with no Support Type). Fail-first (v2.0): the build
REFUSES to guess — a Functional MTI row without PMID, a weak row
without PMID whose experiments are not all CLIP, or a Support Type
the builder does not know all raise ValueError instead of being
silently mapped. (Measured 2026-11-21: v9.0 has 0 weak rows without
PMID and 0 strong rows without PMID; v10.0 has 1,411,229 weak rows
without PMID, all CLIP-only, and 0 strong rows without PMID.)
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
NON_FUNCTIONAL = "Non-Functional MTI"
NON_FUNCTIONAL_WEAK = "Non-Functional MTI (Weak)"
KNOWN_SUPPORT = {STRONG, WEAK, NON_FUNCTIONAL, NON_FUNCTIONAL_WEAK, ""}


def _is_clip_only(experiments):
    """True when every experiment is a CLIP assay (HITS-CLIP, PAR-CLIP,
    CLIP-seq, ...). An empty experiment list is NOT clip-only: without
    any assay named, the category is unknown and the build must fail."""
    return bool(experiments) and all("CLIP" in e.upper() for e in experiments)

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


def _read_rows(path):
    """Yield data rows (as dicts keyed by the v9.0/v10.0 column names).

    Accepts the v9.0 xlsx and the v10.0 CSV. The CSV has a UTF-8 BOM
    (encoding="utf-8-sig") and formats Entrez ID / PMID with ".0"
    (e.g. 17179747.0); both are normalized here so the rest of the
    builder sees the same values it saw from the xlsx.
    """
    if str(path).lower().endswith(".csv"):
        import csv as _csv
        with open(path, encoding="utf-8-sig", newline="") as fh:
            for row in _csv.DictReader(fh):
                yield {k: _norm_num(v) for k, v in row.items()}
        return
    import openpyxl
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    ws = wb[wb.sheetnames[0]]
    rows = ws.iter_rows(values_only=True)
    header = [str(c) for c in next(rows)]
    for r in rows:
        yield dict(zip(header, r))


def _norm_num(v):
    # v10.0 CSV writes numbers as "17179747.0"; the xlsx gives ints.
    if isinstance(v, str) and v.endswith(".0") and v[:-2].isdigit():
        return int(v[:-2])
    return v


def build(xlsx_path, out_path, version="miRTarBase 9.0"):
    # (mirna, gene) -> {"has_strong", "has_weak", "strong_pmids",
    # "weak_pmids", "pmids", "experiments"}. The level is assigned at
    # the END, from the per-support PMID sets: a "strong" claim needs a
    # PMID on a Functional MTI row, a "weak" claim needs one on a
    # Functional MTI (Weak) row. PMIDs from another support level never
    # rescue a level (fail-first, v2.0).
    index = {}
    n_rows = 0
    n_recovered = 0
    n_skipped = 0
    n_nonhuman = 0
    level_counts = Counter()
    for r in _read_rows(xlsx_path):
        n_rows += 1
        # Human only: v9.0 says "Homo sapiens", v10.0 says "hsa".
        # The v9.0 file is human-only, so this is a no-op there.
        sp_mir = (r.get("Species (miRNA)") or "").strip()
        sp_gene = (r.get("Species (Target Gene)") or "").strip()
        if sp_mir or sp_gene:
            if sp_mir not in ("Homo sapiens", "hsa") or \
                    sp_gene not in ("Homo sapiens", "hsa"):
                n_nonhuman += 1
                continue
        mirna = r.get("miRNA")
        mirna = mirna.strip() if isinstance(mirna, str) else None
        gene = _clean_gene(r.get("Target Gene"))
        if gene and hasattr(r.get("Target Gene"), "month"):
            n_recovered += 1
        if not mirna or not gene:
            n_skipped += 1
            continue
        support = r.get("Support Type")
        support = support.strip() if isinstance(support, str) else ""
        if support not in KNOWN_SUPPORT:
            raise ValueError(
                f"Support Type desconocido: {support!r} en {mirna} -> {gene} "
                f"(fila {n_rows} de {xlsx_path}); el constructor no adivina "
                f"niveles — añade el tipo a KNOWN_SUPPORT o corrige el archivo.")
        pmid = r.get("References (PMID)")
        experiments = r.get("Experiments")
        experiments = experiments.strip() if isinstance(experiments, str) else ""
        entry = index.setdefault(mirna, {}).setdefault(
            gene, {"has_strong": False, "has_weak": False,
                   "strong_pmids": set(), "weak_pmids": set(),
                   "pmids": set(), "experiments": set()})
        pmid_int = None
        if isinstance(pmid, int):
            pmid_int = pmid
        elif isinstance(pmid, str) and pmid.strip().isdigit():
            pmid_int = int(pmid.strip())
        if support == STRONG:
            entry["has_strong"] = True
            if pmid_int is not None:
                entry["strong_pmids"].add(pmid_int)
        elif support == WEAK:
            entry["has_weak"] = True
            if pmid_int is not None:
                entry["weak_pmids"].add(pmid_int)
        if pmid_int is not None:
            entry["pmids"].add(pmid_int)
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
            if ev["has_strong"]:
                if not ev["strong_pmids"]:
                    raise ValueError(
                        f"{mirna} -> {gene}: Functional MTI sin PMID en "
                        f"ninguna fila; no se puede llamar 'strong'.")
                level = "strong"
            elif ev["has_weak"]:
                if ev["weak_pmids"]:
                    level = "weak"
                elif _is_clip_only(ev["experiments"]):
                    level = "clip_only"
                else:
                    raise ValueError(
                        f"{mirna} -> {gene}: Functional MTI (Weak) sin PMID "
                        f"y experimentos no todos CLIP "
                        f"({sorted(ev['experiments'])}); categoría desconocida.")
            else:
                level = "none"
            level_counts[level] += 1
            mirnas[mirna][gene] = {
                "level": level,
                "pmids": sorted(ev["pmids"]),
                "experiments": sorted(ev["experiments"]),
            }

    payload = {
        "version": version,
        "source": "https://mirtarbase.cuhk.edu.cn/",
        "built_from": Path(xlsx_path).name,
        "built_at": datetime.now(timezone.utc).isoformat(),
        "rows_read": n_rows,
        "pairs": n_pairs,
        "levels": dict(level_counts),
        "mirnas": mirnas,
    }
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wt", encoding="utf-8") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    print(f"rows: {n_rows}  pairs: {n_pairs}  miRNAs: {len(mirnas)}")
    print(f"levels: {dict(level_counts)}")
    print(f"recovered Excel-mangled gene symbols: {n_recovered}  skipped rows: {n_skipped}  non-human rows: {n_nonhuman}")
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")
    return payload


if __name__ == "__main__":
    xlsx = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_XLSX
    out = (sys.argv[2] if len(sys.argv) > 2 else
           str(Path(__file__).parent / "local_db" / "mirtarbase_official.json.gz"))
    version = "miRTarBase 10.0" if str(xlsx).lower().endswith(".csv") else "miRTarBase 9.0"
    build(xlsx, out, version=version)
