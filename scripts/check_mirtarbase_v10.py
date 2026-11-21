#!/usr/bin/env python3
"""Check the miRTarBase 2025 (v10.0) CSV for the lipid-benchmark pairs.

Reads the official hsa_MTI.csv (UTF-8 BOM, PMID/Entrez formatted
with ".0"), filters to human (Species codes "hsa" on both sides),
and reports: the four ABCA1 PMIDs, the ABCA1 rows, and the evidence
each panel miRNA has for each of the five panel genes.

Usage: python3 scripts/check_mirtarbase_v10.py [csv_path]
Default: the browser-downloaded copy in ~/workspace/browser_downloads/
"""
import sys
from pathlib import Path

import pandas as pd

DEFAULT_CSV = sorted(
    Path.home().glob("workspace/browser_downloads/*/browser-download-*-hsa_MTI.csv")
)
DEFAULT_CSV = str(DEFAULT_CSV[-1]) if DEFAULT_CSV else "hsa_MTI.csv"

PANEL_MIRNAS = ["hsa-miR-33a-5p", "hsa-miR-33b-5p", "hsa-miR-144-3p",
                "hsa-miR-758-3p", "hsa-miR-106b-5p"]
PANEL_GENES = ["ABCA1", "KPNA3", "SCN1A", "SNTB2", "TSC22D2"]
FOUR_PMIDS = {"21885853", "23519695", "23519696", "22119192"}


def _norm_pmid(v):
    s = str(v).strip()
    return s[:-2] if s.endswith(".0") else s


def main():
    csv = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV
    cols = ["miRNA", "Species (miRNA)", "Target Gene",
            "Species (Target Gene)", "Experiments", "Support Type",
            "References (PMID)"]
    df = pd.read_csv(csv, encoding="utf-8-sig", dtype=str, usecols=cols)
    print(f"rows: {len(df)}")
    df["pmid"] = df["References (PMID)"].map(_norm_pmid)
    human = df[(df["Species (miRNA)"] == "hsa")
               & (df["Species (Target Gene)"] == "hsa")].copy()
    print(f"human rows: {len(human)}")
    print("human pairs:",
          human[["miRNA", "Target Gene"]].drop_duplicates().shape[0])
    print("human Support Type counts:")
    print(human["Support Type"].value_counts(dropna=False).to_string())

    hit = human[human["pmid"].isin(FOUR_PMIDS)]
    print(f"\nrows with the 4 PMIDs: {len(hit)}")
    if len(hit):
        print(hit[["miRNA", "Target Gene", "Support Type", "Experiments",
                   "pmid"]].drop_duplicates().to_string(index=False))
    else:
        print("NONE — no human row carries any of the 4 PMIDs")

    abca1 = human[human["Target Gene"] == "ABCA1"]
    print(f"\nABCA1 human rows (deduped): "
          f"{abca1[['miRNA', 'Support Type', 'pmid']].drop_duplicates().shape[0]}")
    print("ABCA1 miRNAs:", sorted(abca1["miRNA"].unique()))
    print(abca1[["miRNA", "Support Type", "pmid"]]
          .drop_duplicates().sort_values(["miRNA", "pmid"])
          .to_string(index=False))

    print("\nPanel evidence in v10.0 (human rows only):")
    for gene in PANEL_GENES:
        for mir in PANEL_MIRNAS:
            sub = human[(human["Target Gene"] == gene)
                        & (human["miRNA"] == mir)]
            if len(sub):
                supports = sorted(sub["Support Type"].dropna().unique())
                pmids = sorted(sub["pmid"].unique())
                print(f"  {mir} -> {gene}: {supports} PMIDs={pmids}")
            else:
                print(f"  {mir} -> {gene}: no rows")


if __name__ == "__main__":
    main()
