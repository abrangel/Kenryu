"""Audit: gene symbols TargetScan has that no local enrichment library
has — and how many HGNC can recover.

Claude (2026-11-21) reported, on v1.56 data: 2,496 of 13,035 TargetScan
genes (19.1%) appear in NO local enrichment library; among them 318
look like C10orf10 (open-reading-frame names) and 246 like KIAA / FAM /
LOC families. Many were renamed (e.g. AARS is now AARS1); others are
clone-based identifiers (AC000003.2) with no approved symbol. Renamed
genes vanish silently from enrichment — and could vanish from the
validated-evidence layer too, if miRTarBase uses newer symbols.

This script reproduces those counts first (against the union of the
4 local GMTs, and against the 3 non-KEGG GMTs, to pin down which
"three libraries" the report meant), then measures the HGNC recovery:
for each absent gene, the prev_symbol / alias_symbol columns of the
HGNC complete set give candidate current symbols; "recoverable" means
exactly one candidate that IS in the libraries.

HGNC terms: all data released under CC0 1.0
(https://www.genenames.org/about/license/); attribution recommended,
not mandatory. The complete set has no release number; the version is
the download date, recorded in the output.

Run from the repo root (the HGNC file is downloaded, not committed):
    python scripts/audit_gene_symbols.py --hgnc /tmp/hgnc_complete_set.txt
Writes docs/validation/gene-symbol-audit.json and prints the verdict.
"""

import argparse
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

GMT_PATHS = {
    "GO": REPO / "local_db" / "gmt_GO_Biological_Process_2023.gmt",
    "Reactome": REPO / "local_db" / "gmt_Reactome_Pathways_2024.gmt",
    "WikiPathways": REPO / "local_db" / "gmt_WikiPathways_2019_Human.gmt",
    "KEGG": REPO / "local_db" / "gmt_KEGG_2021_Human.gmt",
}
OUT = REPO / "docs" / "validation" / "gene-symbol-audit.json"

ORF_RE = re.compile(r"^C\d+orf\d+$")
FAMILY_RE = re.compile(r"^(KIAA|FAM\d*|LOC\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True,
                    help="path to hgnc_complete_set.txt (downloaded, "
                         "not committed)")
    args = ap.parse_args()

    import kenryu_engine as ke

    db = ke.load_local_data()
    universe = set().union(*db.values())
    libs = {name: ke.load_gmt_library(str(p)) for name, p in GMT_PATHS.items()}
    lib_genes = {name: set().union(*terms.values()) if terms else set()
                 for name, terms in libs.items()}
    union4 = set().union(*lib_genes.values())
    union3 = set().union(*[lib_genes[n] for n in
                            ("GO", "Reactome", "WikiPathways")])

    absent4 = universe - union4
    absent3 = universe - union3
    # Which union reproduces Claude's 2,496?
    absent = absent4 if len(absent4) == 2496 else (
        absent3 if len(absent3) == 2496 else absent4)
    union_used = ("union-of-4-GMTs" if absent is absent4 else
                  "union-of-3-non-KEGG-GMTs" if absent is absent3 else
                  "union-of-4-GMTs (neither matches 2,496)")
    print(f"TargetScan universe: {len(universe)}")
    print(f"absent vs union-of-4: {len(absent4)}; "
          f"vs union-of-3 (no KEGG): {len(absent3)}")
    print(f"using: {union_used} -> absent={len(absent)}")

    orf = sorted(g for g in absent if ORF_RE.match(g))
    fam = sorted(g for g in absent if FAMILY_RE.match(g))
    print(f"absent C<digit>orf<digit>: {len(orf)} (reported 318)")
    print(f"absent KIAA/FAM/LOC: {len(fam)} (reported 246)")

    # HGNC mapping (Claude's snippet, verbatim logic).
    import pandas as pd
    h = pd.read_csv(args.hgnc, sep="\t", dtype=str)
    hgnc_date = None
    prev = {}
    for _, r in h.iterrows():
        for col in ("prev_symbol", "alias_symbol"):
            if pd.notna(r[col]):
                for s in r[col].split("|"):
                    prev.setdefault(s, set()).add(r["symbol"])
    print(f"HGNC rows: {len(h)}; symbols with a prev/alias entry: "
          f"{len(prev)}")

    recoverable, ambiguous, unmapped, mapped_not_in_libs = [], [], [], []
    recoverable_kegg_only = 0
    for g in sorted(absent):
        cands = prev.get(g)
        if not cands:
            unmapped.append(g)
            continue
        # Primary: the candidate must be in the SAME libraries the
        # absence was measured against (union3). A candidate found
        # only in KEGG is counted separately: it helps KEGG, not the
        # three libraries above.
        in_libs = sorted(c for c in cands if c in union3)
        in_libs4 = sorted(c for c in cands if c in union4)
        if len(in_libs) == 1:
            recoverable.append((g, in_libs[0]))
        elif len(in_libs) > 1:
            ambiguous.append((g, in_libs))
        elif len(in_libs4) == 1:
            recoverable_kegg_only += 1
            mapped_not_in_libs.append((g, sorted(cands)))
        else:
            mapped_not_in_libs.append((g, sorted(cands)))
    print(f"recoverable (unique candidate in the 3 libraries): "
          f"{len(recoverable)}")
    print(f"recoverable only via KEGG: {recoverable_kegg_only}")
    print(f"ambiguous (several candidates in libraries): {len(ambiguous)}")
    print(f"mapped, candidate not in libraries: {len(mapped_not_in_libs)}")
    print(f"unmapped (no prev/alias entry): {len(unmapped)}")

    # The AARS test case: AARS must reach AARS1, and AARS1 must be in
    # the libraries.
    aars_cands = sorted(prev.get("AARS", set()))
    aars_ok = ("AARS1" in aars_cands and "AARS1" in union4
               and "AARS" in absent)
    print(f"AARS -> {aars_cands}; AARS1 in libraries: "
          f"{'AARS1' in union4}; AARS absent: {'AARS' in absent}; "
          f"test: {'PASS' if aars_ok else 'FAIL'}")

    report = {
        "targetscan_universe": len(universe),
        "libraries": {n: len(v) for n, v in lib_genes.items()},
        "absent_vs_union4": len(absent4),
        "absent_vs_union3_no_kegg": len(absent3),
        "union_used": union_used,
        "absent": len(absent),
        "absent_orf_like": len(orf),
        "absent_kiaa_fam_loc": len(fam),
        "hgnc_rows": len(h),
        "hgnc_source": "https://storage.googleapis.com/"
                       "public-download-files/hgnc/tsv/tsv/"
                       "hgnc_complete_set.txt",
        "hgnc_licence": "CC0 1.0 (https://www.genenames.org/about/license/)",
        "recoverable_unique_in_libraries": len(recoverable),
        "recoverable_only_via_kegg": recoverable_kegg_only,
        "ambiguous_several_in_libraries": len(ambiguous),
        "mapped_candidate_not_in_libraries": len(mapped_not_in_libs),
        "unmapped_no_prev_or_alias": len(unmapped),
        "recoverable_examples": recoverable[:20],
        "ambiguous_examples": ambiguous[:10],
        "aars_test": {"AARS_candidates": aars_cands,
                      "AARS1_in_libraries": "AARS1" in union4,
                      "AARS_absent": "AARS" in absent,
                      "pass": aars_ok},
    }
    OUT.write_text(json.dumps(report, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print(f"-> {OUT.name}")


if __name__ == "__main__":
    main()
