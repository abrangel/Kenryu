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
exactly one candidate that IS in the libraries. POLICY (2026-11-21):
the previous symbol is preferred over an alias — a prev_symbol
candidate wins when it is unique in the libraries; an ambiguous
previous symbol is discarded, never rescued by an alias; aliases
are used only when no previous symbol reaches the libraries.
Ambiguous genes (several candidates at the same level) are NOT
mapped. The complete set has no release number; the version is
the download date, recorded in the output (--download-date).

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
    ap.add_argument("--download-date", default=None,
                    help="date the HGNC file was downloaded (the version; "
                         "HGNC publishes no release number)")
    args = ap.parse_args()

    import kenryu_engine as ke

    # Raw data: the audit measures the pre-mapping baseline, so the
    # HGNC map must NOT be applied here (load_local_data maps by
    # default since v1.59).
    db = ke.load_local_data(mapped=False)
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

    # HGNC mapping. POLICY (2026-11-21): the previous symbol is the
    # authoritative rename and is preferred over an alias. At each
    # level (the three libraries, then KEGG-only), a unique prev
    # candidate wins; an ambiguous prev is discarded, never rescued
    # by an alias; aliases are used only when no previous symbol
    # reaches that level. GUARD (2026-11-21): a gene that is itself
    # an Approved HGNC symbol is NEVER mapped — the string denotes
    # that gene, not the prev/alias target. Without the guard, 45
    # current genes were wrongly renamed (e.g. BCAM -> BCAT2, while
    # TargetScan's BCAM is ENSG00000187244 = HGNC BCAM itself,
    # verified against the Ensembl IDs in the TargetScan source
    # file). Such genes stay as they are: current genes the
    # libraries simply do not hold.
    import pandas as pd
    h = pd.read_csv(args.hgnc, sep="\t", dtype=str)
    hgnc_date = args.download_date
    approved = {r["symbol"] for _, r in h.iterrows()
                if r["status"] == "Approved"}
    prev_map, alias_map = {}, {}
    for _, r in h.iterrows():
        for col, d in (("prev_symbol", prev_map),
                       ("alias_symbol", alias_map)):
            if pd.notna(r[col]):
                for s in r[col].split("|"):
                    d.setdefault(s, set()).add(r["symbol"])
    print(f"HGNC rows: {len(h)}; approved symbols: {len(approved)}; "
          f"symbols with a prev entry: {len(prev_map)}, with an "
          f"alias entry: {len(alias_map)}")

    def _in(cands, lib):
        return sorted(c for c in cands if c in lib)

    recoverable, ambiguous, unmapped, mapped_not_in_libs = [], [], [], []
    recoverable_kegg_only = 0
    recoverable_via_alias = 0
    approved_current = []
    for g in sorted(absent):
        if g in approved:
            approved_current.append(g)
            continue
        # Level 1: the three libraries, previous symbol first.
        in_libs = _in(prev_map.get(g, set()), union3)
        if len(in_libs) == 1:
            recoverable.append((g, in_libs[0]))
            continue
        if len(in_libs) > 1:
            ambiguous.append((g, in_libs))
            continue
        in_libs = _in(alias_map.get(g, set()), union3)
        if len(in_libs) == 1:
            recoverable.append((g, in_libs[0]))
            recoverable_via_alias += 1
            continue
        if len(in_libs) > 1:
            ambiguous.append((g, in_libs))
            continue
        # Level 2: KEGG-only. A candidate found only in KEGG is
        # counted separately: it helps KEGG, not the three libraries.
        in_libs4 = _in(prev_map.get(g, set()), union4)
        if len(in_libs4) == 1:
            recoverable_kegg_only += 1
            mapped_not_in_libs.append((g, sorted(
                prev_map.get(g, set()) | alias_map.get(g, set()))))
            continue
        if len(in_libs4) > 1:
            ambiguous.append((g, in_libs4))
            continue
        in_libs4 = _in(alias_map.get(g, set()), union4)
        if len(in_libs4) == 1:
            recoverable_kegg_only += 1
            recoverable_via_alias += 1
            mapped_not_in_libs.append((g, sorted(
                prev_map.get(g, set()) | alias_map.get(g, set()))))
            continue
        if len(in_libs4) > 1:
            ambiguous.append((g, in_libs4))
            continue
        cands = prev_map.get(g, set()) | alias_map.get(g, set())
        if not cands:
            unmapped.append(g)
        else:
            mapped_not_in_libs.append((g, sorted(cands)))
    print(f"recoverable (unique candidate in the 3 libraries): "
          f"{len(recoverable)} (of which via alias only: "
          f"{recoverable_via_alias})")
    print(f"recoverable only via KEGG: {recoverable_kegg_only}")
    print(f"ambiguous (several candidates in libraries, discarded): "
          f"{len(ambiguous)}")
    print(f"current Approved symbols the libraries lack (never "
          f"mapped): {len(approved_current)}")
    print(f"mapped, candidate not in libraries: {len(mapped_not_in_libs)}")
    print(f"unmapped (no prev/alias entry): {len(unmapped)}")

    # The AARS test case: AARS must reach AARS1, and AARS1 must be in
    # the libraries.
    aars_cands = sorted(prev_map.get("AARS", set())
                        | alias_map.get("AARS", set()))
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
        "hgnc_downloaded": hgnc_date,
        "hgnc_source": "https://storage.googleapis.com/"
                       "public-download-files/hgnc/tsv/tsv/"
                       "hgnc_complete_set.txt",
        "hgnc_licence": "CC0 1.0 (https://www.genenames.org/about/license/)",
        "policy": "prev_symbol preferred over alias_symbol; an "
                  "ambiguous previous symbol is discarded, never "
                  "rescued by an alias; aliases used only when no "
                  "previous symbol reaches the libraries; a gene "
                  "that is itself an Approved HGNC symbol is never "
                  "mapped",
        "recoverable_unique_in_libraries": len(recoverable),
        "recoverable_via_alias_only": recoverable_via_alias,
        "recoverable_only_via_kegg": recoverable_kegg_only,
        "ambiguous_several_in_libraries": len(ambiguous),
        "approved_current_symbol_not_in_libraries": len(approved_current),
        "approved_current_examples": approved_current[:20],
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
