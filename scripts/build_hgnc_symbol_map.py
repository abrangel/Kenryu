"""Build the HGNC symbol map KENRYU applies at load time.

From the HGNC complete set (CC0, downloaded — never committed), map
each TargetScan gene that is absent from the enrichment libraries to
its unique current symbol in the libraries (the audit in
scripts/audit_gene_symbols.py). POLICY (2026-11-21): the previous
symbol is preferred over an alias — a unique prev_symbol candidate
wins; an ambiguous previous symbol is discarded, never rescued by
an alias; aliases are used only when no previous symbol reaches the
libraries. Ambiguous genes are NOT mapped; unmapped genes stay
as they are.

Collisions: several old symbols can map to the same current symbol,
and some targets are already TargetScan genes under their current
name — the same gene under two symbols. Mapping merges them, so the
universe shrinks. No key is itself a target (no chains), so one
pass suffices. COLLISION RULE (2026-11-21): two old symbols may
merge into one current symbol ONLY if they denote the same gene —
the same Ensembl ID. Distinct genes colliding on one symbol (e.g.
paralogs) must never merge; the builder refuses to write the map
instead, so a human decides. Verified for the four groups in the
2026-11-21 map: SEM1 <- {C7orf76, SHFM1} (ENSG00000127922),
WASHC2A <- {FAM21A, FAM21B} (ENSG00000099290; the paralog
WASHC2C is ENSG00000172661, prev symbol FAM21C — no cross-over),
PTPN20 <- {PTPN20A, PTPN20B} (ENSG00000204179), SCX <- {SCXA,
SCXB} (ENSG00000260428). GUARD (2026-11-21): a gene that is itself an
Approved HGNC symbol is NEVER mapped — the string denotes that
gene, not the prev/alias target (verified against the Ensembl IDs
in the TargetScan source file: TargetScan's BCAM is ENSG00000187244
= HGNC BCAM, not BCAT2).

The map is a versioned data file (local_db/hgnc_symbol_map.json).
HGNC publishes no release number: the version is the download date,
recorded in the file and in data_versions.json. Rebuild with:
    python scripts/build_hgnc_symbol_map.py --hgnc /path/to/hgnc_complete_set.txt --download-date 2026-11-21
The script refuses to write if the counts drift from the audit
— that means the HGNC file changed, and the audit must be
re-run first.
"""

import argparse
import json
import sys
from datetime import date
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

GMT_PATHS = {
    "GO": REPO / "local_db" / "gmt_GO_Biological_Process_2023.gmt",
    "Reactome": REPO / "local_db" / "gmt_Reactome_Pathways_2024.gmt",
    "WikiPathways": REPO / "local_db" / "gmt_WikiPathways_2019_Human.gmt",
    "KEGG": REPO / "local_db" / "gmt_KEGG_2021_Human.gmt",
}
OUT = REPO / "local_db" / "hgnc_symbol_map.json"
HGNC_URL = ("https://storage.googleapis.com/public-download-files/"
            "hgnc/tsv/tsv/hgnc_complete_set.txt")
EXPECTED = {"entries": 498, "in_three_libraries": 487,
            "kegg_only": 11, "ambiguous": 3}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--download-date", default=str(date.today()),
                    help="date the HGNC file was downloaded (the version; "
                         "HGNC publishes no release number)")
    args = ap.parse_args()

    import kenryu_engine as ke

    # Raw data: the map is built from unmapped symbols, so the HGNC
    # map must NOT be applied here (load_local_data maps by default
    # since v1.59).
    db = ke.load_local_data(mapped=False)
    universe = set().union(*db.values())
    libs = {n: ke.load_gmt_library(str(p)) for n, p in GMT_PATHS.items()}
    u3 = set().union(*[set().union(*libs[n].values())
                       for n in ("GO", "Reactome", "WikiPathways")])
    u4 = set().union(*[set().union(*v.values()) for v in libs.values()])
    absent = universe - u3

    import pandas as pd
    h = pd.read_csv(args.hgnc, sep="\t", dtype=str)
    approved = {r["symbol"] for _, r in h.iterrows()
                if r["status"] == "Approved"}
    prev_map, alias_map = {}, {}
    for _, r in h.iterrows():
        for col, d in (("prev_symbol", prev_map),
                       ("alias_symbol", alias_map)):
            if pd.notna(r[col]):
                for s in r[col].split("|"):
                    d.setdefault(s, set()).add(r["symbol"])

    def _in(cands, lib):
        return sorted(c for c in cands if c in lib)

    # POLICY (2026-11-21): previous symbol preferred over alias, at
    # each level (the three libraries, then KEGG-only). An ambiguous
    # previous symbol is discarded, never rescued by an alias.
    # GUARD: a gene that is itself an Approved symbol is never
    # mapped — it denotes that gene, not the prev/alias target.
    mapping, ambiguous = {}, {}
    n_kegg_only = 0
    n_approved_skipped = 0
    for g in sorted(absent):
        if g in approved:
            n_approved_skipped += 1
            continue
        in3 = _in(prev_map.get(g, set()), u3)
        if len(in3) == 1:
            mapping[g] = in3[0]
            continue
        if len(in3) > 1:
            ambiguous[g] = in3
            continue
        in3 = _in(alias_map.get(g, set()), u3)
        if len(in3) == 1:
            mapping[g] = in3[0]
            continue
        if len(in3) > 1:
            ambiguous[g] = in3
            continue
        in4 = _in(prev_map.get(g, set()), u4)
        if len(in4) == 1:
            mapping[g] = in4[0]
            n_kegg_only += 1
            continue
        if len(in4) > 1:
            ambiguous[g] = in4
            continue
        in4 = _in(alias_map.get(g, set()), u4)
        if len(in4) == 1:
            mapping[g] = in4[0]
            n_kegg_only += 1
            continue
        if len(in4) > 1:
            ambiguous[g] = in4
        # else: mapped to symbols no library holds — not mapped, and
        # not ambiguous either (audit: 311 of those).

    # Collisions: several old symbols -> one current symbol.
    by_target = {}
    for old, new in mapping.items():
        by_target.setdefault(new, []).append(old)
    collisions = {t: sorted(v) for t, v in by_target.items()
                  if len(v) > 1}
    # COLLISION RULE (2026-11-21): two old symbols may merge into
    # one current symbol ONLY if they denote the same gene — the
    # same Ensembl ID. Distinct genes colliding on one symbol must
    # never merge; refuse to write the map so a human decides.
    ensembl_of_row = {r["symbol"]: r.get("ensembl_gene_id") or ""
                      for _, r in h.iterrows()}
    denoted = {}  # old symbol -> {(approved, ensembl)} via prev/alias
    for _, r in h.iterrows():
        for col in ("prev_symbol", "alias_symbol"):
            if pd.notna(r[col]):
                for s in r[col].split("|"):
                    denoted.setdefault(s, set()).add(
                        (r["symbol"], r.get("ensembl_gene_id") or ""))
    collision_ensembl = {}
    for target, members in collisions.items():
        ids = set()
        for m in members:
            ids |= {e for (_, e) in denoted.get(m, set()) if e}
        target_ensembl = ensembl_of_row.get(target, "")
        same = len(ids) == 1 and target_ensembl in ids
        collision_ensembl[target] = {
            "members": members,
            "ensembl_gene_id": sorted(ids),
            "target_ensembl_gene_id": target_ensembl,
            "same_gene": same,
        }
        if not same:
            raise SystemExit(
                f"COLLISION: {target} <- {members} denote different "
                f"genes (Ensembl {sorted(ids)}, target "
                f"{target_ensembl}) — do not merge; remove them or "
                f"keep them separate.")
    targets_in_universe = sorted(set(mapping.values()) & universe)
    mapped_universe = {mapping.get(g, g) for g in universe}

    counts = {"entries": len(mapping),
              "in_three_libraries": len(mapping) - n_kegg_only,
              "kegg_only": n_kegg_only,
              "ambiguous": len(ambiguous),
              "approved_current_skipped": n_approved_skipped,
              "distinct_targets": len(by_target),
              "collision_groups": len(collisions),
              "merged_entries": len(mapping) - len(by_target),
              "targets_already_in_universe": len(targets_in_universe),
              "universe_before": len(universe),
              "universe_after": len(mapped_universe)}
    print(json.dumps(counts, indent=2))
    for key in ("entries", "in_three_libraries", "kegg_only",
                "ambiguous"):
        if counts[key] != EXPECTED[key]:
            raise SystemExit(
                f"DRIFT: {key} is {counts[key]}, audit said "
                f"{EXPECTED[key]} — re-run audit_gene_symbols.py before "
                f"rebuilding the map.")

    doc = {
        "hgnc_source": HGNC_URL,
        "hgnc_downloaded": args.download_date,
        "hgnc_licence": "CC0 1.0 "
                         "(https://www.genenames.org/about/license/)",
        "hgnc_rows": len(h),
        "scope": "TargetScan genes absent from GO/Reactome/WikiPathways "
                 "with a unique current HGNC symbol in the libraries "
                 "(audit: scripts/audit_gene_symbols.py). POLICY: "
                 "prev_symbol preferred over alias_symbol; ambiguous "
                 "genes are NOT mapped.",
        "counts": counts,
        "collisions": collisions,
        "collision_ensembl": collision_ensembl,
        "targets_already_in_universe": targets_in_universe,
        "ambiguous_not_mapped": ambiguous,
        "mapping": mapping,
    }
    OUT.write_text(json.dumps(doc, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print(f"-> {OUT.name} ({len(mapping)} entries)")


if __name__ == "__main__":
    main()
