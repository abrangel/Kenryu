"""Regenerate the KENRYU baseline for the tool comparison
(docs/validation/comparison-kenryu-baseline.json): for each panel in
docs/validation/comparison-panels.json, the Strict common-gene count
and local-term count from the local engine (offline), plus the
convergence p range (2,000 permutations, tolerances 10/25/50%, seed
from the family labels — deterministic). Re-run after any data change
that moves gene counts (arm fix, HGNC map, source update).

Run from the repo root:
    python scripts/regen_comparison_baseline.py
"""

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

PANELS_PATH = REPO / "docs" / "validation" / "comparison-panels.json"
OUT = REPO / "docs" / "validation" / "comparison-kenryu-baseline.json"


def main():
    import kenryu_engine as ke
    import kenryu_core as kc

    panels_doc = json.loads(PANELS_PATH.read_text(encoding="utf-8"))
    panels = ([{"name": name, "mirnas": mirnas}
               for name, mirnas in panels_doc["reference"].items()]
              + [{"name": r["id"], "mirnas": r["mirnas"]}
                 for r in panels_doc["random"]])
    db = ke.load_local_data()
    ke.load_gmt_libraries()
    out = []
    for p in panels:
        mirnas = p["mirnas"]
        found, genes = ke.panel_core(mirnas, "strict")
        # Local enrichment terms (the method): Fisher + BH on the GMTs.
        enrich = kc.local_enrichment(
            genes, ke.GMT_LIBRARIES, ke.TARGETSCAN_DB)
        sig = [t for t in enrich if t.get("Pval", 1) < 0.05]
        # Convergence p range (API defaults).
        sets = [set(ke.get_targets_local(m, db)) for m in mirnas]
        conv = ke.panel_convergence(
            sets, ke.family_labels_for_targets(mirnas, sets, db), db)
        out.append({
            "panel": p["name"],
            "n_mirnas": len(mirnas),
            "n_genes": len(genes),
            "local_terms": len(enrich),
            "local_sig_terms": len(sig),
            "convergence_p": conv["p_by_tolerance"],
        })
        print(f"{p['name']:8s} genes={len(genes):4d} "
              f"local_terms={len(enrich)} sig={len(sig)} "
              f"p={conv['p_by_tolerance']}")
    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print(f"-> {OUT.name}")


if __name__ == "__main__":
    main()
