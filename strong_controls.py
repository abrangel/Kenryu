"""Block 3: positive controls from STRONG miRTarBase evidence, and
TargetScan sensitivity against known strong targets.

Why: the earlier positive controls (miR-17, let-7, miR-34) used
PREDICTED targets and were weak. These use the strong-evidence
target list of each miRNA from the official miRTarBase 9.0 index
(Functional MTI: luciferase, qRT-PCR, western blot) and run the
calibrated LOCAL enrichment (TargetScan-universe background) on it.

For each miRNA it reports:
  - n_strong / n_weak targets in the official index
  - TargetScan sensitivity: fraction of strong (and weak) targets that
    TargetScan predicts (context: TargetScan only predicts conserved
    sites, so sensitivity below 100% is expected, not a failure)
  - the top 5 local enrichment terms of the strong list, to check
    whether the miRNA's known biology is recovered

Run from the repo root with the test venv:
    ~/workspace/defensa-kenryu/.venv-test/bin/python \\
        ~/workspace/defensa-kenryu/validation-panels/strong_controls.py
"""
import json
import sys

sys.path.insert(0, "/home/hatch/workspace/defensa-kenryu/Kenryu-repo")

import kenryu_engine as ke

CANDIDATES = {
    # miRNA: expected biology (checked by eye against the top terms)
    "hsa-miR-21-5p": "apoptosis / PI3K-Akt / cancer (PTEN, PDCD4)",
    "hsa-miR-33a-5p": "cholesterol / lipid metabolism (ABCA1)",
    "hsa-let-7a-5p": "cell cycle / development (HMGA2, LIN28)",
    "hsa-miR-155-5p": "immune response / inflammation (SOCS1)",
    "hsa-miR-34a-5p": "cell cycle / p53 / apoptosis (BCL2, CDK4/6)",
    "hsa-miR-17-5p": "cell cycle (E2F)",
    "hsa-miR-16-5p": "cell cycle / apoptosis (BCL2)",
    "hsa-miR-29a-3p": "extracellular matrix / collagen",
}

OUT = "/home/hatch/workspace/defensa-kenryu/validation-panels/strong-controls.json"


def main():
    ke.TARGETSCAN_DB = ke.load_local_data()
    db = ke.load_mirtarbase_local()
    ke.load_gmt_libraries()
    universe = set().union(*ke.TARGETSCAN_DB.values())
    print(f"universe: {len(universe)} genes; miRTarBase official: "
          f"{len(db)} miRNAs")

    rows = []
    tot_strong = tot_rec = 0
    for mirna, expected in CANDIDATES.items():
        ev = ke.get_validated_evidence(mirna, db)
        strong = {g for g, v in ev.items() if v.get("level") == "strong"}
        weak = {g for g, v in ev.items() if v.get("level") == "weak"}
        ts = ke.get_targets_local(mirna)
        sens = len(strong & ts) / len(strong) if strong else None
        sens_w = len(weak & ts) / len(weak) if weak else None
        tot_strong += len(strong)
        tot_rec += len(strong & ts)
        terms = ke.local_enrichment(sorted(strong), ke.GMT_LIBRARIES,
                                    background=universe)
        top = [{"term": t["Term"], "source": t["Source"],
                "p_adj": t["Pval"]} for t in terms[:5]]
        rows.append({"mirna": mirna, "expected": expected,
                     "n_strong": len(strong), "n_weak": len(weak),
                     "n_targetscan": len(ts),
                     "sensitivity_strong": round(sens, 3) if sens is not None else None,
                     "sensitivity_weak": round(sens_w, 3) if sens_w is not None else None,
                     "top_terms": top})
        print(f"\n{mirna} — expected: {expected}")
        print(f"  strong {len(strong)}, weak {len(weak)}, "
              f"TargetScan {len(ts)}; sensitivity strong {sens:.3f}, weak {sens_w:.3f}")
        for t in top:
            print(f"    {t['p_adj']:.4f} [{t['source']}] {t['term'][:70]}")

    pooled = tot_rec / tot_strong if tot_strong else None
    print(f"\npooled sensitivity (strong targets): {tot_rec}/{tot_strong} = {pooled:.3f}")
    json.dump({"version": "1.54", "pooled_sensitivity_strong": pooled,
               "controls": rows}, open(OUT, "w"), indent=2, ensure_ascii=False)
    print(f"saved to {OUT}")


if __name__ == "__main__":
    main()
