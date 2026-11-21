"""Audit: does any gene reach level "alta" only thanks to validated
pairs whose arm could NOT be verified?

The v1.58 guard keeps a bare-name fallback pair when the arm cannot be
checked (the bare key is not a TargetScan family, or the seed files are
missing) and logs a warning. Those pairs are experimental evidence, but
their arm is unverified, so this audit asks whether any reference-panel
gene reaches "alta" only because of them.

Method: run the three reference panels (lipid, cancer, control) through
/api/v1/analyze twice, in separate processes (no cache is shared):
  - current: the v1.58 guard (unverifiable fallback pairs are kept).
  - strict:  the guard replaced by a strict variant that keeps a
             fallback pair ONLY when the seed positively verifies it.
Compare the per-gene level in both modes. Expected: no gene is "alta"
in current and lower in strict.

Also counts, per panel miRNA, how many fallback pairs are unverifiable
(and how many of them are strong), from the real local index.

Run from the repo root:
    python scripts/audit_arm_unverified_high.py
Writes docs/validation/arm-unverified-high.json and prints the verdict.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

PANELS = {
    "lipid": ["hsa-miR-33a-5p", "hsa-miR-144-3p", "hsa-miR-758-3p",
              "hsa-miR-106b-5p", "hsa-miR-33b-5p"],
    "cancer": ["hsa-miR-21-5p", "hsa-miR-155-5p", "hsa-let-7a-5p"],
    "control": ["hsa-miR-375", "hsa-miR-29a-3p", "hsa-miR-21-5p"],
}
OUT_DIR = REPO / "docs" / "validation"


def classify_fallback(ke, mirna, fallback_key):
    """Why a bare-name fallback is kept or dropped, replicating the
    v1.58 guard's logic: 'verified' (seed matches), 'wrong-arm'
    (seed absent), 'unverifiable' (no sequence or no family for the
    key), 'files-missing'."""
    seqs, fam = ke._load_arm_seed_maps()
    if seqs is None or fam is None:
        return "files-missing"
    qseq = seqs.get(ke.deep_clean_mirna(mirna))
    fseeds = fam.get(ke.deep_clean_mirna(fallback_key))
    if not qseq or not fseeds:
        return "unverifiable"
    return "verified" if ke.verify_arm_seed(qseq, fseeds) else "wrong-arm"


def run_mode(mode, out_path):
    import kenryu_engine as ke
    from fastapi.testclient import TestClient

    classifications = []

    if mode == "strict":
        def strict_guard(mirna, fallback_key):
            verdict = classify_fallback(ke, mirna, fallback_key)
            classifications.append(
                {"mirna": mirna, "fallback_key": fallback_key,
                 "verdict": verdict})
            return verdict == "verified"
        ke._validated_fallback_arm_ok = strict_guard
    else:
        orig = ke._validated_fallback_arm_ok

        def counting_guard(mirna, fallback_key):
            verdict = classify_fallback(ke, mirna, fallback_key)
            classifications.append(
                {"mirna": mirna, "fallback_key": fallback_key,
                 "verdict": verdict})
            return orig(mirna, fallback_key)
        ke._validated_fallback_arm_ok = counting_guard

    # Unverifiable fallback pairs per panel miRNA, from the real index:
    # a miRNA with no exact hit in the index gets ALL its evidence via
    # the fallback; if that fallback is unverifiable, those are the pairs.
    ke.load_mirtarbase_local()
    db = ke.MIRTARBASE_DB
    unverified_pairs = {}
    for panel, mirnas in PANELS.items():
        for m in mirnas:
            name = ke.clean_mirna_name_cloud(m, "hsa")
            if name in db:
                continue
            base = re.sub(r"-(5p|3p)$", "", name, flags=re.IGNORECASE)
            verdict = classify_fallback(ke, m, base)
            if verdict in ("unverifiable", "files-missing"):
                ev = ke.get_validated_evidence(m, db, "hsa")
                unverified_pairs[m] = {
                    "fallback_key": base, "verdict": verdict,
                    "genes": sorted(ev.keys()),
                    "strong": sorted(g for g, e in ev.items()
                                     if e.get("level") == "strong"),
                }

    levels = {}
    with TestClient(ke.app) as client:
        for panel, mirnas in PANELS.items():
            r = client.post("/api/v1/analyze", json={
                "mirnas": mirnas, "years": 10, "mode": "strict",
                "offline": False})
            assert r.status_code == 200, (panel, r.status_code, r.text[:300])
            genes = r.json()["validation"]["genes"]
            levels[panel] = {g: {"score": v["score"], "level": v["level"]}
                             for g, v in genes.items()}

    result = {"mode": mode, "levels": levels,
              "classifications": classifications,
              "unverified_pairs": unverified_pairs}
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False),
                        encoding="utf-8")
    n_unver = sum(len(v["genes"]) for v in unverified_pairs.values())
    n_strong = sum(len(v["strong"]) for v in unverified_pairs.values())
    print(f"[{mode}] panels={len(levels)} genes="
          f"{sum(len(v) for v in levels.values())} "
          f"unverified-fallback-pairs={n_unver} (strong={n_strong}) "
          f"-> {out_path.name}")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["current", "strict"])
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()
    if args.mode:
        run_mode(args.mode, args.out or
                 OUT_DIR / f"arm-unverified-high-{args.mode}.json")
        return

    cur_path = OUT_DIR / "arm-unverified-high-current.json"
    strict_path = OUT_DIR / "arm-unverified-high-strict.json"
    env = dict(os.environ)
    for mode, path in (("current", cur_path), ("strict", strict_path)):
        subprocess.run([sys.executable, __file__, "--mode", mode,
                        "--out", str(path)], check=True, cwd=REPO, env=env)
    cur = json.loads(cur_path.read_text(encoding="utf-8"))
    strict = json.loads(strict_path.read_text(encoding="utf-8"))

    drops = []  # genes "alta" in current, lower in strict
    for panel in PANELS:
        for gene, cur_v in cur["levels"][panel].items():
            st_v = strict["levels"][panel].get(gene)
            if not st_v:
                continue
            if cur_v["level"] == "alta" and st_v["level"] != "alta":
                drops.append({"panel": panel, "gene": gene,
                              "current": cur_v, "strict": st_v})
    n_unver = sum(len(v["genes"])
                  for v in cur["unverified_pairs"].values())
    n_strong = sum(len(v["strong"])
                   for v in cur["unverified_pairs"].values())
    report = {
        "panels": sorted(PANELS),
        "unverified_fallback_pairs": n_unver,
        "unverified_fallback_strong_pairs": n_strong,
        "unverified_pairs_by_mirna": cur["unverified_pairs"],
        "genes_alta_only_via_unverified": drops,
        "verdict": ("PASS: no gene reaches 'alta' only via unverified "
                    "fallback pairs" if not drops else
                    "FAIL: genes reach 'alta' only via unverified pairs"),
    }
    out = OUT_DIR / "arm-unverified-high.json"
    out.write_text(json.dumps(report, indent=2, ensure_ascii=False),
                   encoding="utf-8")
    print(f"unverified fallback pairs in the reference panels: {n_unver} "
          f"(strong: {n_strong})")
    print(f"genes 'alta' only via unverified pairs: {len(drops)}")
    for d in drops:
        print(f"  {d['panel']}: {d['gene']} {d['current']} -> {d['strict']}")
    print(report["verdict"])
    print(f"-> {out.name}")


if __name__ == "__main__":
    main()
