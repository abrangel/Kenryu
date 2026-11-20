"""Calibrate the LIVE Enrichr route with random gene lists.

The paper's enrichment p-values (thesis and cancer panels) come from
Enrichr's own background, but the random-panel calibration in
docs/validation/ measured the LOCAL route, not Enrichr's. This script
measures the live route directly: for N random gene lists of a given
size, drawn from the TargetScan universe, it submits each list to
Enrichr and counts the fraction with at least one significant term
(adjusted p < 0.05) per library. A calibrated method should give ~5%.

The criterion is EXACTLY the engine's (kenryu_engine.py,
enrich_genes_direct): the first 8 rows of each library, row[6]
(adjusted p) < 0.05. An earlier version of this script counted ALL
rows with row[6] < 0.05, which overstates what the engine publishes
(the engine shows at most 8 terms per library, 24 total). The raw
p-value counts (row[2]) are recorded as a control, to show the
difference between the raw and adjusted criteria.

Verified 2026-11-19 on one random 5-gene list (ARL8B, COL4A2,
IFITM3, MAN2A2, SPATA2): the top GO rows have adjusted p 0.025 with
a ONE-gene overlap (e.g. "Protein Deglycosylation", MAN2A2 only,
raw p 0.0017, adjusted p 0.025), and all 8 top rows share the same
adjusted p. Enrichr's "adjusted p" does not behave like a
Benjamini-Hochberg correction over the whole library at these list
sizes; the numbers below are what Enrichr reports, and the paper
reports them as such, with the date and the library versions.

Request limits: one addList POST + four enrich GETs per list, with a
delay between requests and backoff on 429. Results are saved after
every list, so an interrupted run resumes with --resume instead of
starting over. Default: 100 lists of 5 genes (thesis-like), seed 7.

Run from the repo root with the test venv (takes ~10-15 min for 100
lists; run it in the background):
    ~/workspace/defensa-kenryu/.venv-test/bin/python \\
        calibrate_enrichr.py --panels 100 --size 5 --seed 7 \\
        --out ~/workspace/defensa-kenryu/validation-panels/calibration-enrichr-K5-engine.json
"""
import argparse
import asyncio
import json
import os
import random
import sys
import time

sys.path.insert(0, "/home/hatch/workspace/defensa-kenryu/Kenryu-repo")

import httpx

import kenryu_engine as ke
from kenryu_updater import make_async_client

SOURCES = {
    "GO_Biological_Process_2023": "GO",
    "Reactome_Pathways_2024": "Reactome",
    "WikiPathways_2019_Human": "WikiPathways",
    "KEGG_2021_Human": "KEGG",
}


async def enrich_one_list(client, genes, delay):
    """Submit one gene list to Enrichr; return per-library counts.

    The engine's criterion: first 8 rows, row[6] (adjusted p) < 0.05.
    The raw-p counts (row[2]) and the all-rows counts are recorded as
    controls, so the difference between criteria is visible."""
    payload = {"list": (None, "\n".join(genes)),
               "description": (None, "KENRYU calibration")}
    user_list_id = None
    for attempt in range(3):
        try:
            r = await client.post("https://maayanlab.cloud/Enrichr/addList",
                                  files=payload, timeout=20.0)
            if r.status_code == 200:
                user_list_id = r.json().get("userListId")
                break
            if r.status_code == 429:
                await asyncio.sleep(2 ** attempt * 5)
                continue
            return None
        except Exception:
            await asyncio.sleep(2 ** attempt * 5)
    if not user_list_id:
        return None
    out = {}
    for gs_full, label in SOURCES.items():
        counts = {"adj_top8": 0, "raw_top8": 0, "adj_all": 0, "raw_all": 0}
        for attempt in range(3):
            try:
                r = await client.get(
                    "https://maayanlab.cloud/Enrichr/enrich",
                    params={"userListId": user_list_id,
                            "backgroundType": gs_full},
                    timeout=20.0)
                if r.status_code == 200:
                    rows = r.json().get(gs_full, [])
                    counts["adj_top8"] = sum(1 for row in rows[:8] if row[6] < 0.05)
                    counts["raw_top8"] = sum(1 for row in rows[:8] if row[2] < 0.05)
                    counts["adj_all"] = sum(1 for row in rows if row[6] < 0.05)
                    counts["raw_all"] = sum(1 for row in rows if row[2] < 0.05)
                    break
                if r.status_code == 429:
                    await asyncio.sleep(2 ** attempt * 5)
                    continue
                break
            except Exception:
                await asyncio.sleep(2 ** attempt * 5)
        out[label] = counts
        await asyncio.sleep(delay)
    return out


async def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", type=int, default=100)
    ap.add_argument("--size", type=int, default=5)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--delay", type=float, default=1.0)
    ap.add_argument("--out", required=True)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--population", choices=["targetscan", "gmt-union"],
                    default="targetscan",
                    help="gene population the random lists are drawn from: "
                    "the TargetScan universe (13,035 human genes) or the "
                    "union of the local GMT libraries (17,436 genes, the "
                    "whole-genome proxy)")
    args = ap.parse_args()

    ke.TARGETSCAN_DB = ke.load_local_data()
    ke.load_gmt_libraries()
    if args.population == "gmt-union":
        universe = sorted(set().union(*[set(v) for lib in ke.GMT_LIBRARIES.values() for v in lib.values()]))
    else:
        universe = sorted(set().union(*ke.TARGETSCAN_DB.values()))
    print(f"population: {args.population} | universe: {len(universe)} genes")

    done = []
    if args.resume and os.path.exists(args.out):
        done = json.load(open(args.out)).get("lists", [])
        print(f"resuming: {len(done)} lists already done")
    rng = random.Random(args.seed)
    # Advance the RNG past the already-done lists so resume reproduces
    # the same sequence.
    lists = [sorted(rng.sample(universe, args.size))
             for _ in range(args.panels)]
    todo = lists[len(done):]

    async with make_async_client(timeout=20.0) as client:
        for i, genes in enumerate(todo, start=len(done)):
            res = await enrich_one_list(client, genes, args.delay)
            done.append({"index": i, "genes": genes,
                         "significant": res})
            if res is None:
                print(f"[{i+1}/{args.panels}] FAILED (no userListId)")
            else:
                sig = {k: v["adj_top8"] for k, v in res.items()
                       if v["adj_top8"] > 0}
                print(f"[{i+1}/{args.panels}] {sig or 'no significant term (engine criterion)'}")
            json.dump({"version": "1.54", "route": "live Enrichr",
                       "criterion": ("engine: first 8 rows per library, "
                                     "row[6] adjusted p < 0.05"),
                       "panels": args.panels, "size": args.size,
                       "seed": args.seed, "lists": done},
                      open(args.out, "w"), indent=2)
            await asyncio.sleep(args.delay)

    n = len(done)
    fpr = {}
    for label in SOURCES.values():
        hits = sum(1 for d in done
                   if d["significant"]
                   and d["significant"].get(label, {}).get("adj_top8", 0) > 0)
        fpr[label] = round(hits / n, 4) if n else None
    fpr_any = round(sum(1 for d in done
                        if d["significant"]
                        and any(v.get("adj_top8", 0) > 0
                                for v in d["significant"].values())) / n, 4) if n else None
    # Control: the same rate with the raw p-value (row[2]) instead.
    fpr_raw = {}
    for label in SOURCES.values():
        hits = sum(1 for d in done
                   if d["significant"]
                   and d["significant"].get(label, {}).get("raw_top8", 0) > 0)
        fpr_raw[label] = round(hits / n, 4) if n else None
    print(f"\nfalse-positive rate, engine criterion (n={n}, size={args.size}): {fpr} (any library: {fpr_any})")
    print(f"control, raw p-value (row[2]) top 8: {fpr_raw}")
    summary = {"version": "1.54", "route": "live Enrichr",
               "criterion": ("engine: first 8 rows per library, "
                             "row[6] adjusted p < 0.05"),
               "population": args.population,
               "population_size": len(universe),
               "panels": args.panels, "size": args.size, "seed": args.seed,
               "false_positive_rate": fpr, "false_positive_rate_any": fpr_any,
               "control_raw_p_fpr": fpr_raw, "lists": done}
    json.dump(summary, open(args.out, "w"), indent=2)
    print(f"saved to {args.out}")


if __name__ == "__main__":
    asyncio.run(main())
