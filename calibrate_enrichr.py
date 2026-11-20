"""Calibrate the LIVE Enrichr route with random gene lists.

The paper's enrichment p-values (thesis and cancer panels) come from
Enrichr's own background, but the random-panel calibration in
docs/validation/ measured the LOCAL route, not Enrichr's. This script
measures the live route directly: for N random gene lists of a given
size, drawn from the TargetScan universe, it submits each list to
Enrichr and counts the fraction with at least one significant term
(adjusted p < 0.05) per library. A calibrated method should give ~5%.

Request limits: one addList POST + four enrich GETs per list, with a
delay between requests and backoff on 429. Results are saved after
every list, so an interrupted run resumes with --resume instead of
starting over. Default: 100 lists of 5 genes (thesis-like), seed 7.

Difference from the app's live path: the app keeps the top 8 terms
per library; here ALL rows with adjusted p < 0.05 count, which is the
method's true false-positive rate (slightly more conservative).

Run from the repo root with the test venv (takes ~10-15 min for 100
lists; run it in the background):
    ~/workspace/defensa-kenryu/.venv-test/bin/python \\
        calibrate_enrichr.py --panels 100 --size 5 --seed 7 \\
        --out ~/workspace/defensa-kenryu/validation-panels/calibration-enrichr-K5.json
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
    """Submit one gene list to Enrichr; return {library: n_significant}."""
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
        n_sig = 0
        for attempt in range(3):
            try:
                r = await client.get(
                    "https://maayanlab.cloud/Enrichr/enrich",
                    params={"userListId": user_list_id,
                            "backgroundType": gs_full},
                    timeout=20.0)
                if r.status_code == 200:
                    rows = r.json().get(gs_full, [])
                    n_sig = sum(1 for row in rows if row[6] < 0.05)
                    break
                if r.status_code == 429:
                    await asyncio.sleep(2 ** attempt * 5)
                    continue
                break
            except Exception:
                await asyncio.sleep(2 ** attempt * 5)
        out[label] = n_sig
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
    args = ap.parse_args()

    ke.TARGETSCAN_DB = ke.load_local_data()
    universe = sorted(set().union(*ke.TARGETSCAN_DB.values()))
    print(f"universe: {len(universe)} human genes")

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
                sig = {k: v for k, v in res.items() if v > 0}
                print(f"[{i+1}/{args.panels}] {sig or 'no significant term'}")
            json.dump({"version": "1.54", "route": "live Enrichr",
                       "panels": args.panels, "size": args.size,
                       "seed": args.seed, "lists": done},
                      open(args.out, "w"), indent=2)
            await asyncio.sleep(args.delay)

    n = len(done)
    fpr = {}
    for label in SOURCES.values():
        hits = sum(1 for d in done
                   if d["significant"] and d["significant"].get(label, 0) > 0)
        fpr[label] = round(hits / n, 4) if n else None
    print(f"\nfalse-positive rate (n={n}, size={args.size}): {fpr}")
    summary = {"version": "1.54", "route": "live Enrichr",
               "panels": args.panels, "size": args.size, "seed": args.seed,
               "false_positive_rate": fpr, "lists": done}
    json.dump(summary, open(args.out, "w"), indent=2)
    print(f"saved to {args.out}")


if __name__ == "__main__":
    asyncio.run(main())
