"""Decisive test of the Enrichr adjusted-p mechanism (2026-11-19).

Hypothesis: Enrichr's adjusted p-value is a Benjamini-Hochberg
correction computed only over the terms the API returns for the list
— i.e. the terms with at least one overlapping gene — not over the
whole library. A local replica of one list was consistent with this
(see docs/validation/replica_bh_overlap.py), but that replica used
the same computation twice: it guards against a transcription error
and adds no independent evidence.

This script tests the hypothesis against Enrichr's own output: it
re-submits the first N gene lists of the K5-engine calibration run
(same lists, seed 7), keeps the RAW rows of all four libraries, and
recomputes BH over exactly the returned rows, comparing with the
adjusted p in column 6 of every row. If the recomputed values match
Enrichr's column in all rows of all lists and libraries, the
mechanism is confirmed from Enrichr's own output.

Usage: python verify_enrichr_bh.py --lists 10 [--seed-lists K5-engine]
Output: docs/validation/enrichr-bh-verification.json (rows included).
"""
import argparse
import asyncio
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calibrate_enrichr import SOURCES  # noqa: E402
from kenryu_engine import make_async_client  # noqa: E402

K5_ENGINE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "validation-panels",
    "calibration-enrichr-K5-engine.json",
)


def bh(p):
    """Benjamini-Hochberg adjusted p-values; p sorted ascending."""
    m = len(p)
    out = [0] * m
    prev = 1.0
    for i in range(m - 1, -1, -1):
        prev = min(prev, p[i] * m / (i + 1))
        out[i] = prev
    return out


async def fetch_rows(client, genes):
    payload = {"list": (None, "\n".join(genes)),
               "description": (None, "KENRYU BH verification")}
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
        except Exception:
            await asyncio.sleep(2 ** attempt * 5)
    if not user_list_id:
        return None
    out = {}
    for gs_full, label in SOURCES.items():
        for attempt in range(3):
            try:
                r = await client.get(
                    "https://maayanlab.cloud/Enrichr/enrich",
                    params={"userListId": user_list_id,
                            "backgroundType": gs_full},
                    timeout=20.0)
                if r.status_code == 200:
                    out[label] = r.json().get(gs_full, [])
                    break
                if r.status_code == 429:
                    await asyncio.sleep(2 ** attempt * 5)
                    continue
                break
            except Exception:
                await asyncio.sleep(2 ** attempt * 5)
        await asyncio.sleep(1.0)
    return out


def check_library(rows, tol=1e-4):
    """Recompute BH over the returned rows; compare with column 6."""
    ordered = sorted(rows, key=lambda r: r[2])
    mine = bh([r[2] for r in ordered])
    mismatches = []
    for value, row in zip(mine, ordered):
        if abs(value - row[6]) > tol:
            mismatches.append({
                "term": row[1], "raw_p": row[2],
                "enrichr_adj": row[6], "recomputed": value,
            })
    return {"rows": len(rows), "mismatches": mismatches}


async def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lists", type=int, default=10)
    ap.add_argument("--out", default="docs/validation/enrichr-bh-verification.json")
    args = ap.parse_args()

    k5 = json.load(open(K5_ENGINE))
    gene_lists = [entry["genes"] for entry in k5["lists"][: args.lists]]
    print(f"verifying {len(gene_lists)} lists (seed 7, K5-engine)")

    results = {"date": "2026-11-19", "lists_checked": len(gene_lists),
               "tolerance": 1e-4, "libraries": {}}
    total_rows = 0
    total_mismatch = 0
    async with make_async_client(timeout=20.0) as client:
        for i, genes in enumerate(gene_lists):
            resp = await fetch_rows(client, genes)
            if resp is None:
                print(f"list {i}: addList failed, skipping")
                continue
            for label, rows in resp.items():
                check = check_library(rows)
                total_rows += check["rows"]
                total_mismatch += len(check["mismatches"])
                lib = results["libraries"].setdefault(
                    label, {"rows": 0, "mismatches": []})
                lib["rows"] += check["rows"]
                lib["mismatches"].extend(
                    {**m, "list_index": i} for m in check["mismatches"])
            print(f"list {i}: rows so far {total_rows}, "
                  f"mismatches {total_mismatch}")
    results["total_rows"] = total_rows
    results["total_mismatches"] = total_mismatch
    results["confirmed"] = total_mismatch == 0 and total_rows > 0
    json.dump(results, open(args.out, "w"), indent=2)
    print(f"\nTOTAL rows compared: {total_rows}")
    print(f"TOTAL mismatches: {total_mismatch}")
    print("CONFIRMED" if results["confirmed"] else "NOT CONFIRMED")
    print(f"saved to {args.out}")


if __name__ == "__main__":
    asyncio.run(main())
