"""Replica of the Enrichr adjusted-p pattern with the local GO library.

For the random 5-gene list ARL8B, COL4A2, IFITM3, MAN2A2, SPATA2
(first list of the K5 calibration run, seed 7), submitted to Enrichr
on 2026-11-19 (GO Biological Process 2023), the top 8 terms all had
a single-gene overlap and the identical adjusted p-value 0.02515.

This script recomputes the same list against the repository's local
GO library (5,407 terms, 14,698-gene universe) with a plain Fisher
(hypergeometric) test and Benjamini-Hochberg, two ways:

  * BH over ALL 5,407 terms  -> adjusted p = 1.0 for the top 8.
  * BH over only the 48 terms with at least one overlapping gene
    -> adjusted p = 0.0342, identical across the top 8, same order
    as Enrichr's except the 8th term.

Run from the repository root: python docs/validation/replica_bh_overlap.py
Requires the local GMT libraries (load_gmt_libraries) and scipy.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import numpy as np  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

import kenryu_engine as ke  # noqa: E402

GENES = ["ARL8B", "COL4A2", "IFITM3", "MAN2A2", "SPATA2"]


def benjamini_hochberg(ps):
    ps = np.array(ps)
    m = len(ps)
    order = np.argsort(ps)
    ranked = ps[order]
    adj = ranked * m / (np.arange(m) + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(adj, 1.0)
    return out


def main():
    ke.TARGETSCAN_DB = ke.load_local_data()
    ke.load_gmt_libraries()
    go = ke.GMT_LIBRARIES["GO"]
    universe = set().union(*go.values())
    print(f"GO terms: {len(go)} | universe genes: {len(universe)}")
    gene_set = set(GENES)
    n_total, n_genes = len(universe), len(gene_set)
    rows = []
    for term, tgenes in go.items():
        overlap = len(gene_set & set(tgenes))
        if overlap >= 1:
            k_term = len(set(tgenes) & universe)
            p = hypergeom.sf(overlap - 1, n_total, k_term, n_genes)
            rows.append((term, overlap, k_term, p))
    rows.sort(key=lambda r: r[3])
    print(f"terms with overlap >= 1: {len(rows)}")
    adj_all = benjamini_hochberg(
        [r[3] for r in rows] + [1.0] * (len(go) - len(rows))
    )[: len(rows)]
    adj_overlap = benjamini_hochberg([r[3] for r in rows])
    print("\nTop 8 (raw p, BH over all, BH over overlap):")
    for i, (term, overlap, k_term, p) in enumerate(rows[:8]):
        print(
            f"{i + 1}. {term} | overlap={overlap} K={k_term} | "
            f"raw={p:.6f} | BH-all={adj_all[i]:.4f} | "
            f"BH-overlap={adj_overlap[i]:.4f}"
        )


if __name__ == "__main__":
    main()
