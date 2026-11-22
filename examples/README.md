# Reviewer path — the three reference panels

Run one command, no account, no API key:

```bash
pip install -r requirements.txt
pip install -e .
python examples/verify_examples.py
```

It runs the three reference panels through the real CLI, offline, and
checks the convergent genes (and, when an index is loaded, every
gene's miRTarBase evidence level) against the expected results in
`expected/`. It prints which miRTarBase version each panel used, from
the execution manifest, and exits 0 only if every panel passes.

## The panels

| Panel | miRNAs | Genes (9.0) |
|---|---|---|
| `lipid` — the thesis panel | hsa-miR-33a-5p, hsa-miR-33b-5p, hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p | 5: ABCA1, SNTB2, TSC22D2, KPNA3, SCN1A |
| `cancer` — independent case study | hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p | 19 |
| `control` | hsa-miR-375, hsa-miR-29a-3p, hsa-miR-21-5p | 6 |

## What a reviewer can verify, and what they cannot

The consensus genes come from TargetScan, which ships in this repo
(`data/targetscan_full.txt.zip` + the three GMT libraries), so they
are reproducible from a fresh clone with no extra download. The
evidence tiers (strong / weak / none) come from the miRTarBase
official index, which is runtime data: it is never committed, and the
official site currently returns 404, so a fresh clone cannot download
it.

- **With the 9.0 index loaded** (the version the public service runs):
  genes and evidence levels are checked against `expected/`.
- **With the 10.0 index loaded**: the script refuses to compare
  against the 9.0 expectations and tells you to re-record them with
  `python examples/verify_examples.py --record` (recording needs an
  index; it refuses in reduced mode).
- **Without an index** (reduced mode): the genes are still checked and
  must match; the evidence levels are only checked to be
  legacy/none (Harmonizome GMT snapshot, downloaded on demand as the
  fallback — note that `--offline` does not stop this one download,
  because it is part of loading the local data, not a live layer) or
  all "none" (truly offline, no GMT). The script says which of these
  it is, out loud, for every panel.

That is the honest split: **the genes are reproducible without the
index; the evidence tiers are not**, until the index can be downloaded
again. See `docs/validation/strong-controls.md` ("Reproducibility:
which numbers a stranger can verify") and
`docs/validation/mirtarbase-v10-adoption.md` for the full story.

## Files

- `panels.json` — the three panels (miRNAs, mode, expected file).
- `expected/{lipid,cancer,control}.json` — expected genes + evidence
  levels, recorded with engine 2.0 and the miRTarBase 9.0 index.
- `verify_examples.py` — the verification command (`--record` to
  regenerate the expected files).
