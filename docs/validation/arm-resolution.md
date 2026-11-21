# Arm resolution (-5p/-3p) — audit and fix (v1.58)

## The audit (2026-11-21, on v1.56)

`scripts/audit_arm_resolution.py` tests every mature human miRNA in
miRBase 23.1 that carries an arm (332 names with local TargetScan
targets) and asks whether the TargetScan family that answered models
that arm.

- **33 names from 19 hairpins** resolve to a family whose key carries
  no arm (`docs/validation/arm_swaps.json`, the v1.56 record).
- **In 14 hairpins (28 names) the -5p and -3p arms return exactly the
  same gene set.** With different seeds that is impossible in a seed
  predictor: at least one of the two results is necessarily wrong.
  Examples: miR-107, miR-133b, miR-137, miR-147b, miR-184, miR-217.
- In 5 more hairpins (miR-1193, 1197, 206, 448, 802) only one arm is
  affected.

The name alone cannot tell which arm TargetScan modelled — that is the
limit of the name-based audit, and why the fix uses the seed.

## The fix: verify by seed

TargetScan families are defined by their seed. The fix (v1.58,
`kenryu_core.verify_arm_seed`, `kenryu_engine.check_arm_resolution`)
compares two seeds:

- the **family seed**, from TargetScan's `miR_Family_Info.txt`
  (column "Seed+m8", human rows; `local_db/mir_family_info.txt`,
  downloaded by the updater with TargetScan 8.0);
- the **query seed**, from the query's mature sequence in miRBase 23.1
  (`local_db/mirbase_mature.fa`).

Verified when the family seed occurs in the mature sequence starting
at offset 0-3 (`ARM_SEED_MAX_OFFSET = 3`). The offset window tolerates
the 1-2 nt of 5' annotation drift between the miRBase version
TargetScan 8.0 was built on (≤22.1) and miRBase 23.1:

- miR-144-3p: miRBase 23.1 added 2 nt at the 5' end
  (`ACUACAGUAUAGAUGAUGUACU` vs TargetScan's `UACAGUAUAGAUGAUGUACU`),
  so the family seed `ACAGUAU` sits at offset 3. Still the same arm.
- miR-802-5p: miRBase 23.1 dropped the 5' U, so the family seed
  `CAGUAAC` sits at offset 0. Still the same arm — a naive seed
  comparison (positions 2-8 only) would have wrongly blocked it.

The wrong arm's seed is absent from its sequence entirely, so it is
rejected at any offset. A seed starting at offset 4 or later is not
annotation drift any more and does not verify.

## What changes for users

- **15 arms blocked** (one per hairpin above, plus miR-1193-5p): their
  predicted targets are returned empty, with an explicit warning, and
  the miRNA is listed in the analyze response under `arm_unverified`
  (with the family key and the reason). They are never replaced by the
  other arm's targets. In 12 of the 14 hairpins the -3p arm is the
  modelled one; in miR-217 and miR-346 it is the -5p arm — the seed
  decides, not the convention.
- **A second bug class the seed check catches: 14 letter-drop
  collapses.** The `m_no_letter` fallback drops the a/b/c letter and
  merges distinct paralogs into the wrong family (e.g. hsa-miR-519e-3p
  into the miR-17 family, 1,385 genes; hsa-miR-122b-5p into miR-122;
  hsa-miR-451b-3p/-5p into miR-451). Their seeds do not match, so they
  are blocked too. Total blocked: 29 (`docs/validation/arm_blocked.json`).
- **303 armed names keep their targets.** Every kept result carries
  the family's seed (the audit re-run with the seed criterion reports
  0 swaps; `tests/test_arm_resolution.py` is the CI regression).
- **The thesis panel does not change**: all five lipid miRNAs resolve
  to armed family keys (mir-33-5p, mir-144-3p, mir-758-3p, mir-106-5p)
  and keep their targets. The control panel's bare `hsa-miR-375` is
  not checked (see below).
- **The same swap lived in the miRTarBase lookup.** `get_validated_*`
  fell back from an armed name to the bare name (e.g. hsa-miR-107-5p
  received the -3p arm's 298 validated pairs). The engine now guards
  that fallback the same way (`get_validated_*_checked`); the live
  Harmonizome fetch is guarded too. When the bare key is not a
  TargetScan family at all (268 bare keys in the v9.0 index — viral,
  removed or renamed miRNAs), the seed cannot be checked: the
  experimental pairs are kept, with a warning, because experimental
  evidence is not dropped for an unanswerable question.

## Response fields

`/api/v1/analyze` (both the empty-intersection and the normal path)
now returns:

- `"arm_unverified"`: list of `{"mirna", "reason", "family_key"}` for
  blocked arms (empty list when none).
- `"arm_verification"`: `"active"` when both seed files loaded,
  `"unavailable"` when one is missing. When unavailable the engine
  returns the pre-fix behaviour (targets as before) and logs a warning
  — fail-open, loudly, never silently. The CLI report and the web UI
  show the same warning.

## Limits

- A **bare name** (no -5p/-3p, e.g. `hsa-miR-375`) is not checked: it
  resolves to the arm TargetScan modelled, the historical guide strand.
- **Mouse** (`mmu-`) is not checked: no arm data is loaded for it.
- The check covers the TargetScan family that answered. A query that
  resolves through the fuzzy last-resort step is checked against the
  matched key's seed like any other.
- 4 family keys carry two seeds in `miR_Family_Info.txt`
  (mir-520-3p, mir-487-3p, mir-500b-5p, mir-519-3p: a sub-part string
  shared by two families). The loader logs them; the check accepts any
  of the seeds. None of the four is an armless key, so the arm check is
  unaffected.
