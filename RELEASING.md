# Release policy

- `main` is always deployable: CI (tests + secret scan) must pass before merge.
- Fixes and changes land as normal commits. Each one adds a line under
  "## [Unreleased]" in CHANGELOG.md. No tag per fix.
- A tag/release is created only at a milestone:
  1. results users may already have seen change (scores, enrichment, defaults);
  2. a significant new capability ships;
  3. a version is frozen for citation (paper, Zenodo DOI).
- Versions follow SemVer: MINOR for new capabilities, PATCH for fixes
  grouped into a release. Published tags are never moved or deleted.
- Before tagging: full pytest, secret scan, and the calibration protocol
  (validate_random_panels.py) re-run if enrichment code changed. Record
  data versions and the run date in docs/validation/.
- The running service reports the version AND the short commit hash
  (/health and every analysis JSON), so any result can be traced to the
  exact code that produced it. The hash comes from
  `local_db/build_info.json`, written by the sync step with
  `git rev-parse --short HEAD` (the GitHub hash, never the Space's own
  commit); when the file is absent or malformed the service reports
  "unknown", never an invented value.

Exception: if a change alters results someone may already have used
(for example the "high" level gate or the enrichment background), it
does merit a release with its own note, even if it is small.

## Space sync (lesson from the v1.55 incident, 2026-11-19)

The v1.55 sync uploaded the engine but not `kenryu_core.py`, and the
Space went into "Runtime error": the v1.54 core did not define
`load_mirtarbase_official_with_meta`, which the v1.55 engine
imports. The engine and the core move together. Before every sync:

1. Build the sync package with `git ls-files` or `git archive`,
   never by picking files by hand — a hand-picked list forgets a
   module. Check first which modules changed since the last synced
   version: `git log <prev-tag>..HEAD -- kenryu_core.py
   kenryu_updater.py kenryu_engine.py static/`.
2. Before pushing to the Space, import the engine from the package
   itself: `cd <package> && python -c "import kenryu_engine"`. If
   the import fails there, it will fail on the Space.
3. After pushing, run an automatic `/health` check (version, commit
   hash, `mirtarbase_mode`) plus the two reference panels (lipid and
   cancer) before calling the sync done.
4. A separate test Space costs little and avoids breaking the public
   one. JOSS may open the service during review, so the public Space
   must never be left in "Runtime error" by a sync.
