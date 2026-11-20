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
