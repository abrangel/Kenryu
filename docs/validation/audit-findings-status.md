# Audit findings — status (H-01…H-23)

Claude's independent audit (`~/workspace/user/files/AUDITORIA-KENRYU.md`,
v1.44, 2026-11-19), resolved in rounds 1–3 (v1.45–v1.57) and later
versions. Status verified against the v2.0 tree on 2026-11-22 by Muse,
finding by finding, with executable evidence where the claim allows it.
"Enforced" names the CI guard; where it says *none*, the finding is
closed by inspection and has no regression test — a gap, not a proof.

| ID | Sev | Finding (one line) | Status | Evidence (2026-11-22) | Enforced in CI |
|---|---|---|---|---|---|
| H-01 | Crítica | Update endpoint without auth overwrites server files | **Closed — test-verified** | `KENRYU_ADMIN_TOKEN` required (`kenryu_engine.py:2151-2159`); 403 without token, 403 with wrong token | `tests/test_engine.py:469-483`, `tests/test_updater.py:314` |
| H-02 | Crítica | KEGG + Harmonizome GMTs redistributed under MIT | **Closed — verified in git** | `git ls-files local_db/` has no KEGG/Harmonizome GMT; only GO/Reactome/WikiPathways (redistributable). Snapshot downloads on demand since v1.52 | None (no test guards the absence) |
| H-03 | Crítica | Fuzzy match returns wrong-family targets | **Closed — verified live** | `get_targets_local('hsa-miR-1-5p')` → 0 genes; `hsa-miR-33a-5p` → 545. Fuzzy search stays only as last resort, with a warning | `tests/test_arm_resolution.py` (12 tests, incl. the 332-name audit) |
| H-04 | Crítica | Paper claims "dual evidence" the code does not do | **Closed — verified in text** | No "obligatoria"/"dual" left in `paper.md` or `docs/case-study-cancer-panel.md`; framing is now three evidence levels (`paper.md:34`) | None (documentation claim) |
| H-05 | Alta | Fisher + BH underestimates m (~28× in GO), p_adj < 0.3 gives 40/40 false positives on null lists | **Closed — verified in code + calibration** | `local_enrichment`/`fisher_enrichment` use the TargetScan universe as background; calibration re-run in `docs/validation/` | Enrichment tests in the suite |
| H-06 | Alta | Alphabetical bias: gene_research only on the first 5 alphabetical genes | **Closed — verified** | `prioritize_genes` in core (validation-first order); researched set is the top 5 by evidence, not alphabet | `test_analyze_offline_cancer_panel_researches_prioritized_genes` |
| H-07 | Alta | Frontend sends the calendar year (2026) as `years` | **Closed — verified in code** | `static/script.js:81-85` computes `years` from the selected start year | None cited |
| H-08 | Alta | Reflected/DOM XSS via systematic `innerHTML` | **Closed — verified** | `escapeHtml` used 45×; remaining `innerHTML` are static literals (spinner icons) or integer-built year options | `tests/test_frontend_safety.py` (25 tests) |
| H-09 | Alta | Pinned deps with known CVEs (15 unique) | **Closed — verified executable** | `pip-audit -r requirements.txt` → "No known vulnerabilities found" (2026-11-22) | None (pip-audit is not a CI job) |
| H-10 | Alta | No size/rate/validation limits on /analyze, /compare, /batch | **Closed — partial** | `AnalysisRequest`: miRNAs 1–50, years 1–50, `mode` is a Literal (invalid → 422). No rate limiting | Pydantic validation tests |
| H-11 | Alta | `paper.bib` author/title/year/page errors; Zenodo DOI 404 | **Closed — verified** | `paper.bib` header records the 2026-11-19 corrections; `CITATION.cff` DOI 10.5281/zenodo.20564689 resolves | None (documentation claim) |
| H-12 | Alta | `build_synthesis` asserts causality, cites `[i+1]` with no referent | **Closed — verified in code** | `kenryu_core.py` (v1.47): "sin clinical_evidence no se afirma causalidad ni se cita un índice [i+1]" | None cited |
| H-13 | Media | UTF-8 BOM corrupts gene symbols in `data/hsa-miR-*.txt` | **Closed — verified** | Reads use `encoding='utf-8-sig'` (`kenryu_engine.py:456, 752`) | Suite |
| H-14 | Media | CORS trusts any `*.hf.space`; no security headers | **Closed — verified** | Explicit origins list (localhost, huggingface.co, SPACE_ORIGIN, abrangel.github.io), `allow_credentials=False`, security-headers middleware with CSP | None cited |
| H-15 | Media | Persistent cache: no lock, no atomic write, no size limit, overwritten by pytest | **Closed — verified** | `save_persistent_cache`: `_CACHE_LOCK` + tmp file + `os.replace` + per-section `_CACHE_LIMITS`; `tests/conftest.py` isolates `CACHE_FILE` to a tmp dir | conftest + suite |
| H-16 | Media | miRBase FASTAs unused; no HGNC alias / `C*orf*` normalisation | **Closed — verified** | HGNC map, 498 entries with the Approved-symbol guard (v1.60); collisions verified same-gene by Ensembl; builder refuses to write on a collision mismatch | `tests/test_gene_symbols.py` (5 tests) |
| H-17 | Media | Consensus edge cases: duplicates skew N-1/N-2, invalid `mode` silently falls to n-2, miRNAs without data dropped silently | **Closed — verified** | Dedup in core (`kenryu_core.py:166-168`); `mode` is a Literal (422); `missing_mirnas` reported | Consensus tests in the suite |
| H-18 | Media | Fixed "confidence" values (85/88/92/95) in `gene_details` and the UI | **Closed — verified** | No fixed 85/88/92/95 literals in the engine or the UI (grep clean) | None (a grep-based test could guard it) |
| H-19 | Media | `data_versions.json` stale; `targetscan_full.json.zip` is really a `.txt` | **Closed — verified** | `data_versions.json` updated in v2.0 (local miRTarBase 10.0); versions pinned in the manifest | Manifest/version tests |
| H-20 | Media | pytest fails on Windows (CRLF 161 vs 157 bytes) | **Closed — in code, not on real Windows** | `test_transform_html_fasta_writes_lf_bytes` simulates Windows `write_text`; the transform passes `newline="\n"`. Never run on a real Windows machine | `tests/test_updater.py:451-461` |
| H-21 | Media | YAML/Markdown injection in `exportMarkdown()`; fixed version `v1.38` | **Closed — verified** | Version comes from the backend (`lastData.engine_version`); miRNAs serialised with `JSON.stringify` (`static/script.js:770-795`) | None cited |
| H-22 | Baja | Dockerfile: `chmod 777 local_db`, `build-essential` in the final image, no digest pin | **Closed — verified** | Base image pinned by digest; `build-essential` purged after install; `chmod 755`; non-root user `kenryu` | None (Docker build is not in CI) |
| H-23 | Baja | Minor doc inconsistencies (SECURITY.md placeholder + "NutriKen", old licence in `pyproject.toml`, duplicated `index.html`, "18k+" targets) | **Closed — verified** | `pyproject.toml` licence MIT; `index.html` gone from the root; no "18k" claim in README/paper; `SECURITY.md` clean of NutriKen/placeholder (grep clean) | None (documentation claim) |

## Result

**23 of 23 closed.** 0 open with impact on results. Two are closed with a
stated limitation rather than a full fix:

- **H-10** — size and validation limits exist; there is no rate limit.
  Residual risk is low (the public Space sits behind Hugging Face's
  infrastructure), and it is documented here rather than hidden.
- **H-20** — fixed in code with a simulated-Windows test; never executed
  on a real Windows machine.

## Gaps this table found (no regression test yet)

- H-02: nothing in CI guards that a KEGG/Harmonizome GMT is never
  committed again.
- H-09: `pip-audit` is not a CI job; the 2026-11-22 clean result is a
  one-off check, not a guard.
- H-18: no grep-based test guards against fixed confidence literals
  returning.

## Round 3 pendings (still open, from PAPER-REVIEW §8/§16)

- **Cover letter (May–Jun 2027):** the resubmission letter must state
  the prior desk rejection and what changed. Not writable until the
  evidence exists.
- **Specificity draws:** the app still serves 100-draw values; the
  paper's table uses 10,000 draws (~4 s slower per analysis). Open
  decision, no impact on results.
- **Enrichr libraries milestone (next):** GO → 2026, WikiPathways →
  2024_Human, KEGG → 2026, repeating the calibration. Not started.
- **External-usage evidence (3–5 real researchers):** the author's
  action; still open. This is the resubmission bottleneck, not the
  workflow.
- **Clean-clone install on other machines/OSes:** CI runs Linux only;
  verified on this VM. Windows/macOS never tried.
- **examples/ + one verification command + a small included dataset:**
  does not exist yet (Claude, Week 1). The reviewer path for the
  miRTarBase evidence layer depends on this, because the official site
  returns 404 and the index is not distributed with the repo.
- **Email to the miRTarBase group** asking for the current download
  location and terms: not sent; needs the author's approval.

## Decisions recorded 2026-11-22

- **Space stays on miRTarBase 9.0.** The site returns 404 (no reliable
  source for a runtime download), the v10.0 terms are unverified, and
  the cost is real (20.7 s startup, 1.4 GB RAM). Reconsider when the
  site returns and the licence is legible. The v2.0 code already falls
  back to v9.0 on the Space; nothing changes.
- **Sensitivity denominators are already explicit.** `paper.md:43`
  ties 358/643 to "the official miRTarBase 10.0 index" in the same
  sentence; `docs/validation/strong-controls.md:89-93` contrasts
  358/643 (v10.0) with 404/754 (v9.0). No change needed.
- **"Software design" and "Research impact statement" are NOT
  missing.** Both sections exist in `paper.md:45,48`, substantive,
  since v1.57. The claim that they are missing was true on 2026-11-19
  (PAPER-REVIEW §2) and is outdated now.
