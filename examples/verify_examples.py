"""Reviewer verification command for the three reference panels.

    python examples/verify_examples.py            # verify against expected/
    python examples/verify_examples.py --record   # regenerate expected/ (needs an index)

What it checks, per panel: the CLI runs offline on the panel's miRNAs
and the convergent genes match the expected list exactly (order-free),
and every gene's miRTarBase evidence level matches the expected level.

Which miRTarBase version it uses is printed for every panel, from the
execution manifest — "9.0", "10.0", "legacy" (Harmonizome GMT snapshot)
or "not loaded". The expected levels in expected/ were recorded with
the 9.0 index, the version the public service runs. Without an index
(a fresh clone: the index is runtime data, never committed, and the
site currently returns 404) the engine runs in reduced mode: the
consensus genes are the same — they come from TargetScan, which ships
in the repo — but the evidence tiers are not. With the Harmonizome GMT
snapshot (downloaded on demand as the fallback) a pair present in the
snapshot is "legacy" and everything else is "none"; with no GMT at all
(truly offline) every gene is "none". The script says which of these
it is out loud instead of comparing against the 9.0 levels. That is
the honest reviewer path: the genes are reproducible without the
index; the evidence tiers are not, until the index can be downloaded
again. To verify against the 10.0 index instead, re-record the
expected files with --record while it is loaded.

Exit code 0 only if every panel passes.
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
EXAMPLES = Path(__file__).resolve().parent

PANELS = json.loads((EXAMPLES / "panels.json").read_text())["panels"]


def classify_mirtarbase(version_string):
    """'miRTarBase 9.0 ...' -> '9.0'; '... 10.0 ...' -> '10.0';
    'Harmonizome GMT snapshot (legacy)' -> 'legacy'; 'not loaded' -> None."""
    if not version_string:
        return None
    if version_string.startswith("miRTarBase 9.0"):
        return "9.0"
    if version_string.startswith("miRTarBase 10.0"):
        return "10.0"
    if version_string.startswith("Harmonizome GMT snapshot"):
        return "legacy"
    if version_string == "not loaded":
        return None
    return version_string


def run_cli(panel, out_dir):
    panel_file = out_dir / "panel.txt"
    panel_file.write_text("\n".join(panel["mirnas"]) + "\n", encoding="utf-8")
    cmd = None
    if shutil.which("kenryu"):
        cmd = ["kenryu", str(panel_file), "--offline", "--out", str(out_dir)]
    else:
        cmd = [sys.executable, "-m", "kenryu_cli", str(panel_file),
               "--offline", "--out", str(out_dir)]
    proc = subprocess.run(cmd, cwd=REPO_ROOT, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout[-2000:])
        print(proc.stderr[-2000:])
        raise SystemExit(f"CLI failed for panel {panel['id']} (exit {proc.returncode})")
    return json.loads((out_dir / "analysis.json").read_text())


def check_panel(panel):
    with tempfile.TemporaryDirectory(prefix=f"kenryu-{panel['id']}-") as td:
        result = run_cli(panel, Path(td))
    version = classify_mirtarbase(
        result.get("execution_manifest", {}).get("data_versions", {}).get("mirtarbase")
    )
    genes = sorted(result.get("common_genes", []))
    levels = {g: v.get("evidence_level")
              for g, v in result.get("validation", {}).get("genes", {}).items()}
    expected = json.loads((EXAMPLES / panel["expected"]).read_text())
    problems = []
    if genes != sorted(expected["common_genes"]):
        problems.append(f"genes differ: got {genes}, expected {sorted(expected['common_genes'])}")
    if version is None:
        # Reduced mode, no GMT at all: the genes must still match, but
        # every level is "none" and is NOT compared against the 9.0 file.
        if any(lv != "none" for lv in levels.values()):
            problems.append(f"reduced mode but levels are not all 'none': {levels}")
        note = ("miRTarBase: NOT LOADED — reduced mode: consensus genes verified, "
                "evidence levels all 'none' (no index, no GMT snapshot; the 9.0 "
                f"levels in {panel['expected']} do not apply)")
    elif version == "legacy":
        # Reduced mode, Harmonizome GMT snapshot: a pair in the snapshot
        # is "legacy", the rest "none"; never strong/weak.
        if any(lv not in ("legacy", "none") for lv in levels.values()):
            problems.append(f"legacy mode but a level is not legacy/none: {levels}")
        note = ("miRTarBase: LEGACY (Harmonizome GMT snapshot) — reduced mode: "
                "consensus genes verified, evidence levels legacy/none only "
                f"(the 9.0 levels in {panel['expected']} do not apply)")
    else:
        if version != expected["mirtarbase_version"]:
            problems.append(
                f"expected files are miRTarBase {expected['mirtarbase_version']} "
                f"but the run used {version} — re-record with --record to verify "
                f"against {version}")
        elif levels != expected["evidence"]:
            problems.append(f"evidence levels differ: got {levels}, expected {expected['evidence']}")
        note = f"miRTarBase: {version}; evidence levels checked"
    status = "OK" if not problems else "FAIL"
    print(f"[{panel['id']}] {status}: genes ({len(genes)}); {note}")
    return problems


def record_panel(panel):
    with tempfile.TemporaryDirectory(prefix=f"kenryu-{panel['id']}-") as td:
        result = run_cli(panel, Path(td))
    version = classify_mirtarbase(
        result.get("execution_manifest", {}).get("data_versions", {}).get("mirtarbase")
    )
    if version not in ("9.0", "10.0"):
        raise SystemExit(
            f"cannot record {panel['id']}: miRTarBase is {version or 'not loaded'} — "
            f"recording in reduced mode would bake in legacy/none levels"
        )
    levels = {g: v.get("evidence_level")
              for g, v in result.get("validation", {}).get("genes", {}).items()}
    expected = {
        "panel": panel["id"],
        "name": panel["name"],
        "mirnas": panel["mirnas"],
        "mode": panel["mode"],
        "engine_version": result.get("engine_version"),
        "mirtarbase_version": version,
        "mirtarbase_detail": result.get("execution_manifest", {}).get("data_versions", {}).get("mirtarbase"),
        "common_genes": sorted(result.get("common_genes", [])),
        "evidence": levels,
    }
    path = EXAMPLES / panel["expected"]
    path.parent.mkdir(exist_ok=True)
    path.write_text(json.dumps(expected, indent=2, ensure_ascii=False) + "\n")
    print(f"[{panel['id']}] recorded {path} (miRTarBase {version}, "
          f"{len(expected['common_genes'])} genes)")


def main():
    record = "--record" in sys.argv
    all_problems = []
    for panel in PANELS:
        if record:
            record_panel(panel)
        else:
            all_problems += [f"{panel['id']}: {p}" for p in check_panel(panel)]
    if record:
        print("Recorded. Re-run without --record to verify.")
        return 0
    if all_problems:
        print("\nFAIL:")
        for p in all_problems:
            print(" -", p)
        return 1
    print("\nALL PANELS PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
