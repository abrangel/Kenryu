#!/usr/bin/env python3
"""KENRYU command line: run a miRNA panel without a web server.

The panel runs through the real engine in-process (TestClient), so the
CLI is the same pipeline as the API — local TargetScan, local miRTarBase,
local enrichment when offline — plus a JSON file, a Markdown report and
the PNG plots.

Usage:
    kenryu panel.txt --offline --out results/
    kenryu panel.txt --mode n-1 --years 5

Panel file: one miRNA per line (comments with # and blank lines ignored),
or a single line with comma-separated miRNAs.
"""

import argparse
import base64
import json
import sys
from pathlib import Path


def read_panel(path):
    """Read a panel file: one miRNA per line, or comma-separated."""
    text = Path(path).read_text(encoding="utf-8")
    mirnas = []
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        mirnas.extend(m.strip() for m in line.split(",") if m.strip())
    if not mirnas:
        raise ValueError(f"{path}: no miRNAs found")
    return mirnas


def write_report(result, out_dir):
    """Markdown summary of one analysis (panel, genes, validation,
    enrichment). Pure text, no network."""
    out_dir = Path(out_dir)
    lines = [
        "# KENRYU — panel report",
        "",
        f"Panel: {', '.join(result.get('found_mirnas', []))}",
        f"Genes convergentes: {len(result.get('common_genes', []))}",
        f"Modo offline: {result.get('offline', False)}",
        "",
    ]
    manifest = result.get("execution_manifest", {})
    if result.get("arm_unverified"):
        lines.append("## Brazos no modelados por TargetScan")
        lines.append("")
        for e in result["arm_unverified"]:
            lines.append(f"- {e['mirna']}: {e['reason']} "
                         f"(familia {e.get('family_key')})")
        lines.append("Sus dianas predichas se bloquearon; no se sustituyeron "
                     "por las del otro brazo.")
        lines.append("")
    if manifest:
        params = manifest.get("parameters", {})
        lines.append("## Manifiesto de ejecución")
        lines.append("")
        lines.append(f"Versión: {manifest.get('engine_version')} "
                     f"(commit {manifest.get('commit')})")
        lines.append(f"Parámetros: modo {params.get('mode')}, "
                     f"offline {params.get('offline')}, "
                     f"especie {params.get('species')}")
        lines.append(f"Semilla de convergencia: "
                     f"{manifest.get('seeds', {}).get('convergence')}")
        lines.append(f"Reproducir: `{manifest.get('reproduce', '')}`")
        lines.append("")
    summary = result.get("validation", {}).get("summary", {})
    if summary:
        lines.append(
            "Validación: "
            f"{summary.get('high', 0)} alta, {summary.get('medium', 0)} media, "
            f"{summary.get('low', 0)} baja "
            f"(media {summary.get('mean_score', 0)}). "
            f"Sin validar: {', '.join(summary.get('not_validated', [])) or 'ninguno'}")
        lines.append("")
    lines.append("## Genes")
    lines.append("")
    lines.append("| Gen | Puntaje | Nivel | Capas |")
    lines.append("| --- | --- | --- | --- |")
    for gene, v in result.get("validation", {}).get("genes", {}).items():
        layers = "; ".join(v.get("layers", []))
        lines.append(f"| {gene} | {v.get('score')} | {v.get('level')} | {layers} |")
    lines.append("")
    lines.append("## Enriquecimiento (top 10)")
    lines.append("")
    for e in result.get("enrichment", [])[:10]:
        lines.append(f"- {e.get('Term')} ({e.get('Source')}, p-adj {e.get('Pval'):.4g})")
    lines.append("")
    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="kenryu",
        description="KENRYU: miRNA target convergence and evidence, on your machine.")
    parser.add_argument("panel", help="Panel file (one miRNA per line).")
    parser.add_argument("--mode", default="strict",
                        choices=["strict", "n-1", "n-2"])
    parser.add_argument("--years", type=int, default=10)
    parser.add_argument("--offline", action="store_true",
                        help="No network at all: local data only.")
    parser.add_argument("--out", default="kenryu-output",
                        help="Output directory (default: kenryu-output).")
    args = parser.parse_args(argv)

    mirnas = read_panel(args.panel)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    # The real engine, in-process: no server, no port, same pipeline.
    from fastapi.testclient import TestClient
    import kenryu_engine as ke

    with TestClient(ke.app) as client:
        r = client.post("/api/v1/analyze", json={
            "mirnas": mirnas, "mode": args.mode,
            "years": args.years, "offline": args.offline})
    if r.status_code != 200:
        print(f"kenryu: analysis failed ({r.status_code}): {r.text[:200]}",
              file=sys.stderr)
        return 1
    result = r.json()

    (out_dir / "analysis.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    write_report(result, out_dir)
    for key, name in (("venn_plot", "venn.png"), ("volcano_plot", "volcano.png"),
                      ("ppi_plot", "ppi.png"), ("heatmap_plot", "heatmap.png")):
        if result.get(key):
            (out_dir / name).write_bytes(base64.b64decode(result[key]))

    print(f"kenryu: {len(result.get('common_genes', []))} genes "
          f"from {len(result.get('found_mirnas', []))} miRNAs -> {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
