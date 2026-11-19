#!/usr/bin/env python3
"""
KENRYU — protocolo de validación a nivel de panel (control negativo + positivo).

Pregunta que responde: ¿el pipeline "paneles de miRNAs -> genes convergentes
-> enriquecimiento" inventa señal?

  Control NEGATIVO: paneles aleatorios de K familias de miRNAs (sin relación
  biológica). Si el método está bien calibrado, la fracción de paneles con
  algún término significativo (p ajustado < ALPHA) debería ser pequeña
  (cercana a ALPHA). Se mide con dos universos de fondo:
    A) actual:   unión de los genes de la biblioteca GMT
    B) propuesto: genes con al menos un sitio predicho por TargetScan
       (los genes diana de miRNAs no son una muestra aleatoria del genoma)
  Control POSITIVO: familias con biología conocida; se comprueba si entre
  los términos significativos aparecen palabras clave esperadas. Es una
  comprobación de sensatez, NO una prueba de validez.

Uso (desde la raíz del repo KENRYU, con las dependencias instaladas):
    python validate_random_panels.py --panels 300 --size 5 --mode n-1 --seed 7
Salida: validation_random_panels.json + resumen en pantalla.
"""
import argparse, json, logging, random, sys
from pathlib import Path

logging.disable(logging.CRITICAL)
sys.path.insert(0, ".")
import kenryu_engine as E
import kenryu_core as kc

LIBS = {
    "GO_BP": "local_db/gmt_GO_Biological_Process_2023.gmt",
    "Reactome": "local_db/gmt_Reactome_Pathways_2024.gmt",
    "WikiPathways": "local_db/gmt_WikiPathways_2019_Human.gmt",
}
# Controles positivos: (familias, palabras clave esperadas en los términos)
POSITIVE = {
    "let-7 (RAS / ciclo celular)": (["let-7-5p"], ["cell cycle", "ras", "mapk", "differentiation"]),
    "miR-34 (p53 / apoptosis)": (["mir-34-5p"], ["apoptosis", "p53", "cell cycle", "senescence"]),
    "miR-17 (E2F / ciclo celular)": (["mir-17-5p"], ["cell cycle", "e2f", "tgf", "apoptosis"]),
}


def significant(genes, lib, background, alpha):
    res = kc.fisher_enrichment(set(genes), lib, background=background,
                               top_n=1000, p_adj_max=alpha)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--panels", type=int, default=300)
    ap.add_argument("--size", type=int, default=5)
    ap.add_argument("--mode", default="n-1", choices=["strict", "n-1", "n-2"])
    ap.add_argument("--min-genes", type=int, default=10)
    ap.add_argument("--max-genes", type=int, default=300)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default="validation_random_panels.json")
    a = ap.parse_args()

    db = E.load_local_data()                      # familia -> lista de genes
    fams = sorted(f for f, g in db.items() if len(g) >= 100)
    libs = {n: kc.load_gmt_library(p) for n, p in LIBS.items() if Path(p).exists()}
    ts_universe = set().union(*(set(g) for g in db.values()))
    bg_ts = {n: ts_universe & set().union(*lib.values()) for n, lib in libs.items()}
    rng = random.Random(a.seed)

    used, sizes = 0, []
    hit = {n: {"actual": 0, "targetscan": 0} for n in libs}
    tried = 0
    while used < a.panels and tried < a.panels * 20:
        tried += 1
        panel = rng.sample(fams, a.size)
        common = kc.consensus_genes([set(db[f]) for f in panel], a.mode)
        if not (a.min_genes <= len(common) <= a.max_genes):
            continue
        used += 1
        sizes.append(len(common))
        for n, lib in libs.items():
            if significant(common, lib, None, a.alpha):
                hit[n]["actual"] += 1
            if significant(common, lib, bg_ts[n], a.alpha):
                hit[n]["targetscan"] += 1

    out = {"params": vars(a), "panels_used": used, "panels_tried": tried,
           "genes_per_panel": {"min": min(sizes, default=None),
                               "median": sorted(sizes)[len(sizes) // 2] if sizes else None,
                               "max": max(sizes, default=None)},
           "false_positive_rate": {}, "positive_controls": {}}
    print(f"\nPaneles aleatorios analizados: {used} (de {tried} intentos), "
          f"modo {a.mode}, K={a.size}, alpha={a.alpha}")
    print(f"Genes por panel: {out['genes_per_panel']}")
    print("\nCONTROL NEGATIVO — fracción de paneles aleatorios con >=1 término significativo")
    print(f"{'biblioteca':14} {'fondo actual':>14} {'fondo TargetScan':>18}")
    for n in libs:
        fa = hit[n]["actual"] / used if used else float("nan")
        ft = hit[n]["targetscan"] / used if used else float("nan")
        out["false_positive_rate"][n] = {"actual": round(fa, 3), "targetscan": round(ft, 3)}
        print(f"{n:14} {fa:14.1%} {ft:18.1%}")
    print(f"(esperado si el método está calibrado: cercano a {a.alpha:.0%})")

    print("\nCONTROL POSITIVO — términos significativos con fondo TargetScan")
    for name, (famlist, kws) in POSITIVE.items():
        fl = [f for f in famlist if f in db]
        if not fl:
            print(f"  {name}: familia no encontrada en la base local"); continue
        genes = kc.consensus_genes([set(db[f]) for f in fl], "strict")
        found, terms = [], []
        for n, lib in libs.items():
            for r in significant(genes, lib, bg_ts[n], a.alpha):
                terms.append(f"{n}: {r['Term']} (p_adj={r['Pval']:.2g})")
        low = " | ".join(terms).lower()
        found = [k for k in kws if k in low]
        out["positive_controls"][name] = {"n_genes": len(genes), "n_terms": len(terms),
                                          "keywords_found": found, "top_terms": terms[:8]}
        print(f"  {name}: {len(genes)} genes, {len(terms)} términos significativos, "
              f"palabras clave halladas: {found or 'ninguna'}")
        for t in terms[:3]:
            print(f"      - {t}")
    Path(a.out).write_text(json.dumps(out, ensure_ascii=False, indent=2))
    print(f"\nGuardado: {a.out}")


if __name__ == "__main__":
    main()
