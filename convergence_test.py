#!/usr/bin/env python3
"""
KENRYU — prueba de convergencia de un panel de miRNAs.

Pregunta: ¿cuántos genes en común produce este panel, y cuánto sorprende?
Compara el número observado con el de paneles aleatorios de familias de
miRNAs. Control de tamaño: cada familia observada se reemplaza por una
familia al azar con un número de dianas parecido (±25%), de modo que un
panel de familias "grandes" se compara con familias "grandes".
También informa del tamaño efectivo del panel (familias de semilla
distintas): miR-33a y miR-33b comparten familia y cuentan una sola vez.

Uso (desde la raíz del repo KENRYU):
  python convergence_test.py --mirnas hsa-miR-33a-5p hsa-miR-144-3p hsa-miR-758-3p \
         hsa-miR-106b-5p hsa-miR-33b-5p --mode strict --perms 10000 --seed 7
Salida: resumen en pantalla y convergence_<semilla>.json
"""
import argparse, json, logging, random, sys, statistics
logging.disable(logging.CRITICAL)
sys.path.insert(0, ".")
import kenryu_engine as E
import kenryu_core as kc


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mirnas", nargs="+", required=True)
    ap.add_argument("--mode", default="strict", choices=["strict", "n-1", "n-2"])
    ap.add_argument("--perms", type=int, default=10000)
    ap.add_argument("--tol", type=float, default=0.25, help="tolerancia de tamaño de dianas")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    db = E.load_local_data()
    fam_sets = {f: frozenset(g) for f, g in db.items() if g}
    obs_sets = {m: frozenset(E.get_targets_local(m, db)) for m in a.mirnas}
    missing = [m for m, s in obs_sets.items() if not s]
    if missing:
        sys.exit(f"Sin dianas en TargetScan local para: {missing}")
    distinct = list(dict.fromkeys(obs_sets.values()))      # familias efectivas
    k_eff = len(distinct)
    c_obs = len(kc.consensus_genes([set(s) for s in distinct], a.mode))
    c_obs_raw = len(kc.consensus_genes([set(s) for s in obs_sets.values()], a.mode))

    # Candidatos por familia observada, emparejados por tamaño
    observed_fams = {f for f, s in fam_sets.items() if s in distinct}
    pools = []
    for s in distinct:
        n = len(s)
        pool = [f for f, t in fam_sets.items()
                if f not in observed_fams and (1 - a.tol) * n <= len(t) <= (1 + a.tol) * n]
        pools.append(pool)
    if any(len(p) < 3 for p in pools):
        sys.exit("Pocas familias candidatas para emparejar tamaños: sube --tol")

    rng = random.Random(a.seed)
    null = []
    for _ in range(a.perms):
        chosen, used = [], set()
        for pool in pools:
            f = rng.choice(pool)
            while f in used:
                f = rng.choice(pool)
            used.add(f); chosen.append(set(fam_sets[f]))
        null.append(len(kc.consensus_genes(chosen, a.mode)))
    null.sort()
    ge = sum(1 for c in null if c >= c_obs)
    p = (ge + 1) / (a.perms + 1)
    q = lambda x: null[min(len(null) - 1, int(x * len(null)))]
    res = {"mirnas": a.mirnas, "mode": a.mode, "k_nombres": len(a.mirnas), "k_efectivo": k_eff,
           "tamanos_dianas": [len(s) for s in distinct], "genes_comunes_observados": c_obs,
           "genes_comunes_sin_deduplicar": c_obs_raw, "perms": a.perms, "seed": a.seed,
           "nulo": {"media": round(statistics.mean(null), 3), "mediana": statistics.median(null),
                    "p95": q(0.95), "p99": q(0.99), "max": null[-1]},
           "p_empirico": p, "p_minimo_resoluble": round(1 / (a.perms + 1), 6)}
    print(f"\nPanel: {', '.join(a.mirnas)}  |  modo {a.mode}")
    print(f"Nombres: {len(a.mirnas)}  ->  familias de semilla efectivas: {k_eff}"
          + ("  (¡hay miRNAs de la misma familia!)" if k_eff < len(a.mirnas) else ""))
    print(f"Genes en común observados: {c_obs}"
          + (f"  (sin deduplicar: {c_obs_raw})" if c_obs_raw != c_obs else ""))
    print(f"Nulo ({a.perms} paneles aleatorios, tamaños de dianas ±{a.tol:.0%}): "
          f"media {res['nulo']['media']}, mediana {res['nulo']['mediana']}, "
          f"p95 {res['nulo']['p95']}, p99 {res['nulo']['p99']}, máx {res['nulo']['max']}")
    print(f"p empírico = {p:.4g}  (resolución mínima {res['p_minimo_resoluble']:.1e})")
    out = a.out or f"convergence_{a.seed}.json"
    json.dump(res, open(out, "w"), ensure_ascii=False, indent=2)
    print(f"Guardado: {out}")


if __name__ == "__main__":
    main()
