#!/usr/bin/env python3
"""
KENRYU — auditoría de brazos (-5p/-3p): ¿qué nombres de miRBase se resuelven
en silencio a dianas del brazo OPUESTO?

Para cada miRNA maduro humano de miRBase con brazo, se consultan las dianas
locales (get_targets_local) y se comprueba con la semilla si la familia de
TargetScan que respondió modela de verdad ese brazo (check_arm_resolution):
la semilla de la familia debe aparecer al inicio de la secuencia madura
(offset 0-3, lo que tolera la deriva de anotación 5' entre la miRBase de
TargetScan 8.0 y miRBase 23.1). Un resultado no vacío con estado distinto
de "verified" es un intercambio silencioso de brazo.

v1.56 (criterio por nombre): 33 nombres, 19 horquillas; 14 horquillas
con ambos brazos devolviendo el mismo conjunto. La lista completa de
entonces está en docs/validation/arm_swaps.json.
v1.58 (criterio por semilla): 0 intercambios; 303 nombres conservan
dianas; 29 brazos bloqueados (docs/validation/arm_blocked.json): 15
brazos no modelados más 14 colapsos por letra (el fallback m_no_letter
juntaba paralogos distintos en la familia equivocada).

Uso (raíz del repo):  python scripts/audit_arm_resolution.py
"""
import json
import logging
import re
import sys

logging.disable(logging.CRITICAL)
sys.path.insert(0, ".")
import kenryu_engine as E

db = E.load_local_data()
print(f"arm verification: {E.arm_verification_status()}")

names = []
for line in open("local_db/mirbase_mature.fa"):
    if line.startswith(">hsa-"):
        names.append(line[1:].split()[0])

rows, swaps, blocked = 0, [], []
for n in sorted(set(names)):
    m = re.match(r"^hsa-(?:miR|let)-(\d+[a-z]?(?:-\d+)?)-(5p|3p)$", n, re.I)
    if not m:
        continue
    g = E.get_targets_local(n, db)
    if not g:
        check = E.check_arm_resolution(n, db)
        if check["status"] == "unverified":
            blocked.append({"nombre": n, "brazo_pedido": m.group(2),
                            "familia_TargetScan": check["family_key"],
                            "motivo": check["reason"]})
        continue
    rows += 1
    check = E.check_arm_resolution(n, db)
    if check["status"] != "verified":
        swaps.append({"nombre": n, "brazo_pedido": m.group(2),
                      "familia_TargetScan": check["family_key"],
                      "estado": check["status"], "motivo": check["reason"],
                      "n_genes": len(g)})

print(f"miRNAs maduros humanos con brazo y dianas locales: {rows}")
print(f"Intercambios silenciosos de brazo (criterio por semilla): {len(swaps)}")
for s in swaps[:15]:
    print("  ", s)
print(f"Brazos bloqueados (sin dianas, con aviso): {len(blocked)}")
json.dump(blocked, open("docs/validation/arm_blocked.json", "w"),
          indent=2, ensure_ascii=False)
