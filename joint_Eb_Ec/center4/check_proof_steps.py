# -*- coding: utf-8 -*-
"""Проверка промежуточных шагов доказательства теоремы D на всех найденных конфигурациях (A=10^8):
для каждого класса P: Z != 0, Z ≡ 0 mod prod pi_p^{2e_p}, Pi^2 <= N(Z) <= (N_L a / Pi)^2."""
import json, sys
from circle_lib import *
sys.setrecursionlimit(10000)
A = sys.argv[1] if len(sys.argv) > 1 else '100000000'
# перегенерировать конфигурации (search_twisted.json хранит только статистику) — повторяем поиск в памяти
from collections import defaultdict
from math import isqrt
Ai = int(A)
groups = defaultdict(list)
for be in range(2, isqrt(2 * Ai) + 1):
    b2 = be * be
    for al in range(1 if be % 2 else 2, be, 2):
        s = al * al + b2
        if s > 2 * Ai: break
        groups[s // 2].append((b2 - al * al) // 2)
found = set()
for a, us in groups.items():
    if len(us) < 3: continue
    S = set(us); L = sorted(us)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            x, y = L[i], L[j]
            if x + y in S: found.add((a, (x, y, x + y), (1, 1, -1)))
            if (x + y) % 2 == 0 and (x + y) // 2 in S: found.add((a, (x, y, (x + y) // 2), (1, 1, -2)))
            if (y - x) % 2 == 0 and (y - x) // 2 in S and (y - x) // 2 not in (x, y): found.add((a, ((y - x) // 2, x, y), (2, 1, -1)))
            if x + 2 * y in S: found.add((a, (y, x, x + 2 * y), (2, 1, -1)))
            if y + 2 * x in S: found.add((a, (x, y, y + 2 * x), (2, 1, -1)))
del groups
stats = dict(A=Ai, configs=len(found), classes=0, Z_zero=0, not_divisible=0, lower_fail=0, upper_fail=0)
for a, tr, ns in sorted(found):
    tr = list(tr); ns = list(ns)
    NL = sum(abs(n) for n in ns)
    ok, cls = balance_check(a, tr, NL)
    for sig, ps, Pi, r in cls:
        stats['classes'] += 1
        Z, NZ, Pi2, div = relation_Z(a, tr, ns, ps)
        if Z == (0, 0): stats['Z_zero'] += 1
        if not div: stats['not_divisible'] += 1
        if NZ < Pi2: stats['lower_fail'] += 1
        if NZ * Pi2 > (NL * a) ** 2: stats['upper_fail'] += 1
print(stats)
json.dump(stats, open(f'check_proof_steps_A{Ai}.json', 'w'), indent=1)
