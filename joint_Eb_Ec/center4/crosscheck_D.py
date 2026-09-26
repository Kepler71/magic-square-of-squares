# -*- coding: utf-8 -*-
"""Независимая (без circle_lib) перепроверка теоремы D на конфигурациях до A:
чистота линии u в p: p ∤ gcd(alpha_u, beta_u); ориентация: sigma=+1, если alpha*beta + iota*u ≡ 0 (mod p),
иначе (alpha*beta − iota*u ≡ 0) sigma=−1, где iota^2 ≡ −1 (mod p). Классы — равенство векторов sigma с точностью до знака."""
import sys, json
from math import isqrt, gcd
from collections import defaultdict
from sympy import factorint, sqrt_mod
A = int(sys.argv[1]) if len(sys.argv) > 1 else 10**8
groups = defaultdict(list)
for be in range(2, isqrt(2 * A) + 1):
    b2 = be * be
    for al in range(1 if be % 2 else 2, be, 2):
        s = al * al + b2
        if s > 2 * A: break
        groups[s // 2].append((b2 - al * al) // 2)
found = set()
for a, us in groups.items():
    if len(us) < 3: continue
    S = set(us); L = sorted(us)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            x, y = L[i], L[j]
            if x + y in S: found.add((a, (x, y, x + y), 3))
            if (x + y) % 2 == 0 and (x + y) // 2 in S: found.add((a, (x, y, (x + y) // 2), 4))
            if (y - x) % 2 == 0 and (y - x) // 2 in S and (y - x) // 2 not in (x, y): found.add((a, ((y - x) // 2, x, y), 4))
            if x + 2 * y in S: found.add((a, (y, x, x + 2 * y), 4))
            if y + 2 * x in S: found.add((a, (x, y, y + 2 * x), 4))
del groups
res = dict(A=A, configs=len(found), fails=0, orientation_undefined=0, max_ratio=0.0)
for a, tr, NL in found:
    classes = defaultdict(int)
    for p, e in factorint(a).items():
        if p % 4 != 1: continue
        io = sqrt_mod(p - 1, p)
        sig = []
        pure = True
        for u in tr:
            al, be = isqrt(a - u), isqrt(a + u)
            if gcd(al, be) % p == 0: pure = False; break
            X = al * be
            if (X + io * u) % p == 0: sig.append(1)
            elif (X - io * u) % p == 0: sig.append(-1)
            else: res['orientation_undefined'] += 1; pure = False; break
        if not pure: continue
        if sig[0] < 0: sig = [-s for s in sig]
        key = tuple(sig)
        classes[key] = (classes[key] or 1) * p ** e
    for key, Pi in classes.items():
        r = Pi * Pi / (NL * a)
        res['max_ratio'] = max(res['max_ratio'], r)
        if Pi * Pi > NL * a: res['fails'] += 1
print(res)
json.dump(res, open(f'crosscheck_D_A{A}.json', 'w'), indent=1)
