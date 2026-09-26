# -*- coding: utf-8 -*-
"""Проверка уточнения по простым p = 3 mod 4 (наблюдение из recompute/v3_a51edd4e, t3c/t3e).
Лемма (вывод вручную в REVIEW_LIFT_CLASSES_2026-09-26.md, §5): при p = 3 mod 4 у z нет полюса;
при z в Z_p класс каждой противоположной пары 1 +- mu z определяется вычетом z0 = z mod p через
единичную клетку; [T]_p = c_{r+s}(z0), [L]_p = c_{r-s}(z0); согласованность c_r = c_s, c_{r+s} c_{r-s} = c_r.
(1) I_p как объединение по всем (r0,s0) != (0,0) mod p -- для p = 3,7,11,19,23,31,43.
(2) Минимальный нетривиальный d: бесквадратный, простые = 1 mod 4, d = 1 mod 24, квадрат mod 7.
(3) Три наклона: списки Codex и их фильтр по I_7.
(4) Массово: наклоны 0<r<s<=SMAX, взаимно простые; сколько нетривиальных пар (d_T,d_L) остаётся после
    фильтра по p = 3 mod 4, 7 <= p <= PMAX (верхняя оценка образа).
Все циклы ограничены явно.
"""
import sys, time
from math import gcd
from itertools import product
t0 = time.time()
SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 400

def chi(x, p):
    x %= p
    if x == 0: return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1

def pair_class(w, p):
    a, b = chi(1 + w, p), chi(1 - w, p)
    if a and b and a != b: return None
    return a if a else b

def image(r0, s0, p):
    img = set()
    for z0 in range(p):
        c = [pair_class(mu * z0, p) for mu in (r0, s0, r0 + s0, r0 - s0)]
        if None in c: continue
        if c[0] != c[1] or c[2] * c[3] != c[0]: continue
        img.add((c[2], c[3]))
    return img

print("(1) объединение I_p по всем (r0,s0) mod p:")
for p in (3, 7, 11, 19, 23, 31, 43):
    U = set()
    for r0 in range(p):
        for s0 in range(p):
            if r0 == 0 and s0 == 0: continue
            U |= image(r0, s0, p)
    print("   p=%d: %s" % (p, sorted(U)))

def factor(n):
    n = abs(n); f = {}; d = 2
    while d * d <= n and d < 10**6:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def sqfree(n): return all(e == 1 for e in factor(n).values())
def ok(d): return d > 0 and sqfree(d) and all(q % 4 == 1 for q in factor(d)) and d % 24 == 1

cand = [d for d in range(2, 5000) if ok(d)]
cand7 = [d for d in cand if chi(d, 7) == 1]
print("(2) Codex d>1:", cand[:8], " + квадрат mod 7:", cand7[:8], " минимум:", cand7[0])

def codex_list(N):
    ps = [q for q in factor(N) if q % 4 == 1]
    out = []
    for mask in range(1 << len(ps)):
        d = 1
        for i, q in enumerate(ps):
            if mask >> i & 1: d *= q
        if d % 24 == 1: out.append(d)
    return sorted(out)

P3 = [p for p in range(7, PMAX + 1) if p % 4 == 3 and all(p % q for q in range(2, int(p ** .5) + 1))]
def refine(r, s, primes):
    DT, DL = codex_list(r - s), codex_list(r + s)
    pairs = [(a, b) for a in DT for b in DL]
    left = []
    for (a, b) in pairs:
        good = True
        for p in primes:
            if (chi(a, p), chi(b, p)) not in image(r % p, s % p, p): good = False; break
        if good: left.append((a, b))
    return DT, DL, left

print("(3) три наклона:")
for (r, s) in ((126, 451), (73, 362), (265, 298), (-126, 451), (-73, 362)):
    DT, DL, left7 = refine(r, s, [7])
    print("   %d/%d: D_T=%s D_L=%s -> после I_7: %s" % (r, s, DT, DL, left7))

cnt = dict(slopes=0, codex_nontriv=0, left_nontriv_7=0, left_nontriv_all=0)
examples = []
for s in range(2, SMAX + 1):
    for r in range(1, s):
        if gcd(r, s) != 1: continue
        cnt["slopes"] += 1
        DT, DL = codex_list(r - s), codex_list(r + s)
        if len(DT) * len(DL) == 1: continue
        cnt["codex_nontriv"] += 1
        _, _, l7 = refine(r, s, [7])
        if len(l7) > 1: cnt["left_nontriv_7"] += 1
        _, _, la = refine(r, s, P3)
        if len(la) > 1:
            cnt["left_nontriv_all"] += 1
            if len(examples) < 5: examples.append((r, s, la))
    if s % 100 == 0: print("   s=%d %s t=%.1fs" % (s, cnt, time.time() - t0)); sys.stdout.flush()
print("(4) SMAX=%d PMAX=%d: %s примеры оставшихся: %s" % (SMAX, PMAX, cnt, examples))
print("время %.1f с" % (time.time() - t0))
