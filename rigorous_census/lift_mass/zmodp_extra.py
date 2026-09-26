#!/usr/bin/env python3
"""zmodp_extra.py -- C8: списки Codex перебором делителей (другой код); C9: при каких условиях образ
([T]_p,[L]_p) нетривиален (p = 11, 19, 23, 31, 43; все 27397 наклонов) -- наблюдение, объясняющее,
почему хватает p = 7, 11; C10: явный разбор трёх наклонов Codex и 504 наклонов, где мало p = 7.
Все циклы конечны."""
import json, time
from math import gcd
from zmodp_common import codex_list, local_image_Zp, tl_of, unit_code, primes_3mod4, is_prime
T0 = time.time()
def log(*a): print("[%6.1f c]" % (time.time() - T0), *a, flush=True)
SLOPES = [(r, s) for s in range(2, 301) for r in range(1, s) if gcd(r, s) == 1]

# C8
def sqfree(d):
    k = 2
    while k * k <= d:
        if d % (k * k) == 0: return False
        k += 1
    return True
def brute_list(n):
    n = abs(n)
    return [d for d in range(1, n + 1) if n % d == 0 and d % 24 == 1 and sqfree(d)
            and all(q % 4 == 1 for q in range(2, d + 1) if d % q == 0 and is_prime(q))]
bad = sum(1 for n in range(1, 600) if brute_list(n) != codex_list(n))
log("C8: списки D(n) перебором делителей vs codex_list, n = 1..599: расхождений %d" % bad)
nt = [(r, s) for (r, s) in SLOPES if len(brute_list(s - r)) > 1 or len(brute_list(r + s)) > 1]
log("C8: нетривиальных наклонов по перебору делителей: %d" % len(nt))
cand = sorted(set(d for n in range(1, 600) for d in brute_list(n) if d > 1))
log("C8: все нетривиальные d <= 599:", cand)
qr = lambda d, p: pow(d, (p - 1) // 2, p) == 1
log("C8: из них вычеты mod 7:", [d for d in cand if qr(d, 7)],
    "; из этих вычеты mod 11:", [d for d in cand if qr(d, 7) and qr(d, 11)])

# C9
for p in (11, 19, 23, 31, 43):
    tab = {}
    for (r, s) in SLOPES:
        TL = set(tl_of(v) for v in local_image_Zp(r, s, p))
        tn = any(t for (t, l) in TL); ln = any(l for (t, l) in TL)
        key = ("p|r-s" if (s - r) % p == 0 else "p∤r-s", "p|r+s" if (r + s) % p == 0 else "p∤r+s")
        c = tab.setdefault(key, [0, 0, 0])
        c[0] += 1; c[1] += tn; c[2] += ln
    log("C9 p=%d: {условие: [наклонов, с нетрив.[T], с нетрив.[L]]}:" % p, tab)

# C10
P = primes_3mod4(7, 400)
for (r, s) in [(126, 451), (73, 362), (265, 298)]:
    Dm, Dp = codex_list(s - r), codex_list(r + s)
    out = []
    for a in Dm:
        for b in Dp:
            k = [p for p in P if (unit_code(a, p), unit_code(b, p)) not in set(tl_of(v) for v in local_image_Zp(r, s, p))]
            out.append(((a, b), k[:6]))
    log("C10 %d/%d: D_-=%s D_+=%s; пара -> первые убивающие p:" % (r, s, Dm, Dp), out)
rec = json.load(open("zmodp_mass.json"))
p7 = [x for x in rec["records"] if x["p7_survivors"] != [[1, 1]]]
pairs = {}
for x in p7:
    for pr in x["pairs"]:
        if pr["killers"] and pr["killers"][0] != 7:
            pairs[(pr["dT"], pr["dL"])] = pairs.get((pr["dT"], pr["dL"]), 0) + 1
log("C10: 504 наклона (мало p=7): пары, пережившие p=7 -> число наклонов:", dict(sorted(pairs.items())))
log("C10: первые 10 таких наклонов:", [(x["r"], x["s"], x["p7_survivors"]) for x in p7[:10]])

# C11: от каких простых зависит вывод (совместный фильтр, число наклонов с выжившими кроме (1,1))
NT = [(r, s) for (r, s) in SLOPES if len(codex_list(s - r)) > 1 or len(codex_list(r + s)) > 1]
cacheTL = {}
def TLs(r, s, p):
    k = (r, s, p)
    if k not in cacheTL:
        cacheTL[k] = set(tl_of(v) for v in local_image_Zp(r, s, p))
    return cacheTL[k]
def nsurv(primes):
    bad = 0
    for (r, s) in NT:
        surv = [(a, b) for a in codex_list(s - r) for b in codex_list(r + s)
                if all((unit_code(a, p), unit_code(b, p)) in TLs(r, s, p) for p in primes)]
        bad += surv != [(1, 1)]
    return bad
for name, pr in [("{7}", [7]), ("{11}", [11]), ("{7,11}", [7, 11]), ("P без 7", P[1:]), ("P без 11", [7] + P[2:]),
                 ("P без 7 и 11", P[2:]), ("{19,...,43}", [19, 23, 31, 43]), ("p=1 mod 4: 5,13,17,29,37 (только Z_p-часть)", [5, 13, 17, 29, 37])]:
    log("C11 фильтр по %s: наклонов с выжившими кроме (1,1): %d из %d" % (name, nsurv(pr), len(NT)))
