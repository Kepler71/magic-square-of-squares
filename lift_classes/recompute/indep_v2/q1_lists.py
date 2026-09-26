# -*- coding: utf-8 -*-
"""q1 (независимый пересчёт, не импортирует ни кода Codex, ни r1_*.py соседнего прогона).

По ФОРМУЛИРОВКЕ теоремы Codex:
  D_-(r,s) = {d>0 бесквадратное : p|d => p|(r-s), p = 1 mod 4 ; d = 1 mod 24}   (допустимые d_T)
  D_+(r,s) = то же с (r+s)                                                          (допустимые d_L)
Метод A: Sage factor(m) -> все подмножества простых p=1 mod 4 -> фильтр mod 24.
Метод B: прямой перебор d = 1, 25, 49, ... <= |m| (шаг 24), d | m, бесквадратность и p=1 mod 4 по Sage.
Наклоны: r/s, 2 <= s <= SMAX, -s < r < s, r != 0, gcd(r,s)=1 (r/s и s/r эквивалентны перестановкой).
Все циклы конечны (явные границы SMAX, |m| <= 2*SMAX).
"""
import os, sys, json, time
os.environ.setdefault("DOT_SAGE", "/tmp/claude_indep_v2_sage")
from sage.all import factor, is_squarefree, gcd, prime_divisors
from itertools import combinations

SMAX_B = 300      # метод B (перебор) -- до этой границы
SMAX_A = 1000     # метод A -- расширенная граница
t0 = time.time()

def lists_A(m):
    m = abs(int(m))
    ps = [int(p) for p, e in factor(m)] if m > 1 else []
    ps = [p for p in ps if p % 4 == 1]
    out = []
    for k in range(len(ps) + 1):
        for sub in combinations(ps, k):
            d = 1
            for p in sub:
                d *= p
            if d % 24 == 1:
                out.append(d)
    return sorted(out)

def lists_B(m):
    m = abs(int(m))
    out = []
    for d in range(1, m + 1, 24):          # d = 1 mod 24, d <= |m|
        if m % d:
            continue
        if not is_squarefree(d):
            continue
        if all(int(p) % 4 == 1 for p in prime_divisors(d)):
            out.append(d)
    return out

res = {}
# --- контрольные наклоны (r,s как в формуле f_ij = 1+(i r + j s) z)
codex_table = {(126, 451): ([1], [1, 577]), (73, 362): ([1], [1, 145]), (265, 298): ([1], [1])}
ctl = {}
for (r, s), (cT, cL) in codex_table.items():
    DT_A, DL_A = lists_A(r - s), lists_A(r + s)
    DT_B, DL_B = lists_B(r - s), lists_B(r + s)
    ok = (DT_A == DT_B == cT) and (DL_A == DL_B == cL)
    # наклон -r/s: те же девять клеток, но T и L меняются ролями
    ctl["%d/%d" % (r, s)] = {"r-s": r - s, "fac(r-s)": str(factor(r - s)), "r+s": r + s,
                             "fac(r+s)": str(factor(r + s)), "D_T(A)": DT_A, "D_T(B)": DT_B,
                             "D_L(A)": DL_A, "D_L(B)": DL_B, "codex": [cT, cL], "match": ok,
                             "for -r/s": [lists_A(-r - s), lists_A(-r + s)]}
    print("%d/%d" % (r, s), ctl["%d/%d" % (r, s)], flush=True)
res["control_slopes"] = ctl

# --- все наклоны
stats = {}
mismatchAB = 0
gcd_viol = 0
cor73_viol = 0
maxpairs = (0, None)
both_trivial_examples = []
cntA = {}
for s in range(2, SMAX_A + 1):
    for r in range(-s + 1, s):
        if r == 0 or gcd(r, s) != 1:
            continue
        DT, DL = lists_A(r - s), lists_A(r + s)
        if s <= SMAX_B:
            if DT != lists_B(r - s) or DL != lists_B(r + s):
                mismatchAB += 1
        for a in DT:
            for b in DL:
                if gcd(a, b) != 1:
                    gcd_viol += 1
        if 0 < abs(r - s) < 73 and DT != [1]:
            cor73_viol += 1
        if 0 < abs(r + s) < 73 and DL != [1]:
            cor73_viol += 1
        npairs = len(DT) * len(DL)
        if npairs > maxpairs[0]:
            maxpairs = (npairs, "%d/%d" % (r, s), DT, DL)
        for bound in (300, 600, 1000):
            if s <= bound:
                st = stats.setdefault(bound, {"pos": [0, 0], "all": [0, 0]})
                triv = (DT == [1] and DL == [1])
                st["all"][0] += 1; st["all"][1] += triv
                if r > 0:
                    st["pos"][0] += 1; st["pos"][1] += triv
    if s % 100 == 0:
        print("  s=%d  t=%.1fs" % (s, time.time() - t0), flush=True)

for bound in sorted(stats):
    st = stats[bound]
    for k in ("pos", "all"):
        n, t = st[k]
        st[k] = {"slopes": n, "both_trivial": t, "share": round(t / n, 6)}
    print("s<=%d:" % bound, st, flush=True)
res["stats"] = {str(k): v for k, v in stats.items()}
res["mismatch_A_vs_B(s<=%d)" % SMAX_B] = mismatchAB
res["gcd(d_T,d_L)>1 occurrences"] = gcd_viol
res["corollary73 violations"] = cor73_viol
res["max pairs"] = maxpairs
res["time_s"] = round(time.time() - t0, 1)
print("A/B mismatch:", mismatchAB, " gcd viol:", gcd_viol, " cor73 viol:", cor73_viol, " maxpairs:", maxpairs)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "q1_lists.json"), "w"),
          indent=1, ensure_ascii=False, default=str)
print("done %.1fs" % (time.time() - t0))
