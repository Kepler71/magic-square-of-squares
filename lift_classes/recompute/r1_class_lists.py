# -*- coding: utf-8 -*-
"""r1: независимый пересчёт списков D_-(r,s), D_+(r,s) по ФОРМУЛИРОВКЕ теоремы Codex
(не импортирует код Codex).

D_-(r,s) = {d>0 бесквадратное: p|d => p|(r-s), p = 1 mod 4, d = 1 mod 24}   -- допустимые d_T
D_+(r,s) = то же с (r+s)                                                      -- допустимые d_L

Две реализации:
  A) бесквадратные делители радикала |m| (факторизация пробным делением), затем фильтр;
  B) перебор d = 1..|m| с d | m и проверкой "допустимости" d (своя таблица по решету).
Наклоны: канонические пары 1 <= r < s, gcd(r,s)=1 (s <= SMAX). Для наклона -r/s списки меняются
местами (r-s <-> -(r+s)); свойство "оба тривиальны" от знака не зависит. Обмен r<->s списков не меняет.
"""
import json, sys, time
from math import gcd

SMAX_LIST = [300, 600, 1000, 2000]
SMAX = max(SMAX_LIST)
LIM = 2 * SMAX + 5

# ---- решето: наименьший простой делитель
spf = list(range(LIM + 1))
for i in range(2, int(LIM ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, LIM + 1, i):
            if spf[j] == j:
                spf[j] = i

def factor_sieve(n):
    n = abs(n); out = {}
    while n > 1:
        p = spf[n]; out[p] = out.get(p, 0) + 1; n //= p
    return out

def factor_trial(n):
    n = abs(n); out = {}; p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

def admissible_d(d):
    """d>0 бесквадратное, все простые = 1 mod 4, d = 1 mod 24."""
    if d <= 0 or d % 24 != 1:
        return False
    f = factor_trial(d)
    return all(e == 1 and p % 4 == 1 for p, e in f.items())

ADM = [d for d in range(1, LIM + 1) if admissible_d(d)]

def D_A(m):
    """Реализация A: бесквадратные делители радикала, фильтр теоремы."""
    primes = sorted(factor_trial(m))
    ds = [1]
    for p in primes:
        ds = ds + [d * p for d in ds]
    res = []
    for d in ds:
        f = factor_trial(d)
        if all(p % 4 == 1 for p in f) and d % 24 == 1:
            res.append(d)
    return sorted(res)

def D_B(m):
    """Реализация B: перебор допустимых d <= |m|, делящих m."""
    m = abs(m)
    return [d for d in ADM if d <= m and m % d == 0]

t0 = time.time()
stats = {}
nontriv_examples = []
max_len = (0, None)
mismatch = 0
count_by_bound = {B: [0, 0, 0, 0] for B in SMAX_LIST}  # всего, оба тривиальны, только T нетрив, только L нетрив
both_nontriv = {B: 0 for B in SMAX_LIST}
checked_cor73 = 0; viol_cor73 = 0; viol_gcd = 0
for s in range(2, SMAX + 1):
    if s % 200 == 0:
        print(f"  s={s}  t={time.time()-t0:.1f}s", flush=True)
    for r in range(1, s):
        if gcd(r, s) != 1:
            continue
        m_minus, m_plus = r - s, r + s
        DT = D_B(m_minus); DL = D_B(m_plus)
        if s <= 600:  # вторая реализация на части диапазона (дорогая)
            if DT != D_A(m_minus) or DL != D_A(m_plus):
                mismatch += 1
        # следствие 73
        if 0 < abs(m_minus) < 73:
            checked_cor73 += 1
            if DT != [1]: viol_cor73 += 1
        if 0 < abs(m_plus) < 73:
            checked_cor73 += 1
            if DL != [1]: viol_cor73 += 1
        # попарная взаимная простота всех допустимых пар
        for a in DT:
            for b in DL:
                if gcd(a, b) != 1: viol_gcd += 1
        for B in SMAX_LIST:
            if s <= B:
                c = count_by_bound[B]; c[0] += 1
                tT, tL = (DT == [1]), (DL == [1])
                if tT and tL: c[1] += 1
                elif not tT and tL: c[2] += 1
                elif tT and not tL: c[3] += 1
                else: both_nontriv[B] += 1
        if s <= 300 and (DT != [1] or DL != [1]) and len(nontriv_examples) < 40:
            nontriv_examples.append((r, s, DT, DL))
        L = len(DT) * len(DL)
        if L > max_len[0]:
            max_len = (L, (r, s, DT, DL))

print("время", round(time.time() - t0, 1), "с")
print("расхождений реализаций A/B (s<=600):", mismatch)
print("следствие 73: проверено", checked_cor73, "нарушений", viol_cor73)
print("нарушений gcd(d_T,d_L)=1:", viol_gcd)
for B in SMAX_LIST:
    c = count_by_bound[B]
    print(f"s<={B}: наклонов {c[0]}, оба тривиальны {c[1]} ({c[1]/c[0]:.4%}), только D_- нетрив {c[2]}, "
          f"только D_+ нетрив {c[3]}, оба нетрив {both_nontriv[B]}")
print("макс. число пар классов:", max_len)

# ---- три наклона Codex, обе знаковые конвенции
codex = {"126/451": ([1], [1, 577]), "73/362": ([1], [1, 145]), "265/298": ([1], [1])}
spec = {}
ok_all = True
for key, (cT, cL) in codex.items():
    r, s = map(int, key.split("/"))
    DT, DL = D_A(r - s), D_A(r + s)
    DTb, DLb = D_B(r - s), D_B(r + s)
    same = (DT == cT and DL == cL and DT == DTb and DL == DLb)
    ok_all &= same
    # противоположный знак наклона -r/s: списки меняются местами
    DTn, DLn = D_A(-r - s), D_A(-r + s)
    spec[key] = {"r-s": r - s, "fac(r-s)": factor_trial(r - s), "r+s": r + s, "fac(r+s)": factor_trial(r + s),
                 "D_T": DT, "D_L": DL, "codex_D_T": cT, "codex_D_L": cL, "match": same,
                 "slope_-r/s": {"D_T": DTn, "D_L": DLn}}
    print(key, spec[key])
print("совпадение с таблицей Codex:", ok_all)

out = {"bounds": {B: {"slopes": count_by_bound[B][0], "both_trivial": count_by_bound[B][1],
                      "only_T_nontrivial": count_by_bound[B][2], "only_L_nontrivial": count_by_bound[B][3],
                      "both_nontrivial": both_nontriv[B],
                      "fraction_both_trivial": count_by_bound[B][1] / count_by_bound[B][0]} for B in SMAX_LIST},
       "impl_mismatch_s_le_600": mismatch, "cor73_checked": checked_cor73, "cor73_violations": viol_cor73,
       "gcd_violations": viol_gcd, "max_pairs": max_len, "admissible_d_upto_1200": ADM[:40],
       "first_nontrivial_examples_s_le_300": nontriv_examples, "codex_slopes": spec, "codex_table_match": ok_all}
json.dump(out, open("r1_class_lists.json", "w"), indent=1, default=str)
