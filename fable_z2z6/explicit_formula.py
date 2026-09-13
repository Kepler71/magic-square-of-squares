# -*- coding: utf-8 -*-
"""Явная формула (теорема Fable о 3-изогении) для Sel^phihat и Sel^phi через простые делители шести форм.
Сравнение с сертифицированными расчётами three_isog.sage по сетке."""
import json, glob, sys
from math import gcd
def v(n, p):
    n = abs(n); e = 0
    if n == 0: return 10**9
    while n % p == 0: n //= p; e += 1
    return e
def primes_of(x):
    x = abs(x); out = []; d = 2
    while d*d <= x:
        if x % d == 0:
            out.append(d)
            while x % d == 0: x //= d
        d += 1 if d == 2 else 2
    if x > 1: out.append(x)
    return out
def types(a, b):
    """P1: образ delta_p — всё; P2: образ тривиален."""
    F1 = a*b*(a+b); F2 = (a-b)*(a+2*b)*(2*a+b)
    P1 = set(p for p in primes_of(F1) if p != 2); P2 = set(p for p in primes_of(F2) if p != 2)
    if v(F1, 2) >= 2: P1.add(2)
    else: P2.add(2)   # 2 || ab(a+b) (ровно одно из a,b,a+b чётно)
    assert not (P1 & P2), (a, b, P1, P2)
    return sorted(P1), sorted(P2)
def cubic_index(x, p):
    """p = 1 mod 3: индекс x^{(p-1)/3} в {1, w, w^2} (w — фиксированный первообразный кубический корень из 1)."""
    t = pow(x % p, (p-1)//3, p)
    if t == 1: return 0
    # найдём w: минимальный g с g^{(p-1)/3} != 1
    g = 2
    while pow(g, (p-1)//3, p) == 1: g += 1
    w = pow(g, (p-1)//3, p)
    if t == w: return 1
    assert t == w*w % p
    return 2
def index9(x):
    """x нечётный, 3 ∤ x: дискретный логарифм x mod 9 по основанию 2, взятый mod 3 (кубы = ±1 = <8>)."""
    x %= 9; k = 0; y = 1
    while y != x: y = y*2 % 9; k += 1
    return k % 3
def rank3(rows, ncols):
    """ранг матрицы над GF(3)"""
    rows = [r[:] for r in rows]; rank = 0; col = 0
    for col in range(ncols):
        piv = None
        for i in range(rank, len(rows)):
            if rows[i][col] % 3: piv = i; break
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]
        inv = 1 if rows[rank][col] % 3 == 1 else 2
        rows[rank] = [(x*inv) % 3 for x in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col] % 3:
                f = rows[i][col]
                rows[i] = [(x - f*y) % 3 for x, y in zip(rows[i], rows[rank])]
        rank += 1
    return rank
def explicit_dims(a, b):
    P1, P2 = types(a, b)
    # Sel^phihat: d = prod_{p in P1} p^{e_p}; условия: для q in P2, q=1 mod 3: sum e_p ind_q(p) = 0; если 3 in P2: sum e_p ind9(p) = 0
    conds = []
    for q in P2:
        if q % 3 == 1: conds.append([cubic_index(p, q) for p in P1])
        if q == 3: conds.append([index9(p) for p in P1])
    dim_hat = len(P1) - (rank3(conds, len(P1)) if conds else 0)
    # Sel^phi: характеры с кондуктором на {q in P2: q=1 mod 3} и 9 (если 3 in P2); условие chi(p)=1 для p in P1
    gens = [q for q in P2 if q % 3 == 1] + ([9] if 3 in P2 else [])
    conds2 = []
    for p in P1:
        conds2.append([ (cubic_index(p, q) if q != 9 else index9(p)) for q in gens])
    dim_phi = len(gens) - (rank3(conds2, len(gens)) if (conds2 and gens) else 0)
    return dim_hat, dim_phi, dim_hat + dim_phi - 1, P1, P2
if __name__ == "__main__":
    BMAX = sys.argv[1] if len(sys.argv) > 1 else '60'
    recs = []
    for f in sorted(glob.glob(f'grid/shard_{BMAX}_*.jsonl')):
        for line in open(f): recs.append(json.loads(line))
    bad = 0; n = 0; zeros = 0
    for r in recs:
        if 'isog3' not in r: continue
        n += 1
        dh, dp, bd, P1, P2 = explicit_dims(r['a'], r['b'])
        if (dh, dp) != (r['isog3']['dim_hat'], r['isog3']['dim_phi']):
            bad += 1
            if bad < 10: print("MISMATCH", r['a'], r['b'], (dh, dp), (r['isog3']['dim_hat'], r['isog3']['dim_phi']), P1, P2)
        zeros += (bd == 0)
    print(f"pairs {n}, mismatches {bad}, explicit bound 0: {zeros}")
