# Поиск рациональных точек системы «восемь произведений арифметической сетки — квадраты»
# при общем (b,c) (не фиксируя наклон). Центр 1, клетки 1+ib+jc.
# 4 центральных условия: b, c, b+c, b-c ∈ B = {x: 1-x^2 ∈ Q^2} = {±A/C, ±B/C}.
# WLOG 0<c<b, b+c<1 (симметрии знаков/перестановки; положительность необходима).
# Затем точная проверка всех 8 произведений и выводов теоремы.
# Граница H — гипотенуза примитивной тройки; цикл по b конечный (len(vals)).
import sys, time, json
from math import gcd, isqrt
from fractions import Fraction as F
import numpy as np

H = int(sys.argv[1])
trip = []
for m in range(2, isqrt(H) + 2):
    for n in range(1, m):
        if (m - n) % 2 == 1 and gcd(m, n) == 1:
            C = m*m + n*n
            if C <= H: trip.append((m*m - n*n, 2*m*n, C))
vals = sorted({(A, C) for (a, b2, C) in trip for A in (a, b2)}, key=lambda t: t[0]/t[1])
print('triples', len(trip), 'values', len(vals), flush=True)
An = np.array([v[0] for v in vals], dtype=np.int64)
Cn = np.array([v[1] for v in vals], dtype=np.int64)

def is_sq_arr(M):
    ok = M >= 0
    Mf = np.where(ok, M, 0)
    r = np.round(np.sqrt(Mf.astype(np.float64))).astype(np.int64)
    good = np.zeros_like(ok)
    for d in (-1, 0, 1):
        rr = r + d
        good |= (rr >= 0) & (rr * rr == Mf)
    return ok & good

def issq_frac(x):
    if x <= 0: return False
    a, b = x.numerator, x.denominator
    return isqrt(a)**2 == a and isqrt(b)**2 == b

def sqfree_part(n):
    n = abs(n); d = 1; p = 2
    # n здесь — числитель*знаменатель квадратного класса; факторизация пробным делением с границей
    while p * p <= n and p < 10**7:
        e = 0
        while n % p == 0: n //= p; e += 1
        if e & 1: d *= p
        p += 1
    return d * n  # остаток предполагаем простым (проверяется ниже как флаг)

t0 = time.time(); cand = []; central = 0
for k in range(len(vals)):
    A1, C1 = vals[k]
    # c < b: индексы j<k (vals отсортированы по значению); b+c<1
    A2 = An[:k]; C2 = Cn[:k]
    D = C1 * C2
    Np = A1 * C2 + A2 * C1
    Nm = A1 * C2 - A2 * C1
    mask = (Np < D) & (Nm > 0)
    if not mask.any(): continue
    ok = mask & is_sq_arr((D - Np) * (D + Np)) & is_sq_arr((D - Nm) * (D + Nm))
    for j in np.nonzero(ok)[0]:
        central += 1
        cand.append((A1, C1, int(A2[j]), int(C2[j])))
    if (k + 1) % max(1, len(vals) // 10) == 0:
        print(f'  {k+1}/{len(vals)} central={central} {time.time()-t0:.0f}s', flush=True)
print('points with 4 central conditions:', central, flush=True)

full = []
LINES = [[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] + \
        [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] + \
        [[(t, t) for t in (-1, 0, 1)], [(t, -t) for t in (-1, 0, 1)]]
for (A1, C1, A2, C2) in cand:
    b = F(A1, C1); c = F(A2, C2)
    f = {(i, j): 1 + i*b + j*c for i in (-1, 0, 1) for j in (-1, 0, 1)}
    nsq = sum(issq_frac(f[x]*f[y]*f[w]) for (x, y, w) in LINES)
    if nsq == 8:
        q = b / c; r, s = q.numerator, q.denominator; z = b / r
        T = (1 - r*z)*(1 + (r+s)*z)*(1 - s*z); L = (1 - r*z)*(1 + (r-s)*z)*(1 + s*z)
        dT = sqfree_part(T.numerator * T.denominator); dL = sqfree_part(L.numerator * L.denominator)
        full.append(dict(b=str(b), c=str(c), r=r, s=s, z=str(z), dT=dT, dL=dL,
                         dT_div=((r - s) % dT == 0) if dT > 1 else True,
                         dL_div=((r + s) % dL == 0) if dL > 1 else True,
                         mod24=(dT % 24, dL % 24)))
        print('FULL 8-product point:', full[-1], flush=True)
print('points with all 8 products square:', len(full))
json.dump(dict(H=H, central=central, central_list=[list(x) for x in cand[:200]], full=full),
          open(f'global_search_H{H}.json', 'w'), indent=1)
