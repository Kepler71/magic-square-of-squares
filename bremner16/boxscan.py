#!/usr/bin/env python3
# Независимый прямой поиск в G-семействах: никаких эллиптических кривых и решёток Морделла–Вейля.
# Для (m,n) и p = u/v все девять клеток — целочисленные квадратичные формы от (u,v) с точностью до
# общего множителя; проверяем перебором по прямоугольнику |u| <= U, 1 <= v <= U, gcd(u,v) = 1.
# Это ближайший аналог области поиска самого Бремнера (p+q+num(λ)+den(λ) <= 1000, с. 298).
#   python3 boxscan.py M U [fams]
import sys
from fractions import Fraction
from math import gcd, isqrt
import numpy as np

COEF = [(1, 1, 0), (1, -1, -1), (1, 0, 1), (1, -1, 1), (1, 0, 0), (1, 1, -1), (1, 0, -1), (1, 1, 1), (1, -1, 0)]
GFAM = {'G1': ((1, 7), (3, 5)), 'G2': ((0, 8), (2, 6)), 'G3': ((1, 7), (0, 8))}

M = int(sys.argv[1]) if len(sys.argv) > 1 else 60
U = int(sys.argv[2]) if len(sys.argv) > 2 else 200
FAMS = sys.argv[3].split(',') if len(sys.argv) > 3 else ['G1', 'G2', 'G3']


def forms(gt, m, n):
    """девять клеток как (alpha,beta,gamma): cell = (alpha u^2 + beta uv + gamma v^2)/v^2 (Fraction)"""
    (i, ip), (k, kp) = GFAM[gt]
    m, n = Fraction(m), Fraction(n)
    # P = (mu+nv)/v, Q = (mv-nu)/v, R = (mu-nv)/v, S = (mv+nu)/v
    sq = lambda A, B: (A*A, 2*A*B, B*B)                     # (Au+Bv)^2
    Pf, Qf, Rf = sq(m, n), sq(-n, m), sq(m, -n)
    cf = tuple(Fraction(m*m + n*n, 2)*t for t in (1, 0, 1))
    # a, b из cell_i = P^2, cell_k = R^2
    A1, B1 = COEF[i][1], COEF[i][2]
    A2, B2 = COEF[k][1], COEF[k][2]
    det = Fraction(A1*B2 - A2*B1)
    r1 = tuple(Pf[t] - cf[t] for t in range(3))
    r2 = tuple(Rf[t] - cf[t] for t in range(3))
    a = tuple((B2*r1[t] - B1*r2[t])/det for t in range(3))
    b = tuple((-A2*r1[t] + A1*r2[t])/det for t in range(3))
    out = []
    for (cc, ca, cb) in COEF:
        out.append(tuple(cc*cf[t] + ca*a[t] + cb*b[t] for t in range(3)))
    assert out[i] == Pf and out[k] == Rf
    return out


def scale(fs):
    D = 1
    for f in fs:
        for x in f:
            D = D*x.denominator//gcd(D, x.denominator)
    return D, [tuple(int(x*D) for x in f) for f in fs]


def issq(arr):
    with np.errstate(invalid='ignore'):
        r = np.sqrt(np.maximum(arr, 0).astype(np.float64))
    r = np.rint(r).astype(np.int64)
    return (arr > 0) & (r*r == arr)


mn = [(m, n) for m in range(1, M + 1) for n in list(range(-M, 0)) + list(range(1, M + 1))
      if gcd(m, abs(n)) == 1 and m != abs(n)]
mn.sort(key=lambda t: (max(t[0], abs(t[1])), t[0], t[1]))
print(f"# boxscan: (m,n) {len(mn)} шт., |u|<={U}, 1<=v<={U}, семейства {FAMS}", flush=True)

uu = np.arange(-U, U + 1, dtype=np.int64)
tested = 0
hits = 0
for (m, n) in mn:
    for gt in FAMS:
        fs = forms(gt, m, n)
        D, F = scale(fs)
        (i, ip), (k, kp) = GFAM[gt]
        free = [x for x in range(9) if x not in (i, ip, k, kp)]     # 4 клетки + центр
        for v in range(1, U + 1):
            u = uu[np.gcd(np.abs(uu), v) == 1]
            if len(u) == 0:
                continue
            u2 = u*u; uv = u*v; v2 = v*v
            vals = [D*(F[c][0]*u2 + F[c][1]*uv + F[c][2]*v2) for c in free]
            f0, f1, f2 = issq(vals[0]), issq(vals[1]), issq(vals[2])
            keep = f0 | f1 | f2
            tested += len(u)
            if not keep.any():
                continue
            idx = np.nonzero(keep)[0]
            cnt = (f0[idx].astype(np.int8) + f1[idx] + f2[idx]
                   + issq(vals[3][idx]) + issq(vals[4][idx]))
            for t in np.nonzero(cnt >= 3)[0]:
                j = idx[t]
                uu0 = int(u[j])
                cells = [Fraction(F[c][0]*uu0*uu0 + F[c][1]*uu0*v + F[c][2]*v*v, D*v*v) for c in range(9)]
                nsq = sum(1 for x in cells if x > 0 and isqrt(x.numerator*x.denominator)**2
                          == x.numerator*x.denominator)
                if nsq >= 7 and len(set(cells)) == 9:
                    hits += 1
                    print(f"HIT {gt} (m,n)=({m},{n}) p={uu0}/{v} nsq={nsq} cells={[str(x) for x in cells]}",
                          flush=True)
    if mn.index((m, n)) % 200 == 0:
        print(f"# ... (m,n)=({m},{n}) проверено пар (u,v): {tested}, попаданий {hits}", flush=True)
print(f"# ИТОГ boxscan: проверено {tested} точек (u,v) на семейство-кадр, попаданий {hits}", flush=True)
