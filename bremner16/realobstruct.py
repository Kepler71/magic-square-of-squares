#!/usr/bin/env python3
# Для базовых пар, у которых не нашлось рациональной точки, проверяем СТРОГОЕ вещественное
# препятствие: клетка-квадрат обязана быть >= 0, поэтому если {p: F1(p) >= 0} и {p: F2(p) >= 0}
# не пересекаются, то в этом семействе таких магических квадратов НЕТ ВООБЩЕ (доказано).
#   python3 realobstruct.py M [fams]
import sys
from fractions import Fraction
from math import gcd, sqrt
from itertools import combinations
COEF = [(1, 1, 0), (1, -1, -1), (1, 0, 1), (1, -1, 1), (1, 0, 0), (1, 1, -1), (1, 0, -1), (1, 1, 1), (1, -1, 0)]
GFAM = {'G1': ((1, 7), (3, 5)), 'G2': ((0, 8), (2, 6)), 'G3': ((1, 7), (0, 8))}


def forms(gt, m, n):
    (i, ip), (k, kp) = GFAM[gt]
    m, n = Fraction(m), Fraction(n)
    sq = lambda A, B: (A*A, 2*A*B, B*B)
    Pf, Rf = sq(m, n), sq(m, -n)
    cf = tuple(Fraction(m*m + n*n, 2)*t for t in (1, 0, 1))
    A1, B1 = COEF[i][1], COEF[i][2]
    A2, B2 = COEF[k][1], COEF[k][2]
    det = Fraction(A1*B2 - A2*B1)
    r1 = tuple(Pf[t] - cf[t] for t in range(3))
    r2 = tuple(Rf[t] - cf[t] for t in range(3))
    a = tuple((B2*r1[t] - B1*r2[t])/det for t in range(3))
    b = tuple((-A2*r1[t] + A1*r2[t])/det for t in range(3))
    return [tuple(cc*cf[t] + ca*a[t] + cb*b[t] for t in range(3)) for (cc, ca, cb) in COEF]

M = int(sys.argv[1]) if len(sys.argv) > 1 else 30
FAMS = sys.argv[2].split(',') if len(sys.argv) > 2 else ['G1', 'G2', 'G3']


def nonneg_test_points(quads):
    """точки-свидетели: середины между всеми вещественными корнями плюс концы"""
    rts = []
    for (A, B, C) in quads:
        A, B, C = float(A), float(B), float(C)
        if A == 0:
            if B != 0:
                rts.append(-C/B)
            continue
        d = B*B - 4*A*C
        if d >= 0:
            rts += [(-B - sqrt(d))/(2*A), (-B + sqrt(d))/(2*A)]
    rts.sort()
    pts = [rts[0] - 1.0] if rts else [0.0]
    for i in range(len(rts) - 1):
        pts.append((rts[i] + rts[i+1])/2)
    if rts:
        pts.append(rts[-1] + 1.0)
    pts += rts
    return pts


def compatible(f1, f2):
    """существует ли вещественное p с f1(p) >= 0 и f2(p) >= 0 (плюс проверка p = бесконечность)"""
    for p in nonneg_test_points([f1, f2]):
        v1 = float(f1[0])*p*p + float(f1[1])*p + float(f1[2])
        v2 = float(f2[0])*p*p + float(f2[1])*p + float(f2[2])
        if v1 >= -1e-9 and v2 >= -1e-9:
            return True
    return float(f1[0]) >= 0 and float(f2[0]) >= 0      # p -> бесконечность


mn = [(m, n) for m in range(1, M + 1) for n in list(range(-M, 0)) + list(range(1, M + 1))
      if gcd(m, abs(n)) == 1 and m != abs(n)]
stat = {}
for gt in FAMS:
    (i, ip), (k, kp) = GFAM[gt]
    free = [x for x in range(9) if x not in (i, ip, k, kp)]
    for (c1, c2) in combinations(free, 2):
        bad = tot = 0
        for (m, n) in mn:
            fs = forms(gt, m, n)
            tot += 1
            if not compatible(fs[c1], fs[c2]):
                bad += 1
        stat[(gt, c1, c2)] = (bad, tot)

print("Вещественное препятствие для базовых пар G-семейств "
      f"((m,n) взаимно простые, m,|n| <= {M}, всего {len(mn)} штук):")
print("| семейство | базовая пара | нет вещественных решений (доказано) | всего (m,n) |")
print("|---|---|---|---|")
for (gt, c1, c2), (bad, tot) in sorted(stat.items()):
    print(f"| {gt} | c{c1}/c{c2} | {bad} | {tot} |")
