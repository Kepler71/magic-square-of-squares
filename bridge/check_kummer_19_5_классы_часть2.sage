# -*- coding: utf-8 -*-
# Часть 2: своя точка, независимость точек PARI, ранг без mwrank, перечисление классов.
from sage.all import *
import sys

def hdr(t):
    print("\n" + "="*78); print(t); print("="*78); sys.stdout.flush()

m, n = 19, 5
s = QQ(m**2+n**2)/2
b = s*m**2*n**2
e1, e2, e3 = -b, -s*m**4, -s*n**4
a2 = -(e1+e2+e3); a4 = e1*e2+e1*e3+e2*e3; a6 = -e1*e2*e3
E = EllipticCurve(QQ, [0, a2, 0, a4, a6])
print("E =", E)

hdr("[A] ТОЧКИ, КОТОРЫЕ ВЫДАЛ PARI ПРИ РАЗНЫХ EFFORT — НЕЗАВИСИМЫ ЛИ ОНИ?")
cand = {}
for eff in [0,1,2,3,4,6]:
    try:
        out = pari(E).ellrank(eff)
        print("  effort %s -> %s" % (eff, out))
        for pt in out[3]:
            x = QQ(pt[0]); y = QQ(pt[1])
            P = E(x, y)
            cand[(x,y)] = P
    except Exception as ex:
        print("  effort %s ОШИБКА: %s" % (eff, ex))

pts = list(cand.values())
print("\n  различных точек собрано:", len(pts))
for P in pts:
    print("    P = (%s, %s)   h = %s" % (P[0], P[1], P.height()))

if len(pts) >= 2:
    M = matrix(RR, len(pts), len(pts),
               [[pts[i].height_pairing_matrix if False else pts[i].height() if i==j else 0
                 for j in range(len(pts))] for i in range(len(pts))])
    try:
        HM = E.height_pairing_matrix(pts)
        print("\n  матрица высотного спаривания:")
        print(HM)
        print("  det =", HM.determinant())
        print("  ЕСЛИ det заметно != 0 -> точки НЕЗАВИСИМЫ -> ранг >= %d -> вывод Codex РУШИТСЯ" % len(pts))
    except Exception as ex:
        print("  height_pairing_matrix ОШИБКА:", ex)

hdr("[B] ЯВНЫЕ СООТНОШЕНИЯ МЕЖДУ ЭТИМИ ТОЧКАМИ")
# ищем малые целые соотношения a*P1 + b*P2 + ... = torsion
Tpts = [P for P in E.torsion_points()]
print("  кручение:", Tpts)
if len(pts) >= 2:
    P1 = pts[0]
    for j in range(1, len(pts)):
        Pj = pts[j]
        found = None
        for k in range(-8, 9):
            if k == 0: continue
            for T in Tpts:
                if k*P1 + T == Pj:
                    found = (k, T)
            for T in Tpts:
                if k*P1 + T == -Pj:
                    found = (-k, T)
        print("   P%d = k*P1 + T ?  ->" % (j+1), found)

hdr("[C] НЕЗАВИСИМАЯ ВЕРХНЯЯ ГРАНИЦА РАНГА БЕЗ mwrank/PARI-2-descent")
print("  знак функционального уравнения (root number):", E.root_number())
try:
    ub = E.analytic_rank_upper_bound(max_Delta=2.0, adaptive=True, root_number=-1)
    print("  строгая верхняя граница аналитического ранга:", ub)
except Exception as ex:
    print("  analytic_rank_upper_bound ОШИБКА:", ex)
    try:
        ub = E.analytic_rank_upper_bound()
        print("  (без параметров):", ub)
    except Exception as ex2:
        print("  и без параметров ОШИБКА:", ex2)
print("  Kolyvagin+GZ+модулярность: если аналитический ранг <= 1, то алг. ранг = аналит. ранг.")
print("  аналитический ранг (численно, pari):", E.analytic_rank())

hdr("[D] mwrank С БОЛЬШИМ ВТОРЫМ ПРЕДЕЛОМ")
for lim in [10, 13, 15]:
    try:
        E2 = EllipticCurve(QQ, [0, a2, 0, a4, a6])
        E2.two_descent(second_limit=lim, verbose=False)
        print("  second_limit=%s -> rank_bounds %s" % (lim, (E2.rank_bounds(),)))
    except Exception as ex:
        print("  second_limit=%s ОШИБКА: %s" % (lim, ex))
try:
    print("  E.rank(only_use_mwrank=False) =", E.rank(only_use_mwrank=False))
except Exception as ex:
    print("  ОШИБКА:", ex)

hdr("[E] SELMER И SHA")
print("  selmer_rank (Sage) =", E.selmer_rank())
try:
    print("  rank_bound() =", E.rank_bound())
except Exception as ex:
    print("  rank_bound ОШИБКА:", ex)
try:
    print("  Sha.an() =", E.sha().an())
    print("  Sha.an_numerical() =", E.sha().an_numerical())
except Exception as ex:
    print("  Sha ОШИБКА:", ex)
