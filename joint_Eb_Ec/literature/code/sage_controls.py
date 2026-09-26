# Контроли к литературной записке (Sage 10.9).
# (A) Нужна ли гипотеза P_j ∈ 2E в рациональном варианте? Ищем три точки подгруппы ранга 1
#     {mG+T : 1<=|m|<=M, T ∈ E[2]} с тремя РАЗЛИЧНЫМИ x в арифметической прогрессии на E_N, rank 1.
# (B) Корневое число E_N по N mod 8 (в RESULT_RANK_TWO это «по памяти»).
# (C) Нет ли у y^2=x^3-x рациональной циклической подгруппы порядка 4 (гипотеза теоремы Смита 1702.02325).
import os, itertools, time
from sage.all import EllipticCurve, QQ, ZZ, squarefree_part, is_squarefree

t0 = time.time()
# (C)
E1 = EllipticCurve([-1, 0])
C = E1.isogeny_class()
print("(C) класс изогении y^2=x^3-x:", [c.cremona_label() for c in C.curves])
print("    матрица степеней изогений:\n", C.matrix())
idx = [i for i, c in enumerate(C.curves) if c.is_isomorphic(E1)][0]
print("    индекс E1 в классе:", idx, "; строка степеней:", list(C.matrix().row(idx)))
print("    (циклическая 4-изогения из E1 <=> в строке E1 есть 4)")

# (B)
from collections import defaultdict
rn = defaultdict(set)
for N in range(1, 400):
    if not is_squarefree(N):
        continue
    E = EllipticCurve([-N*N, 0])
    rn[N % 8].add(E.root_number())
print("(B) корневые числа E_N, N<400 бесквадратные, по N mod 8:", dict(sorted((k, sorted(v)) for k, v in rn.items())))

# (A)
rank1 = []
for N in range(1, 120):
    if not is_squarefree(N):
        continue
    E = EllipticCurve([-N*N, 0])
    r = E.rank(proof=False)
    if r == 1:
        rank1.append(N)
print("(A) бесквадратные N<120 с rank E_N = 1 (proof=False):", rank1)
Mmax = 12
found = []
for N in rank1:
    E = EllipticCurve([-N*N, 0])
    G = E.gens(proof=False)[0]
    T = [E(0), E(0, 0), E(N, 0), E(-N, 0)]
    pts = {}
    for m in range(1, Mmax+1):
        mG = m*G
        for t in T:
            P = mG + t
            if P.is_zero():
                continue
            pts.setdefault(P[0], []).append((m, T.index(t)))
    xs = sorted(pts)
    xset = set(xs)
    for i in range(len(xs)):
        for k in range(i+1, len(xs)):
            mid = (xs[i] + xs[k]) / 2
            if mid in xset and mid != xs[i] and mid != xs[k]:
                found.append((N, (pts[xs[i]][0], pts[mid][0], pts[xs[k]][0]), (xs[i], mid, xs[k])))
    # тот же поиск только среди 2E: m чётное, T=0 (т.е. x(2kG))
print(f"(A) найдено x-прогрессий длины 3 в подгруппах ранга 1 (|m|<={Mmax}, с кручением): {len(found)}")
for f in found[:15]:
    N, idxs, x3 = f
    in2E = all(m % 2 == 0 and t == 0 for (m, t) in idxs)
    print("   N=%d  (m,T-индекс)=%s  все в 2E: %s  x=%s" % (N, idxs, in2E, [str(v) for v in x3]))
print("время %.1f c" % (time.time() - t0))
