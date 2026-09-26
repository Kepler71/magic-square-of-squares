# Независимые контроли к утверждениям 7, 8, 9, 10 (Sage 10.9 / PARI).
# ЧИТАЛ до написания: NOTE.md, чужие sage_controls.py/.log, bst_example_check.log (без точек), Br99 с.290–291, BST с.2.
import time, sys
from sage.all import EllipticCurve, QQ, ZZ, PolynomialRing, is_squarefree, pari, gcd
t0 = time.time()
R = PolynomialRing(QQ, 'x'); x = R.gen()

# (C') Утв. 7: у y^2=x^3-x нет рациональной циклической подгруппы порядка 4.
# Подгруппа <P> порядка 4 рациональна <=> x(P) рационально (2P -- рациональная точка 2-кручения).
# x(2P) = (x^2+1)^2 / (4x(x^2-1)) на y^2=x^3-x.
E1 = EllipticCurve([-1, 0])
for xT in (0, 1, -1):
    num = (x**2 + 1)**2 - 4*xT*x*(x**2 - 1)
    rr = num.roots(QQ)
    print(f"(C') половины точки ({xT},0): многочлен {num.factor()}, рациональные корни: {rr}")
# контроль формулы удвоения на случайной точке над числовым полем не нужен: сверим на E1 над Q(i) формально
print("(C') проверка формулы: x(2P) для P=(i, ...) : ", end="")
from sage.all import NumberField
K = NumberField(x**2 + 1, 'i'); i = K.gen()
EK = E1.change_ring(K)
P = EK.lift_x(i); print((2*P).xy() if not (2*P).is_zero() else "O", "(ожидается (0,0))")

# (B') Утв. 8: корневые числа через PARI ellrootno, бесквадратные N < 3000.
bad = []; cnt = 0
for N in range(1, 3000):
    if not is_squarefree(N): continue
    w = int(pari.ellrootno(pari.ellinit([0, 0, 0, -N*N, 0])))
    pred = 1 if N % 8 in (1, 2, 3) else -1
    cnt += 1
    if w != pred: bad.append((N, w))
print(f"(B') бесквадратных N<3000: {cnt}; расхождений с правилом (+1 при 1,2,3; −1 при 5,6,7 mod 8): {len(bad)} {bad[:10]}")

# (E) Утв. 10: точки Br99/BST на E_1254.
n = 1254
pts = [(-528, 26136), (-363, 22869), (-198, 17424)]
on = [yy*yy == xx*(xx*xx - n*n) for xx, yy in pts]
xs = [p[0] for p in pts]
print(f"(E) E_1254: на кривой {on}; x-прогрессия: {xs[1]-xs[0] == xs[2]-xs[1]}; x<0 у всех: {all(v < 0 for v in xs)}")
E = EllipticCurve([-n*n, 0])
Ps = [E(p) for p in pts]
H = E.height_pairing_matrix(Ps)
print(f"(E) det матрицы высот трёх точек = {H.det():.6f} (>0 => независимы)")

# (A') Утв. 9: x-прогрессии в E(Q) ранга 1, бесквадратные N < 200, |m| <= 16, P = mG + T.
M = 16
rank1 = []; found = []; notsat = []
for N in range(1, 200):
    if not is_squarefree(N): continue
    E = EllipticCurve([-N*N, 0])
    r = E.rank(only_use_mwrank=True, proof=False)
    if r != 1: continue
    gens = E.gens(proof=False)
    G = gens[0]
    # насыщенность: G не делится на 2,3,5,7 в E(Q)/E[2]
    sat, idx, _ = E.saturation([G], max_prime=50)
    if idx != 1: notsat.append((N, idx)); G = sat[0]
    rank1.append(N)
    T = [E(0)] + [E(e, 0) for e in (0, N, -N)]
    xs = {}
    for m in range(1, M + 1):
        mG = m*G
        for k, t in enumerate(T):
            Q = mG + t
            xs.setdefault(Q[0], (m, k))
    X = sorted(xs); Xs = set(X)
    for a_ in range(len(X)):
        for c_ in range(a_ + 1, len(X)):
            mid = (X[a_] + X[c_]) / 2
            if mid in Xs:
                found.append((N, xs[X[a_]], xs[mid], xs[X[c_]]))
print(f"(A') бесквадратные N<200 с rank 1: {len(rank1)} шт.; из них N<120: {len([v for v in rank1 if v < 120])}")
print(f"(A') ненасыщенных образующих: {notsat}")
print(f"(A') найдено x-прогрессий длины 3 среди x(mG+T), 1<=m<={M}: {len(found)} {found[:10]}")
print("время %.1f c" % (time.time() - t0))
