# Широкий численный контроль: все бесквадратные N в [N0, N1] ранга >= 1 (образующие mwrank),
# примитивные направления G с коэффициентами |c_i| <= CMAX, M = 2..MMAX:
#  - слабая примитивность (нужная доказательству): l | B_{2M}(G), l не делит B_{2m}(G), m < M;
#  - полная примитивность (формулировка VY): l не делит B_k(G), k < 2M;
#  - нет AP среди x(2mG), m <= MMAX;
#  - отдельно считаем случаи, когда единственный слабо-примитивный простой = 2 (ветка l = 2 доказательства),
#    и проверяем там v_2(x(2MG)) <= -2 (чётность).
from sage.all import *
import sys, itertools, time
N0, N1, MMAX, CMAX = [int(a) for a in sys.argv[1:5]]
t0 = time.time()
tot_dirs = 0; fails = []; aps = []; only2 = []; par_fail = []; curves = []
def prim_part(Bn, prev):
    r = Bn
    for b in prev:
        g = gcd(r, b)
        while g > 1:
            r //= g
            g = gcd(r, b)
    return r
for Nv in range(N0, N1+1):
    if not Integer(Nv).is_squarefree(): continue
    E = EllipticCurve([0,0,0,-Nv**2,0])
    try:
        gens = E.gens()
    except Exception as e:
        print("N=%d: gens не найдены (%s)" % (Nv, e)); continue
    r = len(gens)
    if r == 0: continue
    curves.append((Nv, r))
    for c in itertools.product(range(-CMAX, CMAX+1), repeat=r):
        if all(ci == 0 for ci in c) or gcd(list(c)) != 1: continue
        if next(ci for ci in c if ci != 0) < 0: continue
        if r >= 3 and max(abs(ci) for ci in c) > 1: continue
        G = sum((ci*g for ci, g in zip(c, gens)), E(0))
        tot_dirs += 1
        B = {}; X = {}
        for k in range(1, 2*MMAX+1):
            X[k] = (k*G)[0]; B[k] = X[k].denominator()
        for M in range(2, MMAX+1):
            n = 2*M
            w = prim_part(B[n], [B[2*m] for m in range(1, M)])
            f = prim_part(B[n], [B[k] for k in range(1, n)])
            if w == 1 or f == 1: fails.append((Nv, c, n, w == 1, f == 1))
            if w > 1 and w.prime_to_m_part(2) == 1:
                only2.append((Nv, c, n))
                if X[n].valuation(2) > -2: par_fail.append((Nv, c, n))
        xs = [X[2*m] for m in range(1, MMAX+1)]
        S = set(xs)
        for i in range(MMAX):
            for k in range(i+1, MMAX):
                if (xs[i] + xs[k]) / 2 in S: aps.append((Nv, c, i+1, k+1))
print("кривых ранга>=1:", len(curves), " ранги:", {rr: sum(1 for _, q in curves if q == rr) for rr in set(q for _, q in curves)})
print("направлений:", tot_dirs, " MMAX =", MMAX, " CMAX =", CMAX)
print("провалы примитивности (N, c, n, слабая, полная):", fails)
print("AP среди x(2mG):", aps)
print("случаев, где единственный слабо-примитивный простой = 2:", len(only2), only2[:10])
print("из них нарушений v_2(x) <= -2:", par_fail)
print("время: %.1f c" % (time.time() - t0))
