# Фактор H0 по совместной инволюции (k,u) -> (1/k,1/u). H0 инвариантна только относительно совместной инволюции,
# поэтому образ H0 при (k,u) -> (s,w) = (k+1/k, u+1/u) задаётся M = H0(k,u) * k^8 H0(1/k,u) / (k^8 u^8),
# инвариантным относительно обеих инволюций; H0 -> {M=0} имеет степень 2, т.е. {M=0} ~ H0/iota.
import functools, time
print = functools.partial(print, flush=True)
exec(preparse(open('genus_h0.sage').read().split('# ---------------- контрольные')[0]))   # R, H0, M1, M2
S.<s, w> = QQ[]
# симметрическая редукция: старший моном k^i u^j (i,j >= 0) -> вычесть c (k+1/k)^i (u+1/u)^j
L.<kk, uu> = LaurentPolynomialRing(QQ)
Pn = L(H0(kk, uu)) * L(H0(kk^-1, uu)) * uu^-8
G = S(0)
for _ in range(200):
    if Pn == 0: break
    mons = [(e, cf) for e, cf in zip(Pn.exponents(), Pn.coefficients())]
    (i, j), cf = max(mons, key=lambda m: (m[0][0], m[0][1]))
    assert i >= 0 and j >= 0, "не симметрично"
    G += cf * s^i * w^j
    Pn -= cf * (kk + kk^-1)^i * (uu + uu^-1)^j
assert Pn == 0
check = L(G(kk + kk^-1, uu + uu^-1)) * uu^8 - L(H0(kk, uu)) * L(H0(kk^-1, uu))
print("H0(k,u) H0(1/k,u) == u^8 M(k+1/k, u+1/u):", check == 0)
fl = G.factor(); print("M factorization (bideg, mult):", [(f_.degree(s), f_.degree(w), m_) for f_, m_ in fl])
G = fl[0][0] if len(fl) == 1 else G
print("M bidegree:", G.degree(s), G.degree(w), " irreducible over Q:", len(G.factor()) == 1 and G.factor()[0][1] == 1)
t0 = time.time(); g1 = Curve(G).geometric_genus(); print("M1 (Singular) genus of quotient:", g1, f"({time.time()-t0:.1f}s)")
for p in [101, 103]:
    t0 = time.time()
    print(f"M2 (function field, p={p}):", M2(R(G(k, u)), p), f"({time.time()-t0:.1f}s)")
print("Riemann–Hurwitz: number of fixed points on normalization r = 46 - 4 g' =", 46 - 4*g1)
