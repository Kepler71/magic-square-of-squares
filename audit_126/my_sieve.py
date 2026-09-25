# Claude, 25.09: независимое решето случая A (свои простые, свой код) + точная проверка выживших;
# случай B — собственный вывод многочленов.
from sage.all import *
import time
r, s = 126, 451
S = [s, -s, r, -r, s - r, r - s, s + r, -(s + r)]
E = EllipticCurve([0, 0, 0, -21739732, 23939545056])
G1, G2 = E(-2883, 250305), E(-1610, 234024)
tors = E.torsion_points()
M1 = [[-3439800, 100333], [-2925, 9]]      # x(Q1) = (a z + b)/(c z + d)
a, b, c, d = M1[0][0], M1[0][1], M1[1][0], M1[1][1]
BOX = 72
cands = [(n, m, ti) for n in range(-BOX, BOX + 1) for m in range(-BOX, BOX + 1)
         if n * n + m * m < 5182 for ti in range(len(tors))]
print('кандидатов:', len(cands), flush=True)
def alive_mod_p(p, cands):
    Ep = E.change_ring(GF(p))
    g1, g2 = Ep(G1), Ep(G2); tp = [Ep(t) for t in tors]
    F = GF(p)
    mult1 = {n: n * g1 for n in range(-BOX, BOX + 1)}
    mult2 = {m: m * g2 for m in range(-BOX, BOX + 1)}
    keep = []
    for (n, m, ti) in cands:
        P = mult1[n] + mult2[m] + tp[ti]
        if P.is_zero():
            keep.append((n, m, ti)); continue             # z = M1^{-1}(inf) — не отбрасываем
        x = P[0]
        den = F(c) * x - F(a)                              # z = (b - d x)/(c x - a)
        if den == 0:
            keep.append((n, m, ti)); continue
        zz = (F(b) - F(d) * x) / den
        if all((F(1) + F(cc) * zz).is_square() for cc in S):
            keep.append((n, m, ti))
    return keep
t0 = time.time()
cur = cands
for p in [61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]:
    if E.discriminant() % p == 0 or (a * d - b * c) % p == 0: continue
    cur = alive_mod_p(p, cur)
    print(f'  p={p}: осталось {len(cur)}', flush=True)
    if len(cur) <= 4: break
print('выжившие:', cur, f'({time.time()-t0:.0f}s)', flush=True)
sols = []
for (n, m, ti) in cur:
    P = n * G1 + m * G2 + tors[ti]
    if P.is_zero(): print('  точка O -> z = inf (не конечный параметр)'); continue
    x = P[0]
    if c * x - a == 0: print('  z = inf'); continue
    zz = (b - d * x) / (c * x - a)
    ok = zz != 0 and all((1 + cc * zz) > 0 and (1 + cc * zz).is_square() for cc in S)
    print(f'  ({n},{m},T{ti}) -> z = {zz}: полный квадрат? {ok}')
    if ok: sols.append(str(zz))
print('СЛУЧАЙ A — невырожденных решений:', sols)
# Случай B: Q2 - eps*Q0 = T (кручение) => x(Q0) = x(Q2 - T), x-координаты через z
M0 = [[-2730602, 7490], [-577, 1]]; M2 = [[647388, 7490], [-126, 1]]
Rz = PolynomialRing(QQ, 'z'); zv = Rz.gen()
X0 = (M0[0][0]*zv + M0[0][1]) / (M0[1][0]*zv + M0[1][1])
X2 = (M2[0][0]*zv + M2[0][1]) / (M2[1][0]*zv + M2[1][1])
allz = set()
for t in tors:
    if t.is_zero():
        eq = X0 - X2                                         # Q2 = ±Q0
    else:
        tx = t[0]
        eq = X0 - (tx + (3 * tx**2 - 21739732) / (X2 - tx))   # x(Q2 - T), T = (tx, 0)
    num = eq.numerator()
    rts = [rt for rt, _ in num.roots(QQ)] if num != 0 else ['ТОЖДЕСТВЕННО']
    print(f'  T = {t}: рациональные z = {rts}')
    for rt in rts:
        if rt != 'ТОЖДЕСТВЕННО': allz.add(rt)
bad = [zz for zz in allz if zz != 0 and all((1 + cc * zz) > 0 and QQ(1 + cc * zz).is_square() for cc in S)]
print('СЛУЧАЙ B — невырожденных решений:', bad)
