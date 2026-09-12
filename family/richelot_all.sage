# Все рациональные Richelot-двойственные к C0: D Y^2 = t(t^2-1)(At-C)(Ct-A) (6 рациональных точек ветвления -> 15 разбиений).
# Для каждой: целочисленная примитивная модель (с точностью до квадратов), контроль изогении по многочленам Фробениуса,
# строка Magma RankBounds.
import sys, itertools
b0, h0, n0, D = [ZZ(x) for x in sys.argv[1:5]]
A, C = ZZ((h0^2 + n0^2)/2), ZZ((b0^2 + n0^2)/2)
R = PolynomialRing(QQ, 't'); t = R.gen()
f0 = D*t*(t^2 - 1)*(A*t - C)*(C*t - A)
roots = [QQ(0), oo, QQ(1), QQ(-1), QQ(C)/A, QQ(A)/C]
def quad(r1, r2):
    if r1 is oo: r1, r2 = r2, r1
    return (t - r1) if r2 is oo else (t - r1)*(t - r2)
def matchings(lst):
    if not lst: yield []; return
    for j in range(1, len(lst)):
        for m in matchings(lst[1:j] + lst[j+1:]): yield [(lst[0], lst[j])] + m
def normalize(f):
    # целые коэффициенты, убрать квадратный множитель из содержания
    den = lcm([c.denominator() for c in f.coefficients()])
    g = R(f * den^2)
    cont = gcd([ZZ(c) for c in g.coefficients()])
    sq = prod(p^(2*(e//2)) for p, e in cont.factor()) if cont != 1 else 1
    return R(g / sq)
def frob_ok(fa, fb):
    bad = set((2*fa.discriminant()*fb.discriminant()).numerator().prime_factors()) | set(ZZ(fa.leading_coefficient()).prime_factors()) | set(ZZ(fb.leading_coefficient()).prime_factors())
    return all(HyperellipticCurve(fa.change_ring(GF(p))).frobenius_polynomial() == HyperellipticCurve(fb.change_ring(GF(p))).frobenius_polynomial()
               for p in primes(5, 120) if p not in bad)
lines = []
f0n = normalize(f0)
lines.append(("C0", f0n))
for idx, m in enumerate(matchings(roots)):
    F = [quad(a_, b_) for a_, b_ in m]
    lead = f0 // prod(F); F[0] = lead * F[0]
    co = lambda P: [P[0], P[1], P[2]]
    delta = matrix(QQ, [co(F[0]), co(F[1]), co(F[2])]).det()
    if delta == 0: continue
    G1 = F[1].derivative()*F[2] - F[1]*F[2].derivative()
    G2 = F[2].derivative()*F[0] - F[2]*F[0].derivative()
    G3 = F[0].derivative()*F[1] - F[0]*F[1].derivative()
    f1 = normalize(delta*G1*G2*G3)
    ok = frob_ok(f0n, f1)
    lab = "{" + " ".join(f"({a_},{b_})" for a_, b_ in m) + "}"
    lines.append((f"K{idx} {lab} isog={ok}", f1))
print(f"// section ({b0},{h0},{n0}), D = {D}: {len(lines)} models (C0 + rational Richelot duals)")
for lab, f in lines:
    print(f"// {lab}")
    print(f"f:={f}; lo,hi:=RankBounds(Jacobian(HyperellipticCurve(f))); print \"({b0},{h0},{n0}) {lab.split()[0]}\", lo, hi; delete lo, hi;".replace('t^', 't^'))
