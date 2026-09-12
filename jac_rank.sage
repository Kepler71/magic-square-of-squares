# Jac(C) для C: Y^2 = f(X) = g(X^2), g(t) = (109t-229)(t^2-1).  Jac(C) ~ E+ x E-.
# E+ = C/<(X,Y)->(-X,Y)>:   Y^2 = g(t),            t = X^2
# E- = C/<(X,Y)->(-X,-Y)>:  V^2 = t^3 g(1/t) = (229t-109)(t^2-1),  t = 1/X^2, V = Y/X^3
R.<t> = QQ[]
g = (109*t - 229)*(t^2 - 1)
def ell_from_cubic(c):
    assert c.degree() == 3
    a3, a2, a1, a0 = [c[i] for i in (3, 2, 1, 0)]
    return EllipticCurve([0, a2, 0, a1*a3, a0*a3^2])   # X = a3 x, Y = a3 y
gm = R(t^3 * g(1/t))
print("E- cubic:", gm.factor())
Eplus, Eminus = ell_from_cubic(g).minimal_model(), ell_from_cubic(gm).minimal_model()
E1 = EllipticCurve([0, -556, 0, 73684, 0])
tot = 0
for name, E in [("E+", Eplus), ("E-", Eminus)]:
    print(f"{name}: {E.ainvs()}  cond {E.conductor().factor()}  ~E1: {E.is_isogenous(E1)}  tors {E.torsion_subgroup().invariants()}")
    print(f"    rank_bounds {E.rank_bounds()}  analytic {E.analytic_rank()}")
    r = E.rank(proof=True); tot += r
    print(f"    rank (proof=True) = {r}   gens {E.gens(proof=True)}")
print("rank Jac(C)(Q) = rank E+ + rank E- =", tot)
# контроль разложения: #C(F_p) = p + 1 - a_p(E+) - a_p(E-)
H = HyperellipticCurve((109*t^2 - 229)*(t^4 - 1))
bad = set((2*Eplus.conductor()*Eminus.conductor()*109*229).prime_factors())
checked = 0; ok = True
for p in primes(3, 300):
    if p in bad: continue
    Np = H.change_ring(GF(p)).count_points(1)[0]
    checked += 1
    if Np != p + 1 - Eplus.ap(p) - Eminus.ap(p):
        ok = False; print("mismatch at", p)
print(f"control #C(F_p) = p+1-a_p(E+)-a_p(E-): {ok} on {checked} good primes < 300")
# рациональные точки C из записки должны отображаться в точки E+ и E-
pts = [1, 5, 3/2, 2/3, 1/5, 17/7]
print("known X on C:", [(x, (109*x^2-229)*(x^4-1)).is_square() if False else ((109*x^2-229)*(x^4-1)).is_square() for x in pts])
