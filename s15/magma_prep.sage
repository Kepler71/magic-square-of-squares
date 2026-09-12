import functools
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
D = 5*13*65*65
cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))
a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
Em = E1.global_minimal_model(semi_global=True); iso = E1.isomorphism_to(Em)
Rt2.<T2> = QQ[]
f = T2*(T2^2 - 1)*(109*T2 - 229)*(229*T2 - 109)
t0 = QQ(9)/4; W0 = (f(t0)/D).sqrt()
for s_ in (1, -1):
    for g_ in (sb^3, -sb^3):
        z0 = (t0 - 1)/(t0 + 1); y0 = s_*W0*(1 - z0)^3
        u0 = z0 + beta/z0; v1 = y0*(1 + g_/z0^3)/(u0 - sb)
        X_, Y_ = a3*u0, a3*D*v1
        if Y_^2 == X_^3 + a2*X_^2 + a1*a3*X_ + a0*a3^2:
            G = iso(E1(X_, Y_)); break
print("model:", Em.ainvs())
print("j:", Em.j_invariant())
print("G0:", G.xy())
print("2-torsion x:", [r for r in Em.two_division_polynomial().roots(multiplicities=False)])
print("conductor norm:", Em.conductor().norm().factor())
