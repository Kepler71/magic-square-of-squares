# Явное отображение C_J2 -> E1 над k=Q(sqrt15) и образы рациональных точек; нижняя граница ранга.
import functools, itertools
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
gam_choices = [sb^3, -sb^3]          # gamma^2 = beta^3
Rt2.<T2> = QQ[]
f = T2*(T2^2 - 1)*(109*T2 - 229)*(229*T2 - 109)
def cover_points(d3t, H=60):
    # рациональные точки класса: p,q с d1 r^2 = p^2-q^2, d2 s^2 = p^2+q^2, d3 w^2 = 109p^2-229q^2 (поиск малой высоты)
    d1, d2, d3 = d3t; pts = []
    for q_ in range(1, H):
        for p_ in range(0, H):
            if gcd(p_, q_) != 1: continue
            vals = [(p_^2 - q_^2)/d1, (p_^2 + q_^2)/d2, (109*p_^2 - 229*q_^2)/d3, (229*p_^2 - 109*q_^2)/d3]
            if all(QQ(v).is_square() for v in vals):
                pts.append(QQ(p_)/q_)
    return pts
for d3t, D in [((5, 13, 65), 5*13*65*65), ((6, 26, 39), 6*26*39*39)]:
    print(f"\n===== class {d3t}")
    cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))
    a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
    E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])       # Xw = a3*u, Yw = a3*(D v1)
    Xs = cover_points(d3t)
    print("rational X on augmented cover (height < 60):", Xs)
    pts = []
    for X0 in Xs:
        t0 = X0^2
        W2 = f(t0) / D
        assert W2.is_square()
        for W0 in [W2.sqrt(), -W2.sqrt()]:
            z0 = (t0 - 1)/(t0 + 1); y0 = W0 * (1 - z0)^3        # D y0^2 = F(z0) = c z0 (...)
            assert D*y0^2 == c*z0*(z0^4 + A*z0^2 + B)
            for gam in gam_choices:
                u0 = z0 + beta/z0
                v = y0*(1 + gam/z0^3)
                if u0 == sb: continue
                v1 = v/(u0 - sb)
                if D*v1^2 == c*(u0^2 + A - 2*beta)*(u0 + 2*sb):
                    Pw = E1(a3*u0, a3*D*v1)
                    pts.append(Pw)
    pts = list(set(pts))
    nontors = [P for P in pts if P.order() == oo]
    print(f"images on E1(k): {len(pts)}, of infinite order: {len(nontors)}")
    if nontors:
        hs = [P.height() for P in nontors]
        print("   canonical heights:", [RR(h) for h in hs])
        M = E1.height_pairing_matrix(nontors)
        print("   rank of height pairing matrix (numerical):", M.rank() if M.nrows() < 2 else matrix(RR, M).rank(), "  det:", M.det())
    tors = E1.torsion_subgroup()
    print("   torsion E1(k):", tors.invariants())
