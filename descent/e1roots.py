def E1_roots(A, C, D):
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Frz = Rz.fraction_field()
    Rt = PolynomialRing(QQ, 'T'); T = Rt.gen()
    f = T*(T**2 - 1)*(A*T - C)*(C*T - A)
    Fz = Rz(Frz((1 - z)**6) * Frz(f((1 + z)/(1 - z))))
    cz = Fz.leading_coefficient(); q4 = Rz(Fz/(cz*z)); Ap, Bp = q4[2], q4[0]
    bet = QQ(C - A)/(C + A); assert Bp == bet**2 and Ap == -(1 + bet**2)
    dk = ((C - A)*(C + A)).squarefree_part()
    k = QuadraticField(dk, 'r'); sb = k(bet).sqrt()
    Ru = PolynomialRing(k, 'uu'); uu = Ru.gen()
    cub = Ru(D * cz * (uu**2 + Ap - 2*bet) * (uu + 2*sb))
    a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
    E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3**2])
    rts = [e for e, _ in E1.two_division_polynomial().roots()]
    # целая модель: x -> u**2 x; убрать общий квадратный множитель
    den = lcm([QQ(c).denominator() for t in rts for c in list(t)])
    rts = [t*den**2 for t in rts]
    g = gcd([ZZ(c) for t in rts for c in list(t)])
    s = prod(p**(e//2) for p, e in g.factor()) if g != 1 else 1
    rts = [t/s**2 for t in rts]
    E = EllipticCurve(k, [0, -sum(rts), 0, rts[0]*rts[1] + rts[0]*rts[2] + rts[1]*rts[2], -prod(rts)])
    assert E.is_isomorphic(E1)
    return k, rts
