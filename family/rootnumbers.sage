import functools
print = functools.partial(print, flush=True)
def sqf(n):
    n = ZZ(n); return sign(n) * prod(p^(e % 2) for p, e in n.abs().factor())
secs = [(17,7,13,65), (7,1,5,15), (23,7,17,34), (31,17,25,7), (41,1,29,609), (47,23,37,1295), (49,31,41,41), (73,17,53,265), (71,49,61,671), (89,23,65,910), (79,47,65,455)]
for (b, h, n, D) in secs:
    A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2)
    R = PolynomialRing(QQ, 'z'); z = R.gen(); Fr = R.fraction_field()
    T = PolynomialRing(QQ, 'T').gen()
    f = T*(T^2 - 1)*(A*T - C)*(C*T - A)
    Fz = R(Fr((1 - z)^6) * Fr(f((1 + z)/(1 - z))))
    c = Fz.leading_coefficient(); q4 = R(Fz/(c*z)); Ap = q4[2]
    beta = (C - A)/(C + A); dk = sqf((C - A)*(C + A))
    k = QuadraticField(dk, 'r'); sb = k(beta).sqrt()
    Ru = PolynomialRing(k, 'u'); u = Ru.gen()
    cub = Ru(D*c*(u^2 + Ap - 2*beta)*(u + 2*sb))
    a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
    E = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2]).global_minimal_model(semi_global=True)
    nf = pari.nfinit(pari(f'y^2 - {dk}'))
    Ep = pari.ellinit([pari(str(x.polynomial('y'))) for x in E.ainvs()], nf)
    w = pari.ellrootno(Ep)
    print(f"({b},{h},{n}) k=Q(sqrt({dk})): root number w(E1/k) = {w}  => rank {'odd (1,3,...)' if w == -1 else 'even (>= 2, since >= 1)'}")
