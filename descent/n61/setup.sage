# E1 для сечения (71,49,61) над k = Q(sqrt165): построение как в family/pipeline.sage, сверка с моделью Codex.
A, C, D = 3061, 4381, 671
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Frz = Rz.fraction_field()
Rt = PolynomialRing(QQ, 'T'); T = Rt.gen()
f = T*(T^2 - 1)*(A*T - C)*(C*T - A)
Fz = Rz(Frz((1 - z)^6) * Frz(f((1 + z)/(1 - z))))
cz = Fz.leading_coefficient(); q4 = Rz(Fz/(cz*z)); Ap, Bp = q4[2], q4[0]
bet = QQ(C - A)/(C + A); assert Bp == bet^2 and Ap == -(1 + bet^2)
k.<r> = QuadraticField(165)
sb = k(bet).sqrt(); print("sqrt(beta') =", sb)
Ru.<uu> = PolynomialRing(k)
cub = Ru(D * cz * (uu^2 + Ap - 2*bet) * (uu + 2*sb))
a3, a2, a1, a0 = [cub[j] for j in (3, 2, 1, 0)]
E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
roots = [e for e, _ in E1.two_division_polynomial().roots()]
print("E1 roots of 2-division:", roots)
# модель Codex: y^2 = (x - 163724 r)(x^2 - 2939651^2)  (после x = 244^2 X)
Ec = EllipticCurve(k, [0, -163724*r, 0, -2939651^2, 163724*r*2939651^2])
print("Codex model roots:", [e for e, _ in Ec.two_division_polynomial().roots()])
print("j equal:", E1.j_invariant() == Ec.j_invariant(), " isomorphic over k:", E1.is_isomorphic(Ec))
q = 175015061936
Eo = EllipticCurve(k, [0, -9747472064*r, 0, -q^2, 9747472064*r*q^2])
G = Eo(-173017629136, -784866450593280*r - 4316765478263040)
print("Codex G on original model: ok;  isomorphic Eo ~ Ec:", Eo.is_isomorphic(Ec))
print("disc factor (Ec):", Ec.discriminant().norm().factor())
