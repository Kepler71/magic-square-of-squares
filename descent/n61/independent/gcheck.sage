k.<r> = QuadraticField(165)
q = 175015061936
Eo = EllipticCurve(k, [0, -9747472064*r, 0, -q^2, 9747472064*r*q^2])
G = Eo(-173017629136, -784866450593280*r - 4316765478263040)
print("G torsion?", G.has_finite_order())
print("h(G) ~", G.height(precision=100))
