k.<r> = QuadraticField(165)
e, q = 163724, 2939651
E1 = EllipticCurve(k, [0, -e*r, 0, -q^2, e*r*q^2])            # корни e r, ±q
E0 = EllipticCurve(QQ, [0, -165*e, 0, -165*q^2, 165^2*e*q^2])  # корни 165 e, ±r q
print("E1 = E0 twisted by sqrt165 over k:", E1.is_isomorphic(E0.change_ring(k).quadratic_twist(r)))
E0m = E0.minimal_model()
print("E0 minimal:", E0m.ainvs()); N0 = E0m.conductor(); print("N(E0) =", N0.factor())
print("root number E0:", E0m.root_number(), "; E0 rank bounds (PARI):", pari(E0m.ainvs()).ellinit().ellrank()[:2])
Nk = E1.conductor(); print("conductor E1/k: norm =", Nk.norm().factor())
A = Nk.norm() * k.discriminant()^2
print(f"analytic conductor of L(E1/k,s) (degree 4 over Q): {A} ~ 10^{RR(log(A,10)):.1f}; sqrt ~ 10^{RR(log(A,10)/2):.1f}")
print("root number E1/k:", E1.root_number())
