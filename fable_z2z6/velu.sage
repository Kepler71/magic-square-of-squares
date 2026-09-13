R.<a,b> = QQ[]
K = R.fraction_field()
alpha = a^2+a*b+b^2; beta = a^2*b^2*(a+b)^2
E = EllipticCurve(K, [0, alpha^2, 0, 2*alpha*beta, beta^2])
S.<x> = K[]
phi = E.isogeny(x)      # ядро <T>, T=(0,beta): многочлен ядра x
Ep = phi.codomain()
print("E' ainvs:", [factor(c) if c != 0 else 0 for c in Ep.ainvs()])
print("disc(E) =", factor(E.discriminant()))
print("disc(E') =", factor(Ep.discriminant()))
print("c4(E) =", factor(E.c4()), "\nc4(E') =", factor(Ep.c4()))
print("c6(E) =", factor(E.c6()), "\nc6(E') =", factor(Ep.c6()))
f = phi.formal(prec=3)
print("phi formal:", f)
D3 = Ep.division_polynomial(3)
print("3-div poly E' factors:", [ (factor(g), e) for g,e in D3.factor() ])
# сдвиг E' так, чтобы T' была в x=0: найдём рациональный корень
rts = [r for r,e in D3.roots()]
print("rational roots of 3-div E':", [factor(r) for r in rts])
