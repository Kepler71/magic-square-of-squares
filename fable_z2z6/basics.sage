R.<a,b> = QQ[]
c = a+b
M = b^3*(2*a+b); N = a^3*(a+2*b)
K = R.fraction_field()
E = EllipticCurve(K, [0, M+N, 0, M*N, 0])
x0 = a^2*b^2; y0 = a^2*b^2*(a+b)^2
T = E(x0, y0)
print("T on curve:", T in E, " 3T = O:", 3*T == E(0))
# shifted model y^2 = (x+a^2b^2)(x+b^2c^2)(x+a^2c^2) = x^3 + (alpha x + beta)^2
alpha = a^2+a*b+b^2; beta = a^2*b^2*c^2
S.<x> = K[]
lhs = (x+a^2*b^2)*(x+b^2*c^2)*(x+a^2*c^2); rhs = x^3 + (alpha*x+beta)^2
print("model identity:", lhs == rhs)
# numeric: (a,b)=(2,3)
for (aa,bb) in [(2,3),(1,2),(3,5),(7,2)]:
    Ei = EllipticCurve(QQ,[0, M(aa,bb)+N(aa,bb), 0, M(aa,bb)*N(aa,bb), 0])
    Ti = Ei(aa^2*bb^2, aa^2*bb^2*(aa+bb)^2)
    print((aa,bb), Ei.torsion_order(), Ti.order(), "rank(pari)=", Ei.rank())
    phi = Ei.isogeny(Ti)
    Ep = phi.codomain()
    print("   E' =", Ep.ainvs(), " E'.tors=", Ep.torsion_order(), " E'.rank=", Ep.rank())
    # is E' of form y^2 = x^3 - 3(...)^2 ? check if E' has 3-torsion over Q(sqrt(-3))
    print("   3-division of E' roots:", Ep.division_polynomial(3).roots(QuadraticField(-3,'w')))
