secs = [(7,1,5,15), (17,7,13,65), (23,7,17,34), (31,17,25,7), (41,1,29,609), (47,23,37,1295), (49,31,41,41), (73,17,53,265)]
R = PolynomialRing(QQ, 't'); t = R.gen()
for (b, h, n, D) in secs:
    A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2); S2 = ZZ((A^2 + C^2)/2)
    f0 = D*t*(t^2 - 1)*(A*t - C)*(C*t - A)
    f1 = D*(t^2 - 1)*(t^2 + 1)*(S2*t^2 - 2*A*C*t + S2)
    bad = set((2*D*A*C*(C^2 - A^2)*S2*(S2^2 - (A*C)^2)).prime_factors())
    ok = all(HyperellipticCurve(f0.change_ring(GF(p))).frobenius_polynomial() == HyperellipticCurve(f1.change_ring(GF(p))).frobenius_polynomial()
             for p in primes(3, 300) if p not in bad)
    # также сверка с явной конструкцией Richelot (разбиение {0,oo},{1,-1},{C/A,A/C})
    F1 = (D*A*C)*t; F2 = t^2 - 1; F3 = (A*t - C)*(C*t - A)/(A*C)
    co = lambda P: [P[0], P[1], P[2]]
    delta = matrix(QQ, [co(F1), co(F2), co(F3)]).det()
    G1 = F2.derivative()*F3 - F2*F3.derivative(); G2 = F3.derivative()*F1 - F3*F1.derivative(); G3 = F1.derivative()*F2 - F1*F2.derivative()
    q, r = (delta*G1*G2*G3).quo_rem(f1)
    print(f"({b},{h},{n}) D={D}: Frobenius(C0)=Frobenius(C1) on p<300: {ok};  Richelot dual = f1 * {q.factor() if q.degree()==0 else q}  square: {r == 0 and q.degree() == 0 and QQ(q).is_square()}")
