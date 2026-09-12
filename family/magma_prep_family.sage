import sys
b0, h0, n0, D = [ZZ(x) for x in sys.argv[1:5]]; t0 = QQ(sys.argv[5])
A, C = ZZ((h0^2 + n0^2)/2), ZZ((b0^2 + n0^2)/2)
def sqf(n):
    n = ZZ(n); return sign(n) * prod(pr^(e % 2) for pr, e in n.abs().factor())
R = PolynomialRing(QQ, 't'); t = R.gen()
f0 = D*t*(t^2 - 1)*(A*t - C)*(C*t - A)
F1 = (D*A*C)*t; F2 = t^2 - 1; F3 = (A*t - C)*(C*t - A)/(A*C)
assert F1*F2*F3 == f0
co = lambda P: [P[0], P[1], P[2]]
delta = matrix(QQ, [co(F1), co(F2), co(F3)]).det()
G1 = F2.derivative()*F3 - F2*F3.derivative(); G2 = F3.derivative()*F1 - F3*F1.derivative(); G3 = F1.derivative()*F2 - F1*F2.derivative()
f1 = delta*G1*G2*G3
# убрать квадратный множитель из коэффициента
cf = QQ(f1.leading_coefficient()); num, den = cf.numerator(), cf.denominator()
sq = (num.squarefree_part()*den.squarefree_part())
f1n = R(f1 * sq / cf) * 1
f1n = f1n / f1n.content() if False else f1n
print("C0: Y^2 =", f0)
print("Richelot dual (matching {0,oo},{1,-1},{C/A,A/C}): W^2 =", f1n.factor())
ok = all(HyperellipticCurve(f0.change_ring(GF(p))).frobenius_polynomial() == HyperellipticCurve(f1n.change_ring(GF(p))).frobenius_polynomial()
         for p in primes(29, 200) if (2*D*A*C*(C^2 - A^2)*delta.numerator()*delta.denominator()) % p != 0)
print("isogeny control (Frobenius equal, p in [29,200)):", ok)
