# Verzobio (arXiv:2001.09634v3, Prop. 3.17): Psi_30(X,Y) раскладывается на p1 (deg 16), p2 (deg 32), p3 (deg 32), p4 (deg 64);
# уравнения p1 = +-2^k, p2 = +-2^k решены PARI `thue` без указания флага. По документации PARI при flag=0 результат
# безусловен, в частности, если форма не имеет вещественных корней (перебор) -- проверяем вещественные корни p1, p2.
import time
from sage.all import EllipticCurve, QQ, PolynomialRing, lcm
E = EllipticCurve(QQ, [1, 0])            # a = 1:  y^2 = x^3 + x
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
def p(n):   # psi_n/psi_2 для чётного n, psi_n для нечётного (с точностью до константы)
    return Rx(E.division_polynomial(n, two_torsion_multiplicity=0))
t0 = time.time()
P30 = p(30); L = lcm([p(15), p(10), p(6)])
Q, r = P30.quo_rem(L); assert r == 0
print("deg Psi30(x) =", Q.degree(), " (ожидается 288 = 2*144)", "%.1fs" % (time.time()-t0))
# чётность по x и переход к X = x^2
assert all(c == 0 for i, c in enumerate(Q.list()) if i % 2 == 1)
RX = PolynomialRing(QQ, 'X'); X = RX.gen()
Psi = RX([Q.list()[2*i] for i in range(Q.degree()//2 + 1)])
fac = Psi.factor()
print("разложение Psi30(X,1): степени", [ (f.degree(), e) for f, e in fac ])
for f, e in fac:
    rr = f.real_roots() if hasattr(f, 'real_roots') else None
    nreal = len(f.change_ring(QQ).roots(ring=__import__('sage.all', fromlist=['RealField']).RealField(200)))
    print("  фактор степени", f.degree(), ": вещественных корней", nreal)
# Какой из двух факторов степени 32 -- это p2 Verzobio (коэффициент при X^31 Y равен -4256), а какой -- p3 (-416)?
from sage.all import ZZ, RealField
for f, e in fac:
    g = (f * f.denominator()).change_ring(ZZ); g = g // g.content()
    if g.leading_coefficient() < 0: g = -g
    nreal = len(g.change_ring(QQ).roots(ring=RealField(200)))
    print("  deg", g.degree(), "старший коэф.", g.leading_coefficient(), "коэф. X^(d-1):", g[g.degree()-1], "вещ. корней", nreal)
