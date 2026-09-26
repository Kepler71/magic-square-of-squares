# Безусловный вариант шага n=30 у Verzobio: вместо (p1 или p2) используем (p1 или p3) -- обе формы без вещественных
# корней, PARI thue тогда перебирает решения (безусловно, см. документацию thue/thueinit). Из Psi30 | 2^434 и
# p1*p3 | Psi30 следует v2(p1)+v2(p3) <= 434, значит p1 = +-2^k или p3 = +-2^k с k <= 217.
# Нетривиальными считаются решения, отличные от X=0, Y=0, X=+-Y (как у Verzobio).
import time
from sage.all import EllipticCurve, QQ, ZZ, PolynomialRing, lcm, pari, RealField
E = EllipticCurve(QQ, [1, 0]); Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
p = lambda n: Rx(E.division_polynomial(n, two_torsion_multiplicity=0))
Q, r = p(30).quo_rem(lcm([p(15), p(10), p(6)])); assert r == 0
RX = PolynomialRing(QQ, 'X')
Psi = RX([Q.list()[2*i] for i in range(Q.degree()//2 + 1)])
facs = []
for f, e in Psi.factor():
    g = (f * f.denominator()).change_ring(ZZ); g = g // g.content()
    facs.append(g)
p1 = [g for g in facs if g.degree() == 16][0]
p3 = [g for g in facs if g.degree() == 32 and g[31] == -416][0]
for name, g in (("p1", p1), ("p3", p3)):
    assert len(g.change_ring(QQ).roots(ring=RealField(200))) == 0
    tnf = pari.thueinit(pari(str(g).replace('X', 'x')), 0)
    t0 = time.time(); nontriv = []; total = 0
    for k in range(0, 218):
        for s in (1, -1):
            sols = pari.thue(tnf, s * 2**k)
            for sol in sols:
                X, Y = ZZ(sol[0]), ZZ(sol[1]); total += 1
                if X.gcd(Y) != 1: continue
                if X == 0 or Y == 0 or X == Y or X == -Y: continue
                nontriv.append((k, s, X, Y))
    print(name, "deg", g.degree(), ": всего решений", total, "нетривиальных взаимно простых:", nontriv, "%.1fs" % (time.time()-t0), flush=True)
