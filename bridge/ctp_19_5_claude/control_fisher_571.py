"""Контроль формулы Фишера на его собственном примере 3.4: E = 571a1, y^2+y = x^3-x^2-929x-10595.
Здесь 2-кручение тривиально, L = кубическое ПОЛЕ, поэтому пишется обобщённая версия
той же формулы (та же арифметика H, m, gamma, тот же набор мест, те же символы Гильберта).
Известно независимо (Cremona/LMFDB): rank E(Q) = 0, |Sha| = 4, Sha[2] = (Z/2)^2,
значит Sel^2 = (Z/2)^2 и спаривание Касселса-Тейта НЕВЫРОЖДЕНО: <s1,s2> = -1, <si,si> = +1.
Этот контроль ловит «ложные пропуски» (нули там, где должно быть -1).
"""
from sage.all import *
import random
random.seed(571); set_random_seed(571)

def quartI(g):
    a, b, c, d, e = g; return 12*a*e - 3*b*d + c*c
def quartJ(g):
    a, b, c, d, e = g; return 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3
def hessian(g):
    a, b, c, d, e = g
    return [3*b*b-8*a*c, 4*(b*c-6*a*d), 2*(2*c*c-24*a*e-3*b*d), 4*(c*d-6*b*e), 3*d*d-8*c*e]
def ev(cs, x, z):
    k = len(cs)-1; return sum(cs[i]*x**(k-i)*z**i for i in range(k+1))

E = EllipticCurve([0, -1, 1, -929, -10595])
print('E =', E, 'conductor', E.conductor())
I = E.c4(); J = 2*E.c6()
print('I,J =', I, J, ' (Фишер: 44608, 18842960)')
assert (I, J) == (44608, 18842960)
DISC = 16*(4*I**3 - J**2)/27
print('Disc =', DISC, ' (Фишер: -2338816 = -2^12*571)')
assert DISC == -2338816

Rx = PolynomialRing(QQ, 'X'); Xv = Rx.gen()
fpoly = Xv**3 - 3*I*Xv + J
assert fpoly.is_irreducible()
L = NumberField(fpoly, 'ph'); ph = L.gen()

def z_inv(g):
    a, b, c, _, _ = g; return (4*a*ph + 3*b*b - 8*a*c)/3

def H_form(g):
    h = hessian(g)
    G = [(4*ph*g[j] + h[j])/3 for j in range(5)]
    H = [G[0], G[1]/2, G[2]/6 + QQ(2)/9*(I - ph*ph)]
    Hsq = [sum(H[i]*H[j] for i in range(3) for j in range(3) if i+j == k) for k in range(5)]
    assert Hsq == [G[0]*v for v in G], 'G(1,0)G=H^2 нарушено'
    return H

def val_unit(q, p):
    q = QQ(q); num, den = ZZ(q.numerator()), ZZ(q.denominator()); k = 0
    while num % p == 0: num //= p; k += 1
    while den % p == 0: den //= p; k -= 1
    return k, QQ(num)/QQ(den)
def is_loc_square(q, pl):
    q = QQ(q)
    if q == 0: return False
    if pl == 'real': return q > 0
    p = ZZ(pl); k, u = val_unit(q, p)
    if k % 2: return False
    r = u.numerator()*ZZ(u.denominator()).inverse_mod(8 if p == 2 else p)
    return (r % 8 == 1) if p == 2 else kronecker(r % p, p) == 1
def my_hilbert(a, b, pl):
    a, b = QQ(a), QQ(b)
    if pl == 'real': return -1 if (a < 0 and b < 0) else 1
    p = ZZ(pl); va, ua = val_unit(a, p); vb, ub = val_unit(b, p)
    if p == 2:
        ra = ua.numerator()*ZZ(ua.denominator()).inverse_mod(8) % 8
        rb = ub.numerator()*ZZ(ub.denominator()).inverse_mod(8) % 8
        e = ((ra-1)//2)*((rb-1)//2) + va*((rb*rb-1)//8) + vb*((ra*ra-1)//8)
        return -1 if e % 2 else 1
    ra = ua.numerator()*ZZ(ua.denominator()).inverse_mod(p) % p
    rb = ub.numerator()*ZZ(ub.denominator()).inverse_mod(p) % p
    res = 1
    if (va*vb) % 2 and p % 4 == 3: res = -res
    if vb % 2 and kronecker(ra, p) == -1: res = -res
    if va % 2 and kronecker(rb, p) == -1: res = -res
    return res
def supp(q):
    q = QQ(q)
    if q == 0: return set()
    return set(p for p, _ in factor(q.numerator())) | set(p for p, _ in factor(q.denominator()))
def local_points(g, pl, need=1, avoid=None):
    res = []; p = None if pl == 'real' else ZZ(pl)
    cands = [(ZZ(1), ZZ(0)), (ZZ(0), ZZ(1))]
    for i in range(-120, 121): cands.append((ZZ(i), ZZ(1)))
    for j in range(-120, 121): cands.append((ZZ(1), ZZ(j)))
    if p is not None:
        for k in range(1, 7):
            for i in range(-40, 41):
                cands.append((ZZ(i), p**k)); cands.append((ZZ(i)*p**k, ZZ(1)))
    for x, z in cands:
        if (x, z) == (0, 0): continue
        v = ev(g, QQ(x), QQ(z))
        if v == 0 or not is_loc_square(v, pl): continue
        if avoid is not None and ev(avoid, QQ(x), QQ(z)) == 0: continue
        if (QQ(x), QQ(z)) in res: continue
        res.append((QQ(x), QQ(z)))
        if len(res) >= need: return res
    return res

Rp = PolynomialRing(QQ, 'x')
cover = pari(E).ell2cover()
quarts = []
for c in cover:
    q = Rp(c[0]); g = [QQ(q[4-j]) for j in range(5)]
    lam2 = QQ(quartJ(g)*I)/QQ(J*quartI(g))
    assert QQ(lam2).is_square(), lam2
    g = [v/lam2 for v in g]
    assert quartI(g) == I and quartJ(g) == J
    quarts.append(g)
print('ell2cover дал', len(quarts), 'квартик:')
for g in quarts: print('   ', g, ' z(g) =', z_inv(g))

def pair(gA, gB, gC, reps=2):
    zA, zB, zC = z_inv(gA), z_inv(gB), z_inv(gC)
    pr = zA*zB*zC
    ok, mroot = pr.is_square(root=True)
    assert ok, 'z1z2z3 не квадрат в L'
    H = H_form(gA)
    coef = (mroot/zA)
    gam = [(coef*H[j]).list()[2] for j in range(3)]     # коэффициент при ph^2
    a = gB[0]; assert a != 0
    S = set([2, 3, 5, 7]) | supp(DISC) | supp(a)
    for c in list(gA)+list(gam): S |= supp(c)
    vals = []
    for r in range(reps):
        prod_ = 1; nontriv = []
        for pl in [str(p) for p in sorted(S)] + ['real']:
            pts = local_points(gA, pl, need=r+1, avoid=gam)
            assert len(pts) >= r+1, ('нет локальной точки', pl)
            x, z = pts[r]
            hb = my_hilbert(a, ev(gam, x, z), pl)
            prod_ *= hb
            if hb == -1: nontriv.append(pl)
        vals.append((prod_, nontriv))
    assert len(set(v[0] for v in vals)) == 1, vals
    return vals[0][0], vals[0][1], gam, mroot

# находим тройки: g3 с z(g3) ~ z(g1)z(g2)
import itertools
print()
for i, j in itertools.combinations_with_replacement(range(len(quarts)), 2):
    gA, gB = quarts[i], quarts[j]
    target = z_inv(gA)*z_inv(gB)
    gC = None
    for k, g in enumerate(quarts):
        if (target*z_inv(g)).is_square():
            gC = g; break
    if gC is None:
        print('пара', i, j, ': нет g3 среди квартик покрытия'); continue
    val, nt, gam, mroot = pair(gA, gB, gC)
    print('<g%d,g%d> = %+d   нетривиальные места: %s   gamma = %s' % (i, j, val, nt, gam))
    if (i, j) == (0, 1) or (i, j) == (1, 0):
        print('    m =', mroot)
print()
print('Ожидание (Sha(571a1)[2]=(Z/2)^2, спаривание невырождено):')
print('  диагональ 0 (+1), недиагональ нетривиальна (-1), вклад только вещественного места.')
