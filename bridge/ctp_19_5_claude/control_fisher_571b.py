"""Контроль формулы Фишера на его примере 3.4 (571a1) — версия с построением третьей квартики.
L = Q[Th]/f(Th) — кубическое ПОЛЕ; конструкция 2-накрытия из класса делается общей
(та же, что у меня для полного 2-кручения, но без покомпонентного разложения).
"""
from sage.all import *
import itertools, random
random.seed(5711); set_random_seed(5711)

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
I = E.c4(); J = 2*E.c6(); b2 = E.b2()
DISC = 16*(4*I**3 - J**2)/27
assert (I, J, DISC) == (44608, 18842960, -2338816)
print('571a1: I,J,Disc =', I, J, DISC, ' b2 =', b2)

# модель y^2 = f(x) с теми же c4,c6 (пополнение квадрата)
Rx = PolynomialRing(QQ, 'x'); xg = Rx.gen()
f = xg**3 + QQ(E.b2())/4*xg**2 + QQ(E.b4())/2*xg + QQ(E.b6())/4
assert EllipticCurve([0, f[2], 0, f[1], f[0]]).c4() == I
f2, f1, f0 = f[2], f[1], f[0]
print('f =', f)

# --- арифметика L = Q[Th]/f  на координатных тройках (над любым коммутативным кольцом)
def mulL(u, v):
    c = [0]*5
    for i in range(3):
        for j in range(3):
            c[i+j] += u[i]*v[j]
    # Th^3 = -f2 Th^2 - f1 Th - f0 ;  Th^4 = Th*Th^3
    for d in (4, 3):
        if c[d] == 0: continue
        t = c[d]; c[d] = 0
        c[d-1] += -f2*t; c[d-2] += -f1*t; c[d-3] += -f0*t
        if d == 4:
            # снова свернуть возникший Th^3
            t3 = c[3]; c[3] = 0
            c[2] += -f2*t3; c[1] += -f1*t3; c[0] += -f0*t3
    return [c[0], c[1], c[2]]

K = NumberField(f, 'Th'); Th = K.gen()
def tolist(a): return list(K(a))
def fromlist(c): return c[0] + c[1]*Th + c[2]*Th**2
for _ in range(200):
    u = [QQ(random.randint(-9, 9)) for _ in range(3)]
    v = [QQ(random.randint(-9, 9)) for _ in range(3)]
    assert mulL(u, v) == tolist(fromlist(u)*fromlist(v))
print('OK: умножение в L сверено с NumberField')

Lphi = NumberField(xg**3 - 3*I*xg + J, 'ph'); ph = Lphi.gen()
def z_inv(g):
    a, b, c, _, _ = g; return (4*a*ph + 3*b*b - 8*a*c)/3
def H_form(g):
    h = hessian(g)
    G = [(4*ph*g[j] + h[j])/3 for j in range(5)]
    H = [G[0], G[1]/2, G[2]/6 + QQ(2)/9*(I - ph*ph)]
    Hsq = [sum(H[i]*H[j] for i in range(3) for j in range(3) if i+j == k) for k in range(5)]
    assert Hsq == [G[0]*v for v in G]
    return H
# изоморфизм Q[ph] -> Q[Th]:  ph = -(12 Th + b2)
iso_ph = -(12*Th + b2)
assert (iso_ph**3 - 3*I*iso_ph + J) == 0
def phi_to_theta(a):
    c = list(Lphi(a)); return c[0] + c[1]*iso_ph + c[2]*iso_ph**2

def quartic_from_delta(delta):
    """delta in K=Q[Th];  строим квартику 2-накрытия x-Th = delta t^2/u^2."""
    dl = tolist(delta)
    P = PolynomialRing(QQ, ['t0', 't1', 't2']); t0, t1, t2 = P.gens()
    tt = [t0, t1, t2]
    sq = mulL(tt, tt)
    dt2 = mulL([P(c) for c in dl], sq)
    Q2 = dt2[2]
    con = Conic(P(Q2))
    ok, pt = con.has_rational_point(point=True)
    assert ok, 'коника без точки'
    par = con.parametrization(pt)[0].defining_polynomials()
    Rxz = PolynomialRing(QQ, ['xx', 'zz']); xx, zz = Rxz.gens()
    T = [Rxz(p) for p in par]
    assert Q2(T[0], T[1], T[2]) == 0
    dt2b = mulL([Rxz(c) for c in dl], mulL(T, T))
    assert dt2b[2] == 0
    gpoly = -dt2b[1]
    gc = [QQ(gpoly.coefficient({xx: 4-i, zz: i})) for i in range(5)]
    assert gpoly == sum(gc[i]*xx**(4-i)*zz**i for i in range(5))
    lam2 = QQ(quartJ(gc)*I)/QQ(J*quartI(gc))
    gc = [c/lam2 for c in gc]
    assert quartI(gc) == I and quartJ(gc) == J
    return gc

# --- локальная арифметика
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
    for i in range(-150, 151): cands.append((ZZ(i), ZZ(1)))
    for j in range(-150, 151): cands.append((ZZ(1), ZZ(j)))
    if p is None:
        # вещественное место: берём середины интервалов между вещественными корнями
        pol = Rx([g[4-i] for i in range(5)])
        rts = sorted([r for r, _ in pol.roots(RR)])
        marks = [RR(rts[0]-1)] if rts else [RR(0)]
        for i in range(len(rts)-1): marks.append((rts[i]+rts[i+1])/2)
        if rts: marks.append(RR(rts[-1]+1))
        for mk in marks:
            for den in [1, 2, 4, 8, 16, 64, 256, 1024]:
                cands.append((ZZ((mk*den).round()), ZZ(den)))
        for i in range(-4000, 4001, 7):
            cands.append((ZZ(i), ZZ(1000)))
    if p is not None:
        for k in range(1, 8):
            for i in range(-50, 51):
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

def pair(gA, gB, gC, reps=2):
    zA, zB, zC = z_inv(gA), z_inv(gB), z_inv(gC)
    ok, mroot = (zA*zB*zC).is_square(root=True)
    assert ok, 'z1z2z3 не квадрат'
    H = H_form(gA); coef = mroot/zA
    gam = [(coef*H[j]).list()[2] for j in range(3)]
    a = gB[0]; assert a != 0
    S = set([2, 3, 5, 7]) | supp(DISC) | supp(a)
    for c in list(gA)+list(gam): S |= supp(c)
    out = []
    for r in range(reps):
        pr = 1; nt = []
        for pl in [str(q) for q in sorted(S)] + ['real']:
            pts = local_points(gA, pl, need=r+1, avoid=gam)
            assert len(pts) >= r+1, ('нет локальной точки', pl)
            x, z = pts[r]
            hb = my_hilbert(a, ev(gam, x, z), pl); pr *= hb
            if hb == -1: nt.append(pl)
        out.append((pr, nt))
    assert len(set(o[0] for o in out)) == 1, out
    return out[0][0], out[0][1], gam, mroot

Rp = PolynomialRing(QQ, 'x')
quarts = []
for c in pari(E).ell2cover():
    q = Rp(c[0]); g = [QQ(q[4-j]) for j in range(5)]
    lam2 = QQ(quartJ(g)*I)/QQ(J*quartI(g)); assert QQ(lam2).is_square()
    g = [v/lam2 for v in g]
    assert quartI(g) == I and quartJ(g) == J
    quarts.append(g)
print('базис Sel^2 из ell2cover:', quarts)

# третий нетривиальный класс и тривиальный класс -> строим квартики сами
def build(cls_in_phi):
    dl = phi_to_theta(cls_in_phi)
    g = quartic_from_delta(dl)
    assert (z_inv(g)*cls_in_phi).is_square(), 'z(g) не в нужном классе!'
    return g

z0, z1 = z_inv(quarts[0]), z_inv(quarts[1])
g_sum = build(z0*z1)
g_triv = build(Lphi(1))
print('квартика суммы:', g_sum)
print('квартика тривиального класса:', g_triv)
allq = quarts + [g_sum, g_triv]
names = ['s0', 's1', 's0+s1', '0']

print()
for i in range(4):
    for j in range(4):
        zi, zj = z_inv(allq[i]), z_inv(allq[j])
        gC = None
        for g in allq:
            if (zi*zj*z_inv(g)).is_square(): gC = g; break
        if gC is None:
            print('%-6s %-6s : нет g3' % (names[i], names[j])); continue
        try:
            val, nt, gam, mr = pair(allq[i], allq[j], gC)
            print('<%s,%s> = %+d   нетривиальные места: %s' % (names[i], names[j], val, nt))
        except Exception as ex:
            print('<%s,%s> ОШИБКА %r' % (names[i], names[j], ex))
print()
print('Ожидание: Sel^2(571a1) = Sha[2] = (Z/2)^2, спаривание НЕВЫРОЖДЕНО:')
print('  <s0,s1> = <s1,s0> = -1, диагональ и всё с классом 0 = +1.')
