# -*- coding: utf-8 -*-
"""Независимое ОТКРЫТИЕ сертификата CTP для пары G1 (m,n)=(11,4).

Своя реализация: этап 1 (конструкция квартик из класса delta через одну конику),
этап 2 (формула Фишера, Theorem 3.1) — всё явно, без загрузки ctp_quartic.sage.
Ранг E нигде не используется.
"""
from sage.all import *
import json, itertools
from pathlib import Path

OUT = Path('/home/kep/magicKube/bridge/ctp_cert_11_4')
set_random_seed(114)

# ---------------------------------------------------------------- модель G1
m, n = 11, 4
s = QQ(137) / 2
b = s * m**2 * n**2
e = [-b, -s * m**4, -s * n**4]                 # e1,e2,e3 в порядке отображения C -> E
sh = -sum(4 * ei for ei in e) / 3
assert sh in ZZ
r = [4 * ei + sh for ei in e]                  # корни целочисленной депрессированной модели
assert sum(r) == 0 and len(set(r)) == 3
Rx = PolynomialRing(QQ, 'X'); X = Rx.gen()
f = prod(X - ri for ri in r)
A, B = f[1], f[0]
M = EllipticCurve([0, 0, 0, A, B])
I, J = -48 * A, -1728 * B
assert I == M.c4() and J == 2 * M.c6()
phi = [-12 * ri for ri in r]
assert all(ph**3 - 3 * I * ph + J == 0 for ph in phi)
fp = [prod(r[i] - r[j] for j in range(3) if j != i) for i in range(3)]

DELTA_TARGET = [QQ(1), QQ(274), QQ(274)]       # (1, s, s), s=137/2 ~ 274 mod квадратов

# ------------------------------------------------------- квартики: инварианты
def qI(g):
    a, bb, c, d, ee = g
    return 12 * a * ee - 3 * bb * d + c * c

def qJ(g):
    a, bb, c, d, ee = g
    return 72 * a * c * ee - 27 * a * d * d - 27 * bb * bb * ee + 9 * bb * c * d - 2 * c**3

def qdisc(g):
    a, bb, c, d, ee = g
    return (256 * a**3 * ee**3 - 192 * a**2 * bb * d * ee**2 - 128 * a**2 * c**2 * ee**2
            + 144 * a**2 * c * d**2 * ee - 27 * a**2 * d**4 + 144 * a * bb**2 * c * ee**2
            - 6 * a * bb**2 * d**2 * ee - 80 * a * bb * c**2 * d * ee + 18 * a * bb * c * d**3
            + 16 * a * c**4 * ee - 4 * a * c**3 * d**2 - 27 * bb**4 * ee**2 + 18 * bb**3 * c * d * ee
            - 4 * bb**3 * d**3 - 4 * bb**2 * c**3 * ee + bb**2 * c**2 * d**2)

def zinv(g, ph):
    a, bb, c, d, ee = g
    return (4 * a * ph + 3 * bb * bb - 8 * a * c) / 3

def hess(g):
    a, bb, c, d, ee = g
    return [3 * bb * bb - 8 * a * c, 4 * (bb * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * ee - 3 * bb * d), 4 * (c * d - 6 * bb * ee), 3 * d * d - 8 * c * ee]

def sqcl(q):
    q = QQ(q)
    if q == 0:
        return 0
    num = q.numerator() * q.denominator()
    return (1 if q > 0 else -1) * prod([p for p, k in factor(abs(num)) if k % 2 == 1])

def zclass(g):
    return [sqcl(zinv(g, ph)) for ph in phi]

def act(g, mat):
    """g -> (det)^(-2) * g(p x + q z, u x + v z): собственная эквивалентность, I,J сохраняются."""
    P = PolynomialRing(QQ, 'x,z'); x, z = P.gens()
    G = sum(g[i] * x**(4 - i) * z**i for i in range(5))
    p, q, u, v = mat
    G2 = G.subs({x: p * x + q * z, z: u * x + v * z}) / (p * v - q * u)**2
    return [QQ(G2.coefficient({x: 4 - i, z: i})) for i in range(5)]

def make_z_unit(g):
    if all(zinv(g, ph) != 0 for ph in phi):
        return g, [1, 0, 0, 1]
    for trial in range(500):
        mat = [ZZ.random_element(-4, 5) for _ in range(4)]
        if mat[0] * mat[3] - mat[1] * mat[2] == 0:
            continue
        g2 = act(g, mat)
        if all(zinv(g2, ph) != 0 for ph in phi):
            assert qI(g2) == qI(g) and qJ(g2) == qJ(g)
            return g2, mat
    raise RuntimeError('не удалось сделать z(g) единицей')

# ----------------------------------------- конструкция квартики из класса delta
def quartic_from_delta(d, verbose=True):
    """d = (d1,d2,d3) в Q*, произведение — квадрат. Одна коника над Q, затем нормировка I,J."""
    assert sqcl(d[0] * d[1] * d[2]) == 1, 'норма delta не квадрат'
    coef = [QQ(d[i]) / fp[i] for i in range(3)]
    # убрать квадратные множители: coef_i = sf_i * k_i^2, замена T_i = k_i t_i
    sf, k = [], []
    for c in coef:
        t = sqcl(c)
        sf.append(t)
        k.append((c / t).sqrt())
        assert sf[-1] * k[-1]**2 == c
    C = Conic(QQ, [sf[0], sf[1], sf[2]])
    ok, pt = C.has_rational_point(point=True)
    assert ok, 'коника без рациональной точки'
    par = C.parametrization()[0]
    Tpoly = par.defining_polynomials()
    Pxz = PolynomialRing(QQ, 'x,z'); xx, zz = Pxz.gens()
    T = [Pxz(t) for t in Tpoly]
    assert sum(sf[i] * T[i]**2 for i in range(3)) == 0
    # g = -sum d_i r_i / f'(r_i) * t_i^2,  t_i = T_i / k_i
    G = -sum(QQ(d[i]) * r[i] / fp[i] / k[i]**2 * T[i]**2 for i in range(3))
    g = [QQ(G.coefficient({xx: 4 - i, zz: i})) for i in range(5)]
    assert sum(g[i] * xx**(4 - i) * zz**i for i in range(5)) == G, 'G не квартика'
    Ig, Jg = qI(g), qJ(g)
    assert Ig != 0 or Jg != 0
    mu = QQ(J * Ig / (Jg * I))
    assert mu.is_square(), 'нормировочный множитель не квадрат'
    g = [mu * a for a in g]
    assert qI(g) == I and qJ(g) == J
    g, mat = make_z_unit(g)
    cl = zclass(g)
    assert cl == [sqcl(x) for x in d], ('класс z(g) не совпал с delta', cl, [sqcl(x) for x in d])
    if verbose:
        print('   delta', [str(x) for x in d], '-> g', [str(x) for x in g], 'z-class', cl, flush=True)
    return g

# ------------------------------------------------ формула Фишера, Theorem 3.1
def Hform(g, ph):
    """H = G0 x^2 + (G1/2) xz + (G2/6 + 2/9 (I - phi^2)) z^2, где G = (4 phi g + h)/3."""
    h = hess(g)
    G = [(4 * ph * g[i] + h[i]) / 3 for i in range(5)]
    H = [G[0], G[1] / 2, G[2] / 6 + QQ(2) / 9 * (I - ph * ph)]
    # тождество G(1,0) * G = H^2
    Hsq = [sum(H[i] * H[j] for i in range(3) for j in range(3) if i + j == kk) for kk in range(5)]
    assert Hsq == [G[0] * v for v in G], 'тождество G(1,0)G=H^2 нарушено'
    return H

def gamma_form(g1, g2, g3):
    z1 = [zinv(g1, ph) for ph in phi]
    z2 = [zinv(g2, ph) for ph in phi]
    z3 = [zinv(g3, ph) for ph in phi]
    assert all(x != 0 for x in z1 + z2 + z3)
    prod3 = [z1[i] * z2[i] * z3[i] for i in range(3)]
    assert all(p.is_square() for p in prod3), 'z1z2z3 не квадрат в L'
    mm = [p.sqrt() for p in prod3]
    Hs = [Hform(g1, phi[i]) for i in range(3)]
    # gamma = коэффициент при phi^2 в (m/z1) * H1: интерполяция Лагранжа по трём компонентам
    den = [prod(phi[i] - phi[j] for j in range(3) if j != i) for i in range(3)]
    gam = [sum((mm[i] / z1[i]) * Hs[i][j] / den[i] for i in range(3)) for j in range(3)]
    return gam, z1, z2, z3, mm

def evform(co, x, z):
    d = len(co) - 1
    return sum(co[i] * x**(d - i) * z**i for i in range(len(co)))

def places_for(g1, gam, a):
    S = set([2, 3, 5, 7, 11, 13])
    for q in [qdisc(g1), a]:
        S |= set(p for p, _ in factor(QQ(q)))
    for c in list(g1) + list(gam):
        if c != 0:
            S |= set(p for p, _ in factor(QQ(c).denominator())) if QQ(c).denominator() > 1 else set()
    cont = gcd([QQ(c).numerator() for c in gam if c != 0])
    if cont != 0 and abs(cont) > 1:
        S |= set(p for p, _ in factor(abs(cont)))
    return sorted(S)

def unit_val(a, p):
    a = QQ(a); u, v = ZZ(a.numerator()), ZZ(a.denominator()); k = 0
    assert u != 0
    while u % p == 0:
        u //= p; k += 1
    while v % p == 0:
        v //= p; k -= 1
    return k, QQ(u) / QQ(v)

def is_local_square(q, p):
    """Собственная реализация, без Sage-примитивов."""
    q = QQ(q)
    if q == 0:
        return True
    if p == 'real':
        return q > 0
    p = ZZ(p)
    k, u = unit_val(q, p)
    if k % 2 != 0:
        return False
    num, den = ZZ(u.numerator()), ZZ(u.denominator())
    if p == 2:
        return (num * inverse_mod(den, 8)) % 8 == 1
    return power_mod(int((num * inverse_mod(den, p)) % p), int((p - 1) // 2), int(p)) == 1

def hilbert(a, bq, pl):
    """Символ Гильберта (a,b)_v, своя реализация."""
    a, bq = QQ(a), QQ(bq)
    assert a != 0 and bq != 0
    if pl == 'real':
        return -1 if (a < 0 and bq < 0) else 1
    p = ZZ(pl)
    va, ua = unit_val(a, p)
    vb, ub = unit_val(bq, p)
    if p == 2:
        ra = int((ZZ(ua.numerator()) * inverse_mod(ZZ(ua.denominator()), 8)) % 8)
        rb = int((ZZ(ub.numerator()) * inverse_mod(ZZ(ub.denominator()), 8)) % 8)
        ex = ((ra - 1) // 2) * ((rb - 1) // 2) + va * ((rb * rb - 1) // 8) + vb * ((ra * ra - 1) // 8)
        return -1 if ex % 2 else 1
    ra = int((ZZ(ua.numerator()) * inverse_mod(ZZ(ua.denominator()), p)) % p)
    rb = int((ZZ(ub.numerator()) * inverse_mod(ZZ(ub.denominator()), p)) % p)
    la = power_mod(ra, int((p - 1) // 2), int(p))
    lb = power_mod(rb, int((p - 1) // 2), int(p))
    res = -1 if (va * vb * ((p - 1) // 2)) % 2 else 1
    if la == p - 1 and vb % 2:
        res = -res
    if lb == p - 1 and va % 2:
        res = -res
    return res

def local_points(g, gam, pl, need_gamma=True, bound=60):
    """(x:z) c g(x,z) локальный квадрат и gamma(x,z) != 0."""
    cands = [(QQ(1), QQ(0)), (QQ(0), QQ(1))]
    for zd in range(1, 6):
        for xn in range(-bound, bound + 1):
            if gcd(xn, zd) == 1:
                cands.append((QQ(xn), QQ(zd)))
    for (x, z) in cands:
        if x == 0 and z == 0:
            continue
        q = evform(g, x, z)
        if q == 0:
            continue
        if need_gamma and evform(gam, x, z) == 0:
            continue
        if is_local_square(q, pl):
            return (x, z)
    return None

def pairing(g1, g2, g3, tag=''):
    gam, z1, z2, z3, mm = gamma_form(g1, g2, g3)
    a = g2[0]
    assert a != 0
    pls = places_for(g1, gam, a)
    rows = []
    tot = 1
    for pl in list(pls) + ['real']:
        pt = local_points(g1, gam, pl)
        assert pt is not None, ('нет локальной точки', pl)
        x, z = pt
        gv = evform(gam, x, z)
        hs = hilbert(a, gv, pl)
        tot *= hs
        rows.append({'place': str(pl), 'x': str(x), 'z': str(z),
                     'g': str(evform(g1, x, z)), 'gamma': str(gv), 'hilbert': int(hs)})
    return tot, rows, gam, mm, [z1, z2, z3]

# --------------------------------------------------------------------- прогон
print('M =', M, flush=True)
print('roots =', r, ' I,J =', I, J, flush=True)
print('phi =', phi, flush=True)

rec = {'m': int(m), 'n': int(n), 's': str(s), 'b': str(b),
       'orig_roots': [str(x) for x in e], 'minimal_roots': [str(x) for x in r],
       'A': str(A), 'B': str(B), 'I': str(I), 'J': str(J),
       'phi_roots': [str(x) for x in phi],
       'delta_target': [str(x) for x in DELTA_TARGET], 'runs': []}

print('строю g1 из целевого класса', flush=True)
g1 = quartic_from_delta(DELTA_TARGET)
rec['g1'] = [str(x) for x in g1]
print('g1 =', g1, flush=True)

# кандидаты delta' — классы 2-накрытий PARI (только как источник кандидатов)
cands = []
cov = pari(M).ell2cover()
Rq = PolynomialRing(QQ, 'x'); xq = Rq.gen()
for i, c in enumerate(cov):
    q = Rq(c[0])
    g = [QQ(q[4 - j]) for j in range(5)]
    mu = QQ(J * qI(g) / (qJ(g) * I))
    assert mu.is_square()
    g = [mu * aa for aa in g]
    assert qI(g) == I and qJ(g) == J
    g, _ = make_z_unit(g)
    cands.append(zclass(g))
print('кандидаты из ell2cover:', cands, flush=True)
rec['candidates'] = [[str(x) for x in cc] for cc in cands]

for idx, cl in enumerate(cands):
    try:
        d2 = [QQ(x) for x in cl]
        print('--- кандидат', idx, cl, flush=True)
        g2 = quartic_from_delta(d2)
        d3 = [sqcl(DELTA_TARGET[i] * d2[i]) for i in range(3)]
        print('   delta3 =', d3, flush=True)
        g3 = quartic_from_delta([QQ(x) for x in d3])
        val, rows, gam, mm, zs = pairing(g1, g2, g3)
        print('   PAIRING =', val, ' мест:', len(rows),
              ' нетривиальные:', [x['place'] for x in rows if x['hilbert'] == -1], flush=True)
        rec['runs'].append({'idx': idx, 'delta2': [str(x) for x in d2], 'delta3': [str(x) for x in d3],
                            'g2': [str(x) for x in g2], 'g3': [str(x) for x in g3],
                            'gamma': [str(x) for x in gam], 'm': [str(x) for x in mm],
                            'z': [[str(x) for x in zz] for zz in zs],
                            'pairing': int(val), 'rows': rows})
        (OUT / 'discover.json').write_text(json.dumps(rec, indent=1) + '\n')
        if val == -1:
            print('НАЙДЕНО ненулевое спаривание на кандидате', idx, flush=True)
    except Exception as ex:
        print('   ОШИБКА', repr(ex), flush=True)
        rec['runs'].append({'idx': idx, 'error': repr(ex)})
        (OUT / 'discover.json').write_text(json.dumps(rec, indent=1) + '\n')

(OUT / 'discover.json').write_text(json.dumps(rec, indent=1) + '\n')
print('DONE', flush=True)
