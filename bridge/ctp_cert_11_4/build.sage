# -*- coding: utf-8 -*-
"""Сборка сертификата CTP для пары G1 (m,n)=(11,4), класс delta=(1,274,274).

Самостоятельная реализация (не загружает descent/ctp_quartic.sage).
Ранг E(Q) нигде не используется и не вычисляется.
Выход: certificate.json — читается чистым Python-верификатором.
"""
from sage.all import *
import json
from pathlib import Path

OUT = Path('/home/kep/magicKube/bridge/ctp_cert_11_4')
set_random_seed(114)

# ---------------------------------------------------------------- модель G1
m, n = 11, 4
s = QQ(137) / 2
b = s * m**2 * n**2
e = [-b, -s * m**4, -s * n**4]
sh = -sum(4 * ei for ei in e) / 3
assert sh in ZZ
r = [4 * ei + sh for ei in e]
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
DELTA_TARGET = [QQ(1), QQ(274), QQ(274)]

# ------------------------------------------------------------ инварианты квартик
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
    P = PolynomialRing(QQ, 'x,z'); x, z = P.gens()
    G = sum(g[i] * x**(4 - i) * z**i for i in range(5))
    p, q, u, v = mat
    G2 = G.subs({x: p * x + q * z, z: u * x + v * z}) / (p * v - q * u)**2
    return [QQ(G2.coefficient({x: 4 - i, z: i})) for i in range(5)]

def make_z_unit(g):
    if all(zinv(g, ph) != 0 for ph in phi):
        return g
    for _ in range(2000):
        mat = [ZZ.random_element(-4, 5) for _ in range(4)]
        if mat[0] * mat[3] - mat[1] * mat[2] == 0:
            continue
        g2 = act(g, mat)
        if all(zinv(g2, ph) != 0 for ph in phi):
            assert qI(g2) == qI(g) and qJ(g2) == qJ(g)
            return g2
    raise RuntimeError('не удалось сделать z(g) единицей')

def quartic_from_delta(d):
    assert sqcl(d[0] * d[1] * d[2]) == 1
    coef = [QQ(d[i]) / fp[i] for i in range(3)]
    sf, k = [], []
    for c in coef:
        t = sqcl(c); sf.append(t); k.append((c / t).sqrt())
        assert sf[-1] * k[-1]**2 == c
    C = Conic(QQ, [sf[0], sf[1], sf[2]])
    ok, pt = C.has_rational_point(point=True)
    assert ok
    par = C.parametrization()[0]
    Pxz = PolynomialRing(QQ, 'x,z'); xx, zz = Pxz.gens()
    T = [Pxz(t) for t in par.defining_polynomials()]
    assert sum(sf[i] * T[i]**2 for i in range(3)) == 0
    G = -sum(QQ(d[i]) * r[i] / fp[i] / k[i]**2 * T[i]**2 for i in range(3))
    g = [QQ(G.coefficient({xx: 4 - i, zz: i})) for i in range(5)]
    assert sum(g[i] * xx**(4 - i) * zz**i for i in range(5)) == G
    mu = QQ(J * qI(g) / (qJ(g) * I))
    assert mu.is_square()
    g = [mu * a for a in g]
    assert qI(g) == I and qJ(g) == J
    g = make_z_unit(g)
    assert zclass(g) == [sqcl(x) for x in d]
    return g

# ---------------------------------------------------------- локальная арифметика
def unit_val(a, p):
    a = QQ(a); u, v = ZZ(a.numerator()), ZZ(a.denominator()); k = 0
    assert u != 0
    while u % p == 0:
        u //= p; k += 1
    while v % p == 0:
        v //= p; k -= 1
    return k, QQ(u) / QQ(v)

def is_local_square(q, p):
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
    a, bq = QQ(a), QQ(bq)
    assert a != 0 and bq != 0
    if pl == 'real':
        return -1 if (a < 0 and bq < 0) else 1
    p = ZZ(pl)
    va, ua = unit_val(a, p); vb, ub = unit_val(bq, p)
    if p == 2:
        ra = int((ZZ(ua.numerator()) * inverse_mod(ZZ(ua.denominator()), 8)) % 8)
        rb = int((ZZ(ub.numerator()) * inverse_mod(ZZ(ub.denominator()), 8)) % 8)
        ex = ((ra - 1) // 2) * ((rb - 1) // 2) + va * ((rb * rb - 1) // 8) + vb * ((ra * ra - 1) // 8)
        return -1 if ex % 2 else 1
    ra = int((ZZ(ua.numerator()) * inverse_mod(ZZ(ua.denominator()), p)) % p)
    rb = int((ZZ(ub.numerator()) * inverse_mod(ZZ(ub.denominator()), p)) % p)
    la = power_mod(ra, int((p - 1) // 2), int(p)); lb = power_mod(rb, int((p - 1) // 2), int(p))
    res = -1 if (va * vb * ((p - 1) // 2)) % 2 else 1
    if la == p - 1 and vb % 2:
        res = -res
    if lb == p - 1 and va % 2:
        res = -res
    return res

def evform(co, x, z):
    d = len(co) - 1
    return sum(co[i] * x**(d - i) * z**i for i in range(len(co)))

def witnesses(g, gam, pl, want=1, bound=80):
    cands = [(QQ(1), QQ(0)), (QQ(0), QQ(1))]
    for zd in range(1, 7):
        for xn in range(-bound, bound + 1):
            if gcd(xn, zd) == 1:
                cands.append((QQ(xn), QQ(zd)))
    res = []
    for (x, z) in cands:
        if x == 0 and z == 0:
            continue
        q = evform(g, x, z)
        if q == 0:
            continue
        if gam is not None and evform(gam, x, z) == 0:
            continue
        if is_local_square(q, pl):
            res.append((x, z))
            if len(res) >= want:
                return res
    return res

# ------------------------------------------------ формула Фишера, Theorem 3.1
def Hform(g, ph):
    h = hess(g)
    G = [(4 * ph * g[i] + h[i]) / 3 for i in range(5)]
    H = [G[0], G[1] / 2, G[2] / 6 + QQ(2) / 9 * (I - ph * ph)]
    Hsq = [sum(H[i] * H[j] for i in range(3) for j in range(3) if i + j == kk) for kk in range(5)]
    assert Hsq == [G[0] * v for v in G]
    return H

def gamma_form(g1, g2, g3):
    z1 = [zinv(g1, ph) for ph in phi]
    z2 = [zinv(g2, ph) for ph in phi]
    z3 = [zinv(g3, ph) for ph in phi]
    assert all(x != 0 for x in z1 + z2 + z3)
    prod3 = [z1[i] * z2[i] * z3[i] for i in range(3)]
    assert all(p.is_square() for p in prod3)
    mm = [p.sqrt() for p in prod3]
    Hs = [Hform(g1, phi[i]) for i in range(3)]
    den = [prod(phi[i] - phi[j] for j in range(3) if j != i) for i in range(3)]
    gam = [sum((mm[i] / z1[i]) * Hs[i][j] / den[i] for i in range(3)) for j in range(3)]
    return gam, z1, z2, z3, mm

def places_for(g1, gam, a):
    S = set([2, 3, 5, 7, 11, 13])
    for q in [qdisc(g1), QQ(a), QQ(16) * (4 * I**3 - J**2) / 27]:
        S |= set(ZZ(p) for p, _ in factor(QQ(q)))
    for c in list(g1) + list(gam):
        if c != 0 and QQ(c).denominator() > 1:
            S |= set(ZZ(p) for p, _ in factor(QQ(c).denominator()))
    cont = gcd([QQ(c).numerator() for c in gam if c != 0])
    if cont != 0 and abs(cont) > 1:
        S |= set(ZZ(p) for p, _ in factor(abs(cont)))
    return sorted(S)

def pairing(g1, g2, g3, reps=2):
    gam, z1, z2, z3, mm = gamma_form(g1, g2, g3)
    a = g2[0]
    assert a != 0
    pls = places_for(g1, gam, a)
    runs = [[] for _ in range(reps)]
    tots = [1] * reps
    for pl in list(pls) + ['real']:
        ws = witnesses(g1, gam, pl, want=reps)
        assert len(ws) >= 1, ('нет локальной точки', pl)
        for j in range(reps):
            x, z = ws[min(j, len(ws) - 1)]
            gv = evform(gam, x, z)
            hs = hilbert(a, gv, pl)
            tots[j] *= hs
            runs[j].append({'place': str(pl), 'x': str(x), 'z': str(z),
                            'g': str(evform(g1, x, z)), 'gamma': str(gv), 'hilbert': int(hs)})
    assert len(set(tots)) == 1, ('разные значения при разных свидетелях', tots)
    return tots[0], runs, gam, mm, [z1, z2, z3], pls

def els_witnesses(g):
    S = set([2, 3])
    S |= set(ZZ(p) for p, _ in factor(qdisc(g)))
    for c in g:
        if c != 0 and QQ(c).denominator() > 1:
            S |= set(ZZ(p) for p, _ in factor(QQ(c).denominator()))
    cont = gcd([QQ(c).numerator() for c in g if c != 0])
    if cont != 0 and abs(cont) > 1:
        S |= set(ZZ(p) for p, _ in factor(abs(cont)))
    rows = []
    for pl in sorted(S) + ['real']:
        w = witnesses(g, None, pl, want=1, bound=200)
        assert w, ('квартика не разрешима локально?', pl)
        x, z = w[0]
        rows.append({'place': str(pl), 'x': str(x), 'z': str(z), 'value': str(evform(g, x, z))})
    return rows

# --------------------------------------------------------------------- прогон
print('M =', M, flush=True)
g1 = quartic_from_delta(DELTA_TARGET)
print('g1 =', g1, flush=True)
D2 = [QQ(30), QQ(137), QQ(4110)]                   # выбранный партнёр из Sel^2
D3 = [QQ(sqcl(DELTA_TARGET[i] * D2[i])) for i in range(3)]
print('delta2 =', D2, ' delta3 =', D3, flush=True)
g2 = quartic_from_delta(D2)
g3 = quartic_from_delta(D3)
print('g2 =', g2, flush=True)
print('g3 =', g3, flush=True)

val, runs, gam, mm, zs, pls = pairing(g1, g2, g3, reps=2)
print('ПРЯМОЕ  <g1,g2> =', val, ' места:', [str(p) for p in pls] + ['real'], flush=True)
print('  нетривиальные:', [x['place'] for x in runs[0] if x['hilbert'] == -1], flush=True)

rval, rruns, rgam, rmm, rzs, rpls = pairing(g2, g1, g3, reps=2)
print('ОБРАТНОЕ <g2,g1> =', rval, ' нетривиальные:',
      [x['place'] for x in rruns[0] if x['hilbert'] == -1], flush=True)

# третья проверка: <g1,g3> должно совпадать с <g1,g2> (Fisher 3.2(v))
tval, truns, tgam, tmm, tzs, tpls = pairing(g1, g3, g2, reps=1)
print('КОНТРОЛЬ <g1,g3> =', tval, flush=True)

cert = {
    'pair': '(m,n)=(11,4)', 's': str(s), 'b': str(b),
    'orig_roots': [str(x) for x in e],
    'minimal_roots': [str(x) for x in r],
    'A': str(A), 'B': str(B), 'I': str(I), 'J': str(J),
    'phi_roots': [str(x) for x in phi],
    'delta_target': [str(x) for x in DELTA_TARGET],
    'delta2': [str(x) for x in D2], 'delta3': [str(x) for x in D3],
    'quartics': [[str(x) for x in g] for g in (g1, g2, g3)],
    'quartic_disc': [str(qdisc(g)) for g in (g1, g2, g3)],
    'z_values': [[str(zinv(g, ph)) for ph in phi] for g in (g1, g2, g3)],
    'm_values': [str(x) for x in mm],
    'gamma': [str(x) for x in gam],
    'a': str(g2[0]),
    'places': [str(p) for p in pls] + ['real'],
    'pair_runs': runs,
    'pair_value': int(val),
    'reverse': {'m_values': [str(x) for x in rmm], 'gamma': [str(x) for x in rgam],
                'a': str(g1[0]), 'places': [str(p) for p in rpls] + ['real'],
                'runs': rruns, 'value': int(rval)},
    'control_g1g3': {'value': int(tval), 'gamma': [str(x) for x in tgam],
                     'a': str(g3[0]), 'm_values': [str(x) for x in tmm],
                     'places': [str(p) for p in tpls] + ['real'], 'runs': truns},
    'local_solubility': [els_witnesses(g) for g in (g1, g2, g3)],
}
(OUT / 'certificate.json').write_text(json.dumps(cert, indent=1) + '\n')
print('certificate.json записан', flush=True)
print('DONE', flush=True)
