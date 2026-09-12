# -*- coding: utf-8 -*-
"""НЕЗАВИСИМОЕ ОТКРЫТИЕ (не проверка) сертификата Касселса-Тейта для G1 (m,n)=(11,4).

Реализация с нуля: НЕ загружает descent/ctp_quartic.sage и не читает чужие JSON.
Использует только Sage/PARI (эллиптические кривые, коники, ell2cover).

Выход: candidates.json  —  I, J, phi-корни, тройки квартик (g1,g2,g3).
Всё остальное (z, m, H, gamma, места, символы Гильберта) считает и проверяет
отдельный скрипт на чистом Python.
"""
from sage.all import *
import json
from pathlib import Path

OUT = Path('/home/kep/magicKube/bridge/ctp_11_4_indep')

m, n = 11, 4
s = QQ(137) / 2
b = s * m ** 2 * n ** 2
assert b == 132616
# Порядок корней ЗАДАН постановкой:  e1=-b, e2=-s*m^4, e3=-s*n^4
E_orig = [-b, -s * m ** 4, -s * n ** 4]
print('original roots  e =', E_orig)

# ---- целая модель W: y^2 = f(x),  x = 4X (4 -- квадрат, класс delta сохраняется)
rts0 = [4 * e for e in E_orig]
assert all(r in ZZ for r in rts0), rts0
xx = polygen(QQ, 'x')
f0 = prod(xx - r for r in rts0)
W = EllipticCurve([0, f0[2], 0, f0[1], f0[0]])
print('W =', W)

# ---- минимальная модель, корни в согласованном порядке
M = W.minimal_model()
iso = W.isomorphism_to(M)
rts = [iso(W(r, 0))[0] for r in rts0]
print('M =', M, ' a-inv =', M.ainvs())
print('minimal roots (в исходном порядке) =', rts)
assert M.a1() == 0 and M.a3() == 0, 'нужна модель y^2=f(x)'

# контроль: масштаб между моделями -- квадрат, значит класс delta переносится
sc01 = (rts[0] - rts[1]) / (E_orig[0] - E_orig[1])
sc02 = (rts[0] - rts[2]) / (E_orig[0] - E_orig[2])
sc12 = (rts[1] - rts[2]) / (E_orig[1] - E_orig[2])
assert sc01 == sc02 == sc12, (sc01, sc02, sc12)
assert QQ(sc01).is_square(), sc01
print('масштаб X -> x_M :', sc01, ' квадрат:', QQ(sc01).is_square())

f = prod(xx - r for r in rts)
assert M.a2() == f[2] and M.a4() == f[1] and M.a6() == f[0]

I = M.c4()
J = 2 * M.c6()
b2 = M.b2()
phis = [-(12 * e + b2) for e in rts]
ph = polygen(QQ, 'ph')
assert all(p ** 3 - 3 * I * p + J == 0 for p in phis)
Delta = 16 * (4 * I ** 3 - J ** 2) / 27
print('I =', I)
print('J =', J)
print('phi =', phis)
print('Delta =', Delta, ' = 4096*disc(M)?', Delta == 4096 * M.discriminant())

# ------------------------------------------------------------------ квартики
def quartic_I(g):
    a, bb, c, d, e = g
    return 12 * a * e - 3 * bb * d + c * c

def quartic_J(g):
    a, bb, c, d, e = g
    return 72 * a * c * e - 27 * a * d * d - 27 * bb * bb * e + 9 * bb * c * d - 2 * c ** 3

def z_comps(g):
    a, bb, c, d, e = g
    return [(4 * a * p + 3 * bb * bb - 8 * a * c) / 3 for p in phis]

def sqfree(q):
    q = QQ(q)
    return QQ(q.numerator() * q.denominator()).squarefree_part()

def z_class(g):
    return [sqfree(v) for v in z_comps(g)]

fp = [prod(rts[i] - rts[j] for j in range(3) if j != i) for i in range(3)]   # f'(e_i)
# коэффициент при Theta^1 в базисе Лагранжа L_i
lam = [-(sum(rts[j] for j in range(3) if j != i)) / fp[i] for i in range(3)]
# коэффициент при Theta^0
lam0 = [prod(rts[j] for j in range(3) if j != i) / fp[i] for i in range(3)]

def normalize_invariants(g):
    """g -> mu*g с I(mu g)=I, J(mu g)=J."""
    Ig, Jg = quartic_I(g), quartic_J(g)
    assert Ig != 0 or Jg != 0
    mu = QQ(J * Ig) / QQ(Jg * I)
    g2 = [mu * a for a in g]
    assert quartic_I(g2) == I and quartic_J(g2) == J, (quartic_I(g2), quartic_J(g2))
    return g2, mu

def quartic_from_delta(delta, tag=''):
    """delta = (d1,d2,d3) в компонентах при корнях rts; N(delta) квадрат.
    Конструкция: коника sum_i d_i t_i^2 / f'(e_i) = 0, параметризация t(x,z),
    g(x,z) = -[delta t^2]_1 = -sum_i d_i t_i(x,z)^2 * lam_i.  Затем нормировка I,J."""
    d = [QQ(v) for v in delta]
    assert QQ(prod(d)).is_square(), 'норма не квадрат: %s' % (d,)
    co = [d[i] / fp[i] for i in range(3)]
    D = lcm([QQ(c).denominator() for c in co])
    co = [ZZ(c * D) for c in co]
    k = [ZZ(1)] * 3
    sf = []
    for i in range(3):
        ci = co[i]
        si = ci.squarefree_part()
        ki = ZZ(sqrt(QQ(ci) / si))
        assert si * ki ** 2 == ci
        sf.append(si)
        k[i] = ki
    C = Conic(QQ, [QQ(v) for v in sf])
    assert C.has_rational_point(), 'коника без точки: %s' % (sf,)
    par = C.parametrization()[0]
    R = par.domain().coordinate_ring()
    Tform = [par.defining_polynomials()[i] for i in range(3)]
    # контроль: параметризация лежит на конике
    assert sum(sf[i] * Tform[i] ** 2 for i in range(3)) == 0
    P2 = PolynomialRing(QQ, ['xq', 'zq'])
    xq, zq = P2.gens()
    hom = R.hom([xq, zq], P2)
    t = [hom(Tform[i]) / k[i] for i in range(3)]
    assert sum(co[i] * t[i] ** 2 for i in range(3)) == 0   # = D * [delta t^2]_2
    gpoly = -sum(d[i] * t[i] ** 2 * lam[i] for i in range(3))
    assert gpoly.degree() == 4, gpoly
    g = [QQ(gpoly.coefficient({xq: 4 - j, zq: j})) for j in range(5)]
    g, mu = normalize_invariants(g)
    if any(v == 0 for v in z_comps(g)) or g[0] == 0:
        g = make_z_unit(g)
    zc = z_class(g)
    dc = [sqfree(v) for v in d]
    assert zc == dc, ('класс z(g) != delta', tag, zc, dc)
    return g

def make_z_unit(g):
    """Фишер: можно считать z(g) единицей в L.  Собственно эквивалентная замена GL2."""
    mats = [(1, 0, 0, 1)]
    for k in [1, -1, 2, -2, 3, -3, 5, -5, 7, 11, 13]:
        mats += [(1, k, 0, 1), (1, 0, k, 1), (0, 1, 1, k), (k, 1, 1, 0)]
    for k in [1, -1, 2, -2, 3]:
        for l in [1, -1, 2, -2, 3]:
            mats.append((1, k, l, 1 + k * l + 1))
    mats.append((0, 1, 1, 0))
    for (A, B, Cc, D) in mats:
        if A * D - B * Cc == 0:
            continue
        gg = subst(g, A, B, Cc, D)
        if all(v != 0 for v in z_comps(gg)) and gg[0] != 0:
            assert quartic_I(gg) == I and quartic_J(gg) == J
            return gg
    raise RuntimeError('не удалось сделать z(g) единицей')

def subst(g, A, B, Cc, D):
    """g -> (det)^{-2} * g(A x + B z, C x + D z);  сохраняет I,J при любом GL2."""
    P2 = PolynomialRing(QQ, ['xq', 'zq'])
    xq, zq = P2.gens()
    G = sum(g[j] * xq ** (4 - j) * zq ** j for j in range(5))
    H = G(A * xq + B * zq, Cc * xq + D * zq)
    det = QQ(A * D - B * Cc)
    H = H / det ** 2
    return [QQ(H.coefficient({xq: 4 - j, zq: j})) for j in range(5)]

# ------------------------------------------------------------------ целевой класс
TRIP = [QQ(1), QQ(274), QQ(274)]    # (1, s, s) mod квадратов, порядок (e1,e2,e3)
assert QQ(prod(TRIP)).is_square()
print('\n=== строю g1 для delta =', TRIP, '===')
g1 = quartic_from_delta(TRIP, 'target')
g1 = make_z_unit(g1)
print('g1 =', g1)
print('z(g1) класс =', z_class(g1))
assert z_class(g1) == [sqfree(v) for v in TRIP]

# ------------------------------------------------------------------ кандидаты g2
cover = pari(M).ell2cover()
print('\nell2cover: %d квартик' % len(cover))
Rx = PolynomialRing(QQ, 'x')
cands = []
for i, c in enumerate(cover):
    q = Rx(c[0])
    raw = [QQ(q[4 - j]) for j in range(5)]
    try:
        g2, mu = normalize_invariants(raw)
        g2 = make_z_unit(g2)
    except Exception as ex:
        print('  #%d raw=%s  ПРОПУСК %r' % (i, raw, ex))
        continue
    print('  #%d raw=%s  mu=%s  z-class=%s' % (i, raw, mu, z_class(g2)))
    cands.append((i, str(mu), g2))

# ------------------------------------------------------------------ тройки
records = []
for (i, mus, g2) in cands:
    zz1 = z_comps(g1)
    zz2 = z_comps(g2)
    prod_class = [sqfree(zz1[j] * zz2[j]) for j in range(3)]
    if QQ(prod(prod_class)).is_square():
        try:
            g3 = quartic_from_delta(prod_class, 'sum%d' % i)
            g3 = make_z_unit(g3)
        except Exception as ex:
            print('  #%d g3 FAIL %r' % (i, ex))
            records.append({'i': i, 'mu': mus, 'g2': [str(v) for v in g2], 'g3_error': repr(ex)})
            continue
        records.append({'i': i, 'mu': mus,
                        'g2': [str(v) for v in g2],
                        'g3': [str(v) for v in g3],
                        'z1': [str(v) for v in zz1],
                        'z2': [str(v) for v in zz2],
                        'z3': [str(v) for v in z_comps(g3)]})
        print('  #%d готова тройка, g3 =' % i, g3)
    else:
        print('  #%d норма произведения не квадрат -- пропуск' % i)

# ------------------------------------------------------------------ контроли:
# квартики классов, которые ТОЧНО лежат в образе E(Q)/2E(Q) (2-кручение).
# Спаривание с ними обязано быть нулевым -- контроль на ложные срабатывания.
controls = []
tors = [
    ('T1', [(rts[0] - rts[1]) * (rts[0] - rts[2]), rts[0] - rts[1], rts[0] - rts[2]]),
    ('T2', [rts[1] - rts[0], (rts[1] - rts[0]) * (rts[1] - rts[2]), rts[1] - rts[2]]),
    ('T3', [rts[2] - rts[0], rts[2] - rts[1], (rts[2] - rts[0]) * (rts[2] - rts[1])]),
]
for name, dd in tors:
    dd = [sqfree(v) for v in dd]
    assert QQ(prod(dd)).is_square(), (name, dd)
    try:
        gt = quartic_from_delta(dd, name)
        # партнёр: g3 для суммы [gt] + [g2 из выбранного кандидата]
        controls.append({'name': name, 'delta': [str(v) for v in dd],
                         'g': [str(v) for v in gt], 'z': [str(v) for v in z_comps(gt)]})
        print('контроль %s: delta=%s  g=%s' % (name, dd, gt))
    except Exception as ex:
        print('контроль %s FAIL %r' % (name, ex))

# тривиальный класс
try:
    g0 = quartic_from_delta([QQ(1), QQ(1), QQ(1)], 'triv')
    controls.append({'name': 'trivial', 'delta': ['1', '1', '1'],
                     'g': [str(v) for v in g0], 'z': [str(v) for v in z_comps(g0)]})
    print('контроль trivial: g=%s' % (g0,))
except Exception as ex:
    print('контроль trivial FAIL %r' % (ex,))

# для каждого контрольного класса -- партнёрская тройка с g2 кандидатов
for ctl in controls:
    gt = [QQ(v) for v in ctl['g']]
    ctl['triples'] = []
    for (i, mus, g2) in cands:
        pc = [sqfree(z_comps(gt)[j] * z_comps(g2)[j]) for j in range(3)]
        if not QQ(prod(pc)).is_square():
            continue
        try:
            g3c = quartic_from_delta(pc, '%s+%d' % (ctl['name'], i))
            ctl['triples'].append({'i': i, 'g2': [str(v) for v in g2],
                                   'g3': [str(v) for v in g3c]})
        except Exception as ex:
            pass

rec = {
    'pair': '(m,n)=(11,4)',
    'note': 'НЕЗАВИСИМОЕ построение; ctp_quartic.sage не загружался',
    's': str(s), 'b': str(b),
    'original_roots': [str(v) for v in E_orig],
    'M_ainvs': [str(v) for v in M.ainvs()],
    'minimal_roots': [str(v) for v in rts],
    'model_scale_X_to_xM': str(sc01),
    'I': str(I), 'J': str(J), 'Delta': str(Delta),
    'phi_roots': [str(v) for v in phis],
    'target_trip': [str(v) for v in TRIP],
    'g1': [str(v) for v in g1],
    'z1_class': [str(v) for v in z_class(g1)],
    'candidates': records,
    'controls': controls,
}
(OUT / 'candidates.json').write_text(json.dumps(rec, indent=2) + '\n')
print('\nзаписано', OUT / 'candidates.json')
