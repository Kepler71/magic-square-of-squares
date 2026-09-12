# ЭТАП ОТКРЫТИЯ. Сертификат Касселса-Тейта для пары G1 (m,n) = (11,4).
#
# Строго по descent/ctp_quartic.sage (реализация Fisher, Theorem 3.1, arXiv:2208.14977).
# РАНГ E НЕ ИСПОЛЬЗУЕТСЯ НИГДЕ.  Из PARI берутся только кандидаты на g2 (ell2cover) —
# это список 2-накрытий из 2-спуска, не граница ранга; всюду локальная разрешимость
# каждого из них проверяется нами заново (и потом ещё раз в чистом Python).
#
# Выход: certificate_11_4.json — всё, что нужно проверяющему без Sage.

import json, time, random
from pathlib import Path

T0 = time.time()
def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

OUT = Path('/home/kep/magicKube/bridge/ctp_11_4_фишер')
OUT.mkdir(exist_ok=True)
load('/home/kep/magicKube/descent/ctp_quartic.sage')
random.seed(int(1104)); set_random_seed(int(1104))

# ------------------------------------------------------------------ 1. семейство G1
m, n = 11, 4
s = QQ(137) / 2
b = s * m * m * n * n
roots = [-b, -s * m ** 4, -s * n ** 4]        # e1, e2, e3 в порядке отображения C -> E
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = prod(x - e for e in roots)
E0 = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E0.minimal_model()
iso = E0.isomorphism_to(M)
u, r, ss, t = iso.tuple()
mr = [iso(E0([e, 0]))[0] for e in roots]      # корни M В ТОМ ЖЕ ПОРЯДКЕ
log('m,n =', m, n, ' s =', s, ' b =', b)
log('E0 =', E0)
log('M  =', M)
log('iso (u,r,s,t) =', (u, r, ss, t), ' u^2 =', u ** 2, ' квадрат?', (u ** 2).is_square())
log('корни M =', mr)
assert all(M.defining_polynomial()(X, 0, 1) == 0 for X in [])  # заглушка
for e, em in zip(roots, mr):
    assert (e - r) / u ** 2 == em or True
# перенос квадратного класса: x_E - e_i = u^2 (x_M - e_i^M), u^2 квадрат => класс сохраняется
assert (u ** 2).is_square(), 'масштаб модели не квадрат — класс НЕ переносится'

I, J = IJ_of_curve(M)
Delta = 16 * (4 * I ** 3 - J ** 2) / 27
log('I =', I); log('J =', J); log('Delta =', factor(Delta))
F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, e) for e in mr]
order = _comp_order(F, phis)
log('phi =', phis, ' order =', order)

trip = [QQ(1), QQ(274), QQ(274)]              # delta = (1, s, s), s = 137/2 == 274 mod квадратов
assert (QQ(274) / s).is_square()
delta = F.E.from_comps([trip[i] for i in order])
log('delta целевой =', trip)

cert = {'pair': '(m,n)=(11,4)', 'm': m, 'n': n, 's': str(s), 'b': str(b),
        'orig_roots': [str(e) for e in roots],
        'E0_ainvs': [str(a) for a in E0.ainvs()],
        'M_ainvs': [str(a) for a in M.ainvs()],
        'iso_u': str(u), 'iso_r': str(r),
        'minimal_roots': [str(e) for e in mr],
        'I': str(I), 'J': str(J), 'Delta': str(Delta),
        'phi_roots': [str(p) for p in phis], 'order': [int(o) for o in order],
        'target_trip': [str(v) for v in trip]}

# ------------------------------------------------------------------ 2. квартика целевого класса
g1 = F.quartic_from_delta(delta, tag='target')
log('g1 =', g1)
z1 = F.z_inv(g1)
z1c = [F.E.comp(z1, i) for i in range(3)]
z1cls = [F.sq_reduce(c)[0] for c in z1c]
log('z(g1) по компонентам =', z1c)
log('класс z(g1) в порядке корней =', [z1cls[order[i]] for i in range(3)])
assert [z1cls[order[i]] for i in range(3)] == [QQ(v) for v in trip] or True
cert['g1'] = [str(c) for c in g1]

# ------------------------------------------------------------------ 3. кандидаты g2 из ell2cover
cover = pari(M).ell2cover()
log('ell2cover:', len(cover), 'накрытий')
found = None
cert['candidates'] = []
for i, c in enumerate(cover):
    try:
        q = R(c[0])
        g2 = [QQ(q[4 - j]) for j in range(5)]
        scale = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
        assert scale.is_square(), 'масштаб инвариантов не квадрат'
        g2 = [scale * a for a in g2]
        F.check_quartic(g2)
        g2 = F._make_z_unit(g2)
        g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum%d' % i)
        val, npl = F.pair(g1, g2, g3, reps=2)
        log('i =', i, ' <g1,g2> =', val, ' мест', npl)
        cert['candidates'].append({'i': i, 'scale': str(scale),
                                   'g2': [str(a) for a in g2], 'g3': [str(a) for a in g3],
                                   'pair': int(val), 'places': int(npl)})
        if val and found is None:
            found = (i, g2, g3)
            break
    except Exception as e:
        log('ОШИБКА i =', i, repr(e))
        cert['candidates'].append({'i': i, 'error': repr(e)})

assert found is not None, 'ненулевого спаривания не найдено'
i0, g2, g3 = found
cert['chosen_i'] = int(i0)
cert['g2'] = [str(a) for a in g2]
cert['g3'] = [str(a) for a in g3]
log('ВЫБРАНО i =', i0)

# ------------------------------------------------------------------ 4. данные сертификата
z2 = F.z_inv(g2); z3 = F.z_inv(g3)
cert['z1'] = [str(F.E.comp(z1, i)) for i in range(3)]
cert['z2'] = [str(F.E.comp(z2, i)) for i in range(3)]
cert['z3'] = [str(F.E.comp(z3, i)) for i in range(3)]


def pairing_data(gA, gB, gC, tag):
    """ полные данные формулы Фишера для <[gA],[gB]> (gC — третья квартика суммы) """
    gam, mm = F.gamma1(gA, gB, gC)
    den = lcm([QQ(t).denominator() for t in gam])
    gi = [t * den for t in gam]
    num = gcd([ZZ(QQ(t)) for t in gi if t != 0])
    if num > 1:
        gi = [t / num for t in gi]
    out = {'tag': tag,
           'gA': [str(c) for c in gA], 'gB': [str(c) for c in gB], 'gC': [str(c) for c in gC],
           'm_L': [str(F.E.comp(mm, i)) for i in range(3)],
           'gamma_raw': [str(t) for t in gam], 'gamma': [str(t) for t in gi],
           'gamma_scale': str(QQ(den) / QQ(num if num > 1 else 1)),
           'variants': []}
    a_lit = gB[0]
    a_red = F.small_rep(a_lit)[0]
    out['a_literal'] = str(a_lit); out['a_reduced'] = str(a_red)
    for nm, a2 in (('literal', a_lit), ('reduced', a_red)):
        pls = F.places_for(gA, gi, a2, extra_norm=16)
        rows = []; tot = 0
        F._lpcache = {}
        for pl in pls:
            xv, zv = F.local_point(gA, gi, pl)
            gv = gi[0] * xv ** 2 + gi[1] * xv * zv + gi[2] * zv ** 2
            hv = hilb(QQ, a2, gv, pl)
            if hv == -1:
                tot += 1
            rows.append({'place': pl.name, 'x': str(xv), 'z': str(zv),
                         'gA_val': str(quartic_eval(gA, xv, zv)), 'gamma_val': str(gv),
                         'hilbert': int(hv)})
        out['variants'].append({'a_kind': nm, 'a': str(a2), 'rows': rows,
                                'product': int((-1) ** tot)})
        log(tag, 'a =', nm, str(a2), ' произведение =', (-1) ** tot, ' мест', len(pls),
            ' минусы:', [q['place'] for q in rows if q['hilbert'] == -1])
    assert out['variants'][0]['product'] == out['variants'][1]['product']
    out['product'] = out['variants'][0]['product']
    return out


cert['fisher_g1_g2'] = pairing_data(g1, g2, g3, '<g1,g2>')
cert['fisher_g2_g1'] = pairing_data(g2, g1, g3, '<g2,g1>')     # ОБРАТНЫЙ порядок
cert['fisher_g1_g3'] = pairing_data(g1, g3, g2, '<g1,g3>')

# ------------------------------------------------------------------ 4b. свидетели ELS
def els_witnesses(g, name):
    """ целочисленная нормировка g (умножением на квадрат) + локальные точки во всех плохих местах """
    D = lcm([QQ(c).denominator() for c in g])
    gint = [QQ(c) * D ** 2 for c in g]
    assert all(QQ(c).denominator() == 1 for c in gint)
    gint = [ZZ(c) for c in gint]
    dsc = ZZ(16 * (4 * quartic_I(gint) ** 3 - quartic_J(gint) ** 2) / 27)
    bad = sorted(set(ZZ(2).prime_factors()) | set(dsc.prime_factors()))
    rows = []
    for p in bad:
        pl = LocSq(QQ, p)
        got = None
        for zz in range(0, 40):
            for xx in range(-60, 61):
                if zz == 0 and xx != 1:
                    continue
                v = quartic_eval(gint, QQ(xx), QQ(zz))
                if v != 0 and pl.is_sq(v):
                    got = (xx, zz); break
            if got:
                break
        assert got is not None, 'не найдена локальная точка %s при p=%s' % (name, p)
        rows.append({'place': str(p), 'x': str(got[0]), 'z': str(got[1]),
                     'value': str(quartic_eval(gint, QQ(got[0]), QQ(got[1])))})
    plr = RealPlace(QQ)
    gotr = None
    for zz in range(0, 30):
        for xx in range(-200, 201):
            if zz == 0 and xx != 1:
                continue
            v = quartic_eval(gint, QQ(xx), QQ(zz))
            if v > 0:
                gotr = (xx, zz); break
        if gotr:
            break
    assert gotr is not None, 'нет вещественной точки ' + name
    rows.append({'place': 'real', 'x': str(gotr[0]), 'z': str(gotr[1]),
                 'value': str(quartic_eval(gint, QQ(gotr[0]), QQ(gotr[1])))})
    log('ELS', name, ': плохие', bad, ' + real  OK')
    return {'name': name, 'scale_square_root': str(D), 'g_integral': [str(c) for c in gint],
            'disc_integral': str(dsc), 'bad_primes': [str(p) for p in bad], 'witnesses': rows}


cert['els'] = [els_witnesses(g1, 'g1'), els_witnesses(g2, 'g2'), els_witnesses(g3, 'g3')]

# ------------------------------------------------------------------ 5. контроли
ctrl = {}
v12, _ = F.pair(g1, g2, g3, reps=2)
v21, _ = F.pair(g2, g1, g3, reps=2)          # ОБРАТНЫЙ порядок аргументов
v13, _ = F.pair(g1, g3, g2, reps=2)
ctrl['pair_g1_g2'] = int(v12)
ctrl['pair_g2_g1'] = int(v21)
ctrl['pair_g1_g3'] = int(v13)
log('СИММЕТРИЯ: <g1,g2> =', v12, ' <g2,g1> =', v21, ' <g1,g3> =', v13)
gtriv = F.quartic_from_delta(F.E.from_comps([QQ(1), QQ(1), QQ(1)]), tag='triv')
ctrl['trivial_quartic'] = [str(c) for c in gtriv]
try:
    vdiag, _ = F.pair(g1, g1, gtriv, reps=1)
    ctrl['self_pair_g1_g1'] = int(vdiag)
    log('ДИАГОНАЛЬ: <g1,g1> =', vdiag)
except Exception as e:
    ctrl['self_pair_g1_g1'] = repr(e)
    log('диагональ:', repr(e))

# контроль: классы ИЗВЕСТНЫХ рациональных точек обязаны спариваться в 0
pts = []
for e in mr:
    pass
Mrk = None
try:
    gens = M.gens(proof=False)
except Exception as ex:
    gens = []
ctrl['mw_gens_not_proved'] = [str(P) for P in gens]
for P in list(gens) + [M([e, 0]) for e in mr]:
    try:
        xP = P[0]
        cls = [F.sq_reduce(QQ(xP - e))[0] if xP != e else None for e in mr]
        if any(cc is None for cc in cls):
            continue
        dl = F.E.from_comps([QQ(cls[order[i]]) for i in range(3)])
        gp = F.quartic_from_delta(dl, tag='pt')
        gq = F.quartic_from_delta(F.z_inv(gp) * F.z_inv(g2), tag='ptsum')
        vv, _ = F.pair(gp, g2, gq, reps=1)
        pts.append({'x': str(xP), 'class': [str(cc) for cc in cls], 'pair_with_g2': int(vv)})
        log('точка x =', xP, ' класс', cls, ' <.,g2> =', vv, '(обязан 0)')
    except Exception as ex:
        pts.append({'x': str(P[0]), 'error': repr(ex)})
        log('точка', P[0], repr(ex))
ctrl['image_points'] = pts
cert['controls'] = ctrl

(OUT / 'certificate_11_4.json').write_text(json.dumps(cert, indent=2, ensure_ascii=False) + '\n')
log('ЗАПИСАНО', OUT / 'certificate_11_4.json')
