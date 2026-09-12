"""Этап 3: экспорт сертификата CTP для G1 (19,5) в certificate_19_5.json.
Всё, что попадает в сертификат, пересчитывается здесь заново и потом ещё раз
независимо проверяется скриптом ctp_cert_19_5_ссzero.py на стандартном Python.
"""
from sage.all import *
import json, time, random
from pathlib import Path

src = open('/home/kep/magicKube/bridge/ctp_19_5_claude/pairing_19_5.py').read()
exec(src.split('# ------------------------------------------------ кандидаты g2 из ell2cover')[0])

OUT = Path('/home/kep/magicKube/bridge/ctp_19_5_claude')
random.seed(777); set_random_seed(777)

# ---- квартики: g1 строим сами, g2 берём из ell2cover (кандидат), g3 строим сами
Rp = PolynomialRing(QQ, 'x'); xp = Rp.gen()
cover = pari(M).ell2cover()
q = Rp(cover[3][0])
g2 = [QQ(q[4-j]) for j in range(5)]
lam2 = QQ(quartJ(g2)*I)/QQ(J*quartI(g2))
assert QQ(lam2).is_square()
g2 = [x/lam2 for x in g2]
assert quartI(g2) == I and quartJ(g2) == J
g2 = make_z_unit(g2)
z1 = z_inv(g1); z2 = z_inv(g2)
d3 = [sqfree(z1[k]*z2[k]) for k in range(3)]
g3 = quartic_from_delta(d3, 'sum3')
z3 = z_inv(g3)
print('g1', g1); print('g2', g2); print('g3', g3)
print('class z1', [sqfree(v) for v in z1])
print('class z2', [sqfree(v) for v in z2])
print('class z3', [sqfree(v) for v in z3])
assert [sqfree(v) for v in z1] == [QQ(1), QQ(193), QQ(193)]

# ---- контроль: выбор знака m не влияет на значение
base_val, gam, mm, zs, S, runs = pair_value(g1, g2, g3, signs=(1, 1, 1), reps=2, verbose=True)
print('<g1,g2> =', base_val, 'gamma =', gam, 'm =', mm)
sign_ctrl = []
for sg in [(1, 1, -1), (1, -1, 1), (-1, 1, 1), (-1, -1, -1), (1, -1, -1)]:
    v2, _, _, _, S2, _ = pair_value(g1, g2, g3, signs=sg, reps=1)
    sign_ctrl.append({'signs': list(sg), 'value': int(v2)})
    print('  знаки', sg, '->', v2, flush=True)
assert all(c['value'] == base_val for c in sign_ctrl), 'значение зависит от выбора m!'

# ---- контроль: обратный порядок аргументов <g2,g1>
rev_val, gam_rev, mm_rev, zs_rev, S_rev, runs_rev = pair_value(g2, g1, g3, reps=2, verbose=True)
print('<g2,g1> =', rev_val, 'gamma_rev =', gam_rev)
assert rev_val == base_val, 'НЕСИММЕТРИЯ!'

# ---- контроль: <g1,g1> должно быть 0 (знакопеременность), через тривиальную квартику
try:
    g_triv = quartic_from_delta([QQ(1), QQ(1), QQ(1)], 'triv')
    diag_val, _, _, _, _, _ = pair_value(g1, g1, g_triv, reps=1)
    print('<g1,g1> =', diag_val, ' (g_triv =', g_triv, ')')
except Exception as ex:
    g_triv = None; diag_val = None
    print('диагональ не посчитана:', repr(ex))

# ---- контроль: нормированная gamma даёт тот же продукт
cont = gcd([QQ(c) for c in gam])
gam_norm = [QQ(c)/cont for c in gam]
Sn = place_set(g1, gam_norm, g2[0])
runs_norm = []
for r in range(1):
    row = []
    for pl in [str(p) for p in Sn] + ['real']:
        pts = local_points(g1, pl, need=r+1, avoid=gam_norm)
        x, z = pts[r]
        row.append({'place': pl, 'x': str(x), 'z': str(z),
                    'g': str(ev(g1, x, z)), 'gamma': str(ev(gam_norm, x, z)),
                    'hilbert': int(my_hilbert(g2[0], ev(gam_norm, x, z), pl))})
    runs_norm.append(row)
norm_val = prod(v['hilbert'] for v in runs_norm[0])
print('нормированная gamma =', gam_norm, ' (scale 1/%s)' % cont, '-> продукт', norm_val)
assert norm_val == base_val

# ---- ELS всех трёх квартик
els_places = sorted(set([2, 3]) | supp(DISC) | set().union(*[supp(c) for g in (g1, g2, g3) for c in g if QQ(c).denominator() > 1]) if any(QQ(c).denominator() > 1 for g in (g1, g2, g3) for c in g) else set([2, 3]) | supp(DISC))
print('ELS-места:', els_places)
els = []
for g in (g1, g2, g3):
    row = []
    for pl in [str(p) for p in els_places] + ['real']:
        pts = local_points(g, pl, need=1)
        assert pts, ('ELS не подтверждена', pl, g)
        x, z = pts[0]
        row.append({'place': pl, 'x': str(x), 'z': str(z), 'value': str(ev(g, x, z))})
    els.append(row)
print('ELS подтверждена для всех трёх квартик', flush=True)

# ---- факторизации, нужные верификатору для контроля полноты набора мест
def facdict(v):
    v = QQ(v)
    d = {}
    for p, e in factor(v.numerator()): d[str(p)] = int(e)
    for p, e in factor(v.denominator()): d[str(p)] = -int(e)
    return d

a_lead = g2[0]
cert = {
    'pair': '(m,n) = (19,5)',
    'm': m, 'n': n, 's': str(s), 'b': str(b_par),
    'E_roots': [str(x) for x in E_roots],
    'target_class': [str(x) for x in DELTA],
    'shift': str(sh), 'u2': str(u2),
    'model_A': str(A), 'model_B': str(B),
    'model_roots': [str(x) for x in rr],
    'I': str(I), 'J': str(J), 'discriminant': str(DISC),
    'phi_roots': [str(x) for x in phis],
    'quartics': [[str(c) for c in g] for g in (g1, g2, g3)],
    'z_values': [[str(c) for c in zz_] for zz_ in (z1, z2, z3)],
    'm_values': [str(c) for c in mm],
    'gamma': [str(c) for c in gam],
    'gamma_normalized': [str(c) for c in gam_norm],
    'gamma_norm_scale': str(QQ(1)/cont),
    'a': str(a_lead),
    'factorizations': {'discriminant': facdict(DISC), 'a': facdict(a_lead),
                       'gamma_content': facdict(cont),
                       'gamma_coeffs': [facdict(c) for c in gam],
                       'g1_coeffs': [facdict(c) for c in g1]},
    'pair_runs': runs,
    'pair_runs_gamma_normalized': runs_norm,
    'reverse': {'a': str(g1[0]), 'gamma': [str(c) for c in gam_rev],
                'm_values': [str(c) for c in mm_rev], 'runs': runs_rev},
    'sign_controls': sign_ctrl,
    'diagonal_control': {'g_trivial': None if g_triv is None else [str(c) for c in g_triv],
                         'value': None if diag_val is None else int(diag_val)},
    'local_solubility': els,
    'pairing_value': int(base_val),
    'reverse_pairing_value': int(rev_val),
}
(OUT/'certificate_19_5.json').write_text(json.dumps(cert, indent=2)+'\n')
print('записан certificate_19_5.json', flush=True)
print('ИТОГ: <[g1],[g2]>_CT =', base_val, ' обратный порядок:', rev_val)
