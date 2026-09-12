"""Этап 3': финальный сертификат CTP для G1 (19,5) c ЯВНЫМИ тождествами 2-накрытия.
Добавлено против первой версии: квадратичные формы T_i параметризации коники и тождества
    [delta t^2]_2 = 0,   X_num - r_i*G_raw = d_i*T_i^2,   G_raw = lam^2 * g1,
которые доказывают (проверкой многочленов), что y^2 = g1(x,z) — это ИМЕННО 2-накрытие
C_delta класса delta = (1,193,193), а не просто квартика с теми же инвариантами.
"""
from sage.all import *
import json, time, random
from pathlib import Path

src = open('/home/kep/magicKube/bridge/ctp_19_5_claude/pairing_19_5.py').read()
exec(src.split('# ------------------------------------------------ кандидаты g2 из ell2cover')[0])

OUT = Path('/home/kep/magicKube/bridge/ctp_19_5_claude')
random.seed(777); set_random_seed(777)
Rxz = PolynomialRing(QQ, ['xx', 'zz']); xx, zz = Rxz.gens()

def covering_data(dd, tag=''):
    """как quartic_from_delta, но возвращает ещё T_i, lam2, G_raw, X_num."""
    dd = [QQ(v) for v in dd]
    assert (dd[0]*dd[1]*dd[2]).is_square()
    cc = [dd[i]*Lc2[i] for i in range(3)]
    den = lcm([QQ(c).denominator() for c in cc])
    cc_int = [ZZ(c*den) for c in cc]
    gI = gcd(cc_int); cc_int = [ZZ(c/gI) for c in cc_int]
    cc_red = [ZZ(QQ(c).squarefree_part()) for c in cc_int]
    scal = [QQ(cc_int[i]/cc_red[i]).sqrt() for i in range(3)]
    con = Conic(QQ, cc_red)
    ok, pt = con.has_rational_point(point=True); assert ok
    par = con.parametrization(pt)[0].defining_polynomials()
    T = [Rxz(par[i])/scal[i] for i in range(3)]
    assert sum(cc[i]*T[i]**2 for i in range(3)) == 0
    G_raw = -sum(dd[i]*T[i]**2*Lc1[i] for i in range(3))
    X_num = sum(dd[i]*T[i]**2*Lc0[i] for i in range(3))
    for i in range(3):
        assert X_num - rr[i]*G_raw == dd[i]*T[i]**2, ('тождество накрытия', tag, i)
    gc = [QQ(G_raw.coefficient({xx: 4-i, zz: i})) for i in range(5)]
    lam2 = QQ(quartJ(gc)*I)/QQ(J*quartI(gc))
    assert QQ(lam2).is_square(), ('lam^2 не квадрат!', lam2)
    gfin = [c/lam2 for c in gc]
    assert quartI(gfin) == I and quartJ(gfin) == J
    g2u = make_z_unit(gfin)
    assert g2u == gfin, 'потребовалась GL2-замена; тождества накрытия надо переносить'
    for i in range(3):
        assert sameclass(z_inv(gfin)[i], dd[i])
    return gfin, T, lam2, G_raw, X_num

g1, T1, lam2_1, Graw1, Xnum1 = covering_data(DELTA, 'target')
print('g1 =', g1, ' lam^2 =', lam2_1, ' sqrt =', QQ(lam2_1).sqrt())
print('T_i =', T1)

Rp = PolynomialRing(QQ, 'x'); xp = Rp.gen()
cover = pari(M).ell2cover()
q = Rp(cover[3][0])
g2 = [QQ(q[4-j]) for j in range(5)]
l2 = QQ(quartJ(g2)*I)/QQ(J*quartI(g2)); assert QQ(l2).is_square()
g2 = [v/l2 for v in g2]
assert quartI(g2) == I and quartJ(g2) == J
g2 = make_z_unit(g2)
z1 = z_inv(g1); z2 = z_inv(g2)
d3 = [sqfree(z1[k]*z2[k]) for k in range(3)]
g3, T3, lam2_3, Graw3, Xnum3 = covering_data(d3, 'sum')
z3 = z_inv(g3)
print('g2 =', g2); print('g3 =', g3)
print('классы z:', [sqfree(v) for v in z1], [sqfree(v) for v in z2], [sqfree(v) for v in z3])
assert [sqfree(v) for v in z1] == [QQ(1), QQ(193), QQ(193)]

# ---- значение спаривания + контроли
base_val, gam, mm, zs, S, runs = pair_value(g1, g2, g3, signs=(1, 1, 1), reps=2, verbose=True)
print('<g1,g2> =', base_val)
sign_ctrl = []
for sg in [(1, 1, -1), (1, -1, 1), (-1, 1, 1), (-1, -1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]:
    v2 = pair_value(g1, g2, g3, signs=sg, reps=1)[0]
    sign_ctrl.append({'signs': list(sg), 'value': int(v2)})
assert all(c['value'] == base_val for c in sign_ctrl)
print('контроль знаков m: все 8 вариантов дали', base_val)

rev_val, gam_rev, mm_rev, zs_rev, S_rev, runs_rev = pair_value(g2, g1, g3, reps=2, verbose=True)
print('<g2,g1> =', rev_val)
assert rev_val == base_val

g_triv, _, _, _, _ = covering_data([QQ(1), QQ(1), QQ(1)], 'triv')
diag11 = pair_value(g1, g1, g_triv, reps=1)[0]
diag22 = pair_value(g2, g2, g_triv, reps=1)[0]
print('<g1,g1> =', diag11, '  <g2,g2> =', diag22)
val13 = pair_value(g1, g3, g2, reps=1)[0]
print('<g1,g3> =', val13, '(по Fisher 3.2(v) должно совпасть с <g1,g2>)')
assert val13 == base_val

cont = gcd([QQ(c) for c in gam])
gam_norm = [QQ(c)/cont for c in gam]
Sn = place_set(g1, gam_norm, g2[0])
rowN = []
for pl in [str(p) for p in Sn] + ['real']:
    x, z = local_points(g1, pl, need=1, avoid=gam_norm)[0]
    rowN.append({'place': pl, 'x': str(x), 'z': str(z), 'g': str(ev(g1, x, z)),
                 'gamma': str(ev(gam_norm, x, z)),
                 'hilbert': int(my_hilbert(g2[0], ev(gam_norm, x, z), pl))})
norm_val = prod(v['hilbert'] for v in rowN)
assert norm_val == base_val
print('нормированная gamma', gam_norm, '-> ', norm_val)

# ---- ELS
els_places = sorted(set([2, 3]) | supp(DISC))
els = []
for g in (g1, g2, g3):
    row = []
    for pl in [str(p) for p in els_places] + ['real']:
        pts = local_points(g, pl, need=1)
        assert pts, ('ELS не подтверждена', pl)
        x, z = pts[0]
        row.append({'place': pl, 'x': str(x), 'z': str(z), 'value': str(ev(g, x, z))})
    els.append(row)
print('ELS всех трёх квартик подтверждена на', els_places, '+ real')

def facdict(v):
    v = QQ(v); d = {}
    for p, e in factor(v.numerator()): d[str(p)] = int(e)
    for p, e in factor(v.denominator()): d[str(p)] = -int(e)
    return d
def qform(P):
    return [str(QQ(P.coefficient({xx: 2-i, zz: i}))) for i in range(3)]

cert = {
    'what': 'Cassels-Tate pairing certificate, G1 family, (m,n)=(19,5)',
    'm': m, 'n': n, 's': str(s), 'b': str(b_par),
    'E_roots': [str(v) for v in E_roots],
    'target_class': [str(v) for v in DELTA],
    'shift': str(sh), 'u2': str(u2),
    'model_A': str(A), 'model_B': str(B), 'model_roots': [str(v) for v in rr],
    'I': str(I), 'J': str(J), 'discriminant': str(DISC),
    'phi_roots': [str(v) for v in phis],
    'quartics': [[str(c) for c in g] for g in (g1, g2, g3)],
    'z_values': [[str(c) for c in zv] for zv in (z1, z2, z3)],
    'm_values': [str(c) for c in mm],
    'gamma': [str(c) for c in gam],
    'gamma_normalized': [str(c) for c in gam_norm],
    'a': str(g2[0]),
    'covering_g1': {'T': [qform(v) for v in T1], 'lambda': str(QQ(lam2_1).sqrt()),
                    'D': [str(v) for v in Dl],
                    'G_raw': [str(QQ(Graw1.coefficient({xx: 4-i, zz: i}))) for i in range(5)],
                    'X_num': [str(QQ(Xnum1.coefficient({xx: 4-i, zz: i}))) for i in range(5)]},
    'factorizations': {'discriminant': facdict(DISC), 'a': facdict(g2[0]),
                       'gamma': [facdict(c) for c in gam],
                       'g1': [facdict(c) for c in g1]},
    'pair_runs': runs,
    'pair_run_gamma_normalized': [rowN],
    'reverse': {'a': str(g1[0]), 'gamma': [str(c) for c in gam_rev],
                'm_values': [str(c) for c in mm_rev], 'runs': runs_rev},
    'sign_controls': sign_ctrl,
    'diagonal': {'g_trivial': [str(c) for c in g_triv], 'g1g1': int(diag11), 'g2g2': int(diag22)},
    'g1g3': int(val13),
    'local_solubility': els,
    'pairing_value': int(base_val),
    'reverse_pairing_value': int(rev_val),
}
(OUT/'certificate_19_5.json').write_text(json.dumps(cert, indent=2)+'\n')
import shutil
shutil.copy(str(OUT/'certificate_19_5.json'), '/home/kep/magicKube/bridge/ctp_cert_19_5.json')
print('ИТОГ: <[g1],[g2]>_CT =', base_val, ' обратный:', rev_val, ' диагональ:', diag11, diag22)
