"""Контроли к сертификату (19,5).
(1) Независимость от PARI: g2 пересобирается МОЕЙ конструкцией из класса (193,193,1).
(2) Сквозной контроль соответствия «класс спуска -> квартика -> z(g)»:
    берём настоящие рациональные точки M(Q), считаем их образ (x-r_i),
    строим квартику и ищем на ней рациональную точку (она обязана быть).
"""
from sage.all import *
import json, random
from pathlib import Path

src = open('/home/kep/magicKube/bridge/ctp_19_5_claude/pairing_19_5.py').read()
exec(src.split('# ------------------------------------------------ кандидаты g2 из ell2cover')[0])
OUT = Path('/home/kep/magicKube/bridge/ctp_19_5_claude')
random.seed(31337); set_random_seed(31337)

res = {}

# (1) собственная квартика для класса z2
cls2 = [QQ(193), QQ(193), QQ(1)]
g2b = quartic_from_delta(cls2, 'own_g2')
print('своя g2 для класса (193,193,1):', g2b, flush=True)
z1 = z_inv(g1); z2b = z_inv(g2b)
d3 = [sqfree(z1[k]*z2b[k]) for k in range(3)]
g3b = quartic_from_delta(d3, 'own_g3')
val_b, gam_b, mm_b, zs_b, S_b, runs_b = pair_value(g1, g2b, g3b, reps=2, verbose=True)
print('<g1, своя g2> =', val_b, ' gamma =', gam_b, flush=True)
res['own_g2'] = {'g2': [str(c) for c in g2b], 'g3': [str(c) for c in g3b],
                 'gamma': [str(c) for c in gam_b], 'value': int(val_b),
                 'nontrivial_places': [v['place'] for v in runs_b[0] if v['hilbert'] == -1]}

# ELS своей g2 и g3
els_pl = sorted(set([2, 3]) | supp(DISC))
for nm, g in [('own_g2', g2b), ('own_g3', g3b)]:
    wit = []
    for pl in [str(p) for p in els_pl] + ['real']:
        pts = local_points(g, pl, need=1)
        assert pts, ('ELS не подтверждена', nm, pl)
        wit.append({'place': pl, 'x': str(pts[0][0]), 'z': str(pts[0][1])})
    res[nm + '_ELS'] = wit
print('ELS своих g2,g3 подтверждена', flush=True)

# (2) сквозной контроль: настоящие точки M(Q) -> класс -> квартика -> рациональная точка
def rat_point_on_quartic(g, N=400):
    for den in range(1, 40):
        for num in range(-N, N+1):
            x = QQ(num)/QQ(den)
            v = ev(g, x, QQ(1))
            if v != 0 and QQ(v).is_square():
                return (x, QQ(1), QQ(v).sqrt())
    if QQ(g[0]) != 0 and QQ(g[0]).is_square():
        return (QQ(1), QQ(0), QQ(g[0]).sqrt())
    return None

pts_found = []
# ищем точки M(Q) перебором x (без вызова .rank()/.gens())
cnt = 0
for x0 in range(-2000000, 2000001, 1):
    pass  # слишком медленно; вместо этого используем известные простые точки
# берём точки малой высоты через перебор x в окрестности корней
cand_x = []
for r0 in rr:
    for dx in range(-3000, 3001):
        cand_x.append(r0 + dx)
for x0 in cand_x:
    y2 = x0**3 + A*x0 + B
    if y0_ok := (y2 > 0 and ZZ(y2).is_square()):
        pts_found.append(ZZ(x0))
pts_found = sorted(set(pts_found))
print('точки M(Q) с целым x рядом с корнями:', pts_found[:20], '... всего', len(pts_found), flush=True)

ctrl = []
for x0 in pts_found[:8]:
    comps = [QQ(x0 - ri) for ri in rr]
    if any(c == 0 for c in comps):
        continue
    cl = [sqfree(c) for c in comps]
    assert (cl[0]*cl[1]*cl[2]).is_square()
    gq = quartic_from_delta(cl, 'pt%d' % x0)
    pt = rat_point_on_quartic(gq)
    ctrl.append({'x': str(x0), 'class': [str(c) for c in cl], 'quartic': [str(c) for c in gq],
                 'rational_point': None if pt is None else [str(v) for v in pt]})
    print('x =', x0, 'класс', cl, '-> квартика', gq, '-> точка', pt, flush=True)
res['descent_roundtrip'] = ctrl

# контроль: класс (1,193,193) НЕ встречается среди образов найденных точек
print('классы найденных точек:', [c['class'] for c in ctrl])
(OUT/'controls_19_5.json').write_text(json.dumps(res, indent=2)+'\n')
print('DONE CONTROLS', flush=True)
