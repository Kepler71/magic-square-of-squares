# Fable, 14.09.2026. Для закрытий, где единственное доказательство ранга 0 — L(E,1) != 0 через Sage L_ratio
# (32 шестиклеточных + 5 рода 1), считаю L(E,1) другой программой: PARI ellL1 (с точностью 38 знаков)
# и сравниваю L(E,1)/Omega с точным L_ratio Sage. Omega = период решётки (для Delta>0 — удвоенный).
import sys, json
sys.path.insert(0, '/home/kep/magicKube/fable_census_review')
from common import *
from sage.all import QQ, ZZ, EllipticCurve, pari, RealField
import six_review as six
R = RealField(120)
pari.set_real_precision(38)

items = []
for l in open('/home/kep/magicKube/fable_census_review/kol_review.jsonl'):
    d = json.loads(l)
    E, _ = six.model(*d['abc'], d['kind'])
    items.append(('six', d['slope'], d['abc'], d['kind'], E, QQ(d['L_ratio'])))
for l in open('/home/kep/magicKube/fable_census_review/g1_review.jsonl'):
    d = json.loads(l)
    if d['claude_method'] != 'L_ratio':
        continue
    E = EllipticCurve(QQ, [QQ(t) for t in d['E_ainv']])
    items.append(('g1', d['slope'], d['cells'], '', E, QQ(d['L_ratio'])))
print('кривых:', len(items), flush=True)
bad = 0
for kind, slope, cells, k, E, Lr in items:
    Em = E.minimal_model()
    ep = pari(Em)
    L1 = R(pari.ellL1(ep, 0))
    om = R(ep.omega()[0])
    if Em.discriminant() > 0:
        om = 2 * om          # две вещественные компоненты: Omega = 2 * omega_1
    ratio = L1 / om
    ok = abs(ratio - R(Lr)) < R(1e-20)
    if not ok:
        bad += 1
    print(kind, slope, cells, k, 'L(E,1)=%.12g' % L1, 'L/Omega(PARI)=%.20f' % ratio, 'Sage L_ratio=%s' % Lr, 'совпало' if ok else 'РАСХОЖДЕНИЕ', flush=True)
print('расхождений:', bad)
