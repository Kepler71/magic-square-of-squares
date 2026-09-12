"""Этап ОТКРЫТИЯ сертификата CTP для пары G1 (m,n)=(11,4).

Строит квартику g1 целевого класса delta=(1,274,274), перебирает 2-накрытия PARI
как g2, строит g3 = -(g1+g2) и считает <g1,g2>_CT по теореме 3.1 Фишера.
Ранг E НИГДЕ не используется.
"""
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random

out = Path(__file__).resolve().parent
load(str(out / 'ctp_quartic_snapshot.sage'))
random.seed(114); set_random_seed(114)

# --- исходное семейство G1, пара (11,4) ---
m, n = 11, 4
s = QQ(137) / 2
b = s * m * m * n * n                      # 132616
roots = [-b, -s * m ** 4, -s * n ** 4]     # e1, e2, e3 в порядке отображения C -> E
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = prod(x - e for e in roots)
E = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E.minimal_model()
iso = E.isomorphism_to(M)
mr = [iso(E([e, 0]))[0] for e in roots]    # корни минимальной модели В ТОМ ЖЕ ПОРЯДКЕ

I, J = IJ_of_curve(M)
F = FisherCTP(QQ, I, J, verbose=True)
order = _comp_order(F, [phi_of_root(M, e) for e in mr])
trip = [QQ(1), QQ(274), QQ(274)]           # delta = (1, s, s), s = 137/2 ~ 274 mod квадратов
delta = F.E.from_comps([trip[i] for i in order])

record = {'m': m, 'n': n, 's': str(s), 'b': str(b),
          'E_ainvs': [str(a) for a in E.ainvs()],
          'M_ainvs': [str(a) for a in M.ainvs()],
          'minimal_roots': [str(e) for e in mr],
          'I': str(I), 'J': str(J), 'order': order,
          'delta_target': [str(t) for t in trip],
          'delta': str(delta), 'pairs': []}
print('M =', M, flush=True)
print('minimal roots =', mr, flush=True)
print('I,J =', I, J, flush=True)

BUDGET = 900

try:
    alarm(BUDGET)
    g1 = F.quartic_from_delta(delta, tag='target')
    cancel_alarm()
    record['target_quartic'] = list(map(str, g1))
    print('TARGET', g1, flush=True)
    (out / 'ctp_probe_11_4.json').write_text(json.dumps(record, indent=2) + '\n')

    cover = pari(M).ell2cover()
    print('ell2cover: получено', len(cover), 'накрытий', flush=True)
    for i, c in enumerate(cover):
        try:
            alarm(BUDGET)
            q = R(c[0])
            g2 = [QQ(q[4 - j]) for j in range(5)]
            # согласование инвариантов (квартики PARI могут иметь I/16, J/64)
            scale = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
            assert scale.is_square(), 'масштаб не квадрат'
            g2 = [scale * a for a in g2]
            F.check_quartic(g2)
            g2 = F._make_z_unit(g2)
            g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum' + str(i))
            val, npl = F.pair(g1, g2, g3, reps=2)
            cancel_alarm()
            row = {'i': i, 'g2': list(map(str, g2)), 'g3': list(map(str, g3)),
                   'pair': int(val), 'places': int(npl)}
            record['pairs'].append(row)
            print('PAIR', row, flush=True)
            (out / 'ctp_probe_11_4.json').write_text(json.dumps(record, indent=2) + '\n')
            if val:
                break
        except (Exception, AlarmInterrupt) as e:
            cancel_alarm()
            record['pairs'].append({'i': i, 'error': repr(e)})
            print('ERROR', i, repr(e), flush=True)
            (out / 'ctp_probe_11_4.json').write_text(json.dumps(record, indent=2) + '\n')
except (Exception, AlarmInterrupt) as e:
    cancel_alarm()
    record['error'] = repr(e)
    print('ERROR', repr(e), flush=True)

(out / 'ctp_probe_11_4.json').write_text(json.dumps(record, indent=2) + '\n')
print('DONE', flush=True)
