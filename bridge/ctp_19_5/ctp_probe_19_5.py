"""Поиск сертификата CTP для G1 (m,n)=(19,5).  Запуск: sage -python ... нет — просто `sage ctp_probe_19_5.py`
Использует descent/ctp_quartic.sage (снимок рядом).  Ранг НЕ вычисляется и НЕ используется."""
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random, time

out = Path('/home/kep/magicKube/bridge/ctp_19_5')
load(str(out / 'ctp_quartic_snapshot.sage'))
random.seed(195); set_random_seed(195)

m, n = 19, 5
s = QQ(m ** 2 + n ** 2) / 2            # = 193
b = s * m * m * n * n                  # = 1741825
roots = [-b, -s * m ** 4, -s * n ** 4]
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = prod(x - e for e in roots)
E = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E.minimal_model()
iso = E.isomorphism_to(M)
mr = [iso(E([e, 0]))[0] for e in roots]

I, J = IJ_of_curve(M)
F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, e) for e in mr]
order = _comp_order(F, phis)
trip = [QQ(1), QQ(193), QQ(193)]       # delta = (1, s, s), s = 193 бесквадратно
delta = F.E.from_comps([trip[i] for i in order])

record = {'m': m, 'n': n, 's': str(s), 'b': str(b),
          'orig_roots': [str(e) for e in roots],
          'M_ainvs': [str(a) for a in M.ainvs()],
          'minimal_roots': [str(e) for e in mr],
          'phi_roots': [str(p) for p in phis],
          'I': str(I), 'J': str(J), 'order': order,
          'delta_trip': [str(t) for t in trip],
          'delta': str(delta), 'pairs': []}
print('SETUP', record['M_ainvs'], record['minimal_roots'], flush=True)

t0 = time.time()
try:
    alarm(1800)
    g1 = F.quartic_from_delta(delta, tag='target')
    cancel_alarm()
    record['target_quartic'] = list(map(str, g1))
    record['target_seconds'] = time.time() - t0
    print('TARGET', g1, 'in', time.time() - t0, 's', flush=True)
    (out / 'ctp_probe_19_5.json').write_text(json.dumps(record, indent=2) + '\n')

    cover = pari(M).ell2cover()
    record['n_cover'] = len(cover)
    for i, c in enumerate(cover):
        t1 = time.time()
        try:
            alarm(1800)
            q = R(c[0]); g2 = [QQ(q[4 - j]) for j in range(5)]
            scale = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
            assert scale.is_square(), f'scale {scale} не квадрат'
            g2 = [scale * a for a in g2]
            F.check_quartic(g2)
            g2 = F._make_z_unit(g2)
            g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum' + str(i))
            val, npl = F.pair(g1, g2, g3, reps=2)
            cancel_alarm()
            row = {'i': i, 'scale': str(scale), 'g2': list(map(str, g2)),
                   'g3': list(map(str, g3)), 'pair': int(val), 'places': int(npl),
                   'seconds': time.time() - t1}
            record['pairs'].append(row)
            print('PAIR', row, flush=True)
            (out / 'ctp_probe_19_5.json').write_text(json.dumps(record, indent=2) + '\n')
            if val:
                break
        except (Exception, AlarmInterrupt) as e:
            cancel_alarm()
            record['pairs'].append({'i': i, 'error': repr(e), 'seconds': time.time() - t1})
            print('ERROR', i, repr(e), flush=True)
            (out / 'ctp_probe_19_5.json').write_text(json.dumps(record, indent=2) + '\n')
except (Exception, AlarmInterrupt) as e:
    cancel_alarm()
    record['error'] = repr(e)
    print('ERROR', repr(e), flush=True)
(out / 'ctp_probe_19_5.json').write_text(json.dumps(record, indent=2) + '\n')
print('DONE', flush=True)
