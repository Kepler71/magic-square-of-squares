# -*- coding: utf-8 -*-
"""ОТКРЫТИЕ (не проверка): ищем g2 в Sel^2(E/Q) с <g1,g2>_CT != 0 для пары (m,n)=(11,4).

g1 — квартика целевого класса delta=(1,274,274) (тот, который обязан иметь любая
конечная точка C_(11,4)).  Если найдётся g2 с ненулевым спариванием, то [g1] не в
образе E(Q)/2E(Q) — независимо от ранга.

Здесь используется Sage/PARI и реализация Фишера descent/ctp_quartic.sage.
Результат этого файла — ТОЛЬКО кандидатные квартики; всё доказательство
перепроверяется отдельно скриптом ctp_cert_11_4_ctnonzero.py на чистом Python.
"""
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random, time

OUT = Path('/home/kep/magicKube/bridge')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
random.seed(114114); set_random_seed(114114)

m, n = 11, 4
s = QQ(137) / 2
b = s * m ** 2 * n ** 2
roots_X = [-b, -s * m ** 4, -s * n ** 4]
xv = polygen(QQ, 'x')
f = prod(xv - 4 * e for e in roots_X)            # целая модель x = 4X, 4 — квадрат
E = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E.minimal_model()
iso = E.isomorphism_to(M)
mr = [iso(E([4 * e, 0]))[0] for e in roots_X]    # корни в согласованном порядке
assert M.a1() == 0 and M.a3() == 0
I, J = IJ_of_curve(M)
F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, e) for e in mr]
order = _comp_order(F, phis)
TRIP = [QQ(1), QQ(274), QQ(274)]
delta = F.E.from_comps([TRIP[i] for i in order])

rec = {'m': m, 'n': n, 's': str(s), 'b': str(b), 'M_ainvs': [str(a) for a in M.ainvs()],
       'minimal_roots': [str(e) for e in mr], 'I': str(I), 'J': str(J),
       'phi_roots': [str(p) for p in phis], 'order': order, 'target_trip': [str(t) for t in TRIP],
       'delta': str(delta), 'pairs': [], 'errors': []}
print('I =', I, ' J =', J, ' roots =', mr, flush=True)


def save():
    (OUT / 'ctp_discover_11_4.json').write_text(json.dumps(rec, indent=2) + '\n')


t0 = time.time()
alarm(600)
g1 = F.quartic_from_delta(delta, tag='target')
cancel_alarm()
rec['target_quartic'] = list(map(str, g1))
print('TARGET g1 =', g1, ' (%.1f s)' % (time.time() - t0), flush=True)
save()

cover = pari(M).ell2cover()
print('ell2cover: %d квартик' % len(cover), flush=True)
Rx = PolynomialRing(QQ, 'x')
for i, c in enumerate(cover):
    try:
        alarm(900)
        q = Rx(c[0])
        g2 = [QQ(q[4 - j]) for j in range(5)]
        # согласование инвариантов: I(g)=lam^4 I, J(g)=lam^6 J -> g/lam^2
        sc = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
        assert sc.is_square(), 'масштаб не квадрат: %s' % sc
        g2 = [sc * a for a in g2]
        F.check_quartic(g2)
        g2 = F._make_z_unit(g2)
        g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum%d' % i)
        val, npl = F.pair(g1, g2, g3, reps=2)
        cancel_alarm()
        row = {'i': i, 'scale': str(sc), 'g2': list(map(str, g2)), 'g3': list(map(str, g3)),
               'pair': int(val), 'places': int(npl)}
        rec['pairs'].append(row)
        print('PAIR', i, 'val =', val, 'places =', npl, ' (%.1f s)' % (time.time() - t0), flush=True)
        save()
        if val:
            print('>>> НЕНУЛЕВОЕ СПАРИВАНИЕ на i =', i, flush=True)
            break
    except (Exception, AlarmInterrupt) as e:
        cancel_alarm()
        rec['errors'].append({'i': i, 'error': repr(e)})
        print('ERROR', i, repr(e), flush=True)
        save()
save()
print('DONE %.1f s' % (time.time() - t0), flush=True)
