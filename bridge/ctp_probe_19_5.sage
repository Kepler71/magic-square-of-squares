# Discovery stage for the Cassels-Tate certificate of G1 pair (m,n)=(19,5).
# Uses the project's Fisher implementation (descent/ctp_quartic.sage) ONLY to FIND
# quartics; the resulting certificate is verified afterwards by pure Python.
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random

out = Path('/home/kep/magicKube/bridge')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
random.seed(int(195195)); set_random_seed(int(195195))

m, n = 19, 5
s = QQ(m**2 + n**2) / 2
b = s * m**2 * n**2
roots = [-b, -s * m**4, -s * n**4]
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = prod(x - e for e in roots)
E0 = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E0.minimal_model()
iso = E0.isomorphism_to(M)
mr = [iso(E0([e, 0]))[0] for e in roots]
I, J = IJ_of_curve(M)
print('s,b', s, b, flush=True)
print('M', M.a_invariants(), flush=True)
print('minimal roots', mr, flush=True)
print('I,J', I, J, flush=True)

F = FisherCTP(QQ, I, J, verbose=True)
phis = [phi_of_root(M, e) for e in mr]
print('phis', phis, flush=True)
order = _comp_order(F, phis)
trip = [QQ(1), QQ(s), QQ(s)]
delta = F.E.from_comps([trip[i] for i in order])

record = {'m': m, 'n': n, 's': str(s), 'b': str(b),
          'E0_ainvs': [str(c) for c in E0.a_invariants()],
          'M_ainvs': [str(c) for c in M.a_invariants()],
          'orig_roots': [str(e) for e in roots],
          'minimal_roots': [str(e) for e in mr],
          'phi_roots': [str(p) for p in phis],
          'I': str(I), 'J': str(J), 'order': order,
          'target_triple': [str(t) for t in trip],
          'pairs': []}

def dump():
    (out / 'ctp_probe_19_5.json').write_text(json.dumps(record, indent=2) + '\n')

try:
    alarm(900)
    g1 = F.quartic_from_delta(delta, tag='target')
    cancel_alarm()
    record['target_quartic'] = list(map(str, g1))
    print('TARGET g1 =', g1, flush=True)
    dump()
except (Exception, AlarmInterrupt) as e:
    cancel_alarm()
    record['target_error'] = str(e)
    print('TARGET ERROR', e, flush=True)
    dump()
    raise SystemExit(1)

cover = pari(M).ell2cover()
print('ell2cover count', len(cover), flush=True)
for i, c in enumerate(cover):
    try:
        alarm(1200)
        q = R(c[0]); g2 = [QQ(q[4 - j]) for j in range(5)]
        scale = QQ(J * quartic_I(g2) / (quartic_J(g2) * I))
        assert scale.is_square(), 'scale not a square'
        g2 = [scale * a for a in g2]
        F.check_quartic(g2)
        g2 = F._make_z_unit(g2)
        g3 = F.quartic_from_delta(F.z_inv(g1) * F.z_inv(g2), tag='sum%d' % i)
        val, npl = F.pair(g1, g2, g3, reps=2)
        cancel_alarm()
        row = {'i': i, 'g2': list(map(str, g2)), 'g3': list(map(str, g3)),
               'pair': int(val), 'places': int(npl)}
        record['pairs'].append(row)
        print('PAIR', row['i'], 'value', row['pair'], 'places', row['places'], flush=True)
        dump()
        if val:
            print('NONZERO PAIRING FOUND at i =', i, flush=True)
            break
    except (Exception, AlarmInterrupt) as e:
        cancel_alarm()
        record['pairs'].append({'i': i, 'error': str(e)})
        print('ERROR', i, e, flush=True)
        dump()
dump()
print('DONE', flush=True)
