from sage.all import *
from cysignals.alarm import alarm, cancel_alarm
from cysignals.signals import AlarmInterrupt
from pathlib import Path
import json, random
out = Path('/home/kep/magicKube/bridge/ctp_11_4_v2')
load(str(out/'ctp_quartic_snapshot.sage'))
random.seed(int(114)); set_random_seed(114)

m, n = 11, 4
s = QQ(137)/2
b = s*m*m*n*n
roots = [-b, -s*m**4, -s*n**4]
R = PolynomialRing(QQ,'x'); x = R.gen()
f = prod(x-e for e in roots)
E = EllipticCurve([0, f[2], 0, f[1], f[0]])
M = E.minimal_model()
iso = E.isomorphism_to(M)
mr = [iso(E([e,0]))[0] for e in roots]
I, J = IJ_of_curve(M)
print("M =", M.ainvs(), flush=True)
print("minimal roots =", mr, flush=True)
print("I =", I, "J =", J, flush=True)

F = FisherCTP(QQ, I, J, verbose=True)
order = _comp_order(F, [phi_of_root(M, e) for e in mr])
trip = [QQ(1), QQ(274), QQ(274)]
delta = F.E.from_comps([trip[i] for i in order])
print("order =", order, " delta =", delta, flush=True)

rec = {'m':int(m),'n':int(n),'s':str(s),'b':str(b),'orig_roots':[str(e) for e in roots],
       'M_ainvs':[str(a) for a in M.ainvs()],'minimal_roots':[str(e) for e in mr],
       'I':str(I),'J':str(J),'order':[int(o) for o in order],'delta':str(delta),'trip':[str(t) for t in trip],
       'pairs':[], 'errors':[]}

def dump():
    (out/'probe_v2.json').write_text(json.dumps(rec, indent=2)+'\n')

try:
    alarm(600); g1 = F.quartic_from_delta(delta, tag='target'); cancel_alarm()
except (Exception, AlarmInterrupt) as e:
    cancel_alarm(); rec['fatal'] = 'quartic_from_delta: '+str(e); dump(); print('FATAL', e); raise SystemExit
rec['g1'] = list(map(str, g1))
print("TARGET g1 =", g1, flush=True)
print("z(g1) =", F.z_inv(g1), flush=True)
dump()

cover = pari(M).ell2cover()
print("ell2cover count =", len(cover), flush=True)
for i, c in enumerate(cover):
    try:
        alarm(1200)
        q = R(c[0]); g2 = [QQ(q[4-j]) for j in range(5)]
        scale = QQ(J*quartic_I(g2)/(quartic_J(g2)*I)); assert scale.is_square(), "scale not square"
        g2 = [scale*a for a in g2]
        F.check_quartic(g2)
        g2 = F._make_z_unit(g2)
        g3 = F.quartic_from_delta(F.z_inv(g1)*F.z_inv(g2), tag='sum%d'%i)
        val, npl = F.pair(g1, g2, g3, reps=2)
        cancel_alarm()
        row = {'i':i,'scale':str(scale),'g2':list(map(str,g2)),'g3':list(map(str,g3)),
               'pair':int(val),'places':int(npl)}
        rec['pairs'].append(row); dump()
        print("PAIR", i, "value", val, "places", npl, flush=True)
        if val:
            print("NONZERO FOUND at i =", i, flush=True)
            break
    except (Exception, AlarmInterrupt) as e:
        cancel_alarm(); rec['errors'].append({'i':i,'error':str(e)}); dump()
        print("ERROR", i, e, flush=True)
dump()
print("DONE", flush=True)
