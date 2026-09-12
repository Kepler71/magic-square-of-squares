# -*- coding: utf-8 -*-
# v2: только быстрые ранговые границы (PARI ellrank + eclib rank_bound), без analytic_rank.
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
from itertools import combinations
import json

R = PolynomialRing(QQ, 't'); t = R.gen()
s = QQ(289)/2
forms = {'F0': 225 + 64*t**2, 'F4': s*(1 + t**2), 'F8': 64 + 225*t**2,
         'L': s*(1 + t**2) - 240*t, 'U': s*(1 + t**2) + 240*t}

def ell(f, secs=90):
    R2 = PolynomialRing(QQ, ['uu', 'yy']); uu, yy = R2.gens()
    E = Jacobian(yy*yy - f(uu)).minimal_model()
    res = {'ainvs': str(E.ainvs()), 'conductor': str(E.conductor()),
           'torsion': str(E.torsion_subgroup().invariants())}
    try:
        alarm(secs); rk = E.pari_curve().ellrank()
        res['pari'] = [int(rk[0]), int(rk[1])]
        res['pari_points'] = [str(P) for P in rk[3]]
    except (Exception, AlarmInterrupt) as e:
        res['pari_err'] = type(e).__name__
    finally: cancel_alarm()
    try:
        alarm(secs); res['eclib_bound'] = int(E.rank_bound())
    except (Exception, AlarmInterrupt) as e:
        res['eclib_err'] = type(e).__name__
    finally: cancel_alarm()
    return res

rows = []
for names in combinations(sorted(forms), 3):
    f = prod(forms[n] for n in names)
    rec = {'names': list(names), 'squarefree': bool(gcd(f, f.derivative()).degree() == 0)}
    even_t = (f(-t) == f)
    Rx = PolynomialRing(QQ, 'x'); x = Rx.gen(); FR = FractionField(Rx)
    val = f(FR((1+x)/(1-x))) * FR((1-x)**6)
    fx = Rx(val.numerator()/val.denominator())
    even_x = (fx(-x) == fx)
    rec['even_in_t'] = bool(even_t); rec['even_in_x'] = bool(even_x)
    if even_t or even_x:
        g6 = f if even_t else fx
        Ru = PolynomialRing(QQ, 'u'); u = Ru.gen()
        cs = g6.list(); g3 = Ru([cs[2*i] for i in range(4)])
        rec['coord'] = 't' if even_t else 'x'
        rec['g3'] = str(g3)
        rec['E_plus'] = ell(g3)
        rec['E_minus'] = ell(u*g3)
        def iv(d):
            return d.get('pari', [None, d.get('eclib_bound')])
        a, b = iv(rec['E_plus']), iv(rec['E_minus'])
        rec['rank_Jac'] = [None if a[0] is None or b[0] is None else a[0]+b[0],
                           None if a[1] is None or b[1] is None else a[1]+b[1]]
    else:
        rec['rank_Jac'] = 'нет биэллиптического расщепления в t/x; ранг не вычислен'
    rows.append(rec)
    print(json.dumps(rec, ensure_ascii=False), flush=True)
open('/home/kep/magicKube/bridge/triples_rank_15_8_claude.json', 'w').write(json.dumps(rows, indent=1, ensure_ascii=False))
print("saved")
