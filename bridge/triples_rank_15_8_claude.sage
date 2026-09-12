# -*- coding: utf-8 -*-
# Обзор ВСЕХ десяти тройных произведений пяти форм (15,8): род 2, ранг Jac.
# Цель: найти фактор с rank Jac <= 1 (тогда работает КЛАССИЧЕСКИЙ Шаботи, род 2),
# либо честно зафиксировать, что такого нет.
# Пары (квартики) уже сделаны Codex: bridge/qc40/elliptic_inventory.log — все ранги >= 1.
from sage.all import *
from itertools import combinations
import json

R = PolynomialRing(QQ, 't'); t = R.gen()
s = QQ(289)/2
forms = {'F0': 225 + 64*t**2,
         'F4': s*(1 + t**2),
         'F8': 64 + 225*t**2,
         'L' : s*(1 + t**2) - 240*t,
         'U' : s*(1 + t**2) + 240*t}

def rank_of_quartic_or_cubic(f):
    """Якобиан кривой y^2=f(u) (deg f = 3 или 4) как эллиптическая кривая + ранг."""
    R2 = PolynomialRing(QQ, ['uu','yy']); uu, yy = R2.gens()
    E = Jacobian(yy*yy - f(uu)).minimal_model()
    res = {'ainvs': str(E.ainvs()), 'conductor': str(E.conductor())}
    try:
        rk = E.pari_curve().ellrank()
        res['pari_rank_interval'] = [int(rk[0]), int(rk[1])]
        res['pari_points'] = [str(P) for P in rk[3]]
    except Exception as e:
        res['pari_error'] = repr(e)[:120]
    try:
        res['eclib_rank_bound'] = int(E.rank_bound())
    except Exception as e:
        res['eclib_error'] = repr(e)[:120]
    try:
        res['sage_rank'] = int(E.rank())
    except Exception as e:
        res['sage_rank_error'] = repr(e)[:120]
    try:
        res['analytic_rank'] = int(E.analytic_rank())
    except Exception as e:
        pass
    return res

rows = []
for names in combinations(sorted(forms), 3):
    f = prod(forms[n] for n in names)
    rec = {'names': list(names), 'f': str(f), 'deg': int(f.degree()),
           'squarefree': bool(gcd(f, f.derivative()).degree() == 0)}
    # чётность по t: f(-t) == f(t)?
    even_t = (f(-t) == f)
    # чётность по x, где x=(t-1)/(t+1): подставляем t=(1+x)/(1-x), домножаем на (1-x)^6
    Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
    fx = Rx(((1-x)**6 * f((1+x)/(1-x))).numerator()) if True else None
    try:
        fx = R(0)
        FR = FractionField(Rx)
        val = f(FR((1+x)/(1-x))) * FR((1-x)**6)
        fx = Rx(val.numerator()/val.denominator())
    except Exception as e:
        fx = None
    even_x = (fx is not None) and (fx(-x) == fx)
    rec['even_in_t'] = bool(even_t); rec['even_in_x'] = bool(even_x)
    rec['f_in_x'] = str(fx)
    # разложение биэллиптической кривой
    if even_t or even_x:
        g6 = f if even_t else fx
        Ru = PolynomialRing(QQ, 'u'); u = Ru.gen()
        cs = g6.list()
        g3 = Ru([cs[2*i] for i in range(4)])   # g6(z)=g3(z^2)
        assert g6 == sum(g3.list()[i]*g6.parent().gen()**(2*i) for i in range(4)), "распад не сошёлся"
        rec['coord'] = 't' if even_t else 'x'
        rec['E_plus']  = rank_of_quartic_or_cubic(g3)        # v^2 = g3(u)
        rec['E_minus'] = rank_of_quartic_or_cubic(u*g3)      # v^2 = u*g3(u)
        def lo_hi(d):
            if 'pari_rank_interval' in d: return d['pari_rank_interval']
            return [None, d.get('eclib_rank_bound')]
        a, b = lo_hi(rec['E_plus']), lo_hi(rec['E_minus'])
        rec['rank_Jac'] = [None if a[0] is None or b[0] is None else a[0]+b[0],
                           None if a[1] is None or b[1] is None else a[1]+b[1]]
    else:
        # эмпирика: факторизуется ли L-многочлен над Q при многих p (признак распада Jac)
        C = HyperellipticCurve(f) if f.degree() in (5,6) else None
        splits = []
        for p in primes(5, 120):
            try:
                if f.discriminant() % p == 0: continue
                Cp = HyperellipticCurve(f.change_ring(GF(p)))
                Lp = Cp.frobenius_polynomial()
                splits.append((p, len(Lp.factor()) > 1 or Lp.is_irreducible() is False))
            except Exception:
                pass
        rec['Lpoly_reducible_fraction'] = [int(sum(1 for _, b in splits if b)), int(len(splits))]
        rec['rank_Jac'] = 'не вычислен (кривая не биэллиптична в t/x)'
    rows.append(rec)
    print(json.dumps({k: v for k, v in rec.items() if k not in ('f', 'f_in_x')},
                     ensure_ascii=False), flush=True)

open('/home/kep/magicKube/bridge/triples_rank_15_8_claude.json', 'w').write(
    json.dumps(rows, indent=1, ensure_ascii=False))
print("saved")
