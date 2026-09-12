# -*- coding: utf-8 -*-
# Все подмножества пяти форм (15,8) размера 4 и 5: чётность в t или x, факторы по инволюции,
# род и ранг факторов. Цель — найти фактор с rank Jac <= 1 (классический Шаботи) или rank 0.
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
from itertools import combinations
import json

R = PolynomialRing(QQ, 't'); t = R.gen()
s = QQ(289)/2
forms = {'F0': 225 + 64*t**2, 'F4': s*(1 + t**2), 'F8': 64 + 225*t**2,
         'L': s*(1 + t**2) - 240*t, 'U': s*(1 + t**2) + 240*t}

def ell(f, secs=120):
    R2 = PolynomialRing(QQ, ['uu', 'yy']); uu, yy = R2.gens()
    E = Jacobian(yy*yy - f(uu)).minimal_model()
    res = {'ainvs': str(E.ainvs())}
    try:
        alarm(secs); rk = E.pari_curve().ellrank(); res['pari'] = [int(rk[0]), int(rk[1])]
    except (Exception, AlarmInterrupt) as e: res['pari_err'] = type(e).__name__
    finally: cancel_alarm()
    try:
        alarm(secs); res['eclib_bound'] = int(E.rank_bound())
    except (Exception, AlarmInterrupt) as e: res['eclib_err'] = type(e).__name__
    finally: cancel_alarm()
    return res

def squarefree_scale(g):
    """убрать квадратные множители из старшего множителя/знаменателей (не меняет кривую)"""
    g = g * lcm([c.denominator() for c in g.list()])
    cont = gcd(g.list())
    g = g / cont
    sq = 1
    for pr, e in ZZ(cont).factor():
        if e % 2: sq *= pr
    return (g*sq).change_ring(QQ)

def analyse(g, label, rows):
    """y^2 = g(u): определить род и, если род 2 — попробовать биэллиптическое расщепление."""
    d = g.degree(); genus = (d - 1)//2 if d % 2 else (d - 2)//2
    rec = {'label': label, 'deg': int(d), 'genus': int(genus), 'g': str(g)}
    if genus == 1:
        rec['elliptic'] = ell(g)
    elif genus == 2:
        Ru = g.parent(); u = Ru.gen()
        # чётность прямо / после u -> -u уже учтена; пробуем u -> (u-1)/(u+1) и u -> 1/u
        cands = []
        gg = g if d == 6 else g.parent()(g)      # deg 5 -> нет "u->-u" симметрии
        if d == 6 and g(-u) == g: cands.append(('u', g))
        FR = FractionField(Ru)
        if d in (5, 6):
            val = g(FR((1+u)/(1-u))) * FR((1-u)**6)
            gx = Ru(val.numerator()/val.denominator())
            if gx(-u) == gx: cands.append(('(u-1)/(u+1)', gx))
        if cands:
            coord, g6 = cands[0]
            cs = g6.list() + [0]*(7-len(g6.list()))
            g3 = Ru([cs[0], cs[2], cs[4], cs[6]])
            rec['bielliptic_coord'] = coord
            rec['E_plus'] = ell(g3); rec['E_minus'] = ell(u*g3)
            a = rec['E_plus'].get('pari', [None, rec['E_plus'].get('eclib_bound')])
            b = rec['E_minus'].get('pari', [None, rec['E_minus'].get('eclib_bound')])
            rec['rank_Jac'] = [None if a[0] is None or b[0] is None else a[0]+b[0],
                               None if a[1] is None or b[1] is None else a[1]+b[1]]
        else:
            # эмпирика: расщепляется ли L-многочлен при хороших p
            red = 0; tot = 0
            for p in primes(5, 150):
                try:
                    if ZZ(g.discriminant()) % p == 0: continue
                    Lp = HyperellipticCurve(g.change_ring(GF(p))).frobenius_polynomial()
                    tot += 1; red += 0 if Lp.is_irreducible() else 1
                except Exception: pass
            rec['Lpoly_reducible'] = [int(red), int(tot)]
            rec['rank_Jac'] = 'нужен RankBounds (Magma): нет расщепления в u / (u-1)/(u+1)'
    else:
        rec['rank_Jac'] = 'род > 2, не рассматривается'
    rows.append(rec); print(json.dumps(rec, ensure_ascii=False), flush=True)

rows = []
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen(); FRx = FractionField(Rx)
subsets = [c for k in (4, 5) for c in combinations(sorted(forms), k)]
for names in subsets:
    f = prod(forms[nm] for nm in names)
    d = f.degree()
    even_t = (f(-t) == f)
    val = f(FRx((1+x)/(1-x))) * FRx((1-x)**d)
    fx = Rx(val.numerator()/val.denominator())
    even_x = (fx(-x) == fx)
    print("#", names, "deg", d, "even_t", even_t, "even_x", even_x, flush=True)
    if not (even_t or even_x):
        rows.append({'names': list(names), 'even': False,
                     'note': 'не чётна ни в t, ни в x'}); continue
    g_even = f if even_t else fx
    coord = 't' if even_t else 'x'
    Ru = PolynomialRing(QQ, 'u'); u = Ru.gen()
    cs = g_even.list() + [0]*(d+1-len(g_even.list()))
    gk = Ru([cs[2*i] for i in range(d//2+1)])
    assert sum(gk.list()[i]*g_even.parent().gen()**(2*i) for i in range(len(gk.list()))) == g_even
    analyse(squarefree_scale(gk),    f"{'*'.join(names)} | coord {coord} | y^2=g(u)",   rows)
    analyse(squarefree_scale(u*gk),  f"{'*'.join(names)} | coord {coord} | y^2=u*g(u)", rows)

open('/home/kep/magicKube/bridge/subcovers_15_8_claude.json', 'w').write(json.dumps(rows, indent=1, ensure_ascii=False))
print("saved")
