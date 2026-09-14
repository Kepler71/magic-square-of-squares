# Fable, 14.09.2026. Независимый контроль родов компонент границы через Sage (гиперэллиптические факторы, Кани–Розен)
# и прямой счёт рода кривой в аффинном пространстве для центра (степень 8) и ребра (степень 16).
from sage.all import *
import itertools, json
out = {}
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
classes = {'M': [t, 1+t, 1-t], 'edge a+b+c': [t, 1-t, 2+t, 1+t], 'corner a+b': [t, 1+t, 1-t, 2+t, 2-t]}
# Кани–Розен для (Z/2)^k: g(C) = сумма родов 2^k-1 гиперэллиптических факторов y^2 = произведение подмножества
def hyp_genus(f):
    d = f.degree()
    return (d - 1) // 2 if d % 2 else (d - 2) // 2
for name, cl in classes.items():
    tot = 0
    for r in range(1, len(cl)+1):
        for sub in itertools.combinations(cl, r):
            tot += hyp_genus(prod(sub))
    out['kani_rosen_' + name] = tot
# прямой счёт: кривая в A^{k+1}
def direct(cl):
    k = len(cl)
    A = AffineSpace(QQ, k+1, names=['t'] + ['y%d' % i for i in range(k)])
    g = A.gens(); tt = g[0]
    eqs = [g[i+1]**2 - cl[i](tt) for i in range(k)]
    C = Curve(eqs, A)
    return int(C.genus())
for name in ['M', 'edge a+b+c']:
    try:
        out['direct_' + name] = direct(classes[name])
    except Exception as e:
        out['direct_' + name] = 'error: ' + repr(e)[:200]
json.dump(out, open('/home/kep/magicKube/fable_integral/genus_check.json', 'w'), indent=1, default=str)
print(out)
