# Независимая проверка RESULT_29_WALLS (Claude): для наклона r/s перебираю ВСЕ тройки
# различных a<b<c из {r, s-r, s, s+r} и оба фактора E+ и E-; ищу фактор с доказанным рангом 0,
# у которого ни одна рациональная точка не даёт допустимого z != 0.
import json, sys
from sage.all import *
slopes = "19/60 41/60 38/63 29/65 53/85 29/94 33/97 80/97 63/101 73/105 41/125 100/129 45/137 92/141 14/143 21/143 43/145 23/146 106/149 26/153 17/168 64/171 80/189 101/126 46/147 39/149 41/166 48/175 29/180 3/86 79/110 11/142 48/163".split()
R = PolynomialRing(QQ,'X'); X = R.gen()
def verdict(E_poly, kind, amax):
    # E: y^2 = E_poly(t), кубика; ранг и все точки ранга 0
    E = EllipticCurve(HyperellipticCurve(E_poly).change_ring(QQ)) if False else None
    f = E_poly; lc = f.leading_coefficient()
    # y^2 = lc*(t-t1)(t-t2)(t-t3): приводим к Вейерштрассу t = T/lc, y = Y/lc
    g = (f(X/lc)*lc**2).monic() if lc != 1 else f
    co = g.list()
    Em = EllipticCurve([0, co[2], 0, co[1], co[0]])
    lo, hi = Em.pari_curve().ellrank()[:2]
    if hi != 0:
        return None, (int(lo), int(hi))
    ts = []
    for P in Em.torsion_points():
        if P.is_zero(): ts.append('inf'); continue
        ts.append(P[0]/lc)
    ok = True
    for t in ts:
        if t == 'inf':
            continue  # E+: u=inf не даёт z; E-: x=inf -> z=0
        if kind == '+':
            if 0 < t < 1/amax**2 and t.is_square(): ok = False
        else:
            if t > amax**2 and t.is_square(): ok = False
    return ok, (0, 0, len(ts))
out = {}
for sl in slopes:
    r, s = map(int, sl.split('/'))
    vals = sorted(set(v for v in [r, s-r, s, s+r] if v > 0))
    found = []; ranks = []
    for a, b, c in Combinations(vals, 3):
        up = (1-a**2*X)*(1-b**2*X)*(1-c**2*X)
        um = (X-a**2)*(X-b**2)*(X-c**2)
        for kind, pol in (('+', up), ('-', um)):
            ok, info = verdict(pol, kind, max(a,b,c))
            ranks.append(((a,b,c), kind, info))
            if ok: found.append(((a,b,c), kind))
    out[sl] = dict(closed=bool(found), by=[str(x) for x in found], ranks=[str(x) for x in ranks])
    print(sl, 'ЗАКРЫТ' if found else 'не закрыт', found[:2], flush=True)
json.dump(out, open('walls_check.json','w'), ensure_ascii=False, indent=1)
