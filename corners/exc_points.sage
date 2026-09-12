# Исключительные точки отображения E -> f (Codex, RESULT_..._LLL_SQUARE_CENTERS.md): O и точки с y = 0,
# на которых BaseCurve.f_of_point возвращает None. Предельные значения f:
#   O:             m -> m0                      => f = num(m0)/den(m0)
#   y = 0, N != 0: m -> oo                      => f = f0
#   y = 0, N = 0:  m = m0 - 2q F_y/F_x (касательная; F_x = 0 => f0),  N = 2q(x + c) - d^2/(2q)
# Проход по всем парам и 6 базам: такие f оцениваются тем же score().
import sys, functools
print = functools.partial(print, flush=True)
args = sys.argv[1:]; sys.argv = ['x']
load('/home/kep/magicKube/corners/bremner_type.sage')
Rx = PolynomialRing(QQ, 'X')

def exceptional_fs(bc):
    E = bc.E; a1, a2, a3, a4, a6 = E.ainvs(); q, c, d = bc.qq, bc.cc, bc.dd
    out = []
    m = bc.m0
    if bc.den(m) != 0: out.append(('O', bc.num(m)/bc.den(m)))
    X = Rx.gen()
    for x0, _ in (X^3 + a2*X^2 + a4*X + a6).roots():
        N = 2*q*(x0 + c) - d^2/(2*q)
        if N != 0: out.append((f'({x0},0)', bc.f0)); continue
        Fx = -(3*x0^2 + 2*a2*x0 + a4); Fy = a1*x0 + a3
        if Fx == 0: out.append((f'({x0},0)', bc.f0)); continue
        m = bc.m0 - 2*q*Fy/Fx
        if bc.den(m) != 0: out.append((f'({x0},0)', bc.num(m)/bc.den(m)))
    return out

def exc_pair(b, h):
    C0, cs = conds(b, h)
    res = []; err = 0
    for c1, c2, f0 in [('TL','BR',h),('TL','D',h),('BR','D',h),('BL','TR',b),('BL','D',b),('TR','D',b)]:
        try:
            bc = BaseCurve(cs[c1], cs[c2], QQ(f0))
        except Exception:
            err += 1; continue
        for where, f in exceptional_fs(bc):
            assert (cs[c1][0] + cs[c1][1]*f^2).is_square() and (cs[c2][0] + cs[c2][1]*f^2).is_square(), (b, h, c1, c2, where, f)
            nsq, distinct, flags, cells = score(b, h, f, C0, cs)
            res.append((f'{c1}/{c2}', where, f, nsq, distinct))
    return C0, res, err

if args and args[0] == 'control':
    C0, res, err = exc_pair(31, 17)
    for r in res: print(r)
else:
    for line in open(args[0]):
        b, h = map(int, line.split())
        C0, res, err = exc_pair(b, h)
        n6 = sum(1 for r in res if r[4] and r[3] == 6)
        hits = [r for r in res if r[4] and r[3] >= 7]
        print(f"({b},{h}) C0={C0} {'□' if C0.is_square() else 'non□'} exc={len(res)} err={err} distinct6={n6} hits={hits}")
