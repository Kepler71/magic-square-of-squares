# Fable, 14.09.2026. Независимая перепроверка закрытий этапа «шесть клеток»
# (census_six_cells/census.jsonl — PARI ellrank; kolyvagin.jsonl — L(E,1) != 0).
# Логика (моя): если 1 +- a z, 1 +- b z, 1 +- c z — квадраты, то с u = z^2 квадраты и
#   (1 - a^2 u)(1 - b^2 u)(1 - c^2 u) = y^2           — точка на E+ (u, y),
#   u (1 - a^2 u)(1 - b^2 u)(1 - c^2 u) = (y z)^2     — точка на E- (u, yz).
# Мои модели: E+: Y^2 = (X + b^2c^2)(X + a^2c^2)(X + a^2b^2), u = -X/(abc)^2;
#             E-: Y^2 = (X - a^2)(X - b^2)(X - c^2),         u = 1/X   (X=0 <-> u=oo, X=oo <-> u=0).
# Ранг 0 — L(E,1) != 0 (Sage L_ratio) и/или mwrank (2-спуск со вторым спуском); затем все точки —
# кручение (Sage), их u должны быть квадратами рациональных чисел, z = +-sqrt(u), девять клеток.
import sys, json, random, time
sys.path.insert(0, '/home/kep/magicKube/fable_census_review')
from common import *
from sage.all import QQ, ZZ, EllipticCurve, gcd, sqrt

SRC = '/home/kep/magicKube/census_six_cells/census.jsonl'
SRC_K = '/home/kep/magicKube/census_six_cells/kolyvagin.jsonl'
OUT = sys.argv[3] if len(sys.argv) > 3 else '/home/kep/magicKube/fable_census_review/six_review.jsonl'


def model(a, b, c, kind):
    a, b, c = ZZ(a), ZZ(b), ZZ(c)
    if kind == '+':
        E = EllipticCurve(QQ, [0, (b*c)**2 + (a*c)**2 + (a*b)**2, 0,
                               (b*c)**2*(a*c)**2 + (b*c)**2*(a*b)**2 + (a*c)**2*(a*b)**2,
                               (a*b*c)**4])
        u_of = lambda P: (None if P.is_zero() else -P[0] / (a*b*c)**2)
    else:
        E = EllipticCurve(QQ, [0, -(a**2 + b**2 + c**2), 0,
                               a**2*b**2 + a**2*c**2 + b**2*c**2, -(a*b*c)**2])
        u_of = lambda P: (QQ(0) if P.is_zero() else (None if P[0] == 0 else 1 / P[0]))
    return E, u_of


def selfcheck_model(a, b, c, kind):
    """контроль модели: случайное u с точкой на (E+/E-)-уравнении переносится на модель"""
    import random as rnd
    a, b, c = ZZ(a), ZZ(b), ZZ(c)
    E, u_of = model(a, b, c, kind)
    # проверяем формулу: для X на E+ значение f(u) = (1-a^2u)(1-b^2u)(1-c^2u) должно быть Y^2/(abc)^4
    X = QQ(rnd.randint(-10**6, 10**6)) / rnd.randint(1, 10**3)
    if kind == '+':
        u = -X / (a*b*c)**2
        lhs = (1 - a**2*u)*(1 - b**2*u)*(1 - c**2*u)
        rhs = (X + (b*c)**2)*(X + (a*c)**2)*(X + (a*b)**2) / (a*b*c)**4
    else:
        u = 1 / X
        lhs = u*(1 - a**2*u)*(1 - b**2*u)*(1 - c**2*u)
        rhs = (X - a**2)*(X - b**2)*(X - c**2) / X**4
    return lhs == rhs


def check_one(slope, abc, kind, tL=900, tM=300):
    r, s = map(int, slope.split('/'))
    rec = {'slope': slope, 'abc': list(abc), 'kind': kind}
    rec['cells_ok'] = set(abc) <= set(coeffs(r, s)) and len(set(abc)) == 3
    a, b, c = abc
    rec['model_selfcheck'] = selfcheck_model(a, b, c, kind)
    E, u_of = model(a, b, c, kind)
    rec['conductor'] = str(E.conductor())
    tors = E.torsion_subgroup()
    rec['torsion_struct'] = str(tors.invariants())
    rec.update(prove_rank0(E, tL, tM))
    rec['rank0'] = rec['rank0_by_L'] or rec['rank0_by_mwrank']
    us = []
    for P in tors.points():
        u = u_of(P)
        if u is None:
            continue
        us.append(u)
    us = sorted(set(us))
    rec['u_values'] = [str(u) for u in us]
    zs = []
    for u in us:
        if u >= 0 and QQ(u).is_square():
            w = sqrt(QQ(u))
            zs += [w, -w] if w != 0 else [QQ(0)]
    zs = sorted(set(zs))
    rec['z_values'] = [str(z) for z in zs]
    rec['z_classes'] = [classify(r, s, z) for z in zs]
    rec['alert'] = any(x.startswith('ALERT') for x in rec['z_classes'])
    rec['closed'] = rec['cells_ok'] and rec['model_selfcheck'] and rec['rank0'] and not rec['alert']
    return rec


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sample'
    random.seed(20260914)
    if mode == 'kol':
        rows = [json.loads(l) for l in open(SRC_K) if l.strip()]
        todo = [(d['slope'], d['by']['abc'], d['by']['kind'], 'L_ratio=' + d['by']['L_ratio']) for d in rows if d.get('closed')]
    else:
        rows = [json.loads(l) for l in open(SRC) if l.strip()]
        closed = [d for d in rows if d.get('closed')]
        print('в файле: закрытых', len(closed), flush=True)
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 45
        random.shuffle(closed)
        todo = [(d['slope'], d['by'][0], d['by'][1], 'ellrank ' + json.dumps(d['info'])) for d in closed[:n]]
    done = set()
    try:
        for l in open(OUT):
            d = json.loads(l); done.add((d['slope'], tuple(d['abc']), d['kind']))
    except FileNotFoundError:
        pass
    for i, (slope, abc, kind, meth) in enumerate(todo):
        if (slope, tuple(abc), kind) in done:
            continue
        t0 = time.time()
        rec = check_one(slope, abc, kind)
        rec['claude_method'] = meth
        rec['sec'] = round(time.time() - t0, 1)
        with open(OUT, 'a') as f:
            f.write(json.dumps(rec, ensure_ascii=False) + '\n')
        print('%4d/%d %s %s%s N=%s tors=%s L=%s(%ss) mw=%s z=%s closed=%s %.0fs' % (
            i + 1, len(todo), slope, abc, kind, rec['conductor'], rec['torsion_struct'], rec['L_ratio'], rec['L_time'],
            rec['mwrank_bound'], rec['z_values'], rec['closed'], rec['sec']), flush=True)
