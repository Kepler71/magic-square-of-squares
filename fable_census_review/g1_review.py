# Fable, 14.09.2026. Независимая перепроверка закрытий этапа «род 1» (census_genus1/g1census.jsonl).
# Модель: якобиан кривой y^2 = prod(1 + c_i z) через PARI ellfromeqn (не перенос корня на бесконечность);
# ранг 0 — через L(E,1) != 0 (Sage L_ratio) и/или 2-спуск mwrank; точки кривой — прямой перебор
# hyperellratpoints, полнота — сравнение числа найденных точек с порядком кручения якобиана
# (кривая с рациональной точкой изоморфна своему якобиану, #C(Q) = #E(Q)_tors при ранге 0).
import sys, json, random, time
sys.path.insert(0, '/home/kep/magicKube/fable_census_review')
from common import *
from sage.all import QQ, ZZ, EllipticCurve, pari, gcd

SRC = '/home/kep/magicKube/census_genus1/g1census.jsonl'
OUT = '/home/kep/magicKube/fable_census_review/g1_review.jsonl'


def check_one(slope, cs, H=10**5, tL=900, tM=300):
    r, s = map(int, slope.split('/'))
    allowed = set()
    for c in coeffs(r, s):
        allowed |= {c, -c}
    rec = {'slope': slope, 'cells': cs}
    rec['cells_ok'] = all(c in allowed for c in cs) and len(set(cs)) == len(cs)
    P = R.prod([1 + c * zz for c in cs])
    y = PolynomialRing(QQ, 'y').gen()
    # якобиан по PARI ellfromeqn: f(x,y) = y^2 - P(x)
    x = pari('x'); yp = pari('y')
    Ppari = pari(P).subst('z', x)
    a = pari.ellfromeqn(yp**2 - Ppari)
    E = EllipticCurve(QQ, [QQ(t) for t in a])
    rec['E_ainv'] = [str(t) for t in E.a_invariants()]
    rec['conductor'] = str(E.conductor())
    T = torsion_order(E)
    rec['torsion_order'] = T
    rec['torsion_struct'] = str(E.torsion_subgroup().invariants())
    rec.update(prove_rank0(E, tL, tM))
    rec['rank0'] = rec['rank0_by_L'] or rec['rank0_by_mwrank']
    # точки
    aff, naff, ninf = points_on_y2_P(P, H)
    if naff + ninf < T:
        aff, naff, ninf = points_on_y2_P(P, 10**7)
        rec['H'] = 10**7
    else:
        rec['H'] = H
    rec['n_points_found'] = naff + ninf
    rec['points_complete'] = (naff + ninf == T)
    rec['z_values'] = [str(z) for z in aff]
    rec['z_classes'] = [classify(r, s, z) for z in aff]
    rec['alert'] = any(c.startswith('ALERT') for c in rec['z_classes'])
    rec['closed'] = rec['cells_ok'] and rec['rank0'] and rec['points_complete'] and not rec['alert']
    return rec


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'sample'
    rows = [json.loads(l) for l in open(SRC) if l.strip()]
    closed = [d for d in rows if d.get('closed')]
    ell = [d for d in closed if d['by'][1] == 'ellrank']
    lr = [d for d in closed if d['by'][1] != 'ellrank']
    print('в файле: закрытых', len(closed), 'ellrank', len(ell), 'иных', len(lr), flush=True)
    random.seed(20260914)
    if mode == 'sample':
        n = int(sys.argv[2]) if len(sys.argv) > 2 else 65
        todo = random.sample(ell, n) + lr
    elif mode == 'one':
        todo = [{'slope': sys.argv[2], 'by': [json.loads(sys.argv[3]), '?']}]
    else:
        todo = closed
    done = set()
    try:
        for l in open(OUT):
            d = json.loads(l); done.add((d['slope'], tuple(d['cells'])))
    except FileNotFoundError:
        pass
    for i, d in enumerate(todo):
        key = (d['slope'], tuple(d['by'][0]))
        if key in done:
            continue
        t0 = time.time()
        rec = check_one(d['slope'], d['by'][0])
        rec['claude_method'] = d['by'][1]
        rec['sec'] = round(time.time() - t0, 1)
        with open(OUT, 'a') as f:
            f.write(json.dumps(rec, ensure_ascii=False) + '\n')
        print('%3d/%d %s %s T=%d L=%s(%ss) mw=%s pts=%d/%d closed=%s %.0fs' % (
            i + 1, len(todo), d['slope'], d['by'][0], rec['torsion_order'], rec['L_ratio'], rec['L_time'],
            rec['mwrank_bound'], rec['n_points_found'], rec['torsion_order'], rec['closed'], rec['sec']), flush=True)
