# Fable, 14.09.2026. Отрицательный контроль: квадрат Бремнера–Саллоуса (центр 425^2, 7 квадратов)
# даёт наклон 247/825 и z0 != 0, при котором шесть из восьми клеток — квадраты. Любая кривая рода 1
# из 3–4 таких клеток имеет невырожденную точку с z = z0; метод обязан НЕ закрыть: либо ранг > 0,
# либо (если вдруг ранг 0) точка z0 должна оказаться в списке и дать не «клетка 0».
# Второй контроль: тройки из g1census.jsonl, где сам Claude нашёл нижнюю границу ранга 1 —
# мой L_ratio должен дать 0 (или mwrank >= 1) и не закрыть.
import sys, json, itertools, time
sys.path.insert(0, '/home/kep/magicKube/fable_census_review')
from common import *
from sage.all import QQ, ZZ, EllipticCurve, pari, gcd
import g1_review as g1
import six_review as six

c0 = ZZ(425)**2
sq = [[373**2, 289**2, 565**2], [360721, c0, 23**2], [205**2, 527**2, 222121]]
p, q = ZZ(138600), ZZ(41496)
k = q / p; r, s = k.numerator(), k.denominator()
print('наклон', r, '/', s, 'коэффициенты', coeffs(r, s))
z0 = p / (c0 * s)           # sz0 = p/c0
cells = nine_cells(r, s, z0)
print('девять нормированных клеток при z0:', [str(x) for x in cells])
print('квадраты:', [QQ(x).is_square() for x in cells])
print('сверка с квадратом Бремнера:', sorted(QQ(x) * c0 for x in cells) == sorted(QQ(x) for row in sq for x in row))
print('classify:', classify(r, s, z0))
# какие c дают квадратную клетку 1 + c z0
good = [c for c in [e * v for v in coeffs(r, s) for e in (1, -1)] if QQ(1 + c * z0).is_square()]
print('квадратные клетки 1 + c z0 при c =', good)

print('\n=== контроль рода 1: все тройки и четвёрки квадратных клеток ===')
nbad = 0
for n in (3, 4):
    for cs in itertools.combinations(good, n):
        t = time.time()
        rec = g1.check_one('%d/%d' % (r, s), list(cs), H=10**5, tL=600, tM=120)
        found = str(z0) in rec['z_values']
        print(cs, 'T=%d L=%s mw=%s rank0=%s pts=%d/%d z0 в списке=%s closed=%s %.0fs' % (
            rec['torsion_order'], rec['L_ratio'], rec['mwrank_bound'], rec['rank0'], rec['n_points_found'],
            rec['torsion_order'], found, rec['closed'], time.time() - t), flush=True)
        if rec['closed']:
            nbad += 1
            print('  !!! ЛОЖНОЕ ЗАКРЫТИЕ', rec)
print('ложных закрытий рода 1:', nbad)

print('\n=== контроль рода 1 на тройках Claude с нижней границей ранга >= 1 ===')
rows = [json.loads(l) for l in open(g1.SRC) if l.strip()]
import random
random.seed(1)
cand = []
for d in rows:
    for t in d.get('tried', []):
        if t[1] >= 1:
            cand.append((d['slope'], t[0], t[1], t[2]))
print('троек/четвёрок с нижней границей >= 1 в файле Claude:', len(cand))
nbad = 0
for slope, cs, lo, hi in random.sample(cand, 12):
    rec = g1.check_one(slope, cs, H=10**4, tL=600, tM=120)
    print(slope, cs, 'Claude (%d,%d)' % (lo, hi), 'L=%s mw=%s rank0=%s closed=%s' % (rec['L_ratio'], rec['mwrank_bound'], rec['rank0'], rec['closed']), flush=True)
    if rec['closed']:
        nbad += 1
print('ложных закрытий:', nbad)

print('\n=== контроль шести клеток на тройках Claude, где E+ и E- имели ранг >= 1 ===')
rows = [json.loads(l) for l in open(six.SRC) if l.strip()]
cand = []
for d in rows:
    for abc, kind, lohi in d.get('info', []):
        lo = int(lohi.strip('()').split(',')[0])
        if lo >= 1:
            cand.append((d['slope'], abc, kind, lohi))
print('кандидатов:', len(cand))
nbad = 0
for slope, abc, kind, lohi in random.sample(cand, 12):
    rec = six.check_one(slope, abc, kind, tL=600, tM=120)
    print(slope, abc, kind, 'Claude', lohi, 'L=%s mw=%s rank0=%s closed=%s' % (rec['L_ratio'], rec['mwrank_bound'], rec['rank0'], rec['closed']), flush=True)
    if rec['closed']:
        nbad += 1
print('ложных закрытий:', nbad)
