#!/usr/bin/env python3
# Сводка по логам прогонов: статистика и НОРМАЛИЗАЦИЯ попаданий (приведение к целым, деление на
# наибольший квадратный общий множитель, канонизация по D8) — чтобы отличить новый пример от
# известного Бремнера–Саллоуса и его симметрий/кратных k^2.
import sys, re, glob
from fractions import Fraction
from math import gcd, isqrt
from configs import d8, CELLS

G = d8()
SALLOWS = [373**2, 289**2, 565**2, 360721, 425**2, 23**2, 205**2, 527**2, 222121]


def normalize(cells):
    fr = [Fraction(x) for x in cells]
    L = 1
    for f in fr:
        L = L*f.denominator//gcd(L, f.denominator)
    ints = [int(f*L) for f in fr]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    ints = [x//g for x in ints]
    # убрать наибольший квадратный делитель, общий для всех клеток
    d = 2
    while d*d <= max(abs(x) for x in ints):
        while all(x % (d*d) == 0 for x in ints):
            ints = [x//(d*d) for x in ints]
        d += 1
    return tuple(min(tuple(v[g_[i]] for i in CELLS) for g_ in G for v in [ints]))


def canon(cells):
    fr = [Fraction(x) for x in cells]
    L = 1
    for f in fr:
        L = L*f.denominator//gcd(L, f.denominator)
    ints = [int(f*L) for f in fr]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    if g:
        ints = [x//g for x in ints]
    d = 2
    while d*d <= max(abs(x) for x in ints):
        if all(x % (d*d) == 0 for x in ints):
            ints = [x//(d*d) for x in ints]
        else:
            d += 1
    return min(tuple(ints[g_[i]] for i in CELLS) for g_ in G)


SAL = canon(SALLOWS)
pat = re.compile(r'HIT (\S+) \((?:P,Q|m,n)\)=\(([-\d]+),([-\d]+)\) [pt]=(\S+) nsq=(\d+) (?:pos=(\w+) )?cells=\[(.*)\]')

stats = {'pairs': 0, 'bases': 0, 'nocurve': 0, 'err': 0, 'timeout': 0, 'inc': 0}
hits = {}
for fn in sys.argv[1:]:
    for line in open(fn):
        if line.startswith('=='):
            stats['pairs'] += 1
        elif ' ERROR ' in line:
            stats['err'] += 1
        elif 'TIMEOUT' in line:
            stats['timeout'] += 1
        elif 'не постро' in line or 'нет рацион' in line:
            stats['nocurve'] += 1
        elif re.match(r'^  [FG]\d ', line):
            stats['bases'] += 1
            if ' INC ' in line:
                stats['inc'] += 1
        m = pat.search(line)
        if m:
            fam, A, B, t, nsq, pos, cells = m.groups()
            cl = [Fraction(x.strip().strip("'")) for x in cells.split(',')]
            key = canon(cl)
            hits.setdefault(key, []).append((fam, A, B, t, int(nsq), pos))

print("Статистика логов:", stats)
print(f"Различных (с точностью до D8 и множителя k^2) попаданий с >= 7 квадратами: {len(hits)}")
for k, v in hits.items():
    known = (k == SAL)
    print(f"  {'ИЗВЕСТНЫЙ (Бремнер–Саллоус)' if known else '*** НОВЫЙ ***'}: {list(k)}")
    print(f"     найден {len(v)} раз(а), например {v[0]}, семейства {sorted(set(x[0] for x in v))}")
