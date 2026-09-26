# -*- coding: utf-8 -*-
"""Случайные/массовые данные для теорем B и D: все конфигурации «три прямые с общим центром a»
(a ± u_j — ненулевые квадраты целых, a ≤ A — любое, квадратное или нет) со связью
  тип 1: u_i + u_j = u_k        (n = (1,1,-1), N_L = 3)   — как тройки {b,c,b±c} квадрата;
  тип 2: u_i + u_j = 2 u_k      (n = (1,1,-2), N_L = 4)   — как {b+c,b−c,b};
  тип 3: u_i + 2 u_j = u_k      (n = (1,2,-1), N_L = 4)   — как {b−c,c,b+c} (семилинейные квадраты).
Для каждой: лемма (min t достигается ≥2 раз), ранг <w_u> ≥ 2 (теорема B), неравенство теоремы D.
Отрицательный контроль: тройки линий БЕЗ связи — как часто «неравенство D» нарушалось бы.
Запуск: python3 search_twisted.py A"""
import sys, json, time, random
from math import isqrt
from collections import defaultdict
from circle_lib import *

A = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
t0 = time.time()
groups = defaultdict(list)
# alpha < beta, одной чётности, a = (alpha^2+beta^2)/2 <= A
bmax = isqrt(2 * A)
for be in range(2, bmax + 1):
    b2 = be * be
    start = 2 - (be % 2) if be % 2 == 0 else 1
    for al in range(1 if be % 2 else 2, be, 2):
        s = al * al + b2
        if s > 2 * A:
            break
        groups[s // 2].append((b2 - al * al) // 2)
print('пар:', sum(len(v) for v in groups.values()), 'центров:', len(groups), f'{time.time()-t0:.1f}s', flush=True)

found = []
for a, us in groups.items():
    if len(us) < 3:
        continue
    S = set(us)
    L = sorted(us)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            x, y = L[i], L[j]
            if x + y in S:
                found.append((a, 1, (x, y, x + y), (1, 1, -1)))
            if (x + y) % 2 == 0 and (x + y) // 2 in S:
                found.append((a, 2, (x, y, (x + y) // 2), (1, 1, -2)))
            if (y - x) % 2 == 0 and (y - x) // 2 in S and (y - x) // 2 not in (x, y):
                found.append((a, 3, ((y - x) // 2, x, y), (2, 1, -1)))   # 2*z + x - y = 0
            if x + 2 * y in S:
                found.append((a, 3, (y, x, x + 2 * y), (2, 1, -1)))
            if y + 2 * x in S:
                found.append((a, 3, (x, y, y + 2 * x), (2, 1, -1)))
# уникальность
found = list({(a, t, tuple(sorted(tr)), ns): (a, t, tr, ns) for a, t, tr, ns in found}.values())
print('конфигураций со связью:', len(found), f'{time.time()-t0:.1f}s', flush=True)

stats = dict(A=A, n=len(found), by_type={1: 0, 2: 0, 3: 0}, square_center=[],
             lemma_fail=[], rank_lt2=[], D_fail=[], max_ratio='0', max_ratio_example=None,
             nprimes_hist={}, rank_hist={})
maxr = Fraction(0)
for a, t, tr, ns in found:
    tr = list(tr); ns = list(ns)
    assert sum(n * u for n, u in zip(ns, tr)) == 0, (a, tr, ns)
    NL = sum(abs(n) for n in ns)
    stats['by_type'][t] += 1
    if is_sq(a):
        stats['square_center'].append((a, t, tr))
    ok1, bad = min_twice_ok(a, tr)
    if not ok1:
        stats['lemma_fail'].append((a, tr, bad))
    rk, M = circle_rank(a, tr)
    stats['rank_hist'][rk] = stats['rank_hist'].get(rk, 0) + 1
    if rk < 2:
        stats['rank_lt2'].append((a, tr, M))
    ok2, cls = balance_check(a, tr, NL)
    if not ok2:
        stats['D_fail'].append((a, tr, [(s, ps, Pi, str(r)) for s, ps, Pi, r in cls]))
    for s, ps, Pi, r in cls:
        if r > maxr:
            maxr = r
            stats['max_ratio'] = str(r)
            stats['max_ratio_example'] = dict(a=a, type=t, us=tr, cls=ps, Pi=Pi, factor=str(factorint(a)))
    k = len(split_primes(a))
    stats['nprimes_hist'][k] = stats['nprimes_hist'].get(k, 0) + 1

# отрицательный контроль: тройки без связи, все линии чистые в простом p и одинаково ориентированы
random.seed(1)
neg = dict(triples=0, would_fail=0, examples=[])
cand = [a for a, us in groups.items() if len(us) >= 3]
random.shuffle(cand)
for a in cand[:20000]:
    us = sorted(groups[a])
    tr = random.sample(us, 3)
    x, y, z = sorted(tr)
    if x + y == z or x + y == 2 * z or 2 * x + y == z or x + 2 * y == z or (y - x) == 2 * z or x + z == 2 * y:
        continue
    neg['triples'] += 1
    ok2, cls = balance_check(a, tr, 4)
    if not ok2:
        neg['would_fail'] += 1
        if len(neg['examples']) < 5:
            neg['examples'].append((a, tr, str(factorint(a)), [(s, ps, Pi, str(r)) for s, ps, Pi, r in cls]))
stats['negative_control'] = neg
stats['time_s'] = round(time.time() - t0, 1)
print(json.dumps(stats, indent=1, default=str)[:4000])
json.dump(stats, open(f'search_twisted_A{A}.json', 'w'), indent=1, default=str)
