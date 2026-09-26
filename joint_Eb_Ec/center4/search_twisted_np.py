# -*- coding: utf-8 -*-
"""То же, что search_twisted.py, но блоками по центру a (numpy) — для A до 10^9.
Проверяет: лемма, ранг <w_u> (теорема B, и наблюдение «ранг 3»), неравенство теоремы D.
Запуск: python3 search_twisted_np.py A BLOCK"""
import sys, json, time
import numpy as np
from math import isqrt
from circle_lib import *

A = int(sys.argv[1]); BLOCK = int(sys.argv[2]) if len(sys.argv) > 2 else 10**8
t0 = time.time()
stats = dict(A=A, n=0, by_type={1: 0, 2: 0, 3: 0}, square_center=[], lemma_fail=[], rank_lt2=[],
             rank_hist={}, D_fail=[], max_ratio='0', max_ratio_example=None, nprimes_hist={})
maxr = Fraction(0)
lo = 1
while lo <= A:
    hi = min(A, lo + BLOCK - 1)
    # пары alpha<beta одной чётности, 2*lo <= alpha^2+beta^2 <= 2*hi
    As, Us = [], []
    for be in range(2, isqrt(2 * hi) + 1):
        b2 = be * be
        # alpha^2 >= 2lo - b2, alpha^2 <= 2hi - b2, 1 <= alpha < be, alpha ≡ be (mod 2)
        amin2 = 2 * lo - b2
        amin = 1 if amin2 <= 1 else isqrt(amin2 - 1) + 1
        amax2 = 2 * hi - b2
        if amax2 < 1: continue
        amax = min(be - 1, isqrt(amax2))
        if amin > amax: continue
        if (amin - be) % 2: amin += 1
        if amin > amax: continue
        al = np.arange(amin, amax + 1, 2, dtype=np.int64)
        s = al * al + b2
        As.append(s // 2); Us.append((b2 - al * al) // 2)
    if not As:
        lo = hi + 1; continue
    a_arr = np.concatenate(As); u_arr = np.concatenate(Us)
    del As, Us
    order = np.lexsort((u_arr, a_arr))
    a_arr = a_arr[order]; u_arr = u_arr[order]
    del order
    # границы групп
    brk = np.flatnonzero(np.diff(a_arr)) + 1
    starts = np.concatenate(([0], brk)); ends = np.concatenate((brk, [len(a_arr)]))
    big = np.flatnonzero(ends - starts >= 3)
    for gi in big:
        s0, e0 = starts[gi], ends[gi]
        a = int(a_arr[s0]); L = [int(v) for v in u_arr[s0:e0]]
        S = set(L)
        rels = set()
        for i in range(len(L)):
            x = L[i]
            for j in range(i + 1, len(L)):
                y = L[j]
                if x + y in S: rels.add((1, (x, y, x + y), (1, 1, -1)))
                if (x + y) % 2 == 0 and (x + y) // 2 in S: rels.add((2, (x, y, (x + y) // 2), (1, 1, -2)))
                if (y - x) % 2 == 0 and (y - x) // 2 in S and (y - x) // 2 not in (x, y): rels.add((3, ((y - x) // 2, x, y), (2, 1, -1)))
                if x + 2 * y in S: rels.add((3, (y, x, x + 2 * y), (2, 1, -1)))
                if y + 2 * x in S: rels.add((3, (x, y, y + 2 * x), (2, 1, -1)))
        for t, tr, ns in rels:
            tr = list(tr); ns = list(ns); NL = sum(abs(n) for n in ns)
            stats['n'] += 1; stats['by_type'][t] += 1
            if is_sq(a): stats['square_center'].append((a, t, tr))
            ok1, bad = min_twice_ok(a, tr)
            if not ok1: stats['lemma_fail'].append((a, tr, bad))
            rk, M = circle_rank(a, tr)
            stats['rank_hist'][rk] = stats['rank_hist'].get(rk, 0) + 1
            if rk < 2: stats['rank_lt2'].append((a, tr, M))
            if rk == 2 and len(stats.setdefault('rank2_examples', [])) < 50:
                stats['rank2_examples'].append((a, t, tr, M, str(factorint(a))))
            ok2, cls = balance_check(a, tr, NL)
            if not ok2: stats['D_fail'].append((a, tr))
            for sg, ps, Pi, r in cls:
                if r > maxr:
                    maxr = r; stats['max_ratio'] = str(r)
                    stats['max_ratio_example'] = dict(a=a, type=t, us=tr, cls=ps, Pi=Pi, factor=str(factorint(a)))
            k = len(split_primes(a)); stats['nprimes_hist'][k] = stats['nprimes_hist'].get(k, 0) + 1
    print(f'блок [{lo},{hi}]: пар {len(a_arr)}, найдено всего {stats["n"]}, ранги {stats["rank_hist"]}, {time.time()-t0:.0f}s', flush=True)
    del a_arr, u_arr
    lo = hi + 1
stats['time_s'] = round(time.time() - t0, 1)
json.dump(stats, open(f'search_twisted_np_A{A}.json', 'w'), indent=1, default=str)
print({k: stats[k] for k in ('A', 'n', 'by_type', 'square_center', 'lemma_fail', 'rank_lt2', 'rank_hist', 'D_fail', 'max_ratio', 'max_ratio_example', 'nprimes_hist', 'time_s')})
