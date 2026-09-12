import functools, time
print = functools.partial(print, flush=True)
sys.argv = ['x']
load('corner_search.sage')
secs = [(7,1,5), (17,7,13), (23,7,17), (31,17,25), (41,1,29), (47,23,37), (49,31,41), (73,17,53), (71,49,61), (89,23,65), (79,47,65), (103,7,73), (113,41,85), (97,71,85), (119,41,89), (137,7,97)]
for sec in secs:
    t0 = time.time()
    res = search_section(*sec, N=6)
    if res is None: print(sec, "no base point"); continue
    base, results = res
    summ = []
    hits = {}
    for (i, rk, found) in results:
        summ.append(f"corner{i+1}: rank {rk}")
        if isinstance(found, dict):
            for t1, k in found.items(): hits[t1] = max(hits.get(t1, 0), k)
    k2 = sum(1 for v in hits.values() if v == 2); k3 = sum(1 for v in hits.values() if v == 3); k4 = sum(1 for v in hits.values() if v == 4)
    print(f"{sec} t0={base}: {', '.join(summ)} | t with 2 corners: {k2}, 3 corners: {k3}, 4 corners: {k4}  ({time.time()-t0:.0f}s)")
    for t1, v in hits.items():
        if v >= 3: print(f"    *** {v} corners at t = {t1}")
        elif v == 2 and max(abs(t1.numerator()), t1.denominator()) < 10^30: print(f"    2 corners at t = {t1}")
