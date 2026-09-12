import functools, time
print = functools.partial(print, flush=True)
b, h, n, N = [int(v) for v in sys.argv[1:5]]
sys.argv = ['x']
load('corner_search.sage')
t0 = time.time()
res = search_section(b, h, n, N=N, max_rank=4)
if res is None:
    print(f"({b},{h},{n}) nobase")
else:
    base, results = res
    ranks = [rk for (i, rk, f) in results]
    hits = {}
    for (i, rk, found) in results:
        if isinstance(found, dict):
            for t1, k in found.items(): hits[t1] = max(hits.get(t1, 0), k)
    ks = [hits[t1] for t1 in hits]
    print(f"({b},{h},{n}) ranks={ranks} k2={ks.count(2)} k3={ks.count(3)} k4={ks.count(4)} time={time.time()-t0:.0f}s" + "".join(f"\n   HIT k={v} t={t1}" for t1, v in hits.items()))
