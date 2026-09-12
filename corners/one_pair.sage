import functools, time
print = functools.partial(print, flush=True)
b, h = [int(v) for v in sys.argv[1:3]]
sys.argv = ['x']
load('bremner_type.sage')
t0 = time.time()
try:
    C0, best = search_pair(b, h, N=5, max_rank=4)
    tag = "□" if C0.is_square() else "non□"
    line = f"({b},{h}) C0={C0} {tag} time={time.time()-t0:.0f}s hits={len(best)}"
    for n_, f_, fl in best:
        line += f"\n   HIT {n_} squares: f = {f_} flags={fl}"
    print(line)
except Exception as ex:
    print(f"({b},{h}) ERROR {ex}")
