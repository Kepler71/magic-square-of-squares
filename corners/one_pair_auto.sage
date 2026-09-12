import functools, time
print = functools.partial(print, flush=True)
b, h = [int(v) for v in sys.argv[1:3]]
sys.argv = ['x']
load('bremner_type.sage')
t0 = time.time()
try:
    res = search_pair(b, h, N='auto', max_rank=4)
    if len(res) == 3:
        C0, best, stats = res
    else:
        C0, best = res; stats = None
    tag = "□" if C0.is_square() else "non□"
    line = f"({b},{h}) C0={C0} {tag} time={time.time()-t0:.0f}s hits={len(best)}"
    if stats is not None:
        ok = sum(1 for s_ in stats if s_[0] == 'ok'); err = sum(1 for s_ in stats if s_[0] == 'error'); skip = sum(1 for s_ in stats if s_[0] == 'skip')
        line += f" bases ok={ok} err={err} skip={skip} detail={[(s_[0], s_[1]) for s_ in stats]}"
    for n_, f_, fl in best:
        line += f"\n   HIT {n_} squares: f = {f_} flags={fl}"
    print(line)
except Exception as ex:
    print(f"({b},{h}) ERROR {ex}")
