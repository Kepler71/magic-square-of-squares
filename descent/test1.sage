import functools, time
print = functools.partial(print, flush=True)
load('g2descent.sage')
R.<t> = QQ[]
tests = [
  ("s=1/5 K0 (odd model)  [Magma: 0,1]", 78*t*(t^2 + 3600)*(t^2 + 28561)),
  ("(23,7,17) C0          [Magma: 0,5]", 34*t*(t^2-1)*(169*t-409)*(409*t-169)),
  ("s=1/5 C0              [Magma E1: 1,3]", 65*t*(t^2-1)*(109*t-229)*(229*t-109)),
]
for name, f in tests:
    t0 = time.time()
    try:
        D = Descent(Rx(f))
        d, rb = D.rank_bound()
        print(f"==> {name}: dim Sel^2 = {d}, rank <= {rb}  ({time.time()-t0:.0f}s)\n")
    except Exception as ex:
        import traceback; traceback.print_exc()
        print(f"==> {name}: FAILED {ex}\n")
