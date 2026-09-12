import functools, time
print = functools.partial(print, flush=True)
src = open('h0_quotient.sage').read()
exec(preparse(src.split('t0 = time.time(); g1 = Curve(G)')[0]))
print("H0 on fixed points of iota:", {(a_, b_): H0(a_, b_) == 0 for a_ in (1, -1) for b_ in (1, -1)})
for p in [101]:
    t0 = time.time()
    print(f"M2 genus of quotient (function field, p={p}):", M2(R(G(k, u)), p), f"({time.time()-t0:.1f}s)")
