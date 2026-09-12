# Поиск рациональных точек на C_J2 (обе закрутки) и ранг их образа в E1(k).
import functools, time
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
Rt2.<T2> = QQ[]
f = T2*(T2^2 - 1)*(109*T2 - 229)*(229*T2 - 109)
H = int(sys.argv[1]) if len(sys.argv) > 1 else 10^4
for d3t, D, dsq in [((5, 13, 65), 5*13*65*65, 65), ((6, 26, 39), 6*26*39*39, 39)]:
    # D W^2 = f  <=>  (D W)^2 = D f ; D = dsq * square
    t0 = time.time()
    pts = pari.hyperellratpoints(pari(dsq * f), H)
    tvals = sorted(set(QQ(P[0]) for P in pts), key=lambda x: max(abs(x.numerator()), x.denominator()))
    print(f"class {d3t}: C_J2 points with height <= {H}: {len(pts)} (t-values: {len(tvals)}) ({time.time()-t0:.1f}s)")
    print("   t =", tvals[:40])
