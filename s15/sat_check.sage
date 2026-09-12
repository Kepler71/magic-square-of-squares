import functools
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
k2.<r> = QuadraticField(15)
E = EllipticCurve(k2, [0, 0, 0, -100326681, 98571779592*r])
G = E(-676*r - 8619, 30420*r + 182520)
print("G on E, order:", G.order(), " height:", G.height())
res = E.saturation([G], verbose=True)
print("RESULT:", res[1], res[2])
for pp in [2, 3, 5, 7, 11, 13]:
    print(f"  one_prime={pp}:", E.saturation([G], one_prime=pp)[1])
print("lower height bound (Cremona–Siksek-type):", E.height_function().min_gr(0.1, 5) if hasattr(E, 'height_function') else 'n/a')
