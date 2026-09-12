import functools, time
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
D = 5*13*65*65
cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))
a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2]).global_minimal_model(semi_global=True)
t0 = time.time()
C = E1.isogeny_class()
print("isogeny class size over Q(sqrt15):", len(C), f"({time.time()-t0:.0f}s)")
print("isogeny matrix:\n", C.matrix())
best = 99
for i, E in enumerate(C.curves):
    t0 = time.time()
    try:
        lo, hi, g = E.simon_two_descent()
        best = min(best, hi)
        print(f"  curve {i}: torsion {E.torsion_subgroup().invariants()}  Simon bounds [{lo},{hi}]  ({time.time()-t0:.0f}s)")
    except Exception as ex:
        print(f"  curve {i}: failed {ex}")
print("best upper bound over the isogeny class:", best, " (root number -1 => rank odd)")
