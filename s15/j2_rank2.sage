import functools, time
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
for d3t, D in [((5, 13, 65), 5*13*65*65), ((6, 26, 39), 6*26*39*39)]:
    print(f"\n===== class {d3t}")
    cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))
    a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
    E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2]).global_minimal_model(semi_global=True)
    print("minimal-ish model:", E1.ainvs())
    print("torsion:", E1.torsion_subgroup().invariants(), " conductor norm:", E1.conductor().norm().factor())
    t0 = time.time()
    try:
        lo, hi, gens = E1.simon_two_descent(lim1=8, lim3=200, limtriv=200, maxprob=40, limbigprime=60)
        print(f"Simon (big limits): rank in [{lo}, {hi}], {len(gens)} gens ({time.time()-t0:.0f}s)")
        G = [P for P in gens if P.order() == oo]
        if G:
            M = E1.height_pairing_matrix(G)
            print("   height pairing det of found points:", M.det(), " rank of height matrix:", M.rank())
    except Exception as ex:
        print("Simon failed:", ex)
    # аналитический ранг: L(E1/k, s) через PARI lfun (численно, не доказательство)
    t0 = time.time()
    try:
        nf = pari(k.defining_polynomial()).nfinit()
        Ep = pari.ellinit([pari(str(x.polynomial()).replace('r15', 'y')) if False else pari(x.polynomial()('y')) for x in E1.ainvs()], nf)
        L = pari.lfuncreate(Ep)
        print("analytic order of vanishing at s=1:", pari.lfunorderzero(L), f"({time.time()-t0:.0f}s)")
    except Exception as ex:
        print("PARI lfun failed:", ex)
