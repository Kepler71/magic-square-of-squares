import functools, time
print = functools.partial(print, flush=True)
exec(preparse(open('j2_rank.sage').read().split("for d3t, D in")[0]))
D = 5*13*65*65
cub = Ru(D * c * (uu^2 + A - 2*beta) * (uu + 2*sb))
a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
E1m = E1.global_minimal_model(semi_global=True)
print("model:", E1m.ainvs())
nf = pari.nfinit(pari('y^2 - 15'))
ai = [pari(str(x.polynomial('y'))) for x in E1m.ainvs()]
Ep = pari.ellinit(ai, nf)
t0 = time.time()
gr = pari.ellglobalred(Ep)
N = pari.idealnorm(nf, gr[0])
print("conductor norm:", N, "=", factor(ZZ(N)), f"({time.time()-t0:.1f}s)")
print("root number w(E1/k):", pari.ellrootno(Ep), " (w=-1 => odd rank)")
print("analytic conductor of L(E1/k,s):", ZZ(N) * 60^2)
t0 = time.time()
L = pari.lfuncreate(Ep)
print("lfuncreate ok", f"({time.time()-t0:.1f}s)")
t0 = time.time()
print("order of vanishing at s=1:", pari.lfunorderzero(L), f"({time.time()-t0:.1f}s)")
