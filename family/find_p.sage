import sys
exec(preparse(open('chabauty.sage').read().split("# --- p-адические вложения")[0].replace("assert kronecker(dk, p) == 1 and all(not E1.has_bad_reduction(P) for P in k.primes_above(p)), \"p must split and be good\"", "pass")))
out = []
for pp in primes(11, 400):
    if kronecker(dk, pp) != 1: continue
    Ps = k.primes_above(pp)
    if any(E1.has_bad_reduction(P) for P in Ps): continue
    try:
        Ns = [E1.reduction(P)(G).order() if False else None for P in Ps]
    except Exception:
        Ns = None
    ords = []
    for P in Ps:
        F = k.residue_field(P); Er = E1.reduction(P)
        Gr = Er([F(G[0]), F(G[1])])
        ords.append(Gr.order())
    out.append((lcm(ords), pp, ords))
print(sorted(out)[:10])
