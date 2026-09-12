load('/home/kep/magicKube/family/independent_chabauty/setup.sage')
import sys
b, h, n, D, u0, PMAX = sys.argv[1:7]
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D)); E = S.E; k = S.k
G = S.point_from_u(QQ(u0)); tors = E.torsion_points()
proved = {}
for pr in primes(3, ZZ(PMAX)):
    for P in k.primes_above(pr):
        if E.has_bad_reduction(P): continue
        F = P.residue_field(); Ered = E.reduction(P)
        rp = lambda pt: Ered(0) if pt.is_zero() else Ered([F(pt[0]), F(pt[1])])
        m = Ered.cardinality(); Ab = Ered.abelian_group()
        invs = [ZZ(d) for d in Ab.invariants()]
        for l in ZZ(m).prime_divisors():
            if l == pr or proved.get(l): continue
            if all(not all(ZZ(c) % gcd(l, d) == 0 for c, d in zip(Ab.discrete_log(rp(G) - rp(T)), invs)) for T in tors):
                proved[l] = (pr, P.norm(), m)
                print(f"  l={l}: P|{pr}, N(P)={P.norm()}, |E(F_P)|={m}", flush=True)
print("proved:", sorted(proved))
print("NOT proved up to 200:", [l for l in primes(2, 200) if l not in proved])
