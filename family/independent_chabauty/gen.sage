load('/home/kep/magicKube/family/independent_chabauty/setup.sage')
import sys
b, h, n, D, u0 = sys.argv[1:6]
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D))
E = S.E
k = S.k
G = S.point_from_u(QQ(u0))
print(f"E = {E.ainvs()}")
print(f"G = {G}, u(G) = {S.u_of(G)}")
tors = E.torsion_points()
print(f"torsion {E.torsion_subgroup().invariants()}: {tors}")
print(f"G has finite order? {G.has_finite_order()};  canonical height h(G) = {G.height()}")

# --- (a) точная проверка l-делимости по модулю кручения
for l in [2, 3, 5, 7, 11, 13]:
    bad = [(T, len((G - T).division_points(l))) for T in tors if (G - T).division_points(l)]
    print(f"   exact  l={l}: {'DIVISIBLE ' + str(bad) if bad else 'not divisible'}")

# --- (b) локальная проверка: G - T in l*E(F_P) ?
def redmap(P):
    F = P.residue_field()
    Ered = E.reduction(P)
    def rp(pt):
        if pt.is_zero(): return Ered(0)
        return Ered([F(pt[0]), F(pt[1])])
    return Ered, rp

proved = {}
for pr in primes(3, 300):
    for P in k.primes_above(pr):
        if E.has_bad_reduction(P): continue
        Ered, rp = redmap(P)
        m = Ered.cardinality()
        Ab = Ered.abelian_group()
        invs = [ZZ(d) for d in Ab.invariants()]
        for l in ZZ(m).prime_divisors():
            if l == pr or proved.get(l): continue
            ok = True
            for T in tors:
                Q = rp(G) - rp(T)
                co = Ab.discrete_log(Q)
                if all(ZZ(c) % gcd(l, d) == 0 for c, d in zip(co, invs)):
                    ok = False; break
            if ok:
                proved[l] = (pr, P.norm(), m)
print("l-saturation locally proved for l in:", sorted(proved))
for l in sorted(proved):
    print(f"    l={l}: prime above {proved[l][0]}, N(P)={proved[l][1]}, |E(F_P)|={proved[l][2]}")
