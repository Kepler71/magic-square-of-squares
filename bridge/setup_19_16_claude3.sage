# -*- coding: utf-8 -*-
# (19,16) : базовые факты, пересчитанные с нуля.
import time
m, n = 19, 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
print("s =", s, " b =", b)

# identities check (exact, polynomial)
Rt.<tv> = QQ[]
assert b*tv^2 + b            == (m*n)^2 * (s*(1+tv^2))
assert b*tv^2 + s*m^4        == s*m^2  * (m^2 + n^2*tv^2)
assert b*tv^2 + s*n^4        == s*n^2  * (n^2 + m^2*tv^2)
print("[proved] three identities X-e_i hold identically in t")

def sqcl(q):
    q = QQ(q); assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()
delta_star = (1, sqcl(s), sqcl(s))
print("delta* =", delta_star)
print("s square in Q?", QQ(s).is_square())

# integral model, scale 4 (a square) -> delta classes unchanged
r = [ZZ(-4*b), ZZ(-4*s*m^4), ZZ(-4*s*n^4)]
print("roots of scaled model (e1,e2,e3) =", r)
Rx.<xv> = QQ[]
f = (xv-r[0])*(xv-r[1])*(xv-r[2])
E0 = EllipticCurve([0, f[2], 0, f[1], f[0]])
E  = E0.minimal_model()
iso = E0.isomorphism_to(E)
print("E0 =", E0)
print("E  =", E)
print("iso u,r,s,t =", iso.u, iso.r, iso.s, iso.t)
R = [ZZ(iso.u^(-2)*(ri - iso.r)) for ri in r]
print("roots on E in same order:", R)
for Ri in R: assert E.is_on_curve(Ri,0)
print("u^2 =", iso.u^2, " (square scale -> delta preserved:", (iso.u^2).is_square(), ")")
print("conductor =", E.conductor(), "=", factor(E.conductor()))
print("disc =", factor(E.discriminant()))
print("torsion =", E.torsion_subgroup().invariants())
print("root number =", E.root_number())
t0=time.time(); print("selmer_rank (dim Sel^2) =", E.selmer_rank(), " %.1fs"%(time.time()-t0))
t0=time.time(); print("rank_bound (mwrank) =", E.rank_bound(), " %.1fs"%(time.time()-t0))
print("isogeny class sizes:", [EE.ainvs() for EE in E.isogeny_class().curves])
save(E, "/home/kep/magicKube/bridge/logs_19_16_claude3/E.sobj")
save(R, "/home/kep/magicKube/bridge/logs_19_16_claude3/R.sobj")
