# G1 pair (19,16): search for a generator of E(Q) of infinite rank part.
# Approach: 2-isogeny descent chain + point search on isogenous curves, pull back.
import time
t00 = time.time()

m = 19; n = 16
s = QQ(m^2 + n^2)/2      # 617/2
b = s*m^2*n^2
e1 = -b
e2 = -s*m^4
e3 = -s*n^4
print("m,n =", m, n)
print("s  =", s)
print("b  =", b)
print("e1,e2,e3 =", e1, e2, e3)

R.<X> = QQ[]
f = (X - e1)*(X - e2)*(X - e3)
print("f =", f)
A2 = f[2]; A1 = f[1]; A0 = f[0]
E = EllipticCurve([0, A2, 0, A1, A0])
print("E  =", E)
print("disc(E) =", factor(E.discriminant()))
print("j(E) =", E.j_invariant())
Emin = E.minimal_model()
print("Emin =", Emin)
print("conductor =", factor(Emin.conductor()))
print("torsion =", E.torsion_subgroup().invariants(), E.torsion_points())
print("elapsed", time.time()-t00)
