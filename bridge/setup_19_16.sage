m = 19; n = 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
e1 = -b; e2 = -s*m^4; e3 = -s*n^4
print("s  =", s, " squarefree class:", QQ(s).numerator()*QQ(s).denominator())
print("b  =", b); print("e1 =", e1, " e2 =", e2, " e3 =", e3)
E1,E2,E3 = 4*e1, 4*e2, 4*e3
print("4e:", E1, E2, E3)
R.<x> = QQ[]
f = (x-E1)*(x-E2)*(x-E3)
print("f =", f)
c = f.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
print("E =", E)
Emin = E.minimal_model(); print("Emin =", Emin)
iso = E.isomorphism_to(Emin); print("iso E->Emin:", iso)
print("u of iso:", iso.u if hasattr(iso,'u') else iso)
print("cond =", factor(Emin.conductor()))
print("disc E =", factor(E.discriminant()))
print("torsion =", E.torsion_subgroup().invariants(), [P for P in E.torsion_points()])
save((E,Emin,E1,E2,E3,s,b), '/home/kep/magicKube/bridge/E_19_16.sobj')
