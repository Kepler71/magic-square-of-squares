m = 19; n = 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
E1r = [-b, -s*m^4, -s*n^4]
r = [4*e for e in E1r]
R.<x> = QQ[]
f = (x-r[0])*(x-r[1])*(x-r[2])
c = f.coefficients(sparse=False)     # c0 + c1 x + c2 x^2 + x^3
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
print("E =", E)
print("a-invariants:", E.a_invariants())
print("discriminant factored:", factor(E.discriminant()))
print("conductor =", E.conductor(), "=", factor(E.conductor()))
T = E.torsion_subgroup()
print("torsion invariants:", T.invariants())
print("torsion points:", [P.xy() if P != E(0) else 'O' for P in T.points()])
Em = E.minimal_model()
print("minimal model:", Em, Em.a_invariants())
u,rr,ss,tt = E.isomorphism_to(Em).tuple()
print("iso (u,r,s,t) E->Em:", u, rr, ss, tt)
save(E, '/home/kep/magicKube/bridge/E_19_16.sobj')
