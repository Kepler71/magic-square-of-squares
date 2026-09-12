# Independent setup for (m,n)=(19,16), approach m2 (Heegner / analytic)
m = 19; n = 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
e1 = -b; e2 = -s*m^4; e3 = -s*n^4
print("s =", s, " b =", b)
print("e1,e2,e3 =", e1, e2, e3)
# x = 4X + c, y = 8V
c = -4*(e1+e2+e3)/3
print("c =", c, " integral:", c in ZZ)
r = [4*e1+c, 4*e2+c, 4*e3+c]
print("roots r =", r, " sum =", sum(r))
assert all(ri in ZZ for ri in r)
R.<x> = QQ[]
f = (x-r[0])*(x-r[1])*(x-r[2])
print("f =", f)
A = f.coefficients(sparse=False)[1]; B = f.coefficients(sparse=False)[0]
print("A =", A, "B =", B)
M = EllipticCurve([0,0,0,A,B])
print("M =", M)
print("disc =", factor(M.discriminant()))
print("conductor =", M.conductor(), "=", factor(M.conductor()))
print("torsion =", M.torsion_subgroup().invariants())
Mm = M.minimal_model()
print("minimal model =", Mm)
print("required class (1,s,s) ->", [QQ(1).squarefree_part() if False else 1, s.numerator()*s.denominator(), s.numerator()*s.denominator()])
save(M, "/home/kep/magicKube/bridge/M_19_16_m2.sobj")
