T.<u>=QQ[]; Z.<z>=QQ[]
S.<x>=QQ[]
F0=289*(1+x^2)-322*x; F8=289*(1+x^2)+322*x; L=49+529*x^2; U=529+49*x^2
P=F0*F8*L*U
g=T([P[2*i] for i in range(P.degree()//2+1)])          # P(x) = g(x^2)
h=u*g                                                    # C: Y^2 = u·g(u), род 2
print("h палиндромный:", h.list()[1:]==h.list()[1:][::-1])
# u = (1+z)/(1-z): инволюция u->1/u становится z->-z
num = sum(h[i]*(1+z)^i*(1-z)^(6-i) for i in range(h.degree()+1))   # h((1+z)/(1-z))*(1-z)^6
num = Z(num)
print("чётный по z:", all(num[i]==0 for i in range(1,7,2)), " степень", num.degree())
W.<w>=QQ[]
G = W([num[2*i] for i in range(num.degree()//2+1)])     # num(z) = G(z^2), G кубический
print("G(w) =", G)
def ell_from(poly):
    # y^2 = poly(w), степень 3 или 4 с рациональной точкой -> эллиптическая (через Jacobian квартики/кубики)
    d=poly.degree()
    if d==3:
        a,b,c,dd=poly[3],poly[2],poly[1],poly[0]
        return EllipticCurve(QQ,[0,b,0,a*c,a^2*dd])
    a,b,c,dd,e=[poly[i] for i in (4,3,2,1,0)]
    I=12*a*e-3*b*dd+c^2; J=72*a*c*e+9*b*c*dd-27*a*dd^2-27*e*b^2-2*c^3
    return EllipticCurve(QQ,[-27*I,-27*J])
E1=ell_from(G); E2=ell_from(w*G)
for nm,E in (("E1: y²=G(w)",E1),("E2: y²=w·G(w)",E2)):
    r=pari(E.minimal_model().a_invariants()).ellinit().ellrank()
    print(nm, " ранг ∈", [int(r[0]),int(r[1])], " кручение", E.torsion_subgroup().invariants(), " кондуктор", E.conductor())
