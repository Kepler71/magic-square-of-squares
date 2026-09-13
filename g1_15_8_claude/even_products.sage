import itertools
S.<x> = QQ[]; T.<u> = QQ[]
forms = {'F0':289*(1+x^2)-322*x,'F4':289*(1+x^2),'F8':289*(1+x^2)+322*x,'L':49+529*x^2,'U':529+49*x^2}
def to_u(P):   # чётный многочлен P(x) = g(x^2)
    assert all(P[i]==0 for i in range(1,P.degree()+1,2))
    return T([P[2*i] for i in range(P.degree()//2+1)])
def jac_quartic(f):
    f = T(f); a,b,c,d,e = [f[i] for i in (4,3,2,1,0)]
    I = 12*a*e - 3*b*d + c^2; J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return EllipticCurve(QQ,[-27*I,-27*J])
def ell_of(g):   # y^2 = g(u), deg 3 или 4, с рациональной точкой -> эллиптическая
    if g.degree()==3:
        A,B,C,D = g[3],g[2],g[1],g[0]
        return EllipticCurve(QQ,[0,B,0,A*C,A^2*D])
    return jac_quartic(g)
def rk(E):
    r = pari(E.minimal_model().a_invariants()).ellinit().ellrank(); return (int(r[0]),int(r[1]))
names = list(forms)
for size in (2,3,4,5):
    for combo in itertools.combinations(names,size):
        P = prod(forms[n] for n in combo)
        if any(P[i]!=0 for i in range(1,P.degree()+1,2)): continue
        g = to_u(P); out=[]
        for lab,h in (("y²=g(u)",g),("y²=u·g(u)",u*g)):
            d = h.degree()
            if d in (3,4) and h.discriminant()!=0:
                out.append(f"{lab} [род 1] ранг ∈ {rk(ell_of(h))}")
            else:
                out.append(f"{lab} [степень {d}, род {(d-1)//2}]")
        print("·".join(combo), "→", ";  ".join(out))
