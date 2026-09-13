# Независимая проверка (15,8) другим путём: ранги эллиптических кривых y^2 = Q_i(x) Q_j(x)
# для пяти квадратичных форм после замены t = (1+x)/(1-x). Код свой, без скриптов Codex.
import itertools
R.<t> = QQ[]
F0 = 225 + 64*t^2; F4 = QQ(289)/2*(1+t^2); F8 = 64 + 225*t^2; L = F4 - 240*t; U = F4 + 240*t
S.<x> = QQ[]
def sub(F):   # F((1+x)/(1-x)) * (1-x)^2
    num = sum(F[i]*(1+x)^i*(1-x)^(2-i) for i in range(3))
    return S(num)
forms = {'F0':sub(F0),'F4':sub(F4),'F8':sub(F8),'L':sub(L),'U':sub(U)}
for k,v in forms.items(): print(k, v, " при x=0:", v(0), "квадрат" if QQ(v(0)).is_square() else "НЕ квадрат")
# сверка с формами из записки Codex (тождества, не его код)
exp = {'F0':289*(1+x^2)-322*x,'F4':289*(1+x^2),'F8':289*(1+x^2)+322*x,'L':49+529*x^2,'U':529+49*x^2}
print("совпадение с записанными формами:", all(forms[k]==exp[k] for k in exp))
def jac_quartic(f):
    a,b,c,d,e = [f[i] for i in (4,3,2,1,0)]
    I = 12*a*e - 3*b*d + c^2
    J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return EllipticCurve(QQ,[-27*I,-27*J])
print()
for (n1,f1),(n2,f2) in itertools.combinations(forms.items(),2):
    f = f1*f2
    if f.discriminant()==0: print(f"{n1}*{n2}: вырожденная квартика"); continue
    E = jac_quartic(f).minimal_model()
    r = pari(E.a_invariants()).ellinit().ellrank()
    print(f"{n1}*{n2}: ранг ∈ [{r[0]},{r[1]}], кручение {E.torsion_subgroup().invariants()}, conductor {E.conductor()}")
