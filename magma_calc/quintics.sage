S.<x>=QQ[]; T.<u>=QQ[]
F0=289*(1+x^2)-322*x; F8=289*(1+x^2)+322*x; L=49+529*x^2; U=529+49*x^2; F4=289*(1+x^2)
def to_u(P): return T([P[2*i] for i in range(P.degree()//2+1)])
for name,P in (('F0F4F8L',F0*F4*F8*L),('F0F4F8U',F0*F4*F8*U),('F0F8LU',F0*F8*L*U)):
    g=to_u(P); h=u*g
    # убрать квадратные множители содержания: h = c*h0, c -> бесквадратная часть
    c=h.content() if hasattr(h,'content') else 1
    cont=gcd([ZZ(a) for a in h.list()]); sqp=ZZ(cont).squarefree_part()
    h0=T([ZZ(a)//cont*sqp for a in h.list()])     # y^2 = sqp*(h/cont), эквивалентно над Q
    # доп. упрощение Sage: минимизация по reduce не делаем, передаём как есть
    print(name, "степень", h0.degree(), "коэффициенты:", h0.list())
    open(f"q_{name}.m","w").write(f"H := HyperellipticCurve(Polynomial({h0.list()}));\npts, ok := RationalPointsGenus2(H);\nprint \"POINTS:\", pts;\nprint \"COMPLETE:\", ok;\n")
