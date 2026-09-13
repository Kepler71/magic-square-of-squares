# Claude: символьная проверка THEOREMS_Z2Z6 §1, §5 своим кодом
R.<a,b,w,u> = QQ[]
c=a+b; M=b^3*(2*a+b); N=a^3*(a+2*b)
# (1) M−N
assert M-N == (b-a)*(a+b)^3
# E+: v²=(1−a²u)(1−b²u)(1−c²u) -> X=−a²b²c²u, V=−a²b²c²v, x=X+a²b² -> y²=x(x+M)(x+N)
P=a^2*b^2*c^2
F=R.fraction_field()
X=-P*u; x=X+a^2*b^2
lhs=(P^2)*(1-a^2*u)*(1-b^2*u)*(1-c^2*u)      # V² = P²·v²
rhs=x*(x+M)*(x+N)
print("(1) замена E+ -> E(a,b):", expand(lhs+rhs)==0 or expand(lhs-rhs)==0, "знак:", expand(lhs-rhs)==0)
# Дуйелла–Пераль
S.<vv>=QQ[]
A=lambda t: 37-84*t+102*t^2-36*t^3-3*t^4
Bf=lambda t: 32*(t-1)^3*(t+1)^3*(3*t-5)
tv=F((2*a-b)/(2*a+b)); lam=F(4/(2*a+b))
print("A(v)=λ⁴(M+N):", F(A(tv))==lam^4*(M+N), "  B(v)=λ⁸MN:", F(Bf(tv))==lam^8*M*N)
# §5: серия
aw=8-w^2; bw=2*(2-w^2); Mw=bw^3*(2*aw+bw); Nw=aw^3*(aw+2*bw)
xw=-(w^2-5)*(w^2+4)^2*(5*w^2-16)
yw=27*(w-2)^2*w*(w+2)^2*(w^2-5)*(w^2+4)*(5*w^2-16)
print("точка на E(a(w),b(w)):", expand(yw^2 - xw*(xw+Mw)*(xw+Nw))==0)
# кручение: координаты x точек кручения (7 шт.), совпадения с xw
T=[0,-Mw,-Nw,aw^2*bw^2, aw*bw*(aw*bw+2*(aw+bw)^2), -aw*bw^2*(aw+2*bw), -aw^2*bw*(2*aw+bw)]
Ew=EllipticCurve(FractionField(QQ['w']),[0,Mw+Nw,0,Mw*Nw,0])
roots=set()
for t in T:
    g=QQ['w'](xw-t)
    for r,_ in g.roots(QQ): roots.add(r)
print("рациональные корни совпадений с кручением:", sorted(roots))
# численно: при w=1 порядок бесконечен и x-координаты кручения полные
E1=EllipticCurve([0,Mw(w=1)+Nw(w=1),0,Mw(w=1)*Nw(w=1),0]) if False else EllipticCurve(QQ,[0,QQ(Mw.subs(w=1)+Nw.subs(w=1)),0,QQ(Mw.subs(w=1)*Nw.subs(w=1)),0])
Pt=E1([QQ(xw.subs(w=1)),QQ(yw.subs(w=1))])
print("w=1: a,b =",aw.subs(w=1),bw.subs(w=1)," порядок точки:",Pt.order()," кручение:",E1.torsion_subgroup().invariants(),
      " x кручения ⊆ T(w=1):", set(P[0] for P in E1.torsion_points() if not P.is_zero()) <= set(QQ(t.subs(w=1)) for t in T))
# k_n
for n in range(1,6):
    kn=QQ(4*n^2-2)/(8*n^2-1); print(n,kn, "= b/a при w=1/n:", kn==QQ(bw.subs(w=1/n)/aw.subs(w=1/n)))
