#!/usr/bin/env python3
# Fable, 14.09.2026. Башня X (полумагические) ⊃ W (7 линий) ⊃ V (8 линий) из диагональных квадрик; конформный локус.
from sage.all import *
import itertools, json, time
out={}
R = PolynomialRing(QQ, 'x', 9); x = R.gens()
X = [[x[3*i+j] for j in range(3)] for i in range(3)]
cell = [[X[i][j]**2 for j in range(3)] for i in range(3)]
rows=[sum(cell[i]) for i in range(3)]; cols=[sum(cell[i][j] for i in range(3)) for j in range(3)]
d1=cell[0][0]+cell[1][1]+cell[2][2]; d2=cell[0][2]+cell[1][1]+cell[2][0]
semi=[rows[0]-rows[1], rows[1]-rows[2], cols[0]-cols[1], cols[1]-cols[2]]
IX = R.ideal(semi); IW = R.ideal(semi+[d2-rows[0]]); IV = R.ideal(semi+[d2-rows[0], d1-rows[0]])
for name,I in [('X',IX),('W',IW),('V',IV)]:
    t=time.time(); d=I.dimension()-1; hp=I.hilbert_polynomial(); deg=hp.leading_coefficient()*factorial(hp.degree())
    hs = I.hilbert_series()
    T = PowerSeriesRing(QQ,'t',default_prec=8).gen()
    HS = (hs.numerator()(T))/(hs.denominator()(T))
    coeffs=[HS[k] for k in range(5)]
    print(f"{name}: проективная размерность {d}, степень {deg}, h^0(O(k)) k=0..4: {coeffs}, {time.time()-t:.1f}s")
    out[name]=dict(dim=int(d),deg=int(deg),h0=[int(c) for c in coeffs])
# канонический класс по присоединению: K = O(2m-9), m = число квадрик
print("K_X = O(-1) (Фано), K_W = O(1), K_V = O(3); p_g(V)=h^0(O(3)) =", out['V']['h0'][3], "(Auel–Singer: 111)")
assert out['X']['dim']==4 and out['W']['dim']==3 and out['V']['dim']==2
assert out['V']['h0'][3]==111

# Полумагический квадрат = Лукас(a,b,c) + lam*I + mu*I'
Rl = PolynomialRing(QQ,'a,b,c,lam,mu'); a,b,c,lam,mu = Rl.gens()
Luc = matrix(Rl,[[a-b,a+b+c,a-c],[a+b-c,a,a-b+c],[a+c,a-b-c,a+b]])
Ip = matrix(Rl,[[0,0,1],[0,1,0],[1,0,0]])
G = Luc + lam*identity_matrix(Rl,3) + mu*Ip
rs=[sum(G.row(i)) for i in range(3)]; cs=[sum(G.column(j)) for j in range(3)]
assert all(r==rs[0] for r in rs) and all(cc==rs[0] for cc in cs)
print("Лукас + lam*I + mu*I' полумагичен, S =", rs[0], "; антидиагональ - S =", G[0,2]+G[1,1]+G[2,0]-rs[0], "; диагональ - S =", G[0,0]+G[1,1]+G[2,2]-rs[0])
# линейное пространство полумагических имеет размерность 5 -> отображение (a,b,c,lam,mu) -> 9 клеток инъективно
Mlin = matrix(QQ,[[G[i][j].monomial_coefficient(v) for v in Rl.gens()] for i in range(3) for j in range(3)])
print("ранг линейной параметризации:", Mlin.rank(), "(ожидается 5 = dim полумагических)")
assert Mlin.rank()==5
print("Клетки 7-линейного (mu=0):", (Luc+lam*identity_matrix(Rl,3)).list())

# Конформный локус C = {A : A A^T = S I} внутри X; размерность 3; параметризация Эйлера–Родрига
Rq = PolynomialRing(QQ,'w,xx,y,z'); w,xx,y,z = Rq.gens()
Rm = matrix(Rq,[[w*w+xx*xx-y*y-z*z, 2*(xx*y-w*z), 2*(xx*z+w*y)],[2*(xx*y+w*z), w*w-xx*xx+y*y-z*z, 2*(y*z-w*xx)],[2*(xx*z-w*y), 2*(y*z+w*xx), w*w-xx*xx-y*y+z*z]])
N = w*w+xx*xx+y*y+z*z
assert Rm*Rm.transpose()==N**2*identity_matrix(Rq,3)
# подстановка в уравнения X: все обращаются в нуль
sub = {x[3*i+j]: Rm[i,j] for i in range(3) for j in range(3)}
assert all(f.subs(sub)==0 for f in semi)
print("Эйлер–Родриг: A A^T = N^2 I, A лежит в X; dim C = 3 (образ P^3)")
# Точка Саллоуса со всеми плюсами лежит в X, но не в C
Sal = matrix(ZZ,[[127,46,58],[2,113,94],[74,82,97]])
print("Саллоус (все +): A A^T =", (Sal*Sal.transpose()).list(), " -> не конформна, X != C")
assert Sal*Sal.transpose()!=21609*identity_matrix(3)

# Лемма: рациональная конформная матрица имеет S = квадрат (дискриминант тернарных форм <1,1,1> и <S,S,S>).
# Численный контроль: все целые A с |a_ij|<=B и A A^T = S I -> S квадрат
B=9; cnt=0; Ss=set()
vecs=[v for v in itertools.product(range(-B,B+1),repeat=3) if any(v)]
bynorm={}
for v in vecs: bynorm.setdefault(sum(t*t for t in v),[]).append(v)
for S0,L in bynorm.items():
    for r1 in L:
        for r2 in L:
            if sum(p*q for p,q in zip(r1,r2))!=0: continue
            for r3 in L:
                if sum(p*q for p,q in zip(r1,r3))==0 and sum(p*q for p,q in zip(r2,r3))==0:
                    cnt+=1; Ss.add(S0)
print(f"целые конформные 3x3 с |a_ij|<={B}: {cnt}; значения S: {sorted(Ss)}; все квадраты: {all(ZZ(s).is_square() for s in Ss)}")
assert all(ZZ(s).is_square() for s in Ss)
out['conformal_S_values']=sorted(int(s) for s in Ss)

# V ∩ C над Qbar: кривая в P^3 кватернионов: 3 a22^2 = N^2, 3(a11^2+a33^2) = 2 N^2.
# Над Q(sqrt3): ветвь sqrt3*a22 = N  =>  w^2+y^2 = (2+sqrt3)(x^2+z^2); P=w^2-y^2, Q=x^2-z^2, m=x^2+z^2: P^2+Q^2=(1+sqrt3)^2 m^2.
K = QuadraticField(3,'r3'); r3=K.gen()
Rs = PolynomialRing(K,'s'); s = Rs.gen()
P = (1+r3)*(1-s**2); Q=(1+r3)*2*s; m=(1+s**2)   # конус P^2+Q^2=(1+r3)^2 m^2
f = [ (2+r3)*m+P, (2+r3)*m-P, m+Q, m-Q ]          # 2w^2, 2y^2, 2x^2, 2z^2
roots=[]
for g in f:
    print("  множитель:", g, " неприводим над Q(sqrt3):", g.is_irreducible(), " disc:", g.discriminant())
    roots += g.roots(QQbar) if False else [g.discriminant()]
# различность 8 точек ветвления: попарные результанты ненулевые
res_ok = all(f[i].resultant(f[j])!=0 for i in range(4) for j in range(i+1,4))
print("  попарные результанты ненулевые (8 различных точек ветвления):", res_ok)
# ещё проверить, что ни одна квадратика не имеет кратного корня
assert res_ok and all(g.discriminant()!=0 for g in f)
print("  (Z/2)^4-накрытие P^1 с 8 точками ветвления: 2g-2 = 16*(-2) + 8*8 = 32, g = 17")
out['VcapC_genus']=17
# Прямая проверка над Q: V∩C(Q) = ∅, т.к. 3 a22^2 = N^2 требует sqrt3 в Q. Контроль: перебор кватернионов
bad=0
for q in itertools.product(range(-6,7),repeat=4):
    if not any(q): continue
    Nq = sum(t*t for t in q); a22 = q[0]**2-q[1]**2+q[2]**2-q[3]**2
    if 3*a22*a22==Nq*Nq: bad+=1
print("кватернионы |q_i|<=6 с 3 a22^2 = N^2:", bad)
assert bad==0
json.dump(out, open('/home/kep/magicKube/fable_symmetry/s2_structure.json','w'), indent=1)
print("OK s2")
