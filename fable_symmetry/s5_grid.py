#!/usr/bin/env python3
# Fable, 14.09.2026. Сеточная форма: магический квадрат = 3x3 двумерная АП a+ib+jc из квадратов. Линии совпадений и Ферма.
from sage.all import *
import itertools, json
out={}
R = PolynomialRing(QQ,'a,b,c'); a,b,c = R.gens()
Luc = [[a-b,a+b+c,a-c],[a+b-c,a,a-b+c],[a+c,a-b-c,a+b]]
grid = {(i,j): a+i*b+j*c for i in (-1,0,1) for j in (-1,0,1)}
cells = sorted(sum(Luc,[]), key=str); gv = sorted(grid.values(), key=str)
print("клетки Лукаса = {a+ib+jc}:", cells==gv)
assert cells==gv
# 3-члённые АП среди 9 клеток при общих (b,c): ровно 8 линий сетки; 4-члённых нет
vals=list(grid.items())
def is_ap(seq):
    d=[seq[k+1]-seq[k] for k in range(len(seq)-1)]; return all(x==d[0] for x in d) and d[0]!=0
ap3=set(); ap4=0
for combo in itertools.combinations(vals,3):
    for perm in itertools.permutations(combo):
        if is_ap([p[1] for p in perm]): ap3.add(frozenset(p[0] for p in perm)); break
for combo in itertools.combinations(vals,4):
    for perm in itertools.permutations(combo):
        if is_ap([p[1] for p in perm]): ap4+=1; break
print("3-члённых АП (общие b,c):", len(ap3), " 4-члённых:", ap4)
assert len(ap3)==8 and ap4==0
out['ap3']=8; out['ap4']=0
# На линиях совпадений появляются длинные АП: b=c -> 1-2b,1-b,1,1+b,1+2b (5 членов); b=2c -> 7 членов
for name,sub in [('b=c',{b:c}),('b=2c',{b:2*c}),('c=0',{c:R(0)})]:
    vs=sorted(set(v.subs(sub) for v in grid.values()), key=lambda f: f.subs({a:0,c:1}))
    print(f"  {name}: различных клеток {len(vs)}: {vs}")
# Ферма: четырёх квадратов в АП нет. Численный контроль на кривой x^2+z^2=2y^2, y^2+w^2=2z^2 (высота <= 300): только |x|=|y|=|z|=|w|
H=300; nontriv=0
sq={k*k:k for k in range(0,4*H)}
for x in range(1,H):
    for y in range(x+1,H):
        z2=2*y*y-x*x
        if z2 in sq:
            w2=2*z2-y*y
            if w2 in sq and gcd(x,y)==1: nontriv+=1; print("   !!", x,y,sq[z2],sq[w2])
print("примитивные 4-АП квадратов с x<y<300:", nontriv)
assert nontriv==0
# Нулевая клетка на V(Q) невозможна: клетка 0 => противоположная = 2*центр = 2e^2 должна быть квадратом => 2 квадрат.
print("клетка 0 на V(Q): противоположная клетка 2e^2 — квадрат лишь при e=0; при центре 0 все u_i=0. => V(Q) не содержит точек с нулевой клеткой и ненулевым центром.")
# Быстрый поиск рациональных точек V малой высоты в координатах (b,c) при a=1: b=p/q, c=r/q, |p|,|r|<=N, q<=N: все 8 клеток квадраты?
N=60; found=[]
for q in range(1,N+1):
    for p in range(-q+1,q):
        for r in range(-q+1,q):
            if gcd(gcd(p,r),q)!=1: continue
            if p==0 or r==0: continue
            us=[p,r,p+r,p-r]
            if any(abs(v)>=q for v in us): continue
            if all(ZZ(q*q - v*q).is_square() and ZZ(q*q+v*q).is_square() for v in us):
                found.append((p,r,q))
print(f"точки V(Q) с a=1, b=p/q, c=r/q, q<={N}, bc!=0, все клетки положительны: {found}")
out['small_points_bc_nonzero']=found
# Мод 24: у примитивного целого магического квадрата из квадратов b,c ≡ 0 (mod 24) — из s4 (mod 8, mod 3); контроль на Б–С: b=41496, c=-138600
print("Б–С: b mod 24 =", 41496%24, " c mod 24 =", (-138600)%24, " (обе клетки-неквадраты всё равно дают сравнение? нет — там 7 квадратов; проверка иллюстративна)")
json.dump(out, open('/home/kep/magicKube/fable_symmetry/s5_grid.json','w'), indent=1)
print("OK s5")
