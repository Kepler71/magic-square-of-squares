# Fable, 15.09.2026. Сечения семейств шести клеток над Qbar(a,b) и поле их определения.
# E^-: y^2=(x-a^2)(x-b^2)(x-c^2); E^+: v^2=(1-a^2X)(1-b^2X)(1-c^2X). Ищем x (resp. X) — отношение квадратичных форм от a,b
# с малыми целыми коэффициентами, при котором правая часть — c·h^2, h ∈ Q[a,b]; тогда сечение определено над Q(√c).
from sage.all import *
import itertools, json
R=PolynomialRing(QQ,'a,b'); a,b=R.gens()
def sqclass(f):
    # f — форма над Q; возвращает (True, c) если f = c*h^2 над Q, иначе (False,None)
    if f==0: return (True,0)
    fa=f.factor(); c=fa.unit()
    for g,e in fa:
        if e%2: return (False,None)
    return (True,QQ(c).squarefree_part())
res={}
quads=[]
rng=range(-3,4)
for al,be,ga in itertools.product(rng,repeat=3):
    if (al,be,ga)==(0,0,0): continue
    quads.append(al*a**2+be*a*b+ga*b**2)
for name,c in (('c=a+b',a+b),('c=a-b:(r,s-r,s+r)->(a=s-r? нет) используем тройку (r,s-r,s+r)',None)):
    pass
triples={'E-:c=a+b':(a,b,a+b),
         'E-:(r,s-r,s+r)':(a,b-a,b+a),     # r=a, s=b
         'E-:(s-r,s,s+r)':(b-a,b,b+a),
         'E-:(r,s-r,s)':(a,b-a,b)}          # тоже c=a+b (a+(b-a)=b) — контроль
for name,(A,B,C) in triples.items():
    found=[]
    for q in quads:
        f=(q-A**2)*(q-B**2)*(q-C**2)
        ok,cst=sqclass(f)
        if ok and cst!=0: found.append((str(q),str(cst)))
    # также x = q/ d^2 с d ∈ {a,b,c,a-b}? (нецелые сечения) — x=q/(den^2): f*(den^6) должен быть квадратом
    for den in (a,b,a+b,a-b,2*a+b,a+2*b):
        for q in quads:
            f=(q-A**2*den**2)*(q-B**2*den**2)*(q-C**2*den**2)
            ok,cst=sqclass(f)
            if ok and cst!=0: found.append((str(q)+'/('+str(den)+')^2',str(cst)))
    res[name]=found
    print(name,len(found),"сечений; поля констант:",sorted(set(c for _,c in found)),flush=True)
    for q,c in found[:12]: print("   x =",q,"  над Q(√",c,")")
# E^+: v^2=(1-a^2X)(1-b^2X)(1-c^2X), X=q1/q2
for name,(A,B,C) in (('E+:c=a+b',(a,b,a+b)),('E+:(r,s-r,s+r)',(a,b-a,b+a)),('E+:(s-r,s,s+r)',(b-a,b,b+a))):
    found=[]
    dens=[a**2,b**2,C**2,a*b,a*C,b*C,a**2+b**2,(a+b)**2,(a-b)**2,a**2-b**2,a*b+a*C+b*C, a**2+a*b+b**2]
    for q2 in dens:
        for q1 in quads:
            f=(q2-A**2*q1)*(q2-B**2*q1)*(q2-C**2*q1)*q2   # v^2 q2^3 = ... ; умножаем на q2^4: (q2 v)^2*q2^2 ... проще: X=q1/q2, v^2 = prod(1-λ^2 q1/q2) = prod(q2-λ^2 q1)/q2^3 -> (v q2^2)^2 = q2 * prod
            ok,cst=sqclass(f)
            if ok and cst!=0 and q1!=0: found.append((str(q1)+'/'+str(q2),str(cst)))
    res[name]=found
    print(name,len(found),"сечений; поля констант:",sorted(set(c for _,c in found)),flush=True)
    for q,c in found[:12]: print("   X =",q,"  над Q(√",c,")")
json.dump(res,open('/home/kep/magicKube/fable_bl/sections_search.json','w'),indent=1)
