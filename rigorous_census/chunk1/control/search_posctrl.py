# Claude 26.09: поиск положительного контроля. q(m) = −4m(m²−1)/(1+m²)²: 1±q — квадраты (точки окружности u²+v²=2).
# T = D·{q(m1),q(m2),q(m3)} (целые), z* = 1/D: все 1±λz* — квадраты. Ищем E_T ранга 1 (ellrank [1,1]).
from sage.all import *
import itertools, json, time, sys
sys.path.insert(0,'/home/kep/magicKube/rigorous_census/chunk1')
from dem_rigorous import model, all_sq
def q(m): m=QQ(m); return -4*m*(m**2-1)/(1+m**2)**2
ms=sorted(set(QQ(a)/b for b in range(1,6) for a in range(-12,13) if gcd(a,b)==1 and QQ(a)/b not in (0,1,-1)), key=lambda x:(x.height(),x))
t0=time.time(); found=[]; n=0
for m1,m2,m3 in itertools.combinations(ms[:40],3):
    if time.time()-t0>500: break
    Q=[q(m1),q(m2),q(m3)]
    if len(set(abs(x) for x in Q))<3: continue
    D=lcm([x.denominator() for x in Q]); T=sorted(ZZ(x*D) for x in Q)
    g=gcd(T); T=[x//g for x in T]; zs=QQ(g)/D
    if sorted(-x for x in T)==T: continue
    assert all_sq(zs, T+[-x for x in T])
    try: E1,L,c0=model(T)
    except Exception: continue
    if E1.discriminant()==0: continue
    n+=1
    lo,hi=[int(v) for v in pari(E1.minimal_model()).ellrank()[:2]]
    if lo==hi==1:
        found.append(dict(T=[int(x) for x in T],z=str(zs),m=[str(m1),str(m2),str(m3)])); print('РАНГ 1:',found[-1],flush=True)
        if len(found)>=5: break
print('проверено',n,'наборов, ранг 1:',len(found),'время',round(time.time()-t0),flush=True)
json.dump(found,open('posctrl_candidates.json','w'))
