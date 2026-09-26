# Арифметические следствия теоремы (скептик Claude/Opus, 2026-09-26). Только целочисленная арифметика.
from math import gcd, isqrt
from itertools import combinations
def factor(n):
    n=abs(n); f={}; d=2
    while d*d<=n:
        while n%d==0: f[d]=f.get(d,0)+1; n//=d
        d+=1
    if n>1: f[n]=f.get(n,0)+1
    return f
def sqfree(n): return all(e==1 for e in factor(n).values())
def Dlist(N):  # {d>0 squarefree: p|d => p|N, p=1 mod 4, d=1 mod 24}
    ps=[p for p in factor(N) if p%4==1]
    out=[]
    for k in range(len(ps)+1):
        for sub in combinations(ps,k):
            d=1
            for p in sub: d*=p
            if d%24==1: out.append(d)
    return sorted(out)
# 1) наименьший нетривиальный положительный квадратсвободный d = 1 mod 24 (с условием p=1 mod4 и без)
c1=[d for d in range(2,2000) if d%24==1 and sqfree(d)]
c2=[d for d in c1 if all(p%4==1 for p in factor(d))]
print('квадратсвободные d=1 mod24, d>1 (первые):',c1[:6],'; из них все p=1 mod4:',c2[:6])
# 2) три наклона во всех вариантах знака/порядка
for (a,b) in [(126,451),(73,362),(265,298)]:
    for (r,s) in [(a,b),(b,a),(a,-b),(-a,b)]:
        print(f'r={r:5d} s={s:5d}: r-s={r-s:5d} {factor(r-s)}  r+s={r+s:5d} {factor(r+s)}  D_-={Dlist(r-s)} D_+={Dlist(r+s)}')
# 3) порог 73: для 0<|N|<73 список {1}; при |N|=73 список {1,73} (граница точна для СПИСКА)
assert all(Dlist(N)==[1] for N in range(1,73)); print('0<|N|<73 => D={1}: да;  D(73)=',Dlist(73))
# 4) взаимная простота элементов D_- и D_+ для всех взаимно простых r,s, |r|,|s|<=300, r!=+-s
bad=0; nontriv=0
for r in range(-300,301):
    for s in range(1,301):
        if r==0 or abs(r)==s or gcd(r,s)!=1: continue
        Dm,Dp=Dlist(r-s),Dlist(r+s)
        if len(Dm)>1 or len(Dp)>1: nontriv+=1
        if any(gcd(x,y)>1 for x in Dm for y in Dp): bad+=1
print('пар (r,s) с нетривиальным списком:',nontriv,'; пар с общим простым в D_- и D_+:',bad)
