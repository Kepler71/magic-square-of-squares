# Скептик-2: пересчёт списков D_-(r,s), D_+(r,s) и порога 73 (конечные циклы).
from itertools import combinations
from math import gcd
def pf(n):
    n=abs(n); out=[]; d=2
    while d*d<=n:
        if n%d==0:
            out.append(d)
            while n%d==0: n//=d
        d+=1
    if n>1: out.append(n)
    return out
def D(m):
    ps=[p for p in pf(m) if p%4==1]
    res=[]
    for k in range(len(ps)+1):
        for c in combinations(ps,k):
            d=1
            for p in c: d*=p
            if d%24==1: res.append(d)
    return sorted(res)
for r,s in [(126,451),(73,362),(265,298),(451,126),(-126,451)]:
    print((r,s),'r-s=',r-s,pf(r-s),'D-=',D(r-s),' r+s=',r+s,pf(r+s),'D+=',D(r+s))
sqf=lambda d: all(d%(q*q) for q in range(2,int(d**0.5)+1))
print('smallest nontrivial squarefree d=1 mod 24 with all primes 1 mod 4:',
      min(d for d in range(2,2000) if d%24==1 and sqf(d) and all(p%4==1 for p in pf(d))))
print('smallest nontrivial squarefree d=1 mod 24 (any primes):', min(d for d in range(2,2000) if d%24==1 and sqf(d)))
