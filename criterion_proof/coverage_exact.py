from math import gcd, isqrt, log2
from functools import lru_cache
from selmer_exact import selmer
import json
def sq(x): return x>=0 and isqrt(x)**2==x
@lru_cache(maxsize=None)
def rank0(m,n):
    S,T=selmer(m,n); return len(S)*len(T)<=4
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
PAIRS=lambda r,s: [(r,s),(r,s-r),(s,s-r),(s,r+s),(r,r+s),(s-r,r+s),(s,abs(2*r-s)),(s,2*r+s),(r,2*s-r),(r,2*s+r),(s,2*r),(r,2*s)]
res={}; pyth_cases=[]
for lo,hi in ((2,48),(49,100),(101,200)):
    tot=c=0; open_=[]
    for s in range(lo,hi+1):
        for r in range(1,s):
            if gcd(r,s)!=1 or 2*r==s: continue
            tot+=1; hit=None
            for a,b in PAIRS(r,s):
                if a<=0 or a==b: continue
                m,n=red(a,b)
                if rank0(m,n):
                    if sq(m) and sq(n) and sq(isqrt(m)**2+isqrt(n)**2): pyth_cases.append((r,s,m,n)); hit=hit or 'pyth'
                    else: hit='ok'; break
            if hit: c+=1
            else: open_.append(f"{r}/{s}")
    res[f"{lo}-{hi}"]={'total':tot,'closed':c,'open':open_}
    print(f"наклоны, знаменатели {lo}..{hi}: закрыто (точный Селмер) {c} из {tot} ({100*c/tot:.1f}%)")
print("пифагоровы случаи, требующие поточечного разбора:", sorted(set((r,s) for r,s,*_ in pyth_cases)))
for N in (20,50,100,200):
    pairs=[(m,n) for m in range(2,N+1) for n in range(1,m) if gcd(m,n)==1]
    c=sum(rank0(*red(m*m,n*n)) for m,n in pairs)
    print(f"G1, m,n ≤ {N}: закрыто {c} из {len(pairs)}")
json.dump(res,open('coverage_exact.json','w'),indent=1)
