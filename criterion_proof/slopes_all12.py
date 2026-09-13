from math import gcd, isqrt
from odd_criterion import proves_rank0
import json
def sq(x): return x>=0 and isqrt(x)**2==x
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
def ok(a,b):
    if a<=0 or b<=0 or a==b: return False
    m,n=red(a,b); return proves_rank0(m,n) and not (sq(m) and sq(n))
PAIRS=lambda r,s: [(r,s),(r,s-r),(s,s-r),(s,r+s),(r,r+s),            # кубические
                   (s-r,r+s),(s,abs(2*r-s)),(s,2*r+s),(r,2*s-r),(r,2*s+r),(s,2*r),(r,2*s)]   # квартичные
out={}
for lo,hi in ((2,48),(49,100),(101,200)):
    tot=c=0; open_=[]
    for s in range(lo,hi+1):
        for r in range(1,s):
            if gcd(r,s)!=1 or 2*r==s: continue
            tot+=1
            if any(ok(a,b) for a,b in PAIRS(r,s)): c+=1
            else: open_.append(f"{r}/{s}")
    out[f"{lo}-{hi}"]={'total':tot,'closed':c,'open':open_}
    print(f"знаменатели {lo}..{hi}: наклонов {tot} (без вырожденного 1/2); закрыто доказанно {c} ({100*c/tot:.1f}%)")
json.dump(out,open("slopes_all12.json","w"),indent=1)
