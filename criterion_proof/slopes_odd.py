from math import gcd
from fractions import Fraction as F
from odd_criterion import proves_rank0
import json
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
CUBIC=lambda r,s: {'{r,s}':(r,s),'{r,s-r}':(r,s-r),'{s,s-r}':(s,s-r),'{s,r+s}':(s,r+s),'{r,r+s}':(r,r+s)}
QUART=lambda r,s: {'{s-r,r+s}':(s-r,r+s),'{s,|2r-s|}':(s,abs(2*r-s)),'{s,2r+s}':(s,2*r+s),'{r,2s-r}':(r,2*s-r),'{r,2s+r}':(r,2*s+r),'{s,2r}':(s,2*r),'{r,2s}':(r,2*s)}
res={}; sample=[]
for lo,hi in ((2,48),(49,100),(101,200)):
    tot=c_cub=c_all=0
    for s in range(lo,hi+1):
        for r in range(1,s):
            if gcd(r,s)!=1: continue
            tot+=1
            cub=any(a!=b and a>0 and proves_rank0(*red(a,b)) for a,b in CUBIC(r,s).values())
            qua=any(a!=b and a>0 and proves_rank0(*red(a,b)) for a,b in QUART(r,s).values())
            c_cub+=cub; c_all+=(cub or qua)
            if cub and len(sample)<150 and (r*7+s)%13==0: sample.append((r,s))
    res[f"{lo}-{hi}"]=(tot,c_cub,c_all)
    print(f"знаменатели {lo}..{hi}: наклонов {tot}; закрыто ДОКАЗАННО (кубические пары) {c_cub}; вместе с квартиками {c_all}")
json.dump({'coverage':res,'sample':sample},open("slopes_odd.json","w"))
