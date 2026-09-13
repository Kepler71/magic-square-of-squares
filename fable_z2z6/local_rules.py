import json, glob, sys
from collections import Counter, defaultdict
def v(n,p):
    n=abs(n); e=0
    if n==0: return 99
    while n%p==0: n//=p; e+=1
    return e
recs=[]
for f in sorted(glob.glob('grid/shard_60_*.jsonl')):
    for line in open(f): recs.append(json.loads(line))
forms = lambda a,b: dict(a=a,b=b,amb=a-b,apb=a+b,ap2b=a+2*b,a2pb=2*a+b)
tab = defaultdict(Counter)
for r in recs:
    a,b=r['a'],r['b']; F=forms(a,b)
    for p,el in r['isog3']['imE'].items():
        p=int(p)
        divs = tuple(sorted((k,v(F[k],p)) for k in F if v(F[k],p)>0))
        if p>=5:
            key=(p%3, divs[0][0] if divs else 'none', 'v>=1')
        else:
            key=(p, divs)
        sub = tuple(sorted(tuple(e) for e in el))
        # subgroup description
        n=len(el)
        tab[key][(n, sub if p in (2,3) else n)] += 1
for key in sorted(tab, key=str):
    print(key, dict(tab[key]))
