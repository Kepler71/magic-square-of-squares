import sys, random
sys.path.insert(0,'../')
from selmer_exact import locally_solvable
from sympy import factorint
def primes_of(x): return sorted(factorint(abs(x)).keys())
from math import gcd
from itertools import combinations
random.seed(20260913)
def sqfree_divs(x):
    ps=primes_of(x); out=[]
    for r in range(len(ps)+1):
        for sub in combinations(ps,r):
            d=1
            for p in sub: d*=p
            out.append(d)
    return out
stats={1:[0,0],2:[0,0],3:[0,0],4:[0,0]}; viol=[]
def check(lem,coeffs):
    stats[lem][0]+=1
    if locally_solvable(coeffs,2): stats[lem][1]+=1; viol.append((lem,coeffs))
tries=0
while min(s[0] for s in stats.values())<150 and tries<200000:
    tries+=1
    k=random.randint(3,30)
    n=random.randrange(1,20000,2)
    lem=random.choice([1,2,3,4])
    if lem==1:     # m,n нечётны, m ≡ ±n mod 2^k (k>=3), d ≡ 3,5 mod 8, d | mn, C_d
        m=n+random.choice([1,-1])*(2**k)*random.randint(1,3)
        if m<=0 or gcd(m,n)!=1 or m==n: continue
        m=abs(m); cand=[s*d for d in sqfree_divs(m*n) for s in (1,-1) if (s*d)%8 in (3,5)]
        if not cand: continue
        d=random.choice(cand); A=m*m+n*n; B=m*m*n*n; check(1,[B//d,0,A,0,d])
    if tries%400==0: print('прогресс', {k:v[0] for k,v in stats.items()}, flush=True)
    if lem==2 and False: pass
    elif lem==2:   # m = 4h, n нечётно, d чётно, C_d
        h=random.randrange(1,3000,2); m=4*h
        if gcd(m,n)!=1: continue
        cand=[s*d for d in sqfree_divs(m*n) for s in (1,-1) if d%2==0]
        d=random.choice(cand); A=m*m+n*n; B=m*m*n*n; check(2,[B//d,0,A,0,d])
    elif lem==3:   # m чётно с v2(m)=e≠2 (включая большие e), n нечётно, d ≡ 5 mod 8, d | n²−m², d>0, C'_d
        e=random.choice([1]+list(range(3,31))); m=(2**e)*random.randrange(1,60,2)
        if gcd(m,n)!=1: continue
        D=n*n-m*m; cand=[d for d in sqfree_divs(abs(D)) if d%8==5]
        if not cand: continue
        d=random.choice(cand); A=m*m+n*n; check(3,[D*D//d,0,-2*A,0,d])
    else:          # m,n нечётны, m ≢ ±n mod 8, d чётно, C'_d
        m=random.randrange(1,20000,2)
        if gcd(m,n)!=1 or (m-n)%8==0 or (m+n)%8==0: continue
        D=n*n-m*m; cand=[d for d in sqfree_divs(abs(D)) if d%2==0]
        d=random.choice(cand); A=m*m+n*n; check(4,[D*D//d,0,-2*A,0,d])
for L,(c,s) in stats.items(): print(f"лемма {L}: случаев {c}, разрешимых над Q_2 (нарушений леммы) {s}")
print("нарушения:", viol[:3])
