"""Явный критерий (кандидат в теорему) для E_{m,n}: y^2=x(x+m^2)(x+n^2), gcd(m,n)=1.
S^phi(E)  = { d=±∏P, P⊆primes(mn): (d/q)=1 ∀ q≡1(4), q|n²-m²;  2-адика: если m,n нечётны и m≡±n (8) → d≡±1 (8);
                                                                         если v2(mn)=2 → d нечётно }
S^phihat(E') = { d=∏Q>0, Q⊆primes(n²-m²): все нечётные q|d ≡1 (4); (d/p)=1 ∀ нечётных p|mn;
                 2-адика: d≡5 (8) допустимо только если (m,n нечётны) или v2(mn)=2;
                          d чётно (только при m,n нечётных): d/2≡1 (4) и m≡±n (8) }
rank E ≤ log2|S^phi|+log2|S^phihat|-2; чистый ранг 0 ⟺ |S^phi||S^phihat|=4.
"""
from math import gcd, log2
from itertools import combinations
_LIM=1_000_000
_spf=list(range(_LIM+1))
for _i in range(2,int(_LIM**0.5)+1):
    if _spf[_i]==_i:
        for _j in range(_i*_i,_LIM+1,_i):
            if _spf[_j]==_j: _spf[_j]=_i
def factorint(x):
    x=abs(x); out={}
    if x<=_LIM:
        while x>1:
            p=_spf[x]; out[p]=out.get(p,0)+1; x//=p
        return out
    from sympy import factorint as _f
    return _f(x)
def jacobi_symbol(a,n):
    a%=n; r=1
    while a:
        while a%2==0:
            a//=2
            if n%8 in (3,5): r=-r
        a,n=n,a
        if a%4==3 and n%4==3: r=-r
        a%=n
    return r if n==1 else 0

def legendre(d,p):  # p нечётное простое, d целое (со знаком)
    return jacobi_symbol(d%p,p)

def v2(x):
    e=0
    while x%2==0: x//=2; e+=1
    return e

def S_phi(m,n):
    fac=factorint(m*n); ps=sorted(fac)
    q1=[q for q in factorint(n*n-m*m) if q%4==1]
    both_odd=(m%2==1 and n%2==1); pm8=both_odd and ((m-n)%8==0 or (m+n)%8==0)
    v2mn=v2(m*n) if m*n%2==0 else 0
    out=[]
    for r in range(len(ps)+1):
        for sub in combinations(ps,r):
            d0=1
            for p in sub: d0*=p
            for d in (d0,-d0):
                if any(legendre(d,q)!=1 for q in q1): continue
                if pm8 and d%8 not in (1,7): continue
                if v2mn==2 and d%2==0: continue
                out.append(d)
    return out

def S_phihat(m,n):
    N=n*n-m*m; fac=factorint(abs(N))
    qs=[q for q in sorted(fac) if q==2 or q%4==1]   # только допустимые простые
    podd=[p for p in factorint(m*n) if p!=2]
    both_odd=(m%2==1 and n%2==1); pm8=both_odd and ((m-n)%8==0 or (m+n)%8==0)
    v2mn=v2(m*n) if m*n%2==0 else 0
    out=[]
    for r in range(len(qs)+1):
        for sub in combinations(qs,r):
            d=1
            for q in sub: d*=q
            if any(legendre(d,p)!=1 for p in podd): continue
            if d%2==1:
                if d%8==5 and not (both_odd or v2mn==2): continue
            else:
                if not both_odd: continue   # 2∤N иначе
                if (d//2)%4!=1 or not pm8: continue
            out.append(d)
    return out

def criterion(m,n):
    S=S_phi(m,n); T=S_phihat(m,n)
    return {'S_phi':sorted(S,key=abs),'S_phihat':sorted(T),'rank_upper':int(log2(len(S)*len(T)))-2,'pure_rank0':len(S)*len(T)==4}

if __name__=='__main__':
    import sys; print(criterion(int(sys.argv[1]),int(sys.argv[2])))
