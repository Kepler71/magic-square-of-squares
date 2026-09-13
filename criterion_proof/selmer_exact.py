"""Точные группы Селмера 2-изогении для E_{m,n}: y^2 = x(x+m^2)(x+n^2) — решатель над Q_p для всех p, включая 2.
Свой код. Даёт 2^rank <= |Sel^phi|·|Sel^phihat|/4 как доказательство (ПО) для конкретной пары.
Места: ∞ и простые, делящие 2·m·n·(m²−n²) (остальные — хорошая редукция, точки есть автоматически)."""
from math import comb, gcd
from itertools import combinations
import sys
sys.setrecursionlimit(10000)
class DepthExhausted(RuntimeError): pass
def v(x,p):
    if x==0: return 10**9
    e=0
    while x%p==0: x//=p; e+=1
    return e
def taylor(coeffs,z0,s,p):
    n=len(coeffs); out=[0]*n; ps=p**s
    for i,c in enumerate(coeffs):
        if c==0: continue
        for j in range(i+1):
            out[j]+=c*comb(i,j)*z0**(i-j)*ps**j
    return out
def unit_is_square(u,p,prec):
    """u — p-адическая единица, известна по модулю p^prec. True/False или None (не определено)."""
    if p!=2:
        return (pow(u%p,(p-1)//2,p)==1) if prec>=1 else None
    if prec>=3: return u%8==1
    if prec==2: return False if u%4==3 else None
    return None
def ball(coeffs,z0,s,p,depth=0):
    if depth>400: raise DepthExhausted(f"p={p}: глубина исчерпана — это НЕ ответ «точек нет»")
    g=taylor(coeffs,z0,s,p); c0=g[0]
    if c0==0: return True
    e=v(c0,p); mu=min(v(c,p) for c in g[1:])
    if v(g[1],p)<10**9 and e>2*v(g[1],p): return True            # Гензель: корень в шаре -> точка с w=0
    if mu>e:
        if e%2==1: return False                                    # все значения имеют нечётную оценку e
        r=unit_is_square(c0//p**e,p,mu-e)
        if r is not None: return r
    return any(ball(coeffs,z0+t*p**s,s+1,p,depth+1) for t in range(p))
def locally_solvable(coeffs,p):
    """w^2 = f(z), coeffs = [c0,c1,c2,c3,c4]; есть ли точка над Q_p на гладкой модели."""
    rev=coeffs[::-1]
    return ball(coeffs,0,0,p) or ball(rev,0,1,p)
def primes_of(x):
    x=abs(x); out=[]; d=2
    while d*d<=x:
        if x%d==0:
            out.append(d)
            while x%d==0: x//=d
        d+=1 if d==2 else 2
    if x>1: out.append(x)
    return out
def selmer(m,n):
    a=m*m+n*n; b=m*m*n*n; c=(m*m-n*n)**2
    bad=sorted(set([2]+primes_of(m*n)+primes_of(m*m-n*n)))
    PM=primes_of(m*n); PN=primes_of(m*m-n*n)
    S=[]
    for r in range(len(PM)+1):
        for sub in combinations(PM,r):
            d0=1
            for q in sub: d0*=q
            for d in (d0,-d0):
                Cd=[b//d,0,a,0,d]                                    # w^2 = d z^4 + a z^2 + b/d   (все d вещественно разрешимы)
                if all(locally_solvable(Cd,p) for p in bad): S.append(d)
    T=[]
    for r in range(len(PN)+1):
        for sub in combinations(PN,r):
            d=1
            for q in sub: d*=q
            Cd=[c//d,0,-2*a,0,d]                                     # w^2 = d z^4 - 2a z^2 + c/d, d>0 (d<0 отсекает R)
            if all(locally_solvable(Cd,p) for p in bad): T.append(d)
    return S,T
def rank_bound(m,n):
    S,T=selmer(m,n)
    from math import log2
    return int(round(log2(len(S)*len(T))))-2, S, T
if __name__=="__main__":
    m,n=int(sys.argv[1]),int(sys.argv[2]); print(rank_bound(m,n))
