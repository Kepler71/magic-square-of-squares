# Свой точный решатель: есть ли точка над Q_p (p нечётно) у w^2 = f(z), f ∈ Z[z] степени 4 без кратных корней.
from fractions import Fraction
import random
def v(x,p):
    if x==0: return 10**9
    x=abs(x); e=0
    while x%p==0: x//=p; e+=1
    return e
def is_sq_unit(u,p): return pow(u%p,(p-1)//2,p)==1
def shift(coeffs,z0,s,p):   # g(t) = f(z0 + p^s t), коэффициенты целые
    n=len(coeffs); out=[0]*n
    for i,c in enumerate(coeffs):          # f = sum c_i z^i
        # (z0 + p^s t)^i
        from math import comb
        for j in range(i+1):
            out[j]+= c*comb(i,j)*z0**(i-j)*(p**s)**j
    return out
def ball(coeffs,z0,s,p,depth=0):
    if depth>60: raise RuntimeError("глубина исчерпана — это НЕ ответ «точек нет»")
    g=shift(coeffs,z0,s,p); c0=g[0]
    if c0==0: return True
    e=v(c0,p); mu=min(v(c,p) for c in g[1:])
    if mu>e:
        return (e%2==0) and is_sq_unit(c0//p**e,p)
    # Гензель на корень: простой корень в шаре -> точка с w=0
    if v(g[1],p)<10**9 and e>2*v(g[1],p) and all(v(g[i],p)>=v(g[1],p) for i in range(1,len(g))): return True
    return any(ball(coeffs,z0+t*p**s,s+1,p,depth+1) for t in range(p))
def solvable(coeffs,p):
    rev=coeffs[::-1]              # z^4 f(1/z)
    return ball(coeffs,0,0,p) or ball(rev,0,1,p)
def Cprime(d,m,n):   # d w^2 = ...  ->  w^2 = d z^4 - 2a z^2 + c/d, умножаем на d^2: (dw)^2 = d^3 z^4 - ... ; удобнее: w'^2 = d*(d^2 z^4 - 2ad z^2 + c)
    a=m*m+n*n; c=(m*m-n*n)**2
    return [d*c, 0, -2*a*d*d, 0, d**3]
# Контроль лемм 1 и 2 на серии 1/ℓ
def primes(N):
    s=[1]*(N+1); s[0]=s[1]=0
    for i in range(2,int(N**.5)+1):
        if s[i]:
            for j in range(i*i,N+1,i): s[j]=0
    return [i for i in range(N+1) if s[i]]
def factor(x):
    out=set(); d=2
    while d*d<=x:
        while x%d==0: out.add(d); x//=d
        d+=1
    if x>1: out.add(x)
    return out
good=[l for l in primes(3000) if l>2 and l%8 in (3,5) and all(q%4==3 for q in factor(l*l-1) if q!=2)]
print("простые ℓ ≤ 3000 с условиями:", good)
viol=0; checks=0
for l in good[:12]:
    N=l*l-1; ps=sorted(factor(N))
    from itertools import combinations
    for r in range(1,len(ps)+1):
        for sub in combinations(ps,r):
            d=1
            for q in sub: d*=q
            oddq=[q for q in sub if q!=2]
            witness = oddq[0] if oddq else l          # лемма 1 при нечётном q | d, иначе лемма 2 при p = ℓ
            checks+=1
            if solvable(Cprime(d,1,l),witness): viol+=1; print("  НАРУШЕНИЕ леммы:", l, d, witness)
print(f"проверено (ℓ, d): {checks}, разрешимых там, где лемма запрещает: {viol}")

# ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ: решатель обязан находить точки там, где они есть.
pos=0; tot=0
for l in good[:12]:
    for pr in [q for q in primes(60) if q>2]:
        tot+=1
        if solvable(Cprime(1,1,l),pr): pos+=1          # d = 1: точка есть всегда (образ кручения)
print(f"положительный контроль d=1: найдено точек {pos} из {tot} (должно быть все)")
# и на «хороших» простых при случайном d
hits=0; n2=0
random.seed(3)
for _ in range(200):
    l=random.choice(good[:12]); d=random.choice([2,3,6,5,7,10,11]); pr=random.choice([q for q in primes(200) if q>2 and (l*l-1)*l*d%q!=0])
    n2+=1; hits+=solvable(Cprime(d,1,l),pr)
print(f"хорошие простые (не делят 2ℓ(ℓ²−1)d): разрешимо {hits} из {n2} (должно быть все)")
