"""2-изогенный спуск для E_{m,n}: y^2 = x(x+m^2)(x+n^2), gcd(m,n)=1.
phi: E -> E' = Y^2 = X(X-(m+n)^2)(X-(m-n)^2).
S^phi(E)   = { d | mn (со знаком, бесквадратное): C_d : d w^2 = (d z^2 + m^2)(d z^2 + n^2) разрешимо всюду локально }
S^phihat(E') = { d | m^2-n^2 : C'_d : d w^2 = (d z^2 - (m+n)^2)(d z^2 - (m-n)^2) разрешимо всюду локально }
2^rank = |S^phi_true| * |S^phihat_true| / 4 <= |S^phi| |S^phihat| / 4.
Локальная разрешимость: точная проверка подъёмом Гензеля (классы вычетов по модулю p^k; класс решён,
когда v_p(g) + (1 или 3) <= k). Корни в Q_p (w=0) проверяются отдельно.
"""
from math import gcd
from itertools import combinations
from sympy import factorint

def vp(x,p):
    if x==0: return 10**9
    e=0
    while x%p==0: x//=p; e+=1
    return e

def is_unit_square(u,p):
    # u — p-адическая единица (целое, не делится на p)
    if p==2: return u%8==1
    return pow(u%p,(p-1)//2,p)==1

def poly_eval(c,z):
    r=0
    for a in reversed(c): r=r*z+a
    return r

def zp_square_value(c,p,z0,k,maxk):
    """Существует ли z = z0 mod p^k в Z_p с g(z) квадратом (включая 0)? g задан коэффициентами c (c[i] при z^i).
    Гарантируется отсутствие корней g в Q_p кроме проверенных заранее; maxk — предел глубины."""
    val=poly_eval(c,z0)
    if val==0: return True
    e=vp(val,p); need=1 if p!=2 else 3
    if e+need<=k:
        return e%2==0 and is_unit_square(val//p**e,p)
    if k>=maxk: raise RuntimeError('depth exceeded')
    step=p**k
    for t in range(p):
        if zp_square_value(c,p,z0+t*step,k+1,maxk): return True
    return False

def qp_soluble(c,p,has_qp_root,maxk=60):
    """g(z)=sum c[i] z^i (степень 4). Есть ли z в P^1(Q_p) с g(z) квадратом (проективно: t^4 g(1/t))?"""
    if has_qp_root: return True
    for z0 in range(p):
        if zp_square_value(c,p,z0,1,maxk): return True
    cr=list(reversed(c))     # t^4 g(1/t)
    # t в p Z_p
    if zp_square_value(cr,p,0,1,maxk): return True
    return False

def qp_is_square(a,p):
    """a ≠ 0 рационально (целое): квадрат в Q_p?"""
    if a==0: return True
    e=vp(a,p)
    if e%2: return False
    u=a//p**e
    if u<0 and p!=2:
        return is_unit_square(u%p,p)
    if p==2: return u%8==1
    return is_unit_square(u,p)

def real_soluble(c):
    # g чётная квартика c[4] z^4 + c[2] z^2 + c[0]: есть ли z с g(z) >= 0?
    if c[4]>0: return True
    A,B,C=c[4],c[2],c[0]
    if C>=0: return True
    D=B*B-4*A*C
    if D<0: return False
    # корни по Z=z^2: нужен Z>=0 с g>=0; при A<0 максимум в Z=-B/(2A); если он >=0 и достижим при Z>=0
    # достаточно: существует Z>=0 с A Z^2+B Z+C>=0  <=> (B>0 и D>=0) или (C>=0)
    return B>0

def squarefree_signed_divisors(N):
    ps=[p for p in factorint(abs(N))]
    out=[]
    for r in range(len(ps)+1):
        for sub in combinations(ps,r):
            d=1
            for p in sub: d*=p
            out.append(d); out.append(-d)
    return out

def selmer_phi(m,n):
    """S^phi(E): d | mn."""
    A2,B2=m*m,n*n
    res=[]
    primes=set(factorint(2*m*n*(n*n-m*m)))
    for d in squarefree_signed_divisors(m*n):
        c=[d*A2*B2,0,d*d*(A2+B2),0,d**3]   # g = d*(d z^2+m^2)(d z^2+n^2)
        if not real_soluble(c): continue
        ok=True
        for p in sorted(primes):
            has_root=qp_is_square(-d*A2,p) or qp_is_square(-d*B2,p)   # d z^2 = -m^2 или -n^2  <=> -d m^2 квадрат
            if not qp_soluble(c,p,has_root): ok=False; break
        if ok: res.append(d)
    return res

def selmer_phihat(m,n):
    """S^phihat(E'): d | (m+n)(m-n)."""
    A2,B2=(m+n)**2,(m-n)**2
    res=[]
    primes=set(factorint(2*m*n*(n*n-m*m)))
    for d in squarefree_signed_divisors((m+n)*(m-n)):
        c=[d*A2*B2,0,-d*d*(A2+B2),0,d**3]  # g = d*(d z^2-(m+n)^2)(d z^2-(m-n)^2)
        if not real_soluble(c): continue
        ok=True
        for p in sorted(primes):
            has_root=qp_is_square(d,p)
            if not qp_soluble(c,p,has_root): ok=False; break
        if ok: res.append(d)
    return res

def descent(m,n):
    S=selmer_phi(m,n); T=selmer_phihat(m,n)
    import math
    return {'S_phi':sorted(S,key=abs),'S_phihat':sorted(T,key=abs),
            'rank_upper':int(math.log2(len(S)*len(T)))-2,'pure_rank0':len(S)*len(T)==4}

if __name__=='__main__':
    import sys
    m,n=int(sys.argv[1]),int(sys.argv[2]); print(descent(m,n))
