from fractions import Fraction as F
from math import gcd, isqrt

def v2(n):
    n=abs(n); e=0
    while n%2==0: n//=2; e+=1
    return e

def sq_Q2(x):
    """x rational != 0 -> is x a square in Q2*"""
    if x==0: return True
    if isinstance(x,int): x=F(x)
    a,b=x.numerator,x.denominator
    if a<0: return False if False else _sqQ2(a,b)   # -1 is not a square in Q2, handled below
    return _sqQ2(a,b)

def _sqQ2(a,b):
    if a<0: 
        # -1 is not a square in Q2 ; sign handled by residue mod 8 of odd part (negative odd part never == 1 mod 8 as integer? use mod 8 of a*b)
        pass
    e=v2(a)-v2(b)
    if e%2: return False
    ao=a>>v2(a); bo=b>>v2(b)
    return (ao*bo)%8==1

def isperfsq(x):
    if isinstance(x,int):
        return x>=0 and isqrt(x)**2==x
    return x>=0 and isperfsq(x.numerator) and isperfsq(x.denominator)

def cells(m,n,t):
    """nine cells of G1 as Fractions, t rational"""
    t=F(t); s=F(m*m+n*n,2)
    c0=m*m+n*n*t*t; c8=n*n+m*m*t*t; c4=s*(1+t*t)
    c1=(m*t+n)**2; c3=(m*t-n)**2; c5=(m+n*t)**2; c7=(m-n*t)**2
    c2=c4-2*m*n*t; c6=c4+2*m*n*t
    return [c0,c1,c2,c3,c4,c5,c6,c7,c8]

def magic_ok(c):
    S=3*c[4]
    lines=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
    return all(c[i]+c[j]+c[k]==S for i,j,k in lines)

print("=== 1. MAGICNESS of the nine cells (random rational t) ===")
import random
bad=0
for _ in range(400):
    m=random.randint(1,40); n=random.randint(1,40)
    if gcd(m,n)!=1 or m==n: continue
    t=F(random.randint(-30,30), random.randint(1,30))
    if not magic_ok(cells(m,n,t)): bad+=1
print("violations:",bad)
print("c1,c3,c5,c7 always perfect squares of rationals: ",
      all(isperfsq(cells(11,4,F(a,b))[i]) for a in range(-9,10) for b in range(1,9) for i in (1,3,5,7)))

print()
print("=== 2. t=1 identity  (the one-line refutation) ===")
for (m,n) in [(13,8),(15,8),(19,16),(16,5),(11,4),(7,2),(11,2),(5,3),(15,1)]:
    c=cells(m,n,1)
    print(f"(m,n)=({m},{n}): c0=c4=c8={c[0]}  odd part mod8={(int(c[0])>>v2(int(c[0])))%8}"
          f"  c2={c[2]}=({m-n})^2  c6={c[6]}=({m+n})^2  five-cells-Q2-square:",
          all(sq_Q2(c[i]) for i in (0,2,4,6,8)))
