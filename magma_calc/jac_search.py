# Поиск точек якобиана y^2 = f(x), deg f = 5 (модель Magma для (23,7,17), d3=34): u = x^2+ax+b, v = cx+d, v^2 ≡ f mod u.
from fractions import Fraction as Fr
from math import isqrt
import itertools, sys
f=[0,1040400,0,-97921,0,1156]          # f(x) = 1156x^5 - 97921x^3 + 1040400x  (коэф. по возрастанию)
def fmod(a,b):
    # остаток f mod (x^2+ax+b) = r1 x + r0; через рекуррентность x^k ≡ p_k x + q_k
    p,q=Fr(1),Fr(0)   # x^1
    P=[(Fr(0),Fr(1)),(Fr(1),Fr(0))]     # x^0 = 0x+1, x^1 = 1x+0
    for k in range(2,6):
        pk,qk=P[-1]
        # x^k = x*(pk x + qk) = pk x^2 + qk x = pk(-a x - b) + qk x
        P.append((qk - a*pk, -b*pk))
    r1=sum(f[k]*P[k][0] for k in range(6)); r0=sum(f[k]*P[k][1] for k in range(6))
    return r1,r0
def is_sq(q):
    q=Fr(q)
    if q<0: return None
    n,d=q.numerator,q.denominator
    rn,rd=isqrt(n),isqrt(d)
    return Fr(rn,rd) if rn*rn==n and rd*rd==d else None
found=set()
H=int(sys.argv[1]) if len(sys.argv)>1 else 60
vals=sorted({Fr(p,q) for q in range(1,H+1) for p in range(-H*H//4 if False else -4*H,4*H+1)})
for a in vals:
    for b in vals:
        D=a*a-4*b
        if D==0: continue
        r1,r0=fmod(a,b)
        # (a^2-4b) s^2 + (2 a r1 - 4 r0) s + r1^2 = 0
        A2,B2,C2=D,2*a*r1-4*r0,r1*r1
        disc=B2*B2-4*A2*C2
        sd=is_sq(disc)
        if sd is None: continue
        for s in {(-B2+sd)/(2*A2),(-B2-sd)/(2*A2)}:
            c=is_sq(s)
            if not c or c==0: continue
            d=(r1+a*s)/(2*c)
            found.add((a,b,c,d))
    if len(found)>=6: break
for t in sorted(found,key=lambda t: max(abs(x.numerator)+abs(x.denominator) for x in t))[:6]:
    print("u = x^2 + (%s)x + (%s), v = (%s)x + (%s)"%t)
print("найдено:",len(found))
