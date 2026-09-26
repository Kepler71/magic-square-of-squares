from fractions import Fraction as F
from itertools import product
from math import gcd
import json

LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]

def mul(xs):
    ans=F(1)
    for x in xs: ans*=x
    return ans

def vp(x,p):
    if x==0: raise ValueError('zero')
    a,b=x.numerator,x.denominator
    v=0
    while a%p==0:a//=p;v+=1
    while b%p==0:b//=p;v-=1
    return v,a,b

def is_p_square(x,p):
    if not x:return True
    v,a,b=vp(x,p)
    if v%2:return False
    modulus=8 if p==2 else p
    u=(a*pow(b,-1,modulus))%modulus
    return u==1 if p==2 else pow(u,(p-1)//2,p)==1

def cells(b,c):return [1+i*b+j*c for i in (-1,0,1) for j in (-1,0,1)]
def line_test(cs,p):return all(is_p_square(mul(cs[i] for i in row),p) for row in LINES)

out={'status':'finite exact controls of separately proved theorem','local_2_3':{},'support':{}}
for p in (2,3):
    vals=sorted({F(n,p**k) for k in range(4) for n in range(-45,46)})
    count=passes=0
    for b,c in product(vals,repeat=2):
        cs=cells(b,c)
        if not b or not c or 0 in cs:continue
        count+=1
        if not line_test(cs,p):continue
        passes+=1
        assert all(is_p_square(x,p) for x in cs),(p,b,c)
        assert vp(b,p)[0]>=(3 if p==2 else 1),(p,b,c)
        assert vp(c,p)[0]>=(3 if p==2 else 1),(p,b,c)
    out['local_2_3'][str(p)]={'tested':count,'passes':passes}

for p in (3,5,7,11,13,17,37):
    count=passes=nontrivial=0
    for r in range(1,14):
      for s in range(r+1,16):
        if gcd(r,s)!=1:continue
        for n in (-3,-2,-1,0,1,2):
          for u in range(1,2*p):
            if u%p==0:continue
            z=F(u)*F(p)**n
            cs=cells(r*z,s*z)
            if 0 in cs:continue
            count+=1
            if not line_test(cs,p):continue
            passes+=1
            T=(1-r*z)*(1+(r+s)*z)*(1-s*z)
            L=(1-r*z)*(1+(r-s)*z)*(1+s*z)
            vt=vp(T,p)[0]%2;vl=vp(L,p)[0]%2
            assert not vt or (r-s)%p==0,(p,r,s,z,vt,vl)
            assert not vl or (r+s)%p==0,(p,r,s,z,vt,vl)
            if vt or vl:
                nontrivial+=1
                assert p%4==1
    out['support'][str(p)]={'tested':count,'passes':passes,'nontrivial_valuation':nontrivial}

print(json.dumps(out,indent=2))
with open('support_controls.json','w') as f:json.dump(out,f,indent=2)
