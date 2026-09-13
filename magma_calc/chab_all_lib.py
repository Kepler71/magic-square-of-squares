import json, re, subprocess, time
from fractions import Fraction as Fr
from math import isqrt
def poly_coeffs(d3,A,C):
    # d3 * t (t^2-1)(A t - C)(C t - A), коэффициенты по возрастанию
    import numpy as np
    p=[0,1]                                     # t
    def mul(a,b):
        r=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            for j,y in enumerate(b): r[i+j]+=x*y
        return r
    for q in ([-1,0,1],[-C,A],[-A,C]): p=mul(p,q)
    return [d3*c for c in p]
def fmod(f,a,b):
    P=[(Fr(0),Fr(1)),(Fr(1),Fr(0))]
    for k in range(2,len(f)):
        pk,qk=P[-1]; P.append((qk-a*pk,-b*pk))
    return sum(f[k]*P[k][0] for k in range(len(f))), sum(f[k]*P[k][1] for k in range(len(f)))
def is_sq(q):
    if q<0: return None
    n,d=q.numerator,q.denominator; rn,rd=isqrt(n),isqrt(d)
    return Fr(rn,rd) if rn*rn==n and rd*rd==d else None
def search(f,H,limit=4,deadline=90):
    t0=time.time(); vals=sorted({Fr(p,q) for q in range(1,H+1) for p in range(-4*H,4*H+1)})
    out=[]
    for a in vals:
        if time.time()-t0>deadline: break
        for b in vals:
            D=a*a-4*b
            if D==0: continue
            r1,r0=fmod(f,a,b); A2,B2,C2=D,2*a*r1-4*r0,r1*r1
            sd=is_sq(B2*B2-4*A2*C2)
            if sd is None: continue
            for s in {(-B2+sd)/(2*A2),(-B2-sd)/(2*A2)}:
                c=is_sq(s)
                if c: out.append((a,b,c,(r1+a*s)/(2*c)))
        if len(out)>=limit: break
    return out
def magma(code):
    open('_q.m','w').write(code)
    r=subprocess.run(['curl','-s','-m','150','-A','Mozilla/5.0','--data-urlencode','input@_q.m','https://magma.maths.usyd.edu.au/xml/calculator.xml'],capture_output=True,text=True)
    body=re.sub(r'<[^>]*>','\n',r.stdout)
    return '\n'.join(l for l in body.split('\n') if l.strip() and not re.fullmatch(r'(60|50000|\d{9,}|2\.29-10|[\d.]+|[\d.]+MB)',l.strip()))
