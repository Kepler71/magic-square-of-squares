import json, re, subprocess, time
from fractions import Fraction as Fr
from math import isqrt
secs=json.load(open('sections.json'))
log=open('../family/sections_batch.log',encoding='utf-8').read()
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
res=open('chab_all.jsonl','a')
for s in secs:
    b,h,n=s['bhn']; A,C=s['A'],s['C']
    m=re.search(rf"section \(b,h,n\) = \({b},{h},{n}\).*?\(2\) points with max\(p,q\) < 150: (\{{.*?\}})\n",log,re.S)
    known=m.group(1) if m else '?'
    for d3 in s['d3']:
        f=poly_coeffs(d3,A,C)
        pts=[]; H=8
        while not pts and H<=24:
            pts=search(f,H); H+=8
        if not pts:
            rec={'section':[b,h,n],'d3':d3,'status':'точка на якобиане не найдена (H<=24)'}
        else:
            polyF="+".join(f"({c})*x^{i}" for i,c in enumerate(f) if c)
            elts="\n".join(f"Q{i}:=J![x^2+({a})*x+({bb}), ({c})*x+({dd})];" for i,(a,bb,c,dd) in enumerate(pts))
            names=",".join(f"Q{i}" for i in range(len(pts)))
            code=f"""P<x>:=PolynomialRing(Rationals());
H:=HyperellipticCurve({polyF});
J:=Jacobian(H);
{elts}
inf:=[Q : Q in [{names}] | Order(Q) eq 0];
print "INF:", #inf;
if #inf gt 0 then S:=Chabauty(inf[1]); print "NPTS:", #S; print "TVALS:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; end if;
"""
            out=magma(code)
            rec={'section':[b,h,n],'d3':d3,'known_X_from_pipeline':known,'magma':out}
        res.write(json.dumps(rec,ensure_ascii=False)+'\n'); res.flush()
        print(json.dumps(rec,ensure_ascii=False)[:330],flush=True)
