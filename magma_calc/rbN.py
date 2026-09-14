# лёгкий путь: RankBound(J) и при ранге ≤1 — Chabauty/Chabauty0; кривые из 5 клеток; печатаем все возвращаемые значения Chabauty
import sys, itertools, json, time, re
from fractions import Fraction as F
from chab_all_lib import magma
out=open(sys.argv[1],'a')
for sl in sys.argv[2:]:
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]; closed=False
    for c in vals:
        for a,b in itertools.combinations([v for v in vals if v!=c],2):
            code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+{c}*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2));
J:=Jacobian(C); rb:=RankBound(J); print "RB:", rb;
if rb eq 1 then pp:=Points(C : Bound:=1000); Q:=[q : q in [J!(R-S) : R in pp, S in pp | R ne S] | Order(q) eq 0];
  if #Q gt 0 then S,I:=Chabauty(Q[1]); print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; print "IDX:", I; end if; end if;
if rb eq 0 then S:=Chabauty0(J); print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; end if;"""
            while True:
                t0=time.time(); res=magma(code); dt=time.time()-t0
                if 'too many connections' in res or res.strip()=='':
                    print('лимит — пауза 1800 с',flush=True); time.sleep(1800); continue
                break
            m=re.search(r"Z: \[(.*?)\]",res,re.S); rec=dict(slope=sl,c=c,pair=[a,b],out=res[:500])
            if m:
                zs=[x.strip() for x in m.group(1).split(',')]
                from math import isqrt
                sq=lambda x: x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
                rec['closed']=not any(z!='inf' and F(z)!=0 and all(sq(1+e*F(z)) for e in (s,-s,r,-r,s-r,r-s,s+r,-s-r)) for z in zs)
                if not rec['closed']: rec['ALERT']=zs
            out.write(json.dumps(rec,ensure_ascii=False)+"\n"); out.flush()
            print(sl,c,a,b,res[:90].replace('\n',' | '),'ЗАКРЫТ' if rec.get('closed') else '',flush=True)
            time.sleep(max(0,30-dt))
            if rec.get('closed'): closed=True; break
        if closed: break
print('готово',flush=True)
