# Claude, 14.09: два последних открытых наклона s<=200. Кривые рода 2 из 5 клеток (c; a,b) и 6 клеток (две пары + две одиночные).
# RationalPointsGenus2 с флагом; при лимите калькулятора — пауза 30 мин.
import itertools, json, time
from fractions import Fraction as F
from chab_all_lib import magma

def curves(r,s):
    vals=[r,s-r,s,s+r]
    for c in vals:
        rest=[v for v in vals if v!=c]
        for a,b in itertools.combinations(rest,2):
            yield (f"(1+{c}*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2)", dict(singles=[c],pairs=[a,b]))
    for pa in itertools.combinations(vals,2):
        rest=[v for v in vals if v not in pa]
        for sg in (1,-1):
            yield (f"(1+{rest[0]}*z)*(1+({sg*rest[1]})*z)*(1-{pa[0]}^2*z^2)*(1-{pa[1]}^2*z^2)", dict(singles=[rest[0],sg*rest[1]],pairs=list(pa)))
from math import isqrt
def _sq(x): return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
def degenerate(zs,vals):
    # закрыто, если ни одно конечное z≠0 не делает все восемь нецентральных клеток квадратами
    r,sm,s,sp=vals
    for z in zs:
        if z=='inf' or F(z)==0: continue
        t=F(z)
        if all(_sq(1+e*t) for e in (s,-s,r,-r,sm,-sm,sp,-sp)): return False
    return True
import sys
SL=[tuple(map(int,x.split('/'))) for x in sys.argv[2:]]
out=open(sys.argv[1],'a')
for r,s in SL:
    vals=[r,s-r,s,s+r]; done=False
    for f,meta in curves(r,s):
        code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve({f});
pts,flag:=RationalPointsGenus2(C);
print "PROVEN:", flag; print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in pts];"""
        while True:
            t0=time.time(); res=magma(code); dt=time.time()-t0
            if 'too many connections' in res or res.strip()=='':
                print('лимит/пусто — пауза 1800 с',flush=True); time.sleep(1800); continue
            break
        rec=dict(slope=f"{r}/{s}",curve=f,**meta,out=res[:400])
        if 'PROVEN: true' in res:
            import re
            zs=[x.strip() for x in re.search(r"Z: \[(.*?)\]",res,re.S).group(1).split(',')]
            rec['closed']=degenerate(zs,vals)
            if not rec['closed']: rec['ALERT']='невырожденная точка'
        out.write(json.dumps(rec,ensure_ascii=False)+"\n"); out.flush()
        print(rec['slope'],meta,('ЗАКРЫТ' if rec.get('closed') else ('ALERT' if 'ALERT' in rec else res[:70].replace('\n',' '))),flush=True)
        time.sleep(max(0,30-dt))
        if rec.get('closed'): break
print('готово',flush=True)
