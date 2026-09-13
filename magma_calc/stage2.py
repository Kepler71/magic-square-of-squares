# Claude, 14.09: второй этап переписи. Наклоны, не закрытые ни критерием Z/2×Z/4, ни шестью клетками.
# Для наклона r/s: кривые C: y²=(1+cz)(1−a²z²)(1−b²z²), c — одно значение, a,b — два из трёх остальных (12 кривых).
# Один запрос Magma: RankBound(J); если ≤1 — Chabauty (при 0 — Chabauty0). Закрыт, если полный список точек вырожден.
import itertools, json, os, re, time
from fractions import Fraction as F
from chab_all_lib import magma
opens=[l.strip() for l in open('../census_six_cells/open_after_six.txt') if l.strip()]
skip={'3/86','79/110','11/142','48/163'}
import json as _j
skip|={_j.loads(l)['slope'] for l in open('../census_six_cells/kolyvagin.jsonl') if _j.loads(l)['closed']}
done=set()
if os.path.exists('stage2.jsonl'):
    for l in open('stage2.jsonl'): 
        r=json.loads(l); done.add((r['slope'],r['c'],tuple(r['pair'])))
closed=set()
if os.path.exists('stage2.jsonl'):
    for l in open('stage2.jsonl'):
        r=json.loads(l)
        if r.get('closed'): closed.add(r['slope'])
out=open('stage2.jsonl','a')
def degenerate(zs, vals):
    for z in zs:
        if z=='inf' or z==0: continue
        if not any(z in (F(1,v),F(-1,v)) for v in vals): return False
    return True
for sl in opens:
    if sl in skip or sl in closed: continue
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]
    for c in vals:
        rest=[v for v in vals if v!=c]
        for a,b in itertools.combinations(rest,2):
            if (sl,c,(a,b)) in done: continue
            code=f"""P<z>:=PolynomialRing(Rationals());
C:=HyperellipticCurve((1+{c}*z)*(1-{a}^2*z^2)*(1-{b}^2*z^2));
J:=Jacobian(C); rb:=RankBound(J); print "RB:", rb;
if rb eq 0 then S:=Chabauty0(J); print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; end if;
if rb eq 1 then
  pts:=Points(C : Bound:=1000);
  Q:=[q : q in [J!(R-S) : R in pts, S in pts | R ne S] | Order(q) eq 0];
  if #Q gt 0 then S:=Chabauty(Q[1]); print "Z:", [ (R[3] eq 0) select "inf" else Sprint(R[1]/R[3]) : R in S]; else print "NOPT"; end if;
end if;"""
            while True:
                t0=time.time(); res=magma(code); dt=time.time()-t0
                if 'too many connections' in res or res.strip()=='':
                    print('лимит/пусто — пауза 1800 с',flush=True); time.sleep(1800); continue
                break
            time.sleep(max(0,20-dt))
            rb=re.search(r"RB: (\d+)",res); zs=re.search(r"Z: \[(.*?)\]",res,re.S)
            rec=dict(slope=sl,c=c,pair=[a,b],rb=int(rb.group(1)) if rb else None,sec=round(dt,1))
            if zs:
                Z=[x.strip() for x in zs.group(1).split(',')]
                Zf=[x if x=='inf' else F(x) for x in Z]
                rec['z']=Z; rec['closed']=degenerate(Zf,vals)
                if not rec['closed']: rec['ALERT']='невырожденная точка — проверить девять клеток'
            else:
                rec['raw']=res[:200]
            out.write(json.dumps(rec,ensure_ascii=False)+"\n"); out.flush()
            print(sl,c,a,b,'RB',rec['rb'],'ЗАКРЫТ' if rec.get('closed') else ('ALERT' if 'ALERT' in rec else '-'),rec['sec'],flush=True)
            if rec.get('closed'): break
        else:
            continue
        break
print("готово",flush=True)
