# Claude, 14.09: для открытых наклонов — факторы E± с PARI нижней границей 0 (верхняя >0):
# L(E,1)/Ω ≠ 0 (Sage L_ratio) ⇒ rank 0 (Гросс–Загир–Колывагин); затем полный разбор кручения.
from sage.all import *
import json, multiprocessing as mp
R=PolynomialRing(QQ,'X'); X=R.gen()
def curve(a,b,c,kind):
    pol=(1-QQ(a)**2*X)*(1-QQ(b)**2*X)*(1-QQ(c)**2*X) if kind=='+' else (X-QQ(a)**2)*(X-QQ(b)**2)*(X-QQ(c)**2)
    lc=pol.leading_coefficient(); g=(pol(X/lc)*lc**2).monic() if lc!=1 else pol
    co=g.list(); return EllipticCurve([0,co[2],0,co[1],co[0]]), lc
def work(item):
    sl,cands=item; res=[]
    for (abc,kind) in cands:
        a,b,c=abc; E,lc=curve(a,b,c,kind); amax=QQ(max(abc))
        try:
            Lr=E.minimal_model().lseries().L_ratio()
        except Exception as e:
            res.append(dict(abc=abc,kind=kind,err=str(e)[:80])); continue
        if Lr==0:
            res.append(dict(abc=abc,kind=kind,L_ratio=0)); continue
        bad=False
        for P in E.torsion_points():
            if P.is_zero(): continue
            t=P[0]/lc
            if kind=='+' and 0<t<1/amax**2 and t.is_square(): bad=True
            if kind=='-' and t>amax**2 and t.is_square(): bad=True
        res.append(dict(abc=abc,kind=kind,L_ratio=str(Lr),closed=not bad))
        if not bad: return dict(slope=sl,closed=True,by=res[-1],tried=res)
    return dict(slope=sl,closed=False,tried=res)
if __name__=='__main__':
    items=[]
    for l in open('census.jsonl'):
        r=json.loads(l)
        if r['closed']: continue
        cands=[(x[0],x[1]) for x in r['info'] if x[2].startswith('(0,')]
        items.append((r['slope'],cands))
    out=open('kolyvagin.jsonl','w'); n=0
    with mp.get_context('fork').Pool(12) as pool:
        for res in pool.imap_unordered(work,items):
            out.write(json.dumps(res)+"\n"); out.flush(); n+=1
            if n%20==0: print(n,'/',len(items),flush=True)
    print('готово',flush=True)
