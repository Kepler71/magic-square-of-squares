from sage.all import *
# Claude, 14.09: перепись метода шести клеток по наклонам, не закрытым критерием Z/2×Z/4 (criterion_proof/coverage_exact.json).
# Закрытие наклона: тройка различных a<b<c из {r,s-r,s,s+r}, фактор E+ или E- с PARI ellrank верхней границей 0,
# и ни одна рациональная точка (всё кручение) не даёт допустимого z≠0. Возобновляемо: census.jsonl.
import json, os, sys
from multiprocessing import Pool
cov=json.load(open('/home/kep/magicKube/criterion_proof/coverage_exact.json'))
todo=[(rng,sl) for rng in ['2-48','49-100','101-200'] for sl in cov[rng]['open']]
done=set()
if os.path.exists('census.jsonl'):
    for l in open('census.jsonl'): done.add(json.loads(l)['slope'])
todo=[t for t in todo if t[1] not in done]
def verdict(pol, kind, amax):
    R=pol.parent(); X=R.gen(); lc=pol.leading_coefficient()
    g=(pol(X/lc)*lc**2).monic() if lc!=1 else pol
    co=g.list(); E=EllipticCurve([0,co[2],0,co[1],co[0]])
    lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
    if hi!=0: return None,(lo,hi)
    for P in E.torsion_points():
        if P.is_zero(): continue
        t=P[0]/lc
        if kind=='+' and 0<t<QQ(1)/amax**2 and t.is_square(): return False,(0,0)
        if kind=='-' and t>amax**2 and t.is_square(): return False,(0,0)
    return True,(0,0)
def work(item):
    rng,sl=item; r,s=map(int,sl.split('/'))
    R=PolynomialRing(QQ,'X'); X=R.gen(); a_=None
    vals=sorted(set(v for v in [r,s-r,s,s+r] if v>0)); info=[]
    for a,b,c in Combinations(vals,3):
        for kind,pol in (('+',(1-QQ(a)**2*X)*(1-QQ(b)**2*X)*(1-QQ(c)**2*X)),('-',(X-QQ(a)**2)*(X-QQ(b)**2)*(X-QQ(c)**2))):
            try:
                ok,rk=verdict(pol,kind,QQ(max(a,b,c)))
            except Exception as e:
                ok,rk=None,str(e)[:80]
            info.append([[a,b,c],kind,str(rk)])
            if ok: return dict(range=rng,slope=sl,closed=True,by=[[a,b,c],kind],info=info)
    return dict(range=rng,slope=sl,closed=False,info=info)
if __name__=='__main__':
    out=open('census.jsonl','a'); n=0
    import multiprocessing as mp
    with mp.get_context('fork').Pool(12) as pool:
        for res in pool.imap_unordered(work,todo):
            out.write(json.dumps(res)+"\n"); out.flush(); n+=1
            if n%50==0: print(n,"/",len(todo),flush=True)
    print("готово",flush=True)
