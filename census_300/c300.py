# Claude, 14.09: перепись 200 < s ≤ 300. Для наклона r/s по порядку: кривые рода 1 из 3 клеток, из 4 клеток, шесть клеток E±.
# Закрытие — PARI ellrank верхняя 0 + все точки кручения не дают невырожденных девяти квадратов. L(E,1) — отдельным проходом.
import sys, json, os, itertools
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier, zvals, nine_ok
from sage.all import *
from math import gcd
def e_pm(a,b,c,kind):
    X=PolynomialRing(QQ,'X').gen()
    pol=(1-QQ(a)**2*X)*(1-QQ(b)**2*X)*(1-QQ(c)**2*X) if kind=='+' else (X-QQ(a)**2)*(X-QQ(b)**2)*(X-QQ(c)**2)
    lc=pol.leading_coefficient(); g=(pol(X/lc)*lc**2).monic() if lc!=1 else pol
    co=g.list(); return EllipticCurve([0,co[2],0,co[1],co[0]]), lc
def work(sl):
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]
    C=sorted(set(vals+[-v for v in vals])); lo0=[]
    for k in (3,4):
        for S in itertools.combinations(C,k):
            if tuple(sorted(-c for c in S))<S: continue
            E,L,c0=weier([QQ(c) for c in S])
            if E.discriminant()==0: continue
            lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
            if hi==0:
                zs=zvals(E,L,c0); bad=[str(z) for z in zs if z!=0 and nine_ok(z,r,s)]
                if bad: return dict(slope=sl,ALERT=bad,S=list(S))
                return dict(slope=sl,closed=True,by=['g1',list(S)])
            if lo==0: lo0.append(['g1',list(S)])
    for a,b,c in itertools.combinations(sorted(set(v for v in vals if v>0)),3):
        for kind in '+-':
            E,lc=e_pm(a,b,c,kind); lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
            if hi==0:
                amax=QQ(max(a,b,c)); bad=False
                for P in E.torsion_points():
                    if P.is_zero(): continue
                    t=P[0]/lc
                    if kind=='+' and 0<t<1/amax**2 and t.is_square(): bad=True
                    if kind=='-' and t>amax**2 and t.is_square(): bad=True
                if not bad: return dict(slope=sl,closed=True,by=['six',[a,b,c],kind])
            if lo==0: lo0.append(['six',[a,b,c],kind])
    return dict(slope=sl,closed=False,lo0=lo0)
if __name__=='__main__':
    import multiprocessing as mp
    SMIN,SMAX=int(sys.argv[1]),int(sys.argv[2])
    todo=[f"{r}/{s}" for s in range(SMIN,SMAX+1) for r in range(1,s) if gcd(r,s)==1]
    fn=f'c{SMIN}_{SMAX}.jsonl'; done=set()
    if os.path.exists(fn): done={json.loads(l)['slope'] for l in open(fn)}
    todo=[t for t in todo if t not in done]
    print('к расчёту',len(todo),flush=True)
    out=open(fn,'a'); n=0; nc=0
    with mp.get_context('fork').Pool(10) as pool:
        for res in pool.imap_unordered(work,todo,chunksize=4):
            out.write(json.dumps(res)+"\n"); out.flush(); n+=1; nc+=bool(res.get('closed'))
            if 'ALERT' in res: print('ALERT',res,flush=True)
            if n%500==0: print(n,'/',len(todo),'закрыто',nc,flush=True)
    print('готово',n,'закрыто',nc,flush=True)
