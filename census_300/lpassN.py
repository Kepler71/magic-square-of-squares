import sys, json, itertools
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier, zvals, nine_ok
from c300 import e_pm
from sage.all import *
import multiprocessing as mp
def work(rec):
    sl=rec['slope']; r,s=map(int,sl.split('/'))
    for cand in rec['lo0']:
        try:
            if cand[0]=='g1':
                S=cand[1]; E,L,c0=weier([QQ(c) for c in S])
                if E.minimal_model().lseries().L_ratio()==0: continue
                zs=zvals(E,L,c0); bad=[str(z) for z in zs if z!=0 and nine_ok(z,r,s)]
                if bad: return dict(slope=sl,ALERT=bad)
                return dict(slope=sl,closed=True,by=cand+['L_ratio'])
            else:
                (a,b,c),kind=cand[1],cand[2]; E,lc=e_pm(a,b,c,kind)
                if E.minimal_model().lseries().L_ratio()==0: continue
                amax=QQ(max(a,b,c)); bad=False
                for P in E.torsion_points():
                    if P.is_zero(): continue
                    t=P[0]/lc
                    if kind=='+' and 0<t<1/amax**2 and t.is_square(): bad=True
                    if kind=='-' and t>amax**2 and t.is_square(): bad=True
                if not bad: return dict(slope=sl,closed=True,by=cand+['L_ratio'])
        except Exception as e:
            continue
    return dict(slope=sl,closed=False)
if __name__=='__main__':
    R=[json.loads(l) for l in open(sys.argv[1])]; R=[r for r in R if not r.get('closed')]
    print('открыто',len(R),[r['slope'] for r in R],flush=True)
    with mp.get_context('fork').Pool(12) as pool:
        res=list(pool.imap_unordered(work,R))
    json.dump(res,open(sys.argv[2],'w'))
    print('закрыто L:',[r['slope'] for r in res if r.get('closed')]); print('осталось:',[r['slope'] for r in res if not r.get('closed')]); print('ALERT',[r for r in res if 'ALERT' in r])
