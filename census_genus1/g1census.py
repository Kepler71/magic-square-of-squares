# Claude, 14.09: кривые рода 1 из 3 или 4 клеток 1+cz (c ∈ {±r,±s,±(s−r),±(s+r)}) для открытых наклонов.
# Ранг 0: PARI ellrank верхняя 0, иначе при нижней 0 — L(E,1)≠0 (L_ratio). Тогда все точки — кручение;
# каждая точка даёт z; наклон закрыт, если ни одно z≠0 не делает все девять клеток квадратами.
from sage.all import *
import json, os, itertools, multiprocessing as mp
def sq(x): return x>=0 and x.is_square()
def nine_ok(z,r,s):
    p=s*z; q=r*z
    return all(sq(1+e) for e in (p,-p,q,-q,p-q,q-p,p+q,-p-q))
def weier(S):
    c0=S[0]; others=S[1:]; R=PolynomialRing(QQ,'w'); w=R.gen()
    g=R.prod((c0-c)*w + c*c0 for c in others)
    if len(S)==3: g=c0*w*g
    L=g.leading_coefficient(); X=w
    h=(g(X/L)*L**2).monic()
    co=h.list()
    return EllipticCurve([0,co[2],0,co[1],co[0]]), L, c0
def zvals(E,L,c0):
    zs=[QQ(-1)/c0]                       # O ↔ w=∞
    for P in E.torsion_points():
        if P.is_zero(): continue
        wv=P[0]/L
        if wv==0: continue               # z=∞
        zs.append(QQ(-1)/c0 + 1/wv)
    return zs
def work(sl):
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]
    C=sorted(set([v for v in vals]+[-v for v in vals]))
    tried=[]
    for k in (3,4):
        for S in itertools.combinations(C,k):
            if -S[0] in S and k==3: pass
            # симметрия z→−z: берём подмножество с лексикографически меньшим из S и −S
            neg=tuple(sorted(-c for c in S))
            if neg<S: continue
            try:
                E,L,c0=weier([QQ(c) for c in S])
            except Exception as e:
                continue
            if E.discriminant()==0: continue
            lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
            how=None
            if hi==0: how='ellrank'
            elif lo==0:
                try:
                    if E.minimal_model().lseries().L_ratio()!=0: how='L_ratio'
                except Exception: pass
            tried.append([list(S),lo,hi,how])
            if how:
                zs=zvals(E,L,c0)
                bad=[str(z) for z in zs if z!=0 and nine_ok(z,r,s)]
                if bad: return dict(slope=sl,ALERT=bad,S=list(S))
                return dict(slope=sl,closed=True,by=[list(S),how],ntors=len(zs))
    return dict(slope=sl,closed=False,tried=tried)
if __name__=='__main__':
    opens=[l.strip() for l in open('../census_six_cells/open_after_six.txt') if l.strip()]
    closed=set(json.loads(l)['slope'] for l in open('../census_six_cells/kolyvagin.jsonl') if json.loads(l)['closed'])
    closed|={'3/86','79/110','11/142','48/163'}
    if os.path.exists('../magma_calc/stage2.jsonl'):
        closed|=set(json.loads(l)['slope'] for l in open('../magma_calc/stage2.jsonl') if json.loads(l).get('closed'))
    todo=[s for s in opens if s not in closed]
    print('открыто к проверке:',len(todo),flush=True)
    out=open('g1census.jsonl','w'); n=0
    with mp.get_context('fork').Pool(12) as pool:
        for res in pool.imap_unordered(work,todo):
            out.write(json.dumps(res)+"\n"); out.flush(); n+=1
            if 'ALERT' in res: print('ALERT',res,flush=True)
            if n%10==0: print(n,'/',len(todo),flush=True)
    print('готово',flush=True)
