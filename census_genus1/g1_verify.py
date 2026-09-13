from sage.all import *
import json, random as _r
R=[json.loads(l) for l in open('g1census.jsonl') if json.loads(l).get('closed')]
_r.seed(7); sample=_r.sample(R,12)
def sq(x): return x>=0 and x.is_square()
for rec in sample:
    sl=rec['slope']; S,how=rec['by']; r,s=map(int,sl.split('/'))
    P=PolynomialRing(QQ,'z'); z=P.gen()
    f=P.prod(1+c*z for c in S)
    if len(S)==3:
        L=f.leading_coefficient(); h=(f(z/L)*L**2).monic(); co=h.list()
        E=EllipticCurve([0,co[2],0,co[1],co[0]])
        rk=E.minimal_model().lseries().L_ratio()
        pts=[P_[0]/L for P_ in E.torsion_points() if not P_.is_zero()]
        bad=[t for t in pts if t!=0 and all(sq(1+e) for e in (s*t,-s*t,r*t,-r*t,(s-r)*t,(r-s)*t,(s+r)*t,-(s+r)*t))]
        print(sl,S,how,'L/Ω:',rk,'кручение:',E.torsion_subgroup().invariants(),'z:',sorted(pts),'плохих:',bad,flush=True)
    else:
        print(sl,S,how,'(4 клетки — пропуск прямой модели)')
