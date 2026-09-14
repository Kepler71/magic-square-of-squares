# Семейство E_k: y² = (1−rz)(1−(s−r)z)(1−sz), k=r/s (тройка a, b, a+b; общий ранг над ℚ(k) равен 0).
# Для всех наклонов s ≤ SMAX: корневое число w и PARI ellrank [lo,hi]. Проверка «минималистской» гипотезы / Каца–Сарнака:
# доля w=+1 → 1/2; при w=+1 ранг 0, при w=−1 ранг 1 почти всегда; ранг ≥ 2 — доля, убывающая с ростом высоты.
from sage.all import *
import json, sys, multiprocessing as mp
from math import gcd
X=PolynomialRing(QQ,'z').gen()
def one(rs):
    r,s=rs
    f=(1-r*X)*(1-(s-r)*X)*(1-s*X)
    if f.discriminant()==0: return None
    L=f.leading_coefficient(); h=(f(X/L)*L**2).monic(); co=h.list()
    E=EllipticCurve([0,co[2],0,co[1],co[0]])
    try:
        lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
    except Exception: lo,hi=None,None
    return dict(r=r,s=s,w=int(E.root_number()),lo=lo,hi=hi)
if __name__=='__main__':
    SMAX=int(sys.argv[1])
    items=[(r,s) for s in range(3,SMAX+1) for r in range(1,s) if gcd(r,s)==1 and 2*r!=s]
    out=[]
    with mp.get_context('fork').Pool(8) as pool:
        for x in pool.imap_unordered(one,items,chunksize=8):
            if x: out.append(x)
    json.dump(out,open(f'family_{SMAX}.json','w'))
    print('кривых',len(out))
