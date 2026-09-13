# Проверка на всех 711 наклонах: минимальная модель каждого из 14 хороших классов совпадает с E_{m,n}
# для предсказанной пары (m,n) (после сокращения общего множителя).
import json
T=json.load(open('full_table_2_48.json')); rows=T['rows']
PAIRS={'C|1,k,1+k':lambda r,s:(r,s),'C|1,1-k,-k':lambda r,s:(r,s),'C|1,k,1-k':lambda r,s:(r,s-r),
 'C|1,k,-(1-k)':lambda r,s:(s,s-r),'C|1,-k,-(1+k)':lambda r,s:(s,r+s),'C|1,1+k,-k':lambda r,s:(r,r+s),
 'Q|1+k,1-k,-(1+k),-(1-k)':lambda r,s:(r,s),'Q|1,k,-1,-k':lambda r,s:(s-r,s+r),'Q|k,1-k,-k,-(1-k)':lambda r,s:(s,abs(2*r-s)),
 'Q|k,1+k,-k,-(1+k)':lambda r,s:(s,2*r+s),'Q|1,1-k,-1,-(1-k)':lambda r,s:(r,2*s-r),'Q|1,1+k,-1,-(1+k)':lambda r,s:(r,2*s+r),
 'Q|k,1+k,1-k,-k':lambda r,s:(s,2*r),'Q|1,1+k,-1,-(1-k)':lambda r,s:(r,2*s)}
def Ekey(m,n):
    g=gcd(m,n); m//=g; n//=g
    return ','.join(map(str,EllipticCurve([0,m*m+n*n,0,m*m*n*n,0]).minimal_model().ainvs()))
bad=0; tot=0
for ks,row in rows.items():
    k=QQ(ks); r,s=k.numerator(),k.denominator()
    for cls,f in PAIRS.items():
        key=row[cls]['key']
        if key is None: continue
        m,n=f(r,s); tot+=1
        if Ekey(m,n)!=key: bad+=1; print('MISMATCH',ks,cls,(m,n))
print('проверено',tot,'несовпадений',bad)
