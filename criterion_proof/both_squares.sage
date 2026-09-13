# Пункт 2: наклоны, выпавшие из-за «оба квадрата». Для каждого — строгий разбор конкретной кривой:
# ранг 0 по критерию (доказано), кручение — Sage (строго), все точки кручения переводятся в p и в 8 клеток.
import sys, json
sys.path.insert(0,'.')
from odd_criterion import proves_rank0
from math import isqrt
def sq(x): return x>=0 and isqrt(int(x))**2==int(x)
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
CLASSES = [  # (коэффициенты клеток как функции k, пара (m,n) как функция r,s)
 ('{1,k,1+k}',        lambda k:[1,k,1+k],            lambda r,s:(r,s)),
 ('{k,1-k,1}',        lambda k:[k,1-k,1],            lambda r,s:(r,s-r)),
 ('{1,-(1-k),k}',     lambda k:[1,-(1-k),k],         lambda r,s:(s,s-r)),
 ('{1,-(1+k),-k}',    lambda k:[1,-(1+k),-k],        lambda r,s:(s,r+s)),
 ('{1+k,-k,1}',       lambda k:[1+k,-k,1],           lambda r,s:(r,r+s)),
 ('{±1,±k}',          lambda k:[1,-1,k,-k],          lambda r,s:(s-r,r+s)),
 ('{±k,±(1-k)}',      lambda k:[k,-k,1-k,-(1-k)],    lambda r,s:(s,abs(2*r-s))),
 ('{±k,±(1+k)}',      lambda k:[k,-k,1+k,-(1+k)],    lambda r,s:(s,2*r+s)),
 ('{±1,±(1-k)}',      lambda k:[1,-1,1-k,-(1-k)],    lambda r,s:(r,2*s-r)),
 ('{±1,±(1+k)}',      lambda k:[1,-1,1+k,-(1+k)],    lambda r,s:(r,2*s+r)),
 ('{k,1+k,1-k,-k}',   lambda k:[k,1+k,1-k,-k],       lambda r,s:(s,2*r)),
 ('{1,1+k,-1,-(1-k)}',lambda k:[1,1+k,-1,-(1-k)],    lambda r,s:(r,2*s)),
]
def all_cells(k): return [1,-1,k,-k,1+k,-(1+k),1-k,-(1-k)]
def class_points(cs):
    """все рациональные p (конечные), дающие точки кручения кривой класса; кривая строится как в доказательстве"""
    R.<x>=QQ[]
    if len(cs)==3:
        f=prod(1+c*x for c in cs); A,B,C=f[3],f[2],f[1]
        E=EllipticCurve(QQ,[0,B,0,A*C,A^2]); back=lambda X: X/A
    else:
        c1=cs[0]; r=-1/c1
        g=c1*prod((1-c/c1)*x+c for c in cs[1:]); A,B,C,D=g[3],g[2],g[1],g[0]
        E=EllipticCurve(QQ,[0,B,0,A*C,A^2*D]); back=lambda X: (r + A/X) if X!=0 else None
    ps=[]
    for P in E.torsion_points():
        if P==E(0):
            ps.append(None if len(cs)==3 else -1/cs[0]); continue
        ps.append(back(P[0]))
    return E, ps
def is_nondeg_square_set(k,p):
    cells=[1+c*p for c in all_cells(k)]
    return all(c>0 and QQ(c).is_square() for c in cells) and len(set(cells+[1]))==9
restored=[]; pyth=[]; problems=[]
for s in range(2,201):
    for r in range(1,s):
        if gcd(r,s)!=1 or 2*r==s: continue
        k=QQ(r)/s
        rows=[(nm,cf(k),red(*pr(r,s))) for nm,cf,pr in CLASSES]
        valid=[(nm,cs,(m,n)) for nm,cs,(m,n) in rows if m>0 and m!=n and proves_rank0(m,n)]
        if not valid: continue
        if any(not(sq(m) and sq(n)) for _,_,(m,n) in valid): continue    # уже закрыт общей теоремой
        # все закрывающие пары — «оба квадрата»: строгий разбор первой из них
        nm,cs,(m,n)=valid[0]
        u,v=isqrt(m),isqrt(n); is_p = sq(u*u+v*v)
        E,ps=class_points([QQ(c) for c in cs])
        tors=E.torsion_subgroup().invariants()
        bad=[p for p in ps if p is not None and is_nondeg_square_set(k,p)]
        (pyth if is_p else restored).append((f"{r}/{s}",nm,(m,n),tors))
        if bad: problems.append((f"{r}/{s}",nm,bad))
print(f"выпавших наклонов разобрано: {len(restored)+len(pyth)}")
print(f"  непифагоровы (кручение должно быть Z/2xZ/4): {len(restored)}; кручение у них: {set(t for *_,t in restored)}")
print(f"  пифагоровы u²+v²=□: {len(pyth)}: {pyth[:6]}")
print(f"  точек кручения, дающих невырожденный квадрат из квадратов: {problems if problems else 0}")
json.dump({'nonpyth':[x[0] for x in restored],'pyth':[x[0] for x in pyth]},open('both_squares.json','w'))
