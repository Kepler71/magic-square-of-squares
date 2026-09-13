import sys
sys.path.insert(0,'.')
from selmer_exact import selmer
exec(preparse(open('both_squares.sage').read().split("restored=[]")[0]))
def rank0(m,n):
    S,T=selmer(int(m),int(n)); return len(S)*len(T)<=4
for (r,s) in [(9,16),(9,25),(16,25),(25,119),(25,144),(25,169),(119,144),(144,169)]:
    k=QQ(r)/s
    rows=[(nm,cf(k),red(*pr(r,s))) for nm,cf,pr in CLASSES]
    closing=[(nm,cs,mn) for nm,cs,mn in rows if mn[0]>0 and mn[0]!=mn[1] and rank0(*mn)]
    nonpyth=[c for c in closing if not (sq(c[2][0]) and sq(c[2][1]) and sq(isqrt(c[2][0])**2+isqrt(c[2][1])**2))]
    if nonpyth:
        print(f"{r}/{s}: закрыт и непифагоровой парой {nonpyth[0][2]} ({nonpyth[0][0]}) — поточечный разбор не нужен"); continue
    for nm,cs,mn in closing:
        E,ps=class_points([QQ(c) for c in cs])
        bad=[p for p in ps if p is not None and is_nondeg_square_set(k,p)]
        print(f"{r}/{s}: только пифагорова пара {mn} ({nm}), кручение {E.torsion_subgroup().invariants()}, точек {len(ps)}, невырожденных квадратов: {bad if bad else 0}")
