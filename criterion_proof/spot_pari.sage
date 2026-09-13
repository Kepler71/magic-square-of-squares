import json, random, sys
sys.path.insert(0,'.')
from odd_criterion import proves_rank0
random.seed(int(11)); hits=[]; bad=[]
while len(hits)<300:
    n=random.randint(3,2000); m=random.randint(1,n-1)
    if gcd(m,n)!=1 or not proves_rank0(m,n): continue
    hits.append((m,n))
    E=EllipticCurve(QQ,[0,m^2+n^2,0,m^2*n^2,0]).minimal_model()
    r=pari(E.a_invariants()).ellinit().ellrank()
    if int(r[1])!=0: bad.append((m,n,r[:2]))
print(f"своя сверка с PARI: 300 случайных пар (n ≤ 2000), где нечётный критерий доказывает ранг 0; противоречий: {bad if bad else 0}")
