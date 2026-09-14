# C_k: y² = (1−F_k z)(1−F_{k+1} z)(1−F_{k+2} z) — ранг по PARI ellrank и L(E,1)
from sage.all import *
F=[1,1]
while len(F)<70: F.append(F[-1]+F[-2])
X=PolynomialRing(QQ,'z').gen()
rows=[]
for k in range(1,50):
    f=(1-F[k]*X)*(1-F[k+1]*X)*(1-F[k+2]*X)
    if f.discriminant()==0: continue
    L=f.leading_coefficient(); h=(f(X/L)*L**2).monic(); co=h.list()
    E=EllipticCurve([0,co[2],0,co[1],co[0]])
    lo,hi=[int(t) for t in E.pari_curve().ellrank()[:2]]
    w=E.root_number()
    rows.append((k,F[k],F[k+1],F[k+2],lo,hi,int(w)))
    print(k,(F[k],F[k+1],F[k+2]),'ellrank',[lo,hi],'корн.число',w,E.torsion_subgroup().invariants(),flush=True)
