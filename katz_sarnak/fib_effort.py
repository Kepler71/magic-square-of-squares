from sage.all import *
F=[1,1]
while len(F)<70: F.append(F[-1]+F[-2])
X=PolynomialRing(QQ,'z').gen()
for k in [30,41,43,45,47]:
    f=(1-F[k]*X)*(1-F[k+1]*X)*(1-F[k+2]*X)
    L=f.leading_coefficient(); h=(f(X/L)*L**2).monic(); co=h.list()
    E=EllipticCurve([0,co[2],0,co[1],co[0]])
    res=None
    for eff in (2,4,6):
        res=pari(E).ellrank(eff)[:2]
        if res[0]==res[1]: break
    print(k,'ellrank усилие',eff,res,flush=True)
