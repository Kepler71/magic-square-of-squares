from sage.all import *
T=[QQ(-451),QQ(-247),QQ(-43)]; L=prod(T)
R=PolynomialRing(QQ,'X'); X=R.gen(); g=R.prod(X+L/t for t in T); co=g.list()
E1=EllipticCurve([0,co[2],0,co[1],co[0]]); Emin=E1.minimal_model(); phi=E1.isomorphism_to(Emin); psi=Emin.isomorphism_to(E1)
P0=[p for p in Emin.saturation([Emin(list(map(QQ,p))) for p in pari(Emin).ellrank()[3]])[0] if p.order()==oo][0]
B=Emin.silverman_height_bound(); mx=0; mz=0
def h(q): return log(max(abs(q.numerator()),abs(q.denominator()))) if q!=0 else 0
for n in range(1,16):
    for t in Emin.torsion_points():
        Q=n*P0+t
        if Q.is_zero(): continue
        mx=max(mx,abs(Q.height()-h(Q[0])))
        z=psi(Q)[0]/L; mz=max(mz,abs(h(z)-h(Q[0])))
print('max|ĥ−h(x)| =',float(mx),'при B =',float(B),'; max|h(z)−h(x)| =',float(mz))
