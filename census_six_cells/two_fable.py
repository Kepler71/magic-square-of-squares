from sage.all import *
R=PolynomialRing(QQ,'X'); X=R.gen()
def curve(a,b,c,kind):
    pol=(1-QQ(a)**2*X)*(1-QQ(b)**2*X)*(1-QQ(c)**2*X) if kind=='+' else (X-QQ(a)**2)*(X-QQ(b)**2)*(X-QQ(c)**2)
    lc=pol.leading_coefficient(); g=(pol(X/lc)*lc**2).monic() if lc!=1 else pol
    co=g.list(); return EllipticCurve([0,co[2],0,co[1],co[0]]), lc
for (a,b,c,kind) in [(43,131,174,'+'),(2,195,197,'+'),(2,195,392,'-')]:
    E,lc=curve(a,b,c,kind); Em=E.minimal_model()
    print((a,b,c,kind),'кондуктор',Em.conductor().factor(),'аналит. ранг',Em.analytic_rank(),'кручение',E.torsion_subgroup().invariants(), 'L(E,1)', Em.lseries().L_ratio(), flush=True)
