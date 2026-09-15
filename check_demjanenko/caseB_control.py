# положительный контроль машинерии случая B: берём P = 3P0 и t, K = X(P+t)+X(P); решаем X(P'+t)+X(P') = K тем же способом — P должна найтись
from sage.all import *
T=[QQ(-451),QQ(-247),QQ(-43)]; L=prod(T)
R=PolynomialRing(QQ,'X'); X=R.gen(); g=R.prod(X+L/t for t in T); co=g.list()
E1=EllipticCurve([0,co[2],0,co[1],co[0]])
P0=[p for p in E1.saturation([E1(list(map(QQ,p))) for p in pari(E1).ellrank()[3]])[0] if p.order()==oo][0]
Rx=PolynomialRing(QQ,'xx'); xx=Rx.gen(); fx=co[0]+co[1]*xx+co[2]*xx**2+xx**3
E1R=PolynomialRing(QQ,['x','y']); x,y=E1R.gens()
ok=0; tot=0
for t in E1.torsion_points():
    if t.is_zero(): continue
    for m in (1,2,3):
        P=m*P0; K=(P+t)[0]+P[0]; tot+=1
        xt,yt=t[0],t[1]
        lam_num=y-yt; lam_den=x-xt
        Fnum=lam_num**2-(co[2]+x+xt)*lam_den**2 + x*lam_den**2 - K*lam_den**2
        Fp=PolynomialRing(PolynomialRing(QQ,'x'),'y')(Fnum); c=Fp.list()+[0,0,0]
        A=Rx(dict(PolynomialRing(QQ,'x')(c[0]).dict()))+Rx(dict(PolynomialRing(QQ,'x')(c[2]).dict()))*fx
        Bc=Rx(dict(PolynomialRing(QQ,'x')(c[1]).dict()))
        pol=A if Bc==0 else A**2-Bc**2*fx
        found=any(rt==P[0] for rt,_ in pol.roots(QQ))
        ok+=found
print('контроль случая B: найдено',ok,'из',tot)
