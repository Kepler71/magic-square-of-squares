# Claude, 15.09: независимая реализация проверки Демьяненко–Манина (теорема 7.1 Fable) на кубических множителях.
from sage.all import *
import sys, json
def sq(q): return q>=0 and q.is_square()
def full_ok(z,S): return all(sq(1+l*z) for l in S)
def check(r,s,T,SAFETY=2):
    S=[s,-s,r,-r,s-r,r-s,s+r,-s-r]
    T=[QQ(t) for t in T]; L=prod(T)
    R=PolynomialRing(QQ,'X'); X=R.gen()
    g=R.prod(X+L/t for t in T)                 # Y² = g(X), X = L z, Y = L y
    co=g.list(); E1=EllipticCurve([0,co[2],0,co[1],co[0]])
    lo,hi=[int(v) for v in pari(E1).ellrank()[:2]]
    if (lo,hi)!=(1,1): return dict(T=[int(t) for t in T],skip=f'ранг {lo,hi}')
    Emin=E1.minimal_model(); phi=E1.isomorphism_to(Emin); psi=Emin.isomorphism_to(E1)
    pts=[phi(E1(p)) for p in [] ]
    G=Emin.gens(proof=True) if False else None
    # образующая: точки из ellrank + насыщение
    gp_pts=pari(Emin).ellrank()[3]
    P=[Emin(list(map(QQ,p))) for p in gp_pts]
    P=[p for p in P if p.order()==oo]
    sat,idx,reg=Emin.saturation(P)
    P0=[p for p in sat if p.order()==oo]
    assert len(P0)>=1
    P0=P0[0]
    tors=Emin.torsion_points()
    h0=P0.height()
    B=Emin.silverman_height_bound()
    # z как функция x_min: z = X/L, X = u² x_min + r  (изоморфизм Emin→E1: psi)
    u,rr,ss,tt=psi.tuple()          # E1 = psi(Emin): X = u²x + r
    alpha=u**2/L; beta=rr/L
    # целочисленная матрица z = (a x + b)/(0·x + d)
    den=lcm(alpha.denominator(),beta.denominator())
    a=alpha*den; b=beta*den; d=den
    # обратная: x = (d z − b)/a → матрица [[d,-b],[0,a]]
    c1=log(2*max(abs(a),abs(b),abs(d)))
    C=SAFETY*2*B+2*c1
    M0=floor((C/h0+1)/2)+1
    bad=[]; ncand=0
    for n in range(-M0,M0+1):
        for t in tors:
            Q=n*P0+t
            if Q.is_zero(): continue
            Xq=psi(Q)[0]; z=Xq/L
            ncand+=1
            if z!=0 and (full_ok(z,S) or full_ok(-z,S)): bad.append(str(z))
    # случай B: z(P+t) = −z(P), t ≠ O, на E1 (a1=a3=0): X(P+t) + X(P) = 0  (т.к. z = X/L)
    caseB=[]
    E1R=PolynomialRing(QQ,['x','y']); x,y=E1R.gens()
    f=co[0]+co[1]*x+co[2]*x**2+x**3
    for t in E1.torsion_points():
        if t.is_zero(): continue
        xt,yt=t[0],t[1]
        if yt==0:
            # X(P+t) = xt + (f'(xt))... используем формулу сложения с λ=(y)/(x-xt)
            pass
        lam_num=y-yt; lam_den=x-xt
        x3_num=lam_num**2 - (co[2]+x+xt)*lam_den**2      # x3 = λ² − a2 − x − xt
        Fnum=x3_num + x*lam_den**2                       # X(P+t)+X(P) = 0 ⇒ числитель
        # заменить y² на f
        Fnum=Fnum.reduce([y**2-f]) if False else Fnum
        Fp=PolynomialRing(PolynomialRing(QQ,'x'),'y')(Fnum)
        # y² → f
        c=Fp.list()+[0,0,0]
        A=c[0]+c[2]*f.polynomial(x) if False else None
        Rx=PolynomialRing(QQ,'xx'); xx=Rx.gen()
        fx=co[0]+co[1]*xx+co[2]*xx**2+xx**3
        coeffs=[Rx(ci.subs({}) if hasattr(ci,'subs') else ci) for ci in [c[0],c[1],c[2]]]
        def tox(pol): return Rx(PolynomialRing(QQ,'x')(pol).subs())(xx) if pol!=0 else Rx(0)
        A=Rx(dict(PolynomialRing(QQ,'x')(c[0]).dict())) + Rx(dict(PolynomialRing(QQ,'x')(c[2]).dict()))*fx
        Bc=Rx(dict(PolynomialRing(QQ,'x')(c[1]).dict()))
        if Bc==0: polys=[A]
        else: polys=[A**2-Bc**2*fx]
        for pol in polys:
            if pol==0: caseB.append(('ТОЖДЕСТВЕННЫЙ НОЛЬ',str(t))); continue
            for rt,_ in pol.roots(QQ):
                if rt==xt: continue
                for Pp in E1.lift_x(rt,all=True):
                    Q=Pp+t
                    if Q.is_zero(): continue
                    if Q[0]+Pp[0]==0:
                        z=Pp[0]/L
                        caseB.append(str(z))
                        if z!=0 and full_ok(z,S): bad.append('B:'+str(z))
    return dict(T=[int(t) for t in T],rank=(lo,hi),hP0=float(h0),B=float(B),c1=float(c1),M0=int(M0),
                cand=ncand,caseB=caseB,nondegenerate=bad,tors=len(tors),sat_index=int(idx))
if __name__=='__main__':
    r,s=map(int,sys.argv[1].split('/'))
    T=list(map(int,sys.argv[2].split(',')))
    print(json.dumps(check(r,s,T),ensure_ascii=False))
