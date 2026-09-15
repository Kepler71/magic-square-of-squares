# диагностика контроля (9,22,31), p=7: rho[0] по дискам (код) против Σ_q ω_q(z) (мой omega_direct) в известных рациональных точках
import sys; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
import qc_load; from qc_load import *
from omega_direct import *; from models import fpoly
qc_load.DEBUG_RHO=True
f=fpoly(9,22,31); p=ZZ(7); n=ZZ(20); K=Qp(p,n)
Om,per=omega_direct(f,p,n,verbose=False)
a6,a4,a2,a0=f[6],f[4],f[2],f[0]
E1=EllipticCurve([0,a4,0,a2*a6,a0*a6**2]); E2=EllipticCurve([0,a2,0,a0*a4,a0**2*a6]); E1min=E1.minimal_model(); E2min=E2.minimal_model()
i1=E1.isomorphism_to(E1min); i2=E2.isomorphism_to(E2min)
def omega_at(x):  # x ∈ Q, точка (x, w) с w^2=f(x)
    tot=K(0); parts={}
    for q in sorted(per):
        Kq=Qp(q,60); xq=Kq(x); fx=Kq(f(x))
        y=Kq(0) if fx==0 else fx.sqrt()
        r1=lam_r(E1min,*to_min(i1,Kq(a6)*xq**2,Kq(a6)*y),q) if x!=0 else None
        r2=lam_r(E2min,*to_min(i2,Kq(a0)/xq**2,Kq(a0)*y/xq**3),q)
        c=-r1+r2-2*xq.valuation(); parts[q]=c; tot+=c*K(q).log()
        assert c in per[q], (q,c,per[q])
    return tot,parts
for x in [9,22,31,QQ(31)/2]:
    tot,parts=omega_at(QQ(x)); print('z: x =',x,'Σω =',tot,'коэффициенты',parts)
print('Ω содержит Σω?', [any(tot==w for w in Om) for x in [9,22,31,QQ(31)/2] for tot in [omega_at(QQ(x))[0]]])
rat,other=quadratic_chabauty_bielliptic(f,p,n,Omega=Om,gens=None)
