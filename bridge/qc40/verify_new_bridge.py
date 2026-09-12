from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from pathlib import Path
import json
out=Path(__file__).resolve().parent
R=PolynomialRing(QQ,'x');x=R.gen();K=R.fraction_field();t=(1+x)/(1-x);s=QQ(289)/2
f=(49+529*x*x)*(83521*x**4+63358*x*x+83521)
h=(225+64*t*t)*(64+225*t*t)*(s*(1+t*t)-240*t)
assert K(f)==(1-x)**6*h
assert gcd(f,f.derivative())==1
c=[f[i] for i in [0,2,4,6]];a0,a2,a4,a6=c
# E1: (X,Y)=(a6*x²,a6*y); E2: (X,Y)=(a0/x²,a0*y/x³).
E1=EllipticCurve([0,a4,0,a2*a6,a0*a6*a6]);E2=EllipticCurve([0,a2,0,a4*a0,a6*a0*a0])
assert a6*a6*f==(a6*x*x)**3+a4*(a6*x*x)**2+a2*a6*(a6*x*x)+a0*a6*a6
assert a0*a0*K(f)/x**6==(a0/x**2)**3+a2*(a0/x**2)**2+a4*a0*(a0/x**2)+a6*a0*a0
r={'f':str(f),'coefficients_a0_a2_a4_a6':list(map(str,c)),'identity_verified':True,'smooth_genus':2,'finite_known_point':['0','2023'],'infinity_y_over_x3':'6647','elliptic_factors':[]}
for E in [E1,E2]:
 M=E.minimal_model();v={'raw_ainvs':list(map(str,E.ainvs())),'minimal_ainvs':list(map(str,M.ainvs()))}
 try:
  alarm(40);ec=mwrank_EllipticCurve(list(map(int,M.ainvs())),verbose=False);ec.two_descent(verbose=False)
  v.update(rank=int(ec.rank()),rank_bound=int(ec.rank_bound()),certain=bool(ec.certain()),gens=[list(map(str,P)) for P in ec.gens()])
 except (Exception,AlarmInterrupt) as e:v['error']=type(e).__name__
 finally:cancel_alarm()
 r['elliptic_factors'].append(v);print(v,flush=True)
r['finite_field_checks']=[]
for p in prime_range(5,100):
 if f.discriminant()%p==0:continue
 fp=GF(p);chi=lambda v:0 if v==0 else (1 if v.is_square() else -1)
 count=p+1+sum(chi(fp(f(k))) for k in range(p))+chi(fp(a6))
 expect=p+1-E1.change_ring(fp).trace_of_frobenius()-E2.change_ring(fp).trace_of_frobenius()
 assert count==expect
 r['finite_field_checks'].append({'p':int(p),'points':int(count),'ordinary':all(E.change_ring(fp).trace_of_frobenius()%p!=0 for E in [E1,E2])})
assert f(0)==2023**2 and a6==6647**2
(out/'new_bridge_verified.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)
