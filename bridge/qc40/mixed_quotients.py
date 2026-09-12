from sage.all import *
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from pathlib import Path
import json
out=Path(__file__).resolve().parent
R=PolynomialRing(QQ,names=('u','v'));u,v=R.gens();a=QQ(120);s=QQ(289)/2;d=QQ(25921);c=QQ(54721)
qs={
 'H_F0F8L_plus':(u+2)*(a*a*u*u+d)*(s*u-2*a),
 'H_F0F8L_minus':(u-2)*(a*a*u*u+d)*(s*u-2*a),
 'K_even':(a*a*u*u+c*u+a*a)*(s*s*u*u+(2*s*s-4*a*a)*u+s*s),
 'K_reciprocal':(a*a*u*u+d)*(s*s*u*u-4*a*a),
 'K_antireciprocal':(a*a*u*u+289**2)*(s*s*u*u+d)}
rows=[]
for name,f in qs.items():
 E=Jacobian(v*v-f).minimal_model();r={'name':name,'quartic':str(f),'ainvs':list(map(str,E.ainvs())),'j':str(E.j_invariant())}
 try:
  alarm(20);rk=E.pari_curve().ellrank();r['rank_interval']=[int(rk[0]),int(rk[1])];r['points']=[list(map(str,P)) for P in rk[3]]
 except (Exception,AlarmInterrupt) as e:r['error']=type(e).__name__
 finally:cancel_alarm()
 rows.append(r);print(json.dumps(r),flush=True);(out/'mixed_quotients.json').write_text(json.dumps(rows,indent=2))
