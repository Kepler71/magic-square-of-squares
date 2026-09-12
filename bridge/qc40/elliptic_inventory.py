from sage.all import *
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from pathlib import Path
from itertools import combinations
import json,time
out=Path(__file__).resolve().parent
R=PolynomialRing(QQ,names=('t','y'));t,y=R.gens();s=QQ(289)/2
forms={'F0':225+64*t*t,'F4':s*(1+t*t),'F8':64+225*t*t,'L':s*(1+t*t)-240*t,'U':s*(1+t*t)+240*t}
rows=[]
for names in combinations(forms,2):
 E=Jacobian(y*y-prod(forms[n] for n in names)).minimal_model()
 rec={'names':names,'ainvs':list(map(str,E.ainvs())),'j':str(E.j_invariant())}
 try:
  alarm(12);rk=E.pari_curve().ellrank();rec['rank_interval']=[int(rk[0]),int(rk[1])];rec['points']=[list(map(str,P)) for P in rk[3]]
 except (Exception,AlarmInterrupt) as e:rec['error']=type(e).__name__
 finally:cancel_alarm()
 rows.append(rec);print(json.dumps(rec),flush=True);(out/'elliptic_inventory.json').write_text(json.dumps(rows,indent=2))
