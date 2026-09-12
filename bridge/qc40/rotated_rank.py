from sage.all import *
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from pathlib import Path
import json
m,n=23,7;s=QQ(m*m+n*n)/2;b=s*m*m*n*n;a=QQ(m*m)/n**2+1+QQ(n*n)/m**2
E=EllipticCurve([0,b*a,0,b*b*a,b**3]).minimal_model();r={'pair':[m,n],'ainvs':list(map(str,E.ainvs()))}
try:
 alarm(15);rr=E.pari_curve().ellrank();r['pari_rank_interval']=[int(rr[0]),int(rr[1])];r['points']=[list(map(str,P)) for P in rr[3]]
except (Exception,AlarmInterrupt) as e:r['error']=type(e).__name__
finally:cancel_alarm()
print(json.dumps(r),flush=True);(Path(__file__).parent/'rotated_rank.json').write_text(json.dumps(r,indent=2))
