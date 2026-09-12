from sage.all import *
from pathlib import Path
import json
out=Path(__file__).resolve().parent
r=json.loads((out/'new_bridge_verified.json').read_text());checks=[]
for row in r['elliptic_factors']:
 E=EllipticCurve(list(map(QQ,row['raw_ainvs'])));M=E.minimal_model();iso=E.isomorphism_to(M);delta=E.discriminant()/M.discriminant()
 checks.append({'delta_ratio':str(delta),'isomorphism':str(iso),'unit_scaling_primes':{str(p):int(delta.valuation(p)) for p in [3,5,13,37,1213,9397]}})
 assert all(delta.valuation(p)==0 for p in [3,5,13,37,1213,9397])
(out/'zero_height_model_checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(checks)
