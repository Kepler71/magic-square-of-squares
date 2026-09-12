from sage.all import *
from pathlib import Path
import json
out=Path(__file__).resolve().parent;r=json.loads((out/'new_bridge_verified.json').read_text());rows=[]
for row in r['elliptic_factors']:
 E=EllipticCurve(list(map(QQ,row['raw_ainvs'])))
 rec={str(p):{'vDelta':int(E.discriminant().valuation(p)),'vc4':int(E.c4().valuation(p)),'kodaira':str(E.kodaira_symbol(p)),'tamagawa':int(E.tamagawa_number(p))} for p in [2,7,17,23]}
 rows.append(rec)
print(rows);(out/'local_model_data.json').write_text(json.dumps(rows,indent=2))
