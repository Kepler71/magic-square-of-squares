# Мелкие числовые утверждения: контроль E_34 (A,B,A+4B) из INDEX_AND_LOCAL_LIMIT §2; числа аудита ранга 2.
from sage.all import *
import json
from pathlib import Path
OUT = Path(__file__).resolve().parent
res = {}
E = EllipticCurve([-34**2, 0]); G = E(-16, 120); H = E(-2, 48)
res['E34_rank'] = int(E.rank()); res['E34_gens_saturated'] = str(E.saturation([G, H])[1])
A, B = 2*G, 2*H
tri = [A, B, A + 4*B]
res['E34_control'] = [dict(v2=int(P[0].valuation(2)), v3=int(P[0].valuation(3)), rootden=int(ZZ(P[0].denominator()).sqrt())) for P in tri]
res['E34_common_rootden'] = int(gcd([ZZ(P[0].denominator()).sqrt() for P in tri]))
xs = [P[0] for P in tri]
res['E34_control_is_AP_any_order'] = any(xs[(k+1)%3] + xs[(k+2)%3] == 2*xs[k] for k in range(3))
d = json.loads(Path('/home/kep/magicKube/joint_Eb_Ec/lattice2d/equal_area_rank_M1500.json').read_text())
rs = [r['det_ratio'] for r in d['R1_triples']]
res['equal_area'] = dict(summary=d['summary'] if not isinstance(d['summary'], dict) else {k: v for k, v in d['summary'].items() if k != 'dependent'},
                         n_R1=len(rs), n_shadow=sum(r['ap_shadow'] for r in d['R1_triples']), det_ratio_min=min(rs), det_ratio_max=max(rs))
(OUT / 'd_misc.json').write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str))
