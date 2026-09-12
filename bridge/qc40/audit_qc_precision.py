from sage.all import *
from pathlib import Path
import json,sys
out=Path(__file__).resolve().parent
sys.path.insert(0,str(out/'g1_lift_pipeline'))
from nine_cells import validate
K=Qp(11,40)
def point(s):
 pieces=s[1:-1].split(' : ')
 x,y,z=[K(sage_eval(v,locals={'O':O})) for v in pieces]
 if z!=1:x,y=x/z,y/z**3
 if x.valuation()<0:return 'infinity_chart',1/x,y/x**3
 return 'finite',x,y
c=[4092529,47287151,37608911,44182609]
def signature(s):
 chart,x,y=point(s)
 coeff=c if chart=='finite' else c[::-1]
 residual=y*y-sum(coeff[i]*x**(2*i) for i in range(4))
 assert residual.valuation()>=10,(s,residual)
 assert x.valuation()>0
 return chart,int(x.lift()%11**8),int(y.lift()%11**8)
rs=[]
for name in ['qc_mixed_15_8_p11_n20_twelve.json','qc_mixed_15_8_p11_n30_twelve_check.json']:
 r=json.loads((out/name).read_text());assert r['status'].startswith('completed')
 sigs=[sorted(signature(s) for s in row) for row in r['extra_points_by_omega']]
 rat={s for row in r['rational_points_by_omega'] for s in row}
 assert rat=={'(0 : -2023 : 1)','(0 : 2023 : 1)','(1 : -6647 : 0)','(1 : 6647 : 0)'}
 assert len({s for row in sigs for s in row})==40
 rs.append(sigs)
assert rs[0]==rs[1]
lifts=[validate(15,8,t) for t in [-1,1]]
assert all(r['all_squares'] and not r['distinct'] and r['all_eight_sums_equal'] for r in lifts)
res={'precision_runs':[20,30],'extra_padic_points':40,'matched_modulus':'11^8','all_candidates_on_curve_to_at_least_precision':10,'known_rational_points':4,'all_known_lifts_degenerate':True,'complete_global_classification':False,'warning':'Precision stability and curve equations checked; not a full independent audit of QC zero isolation.','lifts':lifts}
(out/'qc_precision_audit.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='lifts'},indent=2))
