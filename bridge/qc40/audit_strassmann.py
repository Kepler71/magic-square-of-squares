from sage.all import *
from pathlib import Path
import json
out=Path(__file__).resolve().parent;K=Qp(11,40);Rp=PolynomialRing(GF(11),'z');z=Rp.gen()
fn=json.loads((out/'qc_mixed_15_8_p11_n30_twelve_check_functions.json').read_text())
oms=[r7*K(7).log()+r17*K(17).log()+r23*K(23).log() for r7 in [QQ(-2),QQ(-4)/3] for r17 in [QQ(-2),QQ(0),QQ(2)] for r23 in [QQ(2),QQ(4)/3]]
rows=[]
results=json.loads((out/'qc_mixed_15_8_p11_n30_twelve_check.json').read_text())
for fi,d in enumerate(fn):
 cs=[K(sage_eval(s,locals={'O':O})) for s in d['coefficients']]
 for i,om in enumerate(oms):
  ds=cs.copy();ds[0]-=om
  minv=min(v.valuation() for v in ds);tail=d['tail_bound'];assert minv<tail
  js=[j for j,v in enumerate(ds) if v.valuation()==minv]
  poly=Rp([GF(11)(v/(11**minv)) for v in ds]);roots=poly.roots();simple=all(m==1 for a,m in roots)
  base=(fi==0 and i==10) or(fi==1 and i==1)
  depth=0
  while not base and roots==[(GF(11)(0),2)] and max(js)==2:
   depth+=1
   assert depth<10
   ds=[v*(11**j) for j,v in enumerate(ds)]
   minv=min(v.valuation() for v in ds);assert minv<tail
   js=[j for j,v in enumerate(ds) if v.valuation()==minv]
   poly=Rp([GF(11)(v/(11**minv)) for v in ds]);roots=poly.roots();simple=all(m==1 for a,m in roots)
  certified=not roots or simple or (base and roots==[(GF(11)(0),2)] and max(js)==2)
  expected=sum(1 for ss in results['extra_points_by_omega'][i] if (K(sage_eval(ss[1:-1].split(' : ')[0],locals={'O':O})).valuation()<0)==(fi==0))//2+(1 if base else 0)
  assert expected==(1 if base else len(roots)),(fi,i,expected,roots)
  rows.append({'orbit':fi+1,'omega':i,'min_v':int(minv),'tail_bound':int(tail),'strassmann_bound':max(js),'roots_mod11':[(int(a),int(m)) for a,m in roots],'base_known_double':base,'certified':certified,'root_count':len(roots) if not base else 1,'refinement_depth':depth})
print(json.dumps(rows,indent=2));(out/'strassmann_audit.json').write_text(json.dumps(rows,indent=2)+'\n')
assert all(r['certified'] for r in rows)
