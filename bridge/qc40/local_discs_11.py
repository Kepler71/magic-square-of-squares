from pathlib import Path
import json
p=11;qr={i*i%p for i in range(p)}
rows=[]
for a,b in [(x,1) for x in range(p)]+[(1,0)]:
 c=289*(a*a+b*b)
 vals=[c-322*a*b,c,c+322*a*b,49*b*b+529*a*a,529*b*b+49*a*a]
 rows.append({'a':a,'b':b,'values':[v%p for v in vals],'all_five_square':all(v%p in qr for v in vals)})
assert [(r['a'],r['b']) for r in rows if r['all_five_square']]==[(0,1),(1,0)]
res={'p':p,'rows':rows,'allowed_projective_x':['0','infinity'],'points_on_D':['(0,1)','(0,-1)','infinity: y/x^3=3','infinity: y/x^3=-3'],'zeros_allowed':True}
(Path(__file__).parent/'local_discs_11.json').write_text(json.dumps(res,indent=2)+'\n');print(res['allowed_projective_x'])
