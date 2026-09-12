from pathlib import Path
import json
A=[4092529,47287151,37608911,44182609];a0,a2,a4,a6=A
rows=[]
for p in [3,5,13,37,1213,9397]:
 qr={i*i%p for i in range(p)};bad=[];eligible=[];roots=[]
 for x in range(p):
  c=289*(1+x*x);vals=[c-322*x,c,c+322*x,49+529*x*x,529+49*x*x]
  if all(v%p in qr for v in vals):eligible.append(x)
  fx=(a6*x**6+a4*x**4+a2*x*x+a0)%p
  if fx:continue
  roots.append(x)
  if not all(v%p in qr for v in vals):continue
  xx1=a6*x*x%p
  sing1=(3*xx1**2+2*a4*xx1+a2*a6)%p==0
  sing2=False
  if x:
   xx2=a0*pow(x*x,-1,p)%p
   sing2=(3*xx2**2+2*a2*xx2+a4*a0)%p==0
  if sing1 or sing2:bad.append({'x':x,'E1_singular':sing1,'E2_singular':sing2,'values':[v%p for v in vals]})
 rows.append({'p':p,'a0_a6_units':a0%p!=0 and a6%p!=0,'D_branch_roots':roots,'full_five_square_residue_count':len(eligible),'possible_singular_images':bad})
 print(p,len(eligible),len(bad),flush=True)
(Path(__file__).parent/'unramified_height_primes.json').write_text(json.dumps(rows,indent=2)+'\n')
