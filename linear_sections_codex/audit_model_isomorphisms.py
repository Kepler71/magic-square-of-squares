from sage.all import *
from pathlib import Path
import json,argparse
out=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--inputs',nargs='+',default=['extended_slopes_7_12.json','extended_slopes_13_24.json','four_cell_followup.json'])
parser.add_argument('--output',default='model_isomorphisms_checked.json')
args=parser.parse_args()
R=PolynomialRing(QQ,names=('x','y'));x,y=R.gens()
def equation(ai,x,y):
 a1,a2,a3,a4,a6=ai
 return y*y+a1*x*y+a3*y-x**3-a2*x*x-a4*x-a6
checks=[]
for name in args.inputs:
 data=json.loads((out/name).read_text())
 for row in data['slopes']:
  if row['status']!='excluded':continue
  w=row['tried'][row['witness_index']]
  a=list(map(QQ,w['certificate']['raw_ainvs']));b=list(map(QQ,w['curve_key'].split(',')))
  E=EllipticCurve(a);M=EllipticCurve(b)
  u,r,s,t=E.isomorphism_to(M).tuple()
  assert u!=0
  assert equation(a,u*u*x+r,u**3*y+s*u*u*x+t)==u**6*equation(b,x,y)
  for ff in w['certificate']['torsion_bound']['finite_fields']:
   assert M.discriminant()%ff['prime']!=0
  checks.append({'file':name,'k':row['k'],'u_r_s_t':list(map(str,[u,r,s,t]))})
(out/args.output).write_text(json.dumps({'all_passed':True,'checks':checks},indent=2)+'\n')
print('Exact change-of-variable identities and good reduction checked:',len(checks))
