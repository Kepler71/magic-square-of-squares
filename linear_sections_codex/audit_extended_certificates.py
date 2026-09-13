"""Independent arithmetic audit using only Python Fraction, no Sage elliptic code.
Does not recompute rank; checks polynomials, all displayed torsion, group closure,
finite-field counts, recovered parameters, and nine-cell rejection.
"""
from fractions import Fraction as Q
from pathlib import Path
from math import gcd,isqrt
import json,argparse

out=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--inputs',nargs='+',default=['extended_slopes_7_12.json','extended_slopes_13_24.json','four_cell_followup.json'])
parser.add_argument('--output',default='extended_certificates_audit.json')
args=parser.parse_args()
def mulpoly(a,b):
 c=[Q(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
def add(P,T,a2,a4):
 if P is None:return T
 if T is None:return P
 x,y=P;u,v=T
 if x==u and y==-v:return None
 slope=(3*x*x+2*a2*x+a4)/(2*y) if P==T else (v-y)/(u-x)
 xx=slope*slope-a2-x-u
 return xx,-y+slope*(x-xx)
def square(x):return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator
results=[]
for name in args.inputs:
 data=json.loads((out/name).read_text())
 for row in data['slopes']:
  if row['status']!='excluded':continue
  tr=row['tried'][row['witness_index']];cert=tr['certificate'];k=Q(row['k'])
  if tr['eclib']!={'rank':0,'upper_bound':0,'certain':True}:
   assert cert.get('rank_proof')=='eclib_selmer_plus_fisher_pairing'
   pa=json.loads((out/'pairing_40_43_audit.json').read_text())
   pp=json.loads((out/'pairing_40_43.json').read_text())
   assert pa['all_passed'] and pa['rank_upper_bound']==0
   assert ','.join(pp['curve_ainvs'])==tr['curve_key']
  assert data['rank_cache'][tr['curve_key']]['rank_bounds']==[0,0]
  cs=list(map(Q,tr['coefficients']));A=Q(cert['A'])
  coeffs={Q(1),k,1+k,1-k,Q(-1),-k,-1-k,k-1}
  assert len(set(cs))==len(cs) and set(cs)<=coeffs
  poly=[Q(1)]
  if len(cs)==3:
   for c in cs:poly=mulpoly(poly,[Q(1),c])
   recover=lambda x:x/A
   extra=set()
  else:
   assert len(cs)==4
   a=cs[0];assert Q(cert['root'])==-1/a
   poly=[a]
   for c in cs[1:]:poly=mulpoly(poly,[c,1-c/a])
   recover=lambda x:-1/a+A/x
   extra={-1/a}
  assert poly[3]==A
  ai=list(map(Q,cert['raw_ainvs']))
  assert ai==[Q(0),poly[2],Q(0),A*poly[1],A*A*poly[0]]
  points=[]
  for raw in cert['torsion_points']:
   x,y,z=map(Q,raw)
   if z==0:
    assert x==0 and y!=0;points.append(None);continue
   assert z==1 and y*y==x**3+ai[1]*x*x+ai[3]*x+ai[4]
   points.append((x,y))
  n=len(points);assert n==cert['torsion_order'] and len(set(points))==n and None in points
  group=set(points)
  for P in points:
   for T in points:assert add(P,T,ai[1],ai[3]) in group
   acc=None
   for _ in range(n):acc=add(acc,P,ai[1],ai[3])
   assert acc is None
  ma=list(map(int,tr['curve_key'].split(',')));bound=0
  for ff in cert['torsion_bound']['finite_fields']:
   ell=ff['prime'];a1,a2,a3,a4,a6=ma
   count=ell+1
   for x in range(ell):
    disc=((a1*x+a3)**2+4*(x**3+a2*x*x+a4*x+a6))%ell
    if disc:count+=1 if pow(disc,(ell-1)//2,ell)==1 else -1
   assert count==ff['cardinality'];bound=gcd(bound,count)
  assert bound==n==cert['torsion_bound']['bound']
  ps=set(extra)
  for P in points:
   if P is not None and (len(cs)==3 or P[0]!=0):ps.add(recover(P[0]))
  assert ps=={Q(t['p']) for t in cert['all_p_tests']}
  for test in cert['all_p_tests']:
   p=Q(test['p']);q=k*p
   cells=[1+p,1-p-q,1+q,1-p+q,Q(1),1+p-q,1-q,1+p+q,1-p]
   assert cells==list(map(Q,test['cells']))
   positive=min(cells)>0;distinct=len(set(cells))==9;inds=[i for i,c in enumerate(cells) if square(c)]
   assert positive==test['positive'] and distinct==test['distinct'] and inds==test['square_indices']
   assert not(positive and distinct and len(inds)==9)
 results.append({'file':name,'excluded_certificates_checked':sum(r['status']=='excluded' for r in data['slopes'])})
report={'checks':results,'all_passed':True,'checks_do_not_include':'Independent recomputation of PARI/eclib rank; minimal-model isomorphisms are checked separately.'}
(out/args.output).write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
