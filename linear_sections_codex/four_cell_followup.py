"""Genus-one quartic-product filters for slopes unresolved by all 56 cubic products."""
from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from itertools import combinations
from pathlib import Path
import json,time,argparse

out=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--input',default='extended_slopes_13_24.json')
parser.add_argument('--output',default='four_cell_followup.json')
args=parser.parse_args()
prior=json.loads((out/args.input).read_text())
todo=[QQ(r['k']) for r in prior['slopes'] if r['status']=='unresolved']
R=PolynomialRing(QQ,'x');x=R.gen();K=R.fraction_field()
cache={};rows=[];start=time.monotonic()
def save():
 (out/args.output).write_text(json.dumps({'slopes':rows,'rank_cache':cache,
  'elapsed_seconds':float(time.monotonic()-start),'scope':'Four-cell product quotients; all rational points enumerated only for rank zero.'},indent=2)+'\n')
def check(p,k):
 q=k*p;v=[1+p,1-p-q,1+q,1-p+q,QQ(1),1+p-q,1-q,1+p+q,1-p]
 return {'p':str(p),'cells':list(map(str,v)),'positive':bool(min(v)>0),'distinct':len(set(v))==9,
         'square_indices':[j for j,z in enumerate(v) if z>=0 and z.is_square()]}
for k in todo:
 co=[QQ(1),k,1+k,1-k,QQ(-1),-k,-1-k,k-1]
 row={'k':str(k),'status':'unresolved','tried':[]};rows.append(row)
 for cs in combinations(co,4):
  a=cs[0];root=-1/a
  g=a*prod((1-c/a)*x+c for c in cs[1:]);A=g[3]
  assert K(x**4*prod(1+c*(root+1/x) for c in cs))==g
  E=EllipticCurve([0,g[2],0,A*g[1],A*A*g[0]]);M=E.minimal_model()
  key=','.join(map(str,M.ainvs()))
  if key not in cache:
   try:
    alarm(4);pr=M.pari_curve().ellrank()
    cache[key]={'rank_bounds':[int(pr[0]),int(pr[1])],'minimal_ainvs':list(map(str,M.ainvs()))}
   except (Exception,AlarmInterrupt) as e:cache[key]={'error':type(e).__name__}
   finally:cancel_alarm()
  trial={'coefficients':list(map(str,cs)),'curve_key':key};row['tried'].append(trial)
  if cache[key].get('rank_bounds')!=[0,0]:continue
  try:
   alarm(15);ec=mwrank_EllipticCurve(list(map(int,M.ainvs())),verbose=False);ec.two_descent(verbose=False)
   trial['eclib']={'rank':int(ec.rank()),'upper_bound':int(ec.rank_bound()),'certain':bool(ec.certain())}
  finally:cancel_alarm()
  if trial['eclib']!={'rank':0,'upper_bound':0,'certain':True}:continue
  pts=E.torsion_points();order=len(pts)
  assert all(order*T==E(0) for T in pts)
  bound=ZZ(0);fields=[]
  for ell in prime_range(5,300):
   if M.discriminant()%ell==0:continue
   count=M.change_ring(GF(ell)).cardinality();bound=gcd(bound,count)
   fields.append({'prime':int(ell),'cardinality':int(count)})
   if bound==order:break
  if bound!=order:trial['incomplete_torsion_bound']=int(bound);continue
  # O corresponds to p=root: at least one required cell zero. X=0 gives p=infinity.
  ps=sorted(set([root]+[root+A/T[0] for T in pts if not T.is_zero() and T[0]!=0]))
  tests=[check(p,k) for p in ps]
  full=[v for v in tests if v['positive'] and v['distinct'] and len(v['square_indices'])==9]
  trial['certificate']={'raw_ainvs':list(map(str,E.ainvs())),'A':str(A),'root':str(root),
        'torsion_points':[[str(z) for z in T] for T in pts],'torsion_order':order,
        'torsion_bound':{'bound':int(bound),'finite_fields':fields},'all_p_tests':tests,
        'exceptional_points':'O gives p=root with a zero cell; X=0 gives p=infinity, outside the normalized positive square.'}
  row['status']='full_candidate_found' if full else 'excluded';row['witness_index']=len(row['tried'])-1
  if full:row['full_candidates']=full
  break
 save();print({'k':str(k),'status':row['status'],'quartics_tried':len(row['tried'])},flush=True)
save()
print({'remaining':[r['k'] for r in rows if r['status']=='unresolved'],'models':len(cache)},flush=True)
