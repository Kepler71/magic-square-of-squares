"""Screen slopes of reduced denominator 7..12; certify exclusions, retain failures.
Necessary condition: product of any three noncentral cells is a rational square.
For rank zero, enumerate every rational point and test all nine cells exactly.
"""
from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from itertools import combinations
from pathlib import Path
import json,time,argparse

out=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--min-den',type=int,default=7)
parser.add_argument('--max-den',type=int,default=12)
args=parser.parse_args()
assert 7<=args.min_den<=args.max_den
R=PolynomialRing(QQ,'p');p=R.gen()
cache={};ec_cache={};rows=[];start=time.monotonic()
def save():
 data={'denominators':[args.min_den,args.max_den],'slopes':rows,'rank_cache':cache,
       'elapsed_seconds':float(round(time.monotonic()-start,3)),
       'scope':'Each certified exclusion is for all rational p, not a height-bounded search. Unresolved slopes are not solutions.'}
 (out/f'extended_slopes_{args.min_den}_{args.max_den}.json').write_text(json.dumps(data,indent=2)+'\n')
def cells(v,k):
 q=k*v
 return [1+v,1-v-q,1+q,1-v+q,QQ(1),1+v-q,1-q,1+v+q,1-v]
def validate(v,k):
 c=cells(v,k)
 return {'p':str(v),'cells':list(map(str,c)), 'positive':bool(min(c)>0),
         'distinct':len(set(c))==9,'square_indices':[i for i,z in enumerate(c) if z>=0 and z.is_square()]}

slopes=sorted(set(QQ(a)/b for b in range(args.min_den,args.max_den+1) for a in range(1,b) if gcd(a,b)==1))
for k in slopes:
 co=[QQ(1),k,1+k,1-k,QQ(-1),-k,-1-k,k-1]
 preferred=[(co[0],co[1],co[2]),(co[0],co[1],co[3]),(co[1],co[3],co[2])]
 triples=[];seen=set()
 for tr in preferred+list(combinations(co,3)):
  key=tuple(sorted(tr))
  if key not in seen:triples.append(tr);seen.add(key)
 assert len(triples)==56
 row={'k':str(k),'status':'unresolved','tried':[]};rows.append(row)
 for tr in triples:
  f=prod(1+c*p for c in tr);A=f[3]
  E=EllipticCurve([0,f[2],0,A*f[1],A*A]);M=E.minimal_model()
  assert A*A*f==(A*p)**3+f[2]*(A*p)**2+A*f[1]*(A*p)+A*A
  key=','.join(map(str,M.ainvs()))
  if key not in cache:
   try:
    alarm(4);r=M.pari_curve().ellrank()
    cache[key]={'rank_bounds':[int(r[0]),int(r[1])],'minimal_ainvs':list(map(str,M.ainvs()))}
   except (Exception,AlarmInterrupt) as err:cache[key]={'error':type(err).__name__}
   finally:cancel_alarm()
  trial={'coefficients':list(map(str,tr)),'curve_key':key}
  row['tried'].append(trial)
  if cache[key].get('rank_bounds')!=[0,0]:continue
  if key not in ec_cache:
   try:
    alarm(15);ec=mwrank_EllipticCurve(list(map(int,M.ainvs())),verbose=False);ec.two_descent(verbose=False)
    ec_cache[key]={'rank':int(ec.rank()),'upper_bound':int(ec.rank_bound()),'certain':bool(ec.certain())}
   except (Exception,AlarmInterrupt) as err:ec_cache[key]={'error':type(err).__name__}
   finally:cancel_alarm()
  trial['eclib']=ec_cache[key]
  if trial['eclib']!={'rank':0,'upper_bound':0,'certain':True}:continue
  pts=E.torsion_points();order=len(pts)
  assert all(order*T==E(0) for T in pts)
  bound=ZZ(0);checks=[]
  for ell in prime_range(5,300):
   if M.discriminant()%ell==0:continue
   count=M.change_ring(GF(ell)).cardinality();bound=gcd(bound,count)
   checks.append({'prime':int(ell),'cardinality':int(count)})
   if bound==order:break
  if bound!=order:
   trial['uncertified_torsion_bound']={'found':order,'bound':int(bound)};continue
  ps=sorted(set(T[0]/A for T in pts if not T.is_zero()))
  tested=[validate(v,k) for v in ps]
  full=[v for v in tested if v['positive'] and v['distinct'] and len(v['square_indices'])==9]
  trial['certificate']={'raw_ainvs':list(map(str,E.ainvs())),'A':str(A),
       'torsion_order':order,'torsion_points':[[str(z) for z in T] for T in pts],
       'torsion_bound':{'bound':int(bound),'finite_fields':checks},'all_p_tests':tested}
  if full:
   row['status']='full_candidate_found';row['full_candidates']=full
  else:row['status']='excluded';row['witness_index']=len(row['tried'])-1
  break
 save()
 print({'k':str(k),'status':row['status'],'triples_tried':len(row['tried'])},flush=True)
save()
print({'slopes':len(rows),'excluded':sum(r['status']=='excluded' for r in rows),
       'unresolved':[r['k'] for r in rows if r['status']=='unresolved'],'models':len(cache)},flush=True)
