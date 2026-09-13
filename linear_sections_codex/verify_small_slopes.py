"""Exact algebra, rank bounds from PARI and eclib, complete torsion checks.
Run with: DOT_SAGE=/tmp/codex-linear-sections-sage python3 verify_small_slopes.py
No analytic-rank/BSD assumptions, point-height search, or floating-point root tests.
"""
from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from cysignals.alarm import alarm,cancel_alarm,AlarmInterrupt
from itertools import permutations
from fractions import Fraction
from pathlib import Path
import json

out=Path(__file__).resolve().parent
offsets=[(1,0),(-1,-1),(0,1),(-1,1),(0,0),(1,-1),(0,-1),(1,1),(-1,0)]
ap={}
for ids in permutations(range(9),4):
 eq=[tuple(offsets[ids[j]][i]-2*offsets[ids[j+1]][i]+offsets[ids[j+2]][i] for i in range(2)) for j in (0,1)]
 nz=[(a,b) for a,b in eq if a or b]
 assert nz # No four distinct coefficient vectors are equally spaced on the 3x3 grid.
 if any(b==0 for a,b in nz):continue
 slopes={QQ(-a)/b for a,b in nz}
 if len(slopes)!=1:continue
 k=slopes.pop()
 if not 0<k<1 or k==QQ(1)/2:continue
 vals=[QQ(offsets[i][0])+offsets[i][1]*k for i in ids]
 if len(set(vals))!=4:continue
 assert vals[1]-vals[0]==vals[2]-vals[1]==vals[3]-vals[2]
 ap.setdefault(str(k),{'cell_indices_zero_based':list(ids),'offset_coefficients':list(map(str,vals))})
assert set(ap)=={'1/4','1/3','2/3','3/4'}

P=PolynomialRing(QQ,'p');p=P.gen()
cases=[(QQ(1)/5,[1,QQ(1)/5,QQ(6)/5]),(QQ(2)/5,[1,QQ(2)/5,QQ(3)/5]),
       (QQ(3)/5,[1,QQ(3)/5,QQ(8)/5]),(QQ(4)/5,[1,QQ(4)/5,QQ(9)/5]),
       (QQ(1)/6,[1,QQ(1)/6,QQ(7)/6]),(QQ(5)/6,[1,QQ(5)/6,QQ(11)/6])]
reports=[]
for k,cs in cases:
 cs=list(map(QQ,cs))
 allcoeffs=[QQ(a)+k*b for a,b in offsets]
 assert all(v in allcoeffs for v in cs)
 f=prod(1+a*p for a in cs);A=f[3]
 E=EllipticCurve([0,f[2],0,A*f[1],A*A])
 assert A*A*f==(A*p)**3+f[2]*(A*p)**2+A*f[1]*(A*p)+A*A
 M=E.minimal_model()
 row={'k':str(k),'selected_coefficients':list(map(str,cs)),
      'selected_cell_indices':[allcoeffs.index(v) for v in cs],
      'cubic':str(f),'A':str(A),'raw_ainvs':list(map(str,E.ainvs())),
      'minimal_ainvs':list(map(str,M.ainvs()))}
 try:
  alarm(20)
  pr=M.pari_curve().ellrank()
  row['pari_rank_bounds']=[int(pr[0]),int(pr[1])]
  ec=mwrank_EllipticCurve(list(map(int,M.ainvs())),verbose=False)
  ec.two_descent(verbose=False)
  row['eclib']={'rank':int(ec.rank()),'upper_bound':int(ec.rank_bound()),'certain':bool(ec.certain())}
 finally:cancel_alarm()
 assert row['pari_rank_bounds']==[0,0]
 assert row['eclib']=={'rank':0,'upper_bound':0,'certain':True}
 points=E.torsion_points()
 assert len(points)==8
 assert all(8*T==E(0) for T in points)
 row['torsion_invariants']=list(map(int,E.torsion_subgroup().invariants()))
 row['raw_torsion_points']=[list(map(str,T)) for T in points]
 ps=sorted(set(T[0]/A for T in points if not T.is_zero()))
 row['all_affine_rational_p']=list(map(str,ps))
 assert all(v<=0 for v in ps)
 # Independent finite-field cardinality upper bound for rational torsion.
 cardinal_gcd=ZZ(0);checks=[]
 for ell in prime_range(5,200):
  if M.discriminant()%ell==0:continue
  card=M.change_ring(GF(ell)).cardinality()
  cardinal_gcd=gcd(cardinal_gcd,card)
  checks.append({'prime':int(ell),'cardinality':int(card)})
  if cardinal_gcd==8:break
 assert cardinal_gcd==8
 row['torsion_cardinality_bound']={'gcd':int(cardinal_gcd),'finite_fields':checks}
 row['positive_p_exists']=False
 reports.append(row)
 print({'k':str(k),'ranks':row['pari_rank_bounds'],'eclib':row['eclib'],'all_p':row['all_affine_rational_p']},flush=True)

slopes=sorted(set(QQ(a)/b for b in range(2,7) for a in range(1,b)))
assert len(slopes)==11
assert set(map(str,slopes))==set(ap)|{r['k'] for r in reports}|{'1/2'}
result={'scope':'All rational slopes k=q/p in (0,1) with reduced denominator <=6; no global slope bound.',
        'all_slopes':list(map(str,slopes)),'four_AP_exclusions':ap,'repeated_cell_slope':'1/2',
        'rank_zero_exclusions':reports,'excluded_slopes':11,
        'meaning':'For positive distinct magic squares, use D8 symmetry to arrange 0<q<p. Every listed slope is impossible for all rational p>0, without a parameter-height bound.',
        'theorem_source':'https://kconrad.math.uconn.edu/blurbs/ugradnumthy/4squarearithprog.pdf'}
(out/'small_slopes_verified.json').write_text(json.dumps(result,indent=2)+'\n')
print('All 11 slopes excluded; output written.',flush=True)
