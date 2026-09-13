from sage.all import *
from pathlib import Path
import json
out=Path(__file__).resolve().parent
proof=json.loads((out/'pairing_40_43.json').read_text())
audit=json.loads((out/'pairing_40_43_audit.json').read_text())
assert audit['all_passed'] and audit['rank_upper_bound']==0
prior=json.loads((out/'four_cell_25_48.json').read_text())
row=next(r for r in prior['slopes'] if r['k']=='40/43')
trial=next(t for t in row['tried'] if prior['rank_cache'][t['curve_key']].get('rank_bounds')==[0,0])
assert trial['curve_key']==','.join(proof['curve_ainvs'])
k=QQ(40)/43;cs=list(map(QQ,trial['coefficients']))
R=PolynomialRing(QQ,'x');x=R.gen();K=R.fraction_field()
a=cs[0];root=-1/a;g=a*prod((1-c/a)*x+c for c in cs[1:]);A=g[3]
assert K(x**4*prod(1+c*(root+1/x) for c in cs))==g
E=EllipticCurve([0,g[2],0,A*g[1],A*A*g[0]]);M=E.minimal_model()
assert list(map(str,M.ainvs()))==proof['curve_ainvs']
tors=E.torsion_points();order=len(tors)
assert order==8 and sum(2*T==E(0) for T in tors)==4
assert all(8*T==E(0) for T in tors)
fields=[];bound=ZZ(0)
for ell in prime_range(5,300):
 if M.discriminant()%ell==0:continue
 card=M.change_ring(GF(ell)).cardinality();bound=gcd(bound,card)
 fields.append({'prime':int(ell),'cardinality':int(card)})
 if bound==8:break
assert bound==8
ps=sorted(set([root]+[root+A/T[0] for T in tors if not T.is_zero() and T[0]!=0]))
assert ps==[QQ(-43)/40,QQ(-43)/83,QQ(0),QQ(43)/83,QQ(43)/40]
tests=[]
for p in ps:
 q=k*p;v=[1+p,1-p-q,1+q,1-p+q,QQ(1),1+p-q,1-q,1+p+q,1-p]
 test={'p':str(p),'cells':list(map(str,v)),'positive':bool(min(v)>0),'distinct':len(set(v))==9,
       'square_indices':[i for i,c in enumerate(v) if c>=0 and c.is_square()]}
 assert not(test['positive'] and test['distinct'] and len(test['square_indices'])==9)
 tests.append(test)
trial['certificate']={'raw_ainvs':list(map(str,E.ainvs())),'A':str(A),'root':str(root),'torsion_order':order,
       'torsion_points':[[str(c) for c in T] for T in tors],
       'torsion_bound':{'bound':8,'finite_fields':fields},'all_p_tests':tests,
       'rank_proof':'eclib_selmer_plus_fisher_pairing',
       'rank_proof_files':['pairing_40_43.json','pairing_40_43_audit.json'],
       'exceptional_points':'O gives zero cell p=root; X=0 gives p=infinity, outside normalized finite square.'}
result={'slopes':[{'k':'40/43','status':'excluded','tried':[trial],'witness_index':0}],
        'rank_cache':{trial['curve_key']:prior['rank_cache'][trial['curve_key']]},
        'note':'eclib rank upper bound remains2; its Selmer dimension4 plus independent Fisher pairing rank>=2 and rational2-torsion dimension2 prove rank0.'}
(out/'exception_40_43_certified.json').write_text(json.dumps(result,indent=2)+'\n')
print('40/43 excluded; all rational p:',ps)
