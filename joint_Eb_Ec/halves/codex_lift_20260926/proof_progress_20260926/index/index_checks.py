import json
from sage.all import *
E=EllipticCurve([-34**2,0]);G=E(-16,120);H=E(-2,48)
A=2*G;B=2*H;P=[A,B,A+4*B]
M=matrix(ZZ,[[1,0,1],[0,1,4]])
D=lambda Q:ZZ(Q[0].denominator().sqrt())
out={'curve_N':34,'basis_G':[str(t) for t in G], 'basis_H':[str(t) for t in H], 'point_vectors_in_2basis':[[1,0],[0,1],[1,4]],'reduced_index':int(gcd([ZZ(M.matrix_from_columns([i,j]).det()) for i,j in [(0,1),(0,2),(1,2)]])), 'torsionfree_index':4,'denominators':[str(D(Q)) for Q in P],'common_denominator':str(gcd([D(Q) for Q in P])), 'x_valuations_2':[int(Q[0].valuation(2)) for Q in P], 'x_valuations_3':[int(Q[0].valuation(3)) for Q in P], 'is_AP':bool(P[0][0]+P[2][0]==2*P[1][0])}
assert out['reduced_index']==1 and len(set(out['x_valuations_2']))==len(set(out['x_valuations_3']))==1 and out['common_denominator']=='12' and not out['is_AP']
print(json.dumps(out,indent=2));open('index_checks.json','w').write(json.dumps(out,indent=2))
