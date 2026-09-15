from sage.all import *
import sys, time
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier
r,s=204,247
lam=[s,-s,r,-r,s-r,r-s,s+r,-s-r]
T=[QQ(s),QQ(r),QQ(s+r)]
t0=time.time()
E,L,c0=weier(T); print(E, L, c0)
rk=E.pari_curve().ellrank(); print('ellrank',rk, time.time()-t0)
pts=[E(list(p)) for p in rk[3]]
print(pts)
t0=time.time()
sat=E.saturation(pts, max_prime=3); print('sat 2,3', sat, time.time()-t0)
t0=time.time()
sat2=E.saturation(pts); print('sat full', sat2, time.time()-t0)
print('torsion', E.torsion_subgroup().invariants(), [P for P in E.torsion_points()][:8])
# z=0 point: w = 1/(0+1/c0)=c0, X=L*c0, Y^2=h(X)
X=L*c0; Y2=E.defining_polynomial()  # check
P0=E.lift_x(X); print('P0',P0, P0.height())
# discrete log test
l=1009
El=E.change_ring(GF(l)); G=El.abelian_group(); print(G.invariants(), G.gens())
Q=El(sat[0][0]); t0=time.time(); print(G.discrete_log(Q), time.time()-t0)
