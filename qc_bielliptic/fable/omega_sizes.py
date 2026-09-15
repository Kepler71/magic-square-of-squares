import sys; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
from models import fpoly,E12
a,b,c,p=map(int,sys.argv[1:5]); p=ZZ(p); f=fpoly(a,b,c); E1,E2=E12(f); E1min=E1.minimal_model(); E2min=E2.minimal_model()
K=Qp(p,20); H=HyperellipticCurve(f)
bp1,W1=local_heights_at_bad_primes_new(E1min,E1,K); bp2,W2=local_heights_at_bad_primes_new(E2min,E2,K)
print('E1 bad',[(q,len(w)) for q,w in zip(bp1,W1)]); print('E2 bad',[(q,len(w)) for q,w in zip(bp2,W2)])
print('E1 delta',factor(E1.discriminant()/E1min.discriminant()),'E2 delta',factor(E2.discriminant()/E2min.discriminant()))
print('a0',factor(f[0]))
for q in sorted(set(bp1+bp2+ZZ(f[0]).prime_factors())):
    print(q,'pot good',has_potential_good_reduction(H,q))
