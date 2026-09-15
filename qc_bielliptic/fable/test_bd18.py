# Контроль 2: пример 8.3 Balakrishnan–Dogra, X: y^2 = x^6 - 2x^4 - x^2 + 1. У них p=3; здесь p — первое простое >=5
# хорошей ординарной редукции для обеих кривых. Ожидание (теорема 8.6 BD18): X(Q) = {(0,±1),(±3/2,±1/8),inf±}.
import sys, time; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
R=PolynomialRing(QQ,'x'); x=R.gen()
f = x**6 - 2*x**4 - x**2 + 1
E1=EllipticCurve([0,f[4],0,f[2]*f[6],f[0]*f[6]**2]); E2=EllipticCurve([0,f[2],0,f[0]*f[4],f[0]**2*f[6]])
print('ранги',E1.rank(),E2.rank(), 'кондукторы',E1.conductor(),E2.conductor())
for p in primes(5,60):
    if E1.has_good_reduction(p) and E2.has_good_reduction(p) and E1.is_ordinary(p) and E2.is_ordinary(p):
        print('p =',p); break
t0=time.time()
rat,other = quadratic_chabauty_bielliptic(f,ZZ(p),ZZ(20))
print('время',round(time.time()-t0,1),'с')
print('рациональные',rat); print('лишние',len(other)); print(other)
