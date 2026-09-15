# Контроль 1: пример из README/docstring Bianchi–Padurariu (LMFDB 99856.b.99856.1), p = 5, n = 20.
# Ожидание (docstring): 7 рациональных точек (0:1:0),(0,±2),(±3/4,±43/64) и 8 лишних p-адических точек.
import sys, time; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
R=PolynomialRing(QQ,'x'); x=R.gen()
f = x**6 + 22*x**4 - 19*x**2 + 4
t0=time.time()
rat,other = quadratic_chabauty_bielliptic(f,ZZ(5),ZZ(20))
print('время',round(time.time()-t0,1),'с')
print('рациональные',rat); print('лишние',len(other)); print(other)
