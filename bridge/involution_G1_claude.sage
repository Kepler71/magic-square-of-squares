# -*- coding: utf-8 -*-
# ПРОВЕРКА ИНВОЛЮЦИИ СЕМЕЙСТВА G1:  (m,n,t) -> (m+n, m-n, x=(t-1)/(t+1))
from sage.all import *
Rmn = PolynomialRing(QQ, ['m','n']); m, n = Rmn.gens()
Rx = PolynomialRing(Rmn, 'x'); x = Rx.gen()
FR = FractionField(Rx)
s = (m*m + n*n)/2
tt = FR((1+x)/(1-x))
F0 = m*m + n*n*tt**2
F4 = s*(1 + tt**2)
F8 = n*n + m*m*tt**2
L  = F4 - 2*m*n*tt
U  = F4 + 2*m*n*tt
M, N = m+n, m-n
S = (M*M + N*N)/2
F0p = M*M + N*N*x**2
F4p = S*(1 + x**2)
F8p = N*N + M*M*x**2
Lp  = F4p - 2*M*N*x
Up  = F4p + 2*M*N*x
sc = FR((1-x)**2)
print("F0*(1-x)^2 == L'  :", bool(F0*sc == Lp))
print("F4*(1-x)^2 == F4' :", bool(F4*sc == F4p))
print("F8*(1-x)^2 == U'  :", bool(F8*sc == Up))
print("L *(1-x)^2 == F8' :", bool(L*sc  == F8p))
print("U *(1-x)^2 == F0' :", bool(U*sc  == F0p))
print("S == m^2+n^2      :", bool(S == m*m+n*n))
print("2MN == 2(m^2-n^2) :", bool(2*M*N == 2*(m*m-n*n)))
