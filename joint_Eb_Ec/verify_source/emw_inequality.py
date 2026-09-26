# Проверка финального неравенства EMW 2006 (Theorem 2.2, чётные n):
#   n^2 (3/4 - rho(n)) <= 4 (0.621 eta(n) + 1.216 omega(n) + 9.0776)
# должно НАРУШАТЬСЯ для всех чётных n >= 12 (тогда B_n имеет примитивный делитель).
# Также пересчитываем константу 9.0776 из (16), (17): 5 + 1.5*2.0244 + 1.041/log 5.
from math import log
def primes_of(n):
    ps=[];d=2
    while d*d<=n:
        if n%d==0:
            ps.append(d)
            while n%d==0: n//=d
        d+=1
    if n>1: ps.append(n)
    return ps
const = 5 + 1.5*(log(26)/log(5)) + 1.041/log(5)
print("пересчёт константы: %.4f (у EMW 9.0776, их значение не меньше)" % const)
print("log(26)/log(5) = %.5f (у EMW 2.0244)" % (log(26)/log(5)))
bad=[]
LIM=2*10**6
for n in range(2,LIM+1,2):
    ps=primes_of(n)
    rho=sum(1.0/p**2 for p in ps); eta=2*sum(log(p) for p in ps); om=len(ps)
    lhs=n*n*(0.75-rho); rhs=4*(0.621*eta+1.216*om+9.0776)
    if lhs<=rhs: bad.append(n)
print("чётные n<=%d, где неравенство выполнено (примитивность не гарантирована):"%LIM, bad)
# при n > LIM: lhs >= n^2*(0.75-0.4523) > 0.29 n^2, rhs <= 4(0.621*2*log n*... ) — грубо
import math
n=LIM; print("хвост: 0.29*n^2 = %.3g против 4*(1.242*log n*1.443*log n/ log 2 ...) << "%(0.29*n*n))
