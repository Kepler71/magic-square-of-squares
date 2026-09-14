# Fable, 14.09.2026. «Симметричная стена» спуска: семейства пар (m,n), у которых одна из сторон
# (делители mn или делители n^2-m^2) ограничена по построению, а другая — нет; убивать классы на неограниченной
# стороне нечем. Точные группы Селмера 2-изогении — решатель Claude (оракул), ранг — PARI ellrank.
import sys
sys.path.insert(0,'/home/kep/magicKube/criterion_proof')
from selmer_exact import selmer
from sage.all import *
from math import log2
def row(tag,m,n):
    S,T=selmer(m,n); b=int(round(log2(len(S)*len(T))))-2
    E=EllipticCurve([0,m*m+n*n,0,m*m*n*n,0])
    try: lo,hi=E.rank_bounds()
    except Exception: lo,hi=('?','?')
    print(f'{tag:28s} (m,n)=({m},{n}) |Sphi|={len(S):3d} |Sphihat|={len(T):2d} граница={b}  PARI ранг [{lo},{hi}]  omega(mn)={len(prime_divisors(m*n))} t={sum(1 for q in prime_divisors(n*n-m*m) if q%4==1)}',flush=True)
print('--- n-m=1, n+m=3^j: n^2-m^2=3^j, Sphihat={1} принудительно; Sphi — все делители mn ---')
for j in range(2,9):
    row('(3^j-1)/2,(3^j+1)/2',(3**j-1)//2,(3**j+1)//2)
print('--- (2^j-1, 2^j+1): n^2-m^2=2^(j+2), ни одного q=1(4); Sphi — все делители 4^j-1 ---')
for j in range(2,9):
    row('(2^j-1,2^j+1)',2**j-1,2**j+1)
print('--- n-m=1, n+m=ell простое 5 mod 8: Sphihat={1} принудительно; Sphi — делители (ell^2-1)/4 ---')
for l in [5,13,29,37,53,61,101,109,149,157,173,181,197,229]:
    row('((ell-1)/2,(ell+1)/2)',(l-1)//2,(l+1)//2)
print('--- (1,ell), ell с условием теоремы 1/ell (контроль: ранг 0) ---')
for l in [13,37,43,197,277]:
    row('(1,ell)',1,l)
