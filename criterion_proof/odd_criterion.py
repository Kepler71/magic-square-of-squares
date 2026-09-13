"""Критерий ранга 0 для E_{m,n}: y^2 = x(x+m^2)(x+n^2), ТОЛЬКО из доказанных лемм в нечётных простых.
Свой код (не использует criterion.py Fable).
2^rank = |α(E)|·|α'(E')|/4 <= |Sφ|·|Sφ'|/4, где Sφ, Sφ' — кандидаты, НЕ отсечённые леммами:
  Лемма 3 (C_d, d | mn со знаком):  q ≡ 1 (4), q | n^2-m^2, (d/q) = -1  => нет точек над Q_q
  Вещественное место (C'_d):        d < 0                               => нет точек над R
  Лемма 1 (C'_d, d | n^2-m^2):      нечётное q | d, q ≡ 3 (4)             => нет точек над Q_q
  Лемма 2 (C'_d):                    нечётное p | mn, (d/p) = -1          => нет точек над Q_p
Простое 2 не используется вовсе: всё, что оно могло бы отсечь, остаётся в кандидатах (граница верхняя)."""
from math import gcd, log2
from itertools import combinations
def primes_of(x):
    x=abs(x); out=[]; d=2
    while d*d<=x:
        if x%d==0:
            out.append(d)
            while x%d==0: x//=d
        d+=1 if d==2 else 2
    if x>1: out.append(x)
    return out
def legendre(a,p):
    a%=p
    if a==0: return 0
    return 1 if pow(a,(p-1)//2,p)==1 else -1
def odd_bound(m,n):
    N=n*n-m*m; PM=primes_of(m*n); PN=primes_of(N)
    q1=[q for q in PN if q%4==1]
    podd=[p for p in PM if p!=2]
    Sphi=0
    for r in range(len(PM)+1):
        for sub in combinations(PM,r):
            d0=1
            for p in sub: d0*=p
            for d in (d0,-d0):
                if not any(legendre(d,q)==-1 for q in q1): Sphi+=1
    Sphi_hat=0
    for r in range(len(PN)+1):
        for sub in combinations(PN,r):
            d=1
            for q in sub: d*=q
            if any(q%4==3 for q in sub if q!=2): continue
            if any(legendre(d,p)==-1 for p in podd): continue
            Sphi_hat+=1
    return Sphi, Sphi_hat, int(log2(Sphi*Sphi_hat))-2
def proves_rank0(m,n):
    a,b,_=odd_bound(m,n); return a*b<=4
if __name__=="__main__":
    import json,sys
    rows=json.load(open("../linear_sections_codex/fable/grid_rows.json"))
    proven=0; false=[]; r0=0
    for x in rows:
        m,n=x['m'],x['n']
        if x['r2']==0: r0+=1
        if proves_rank0(m,n):
            proven+=1
            if x['r2']!=0: false.append((m,n))
    print(f"сетка {len(rows)} пар: ранг 0 по PARI {r0}; доказан нечётным критерием {proven}; противоречий с PARI {len(false)} {false[:5]}")
