"""Fable, 14.09.2026. Контроль Предложения D (нижняя граница на Sel^{phihat}(E'_{m,n})):
dim_2 Sel^{phihat}(E') >= t - #{p | mn, p >= 5} - 2*[3 | mn] - 3 - 2*[3 | n^2-m^2] - 2*[5 | n^2-m^2],
t = число простых q = 1 (mod 4), делящих n^2 - m^2.
Точная группа Sel^{phihat} берётся решателем Claude criterion_proof/selmer_exact.py (множество T) — чужой код, использован
как оракул; моя реализация только у границы. Проверка: нарушений быть не должно.
"""
import sys, json
sys.path.insert(0,'/home/kep/magicKube/criterion_proof')
from selmer_exact import selmer, primes_of
from math import gcd, log2
import multiprocessing as mp

def lower_bound(m,n):
    Pmn=primes_of(m*n); Pd=primes_of(n*n-m*m)
    t=sum(1 for q in Pd if q%4==1)
    lb=t-sum(1 for p in Pmn if p>=5)-2*(3 in Pmn)-3-2*(3 in Pd)-2*(5 in Pd)
    return t,lb
def work(pair):
    m,n=pair
    S,T=selmer(m,n)
    t,lb=lower_bound(m,n)
    dimT=int(round(log2(len(T))))
    return dict(m=m,n=n,t=t,lb=lb,dimT=dimT,dimS=int(round(log2(len(S)))))
if __name__=='__main__':
    N=int(sys.argv[1]) if len(sys.argv)>1 else 150
    pairs=[(m,n) for n in range(2,N+1) for m in range(1,n) if gcd(m,n)==1]
    with mp.get_context('fork').Pool(6) as pool:
        res=pool.map(work,pairs,chunksize=50)
    viol=[r for r in res if r['dimT']<r['lb']]
    print('пар:',len(res),' нарушений dimT < LB:',len(viol), viol[:5])
    nontriv=[r for r in res if r['lb']>=1]
    print('пар, где LB >= 1 (граница нетривиальна):',len(nontriv))
    print('примеры:',[(r['m'],r['n'],r['t'],r['lb'],r['dimT']) for r in nontriv[:12]])
    # насколько граница далека от истины среди пар с большим t
    big=[r for r in res if r['t']>=3]
    print('пар с t>=3:',len(big),' средний dimT:',sum(r['dimT'] for r in big)/max(1,len(big)),' средний LB:',sum(r['lb'] for r in big)/max(1,len(big)))
    json.dump(res,open('/home/kep/magicKube/fable_infinity/selmer_lower.json','w'))
