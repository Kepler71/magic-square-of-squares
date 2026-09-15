# Кратности классов изогении эллиптических множителей (3–4 клетки) якобиана кривой наклона X_k над ℚ.
# Демьяненко–Манин: m независимых морфизмов X_k → E (с точностью до изогении) и rank E(ℚ) < m ⇒ X_k(ℚ) конечно (эффективно).
from sage.all import *
import sys, itertools, json
sys.path.insert(0,'/home/kep/magicKube/census_genus1')
from g1census import weier
from collections import defaultdict
def run(sl, P=range(3,200)):
    r,s=map(int,sl.split('/')); vals=[r,s-r,s,s+r]
    C=sorted(set(vals+[-v for v in vals]))
    curves=[]
    for k in (3,4):
        for S in itertools.combinations(C,k):
            E,L,c0=weier([QQ(c) for c in S])
            if E.discriminant()==0: continue
            curves.append((S,E.minimal_model()))
    # подпись a_p по хорошим простым (одинакова у изогенных кривых)
    groups=defaultdict(list)
    for S,E in curves:
        N=E.conductor()
        sig=tuple(E.ap(p) if N%p else None for p in primes(3,200))
        groups[sig].append((S,E))
    out=[]
    for sig,lst in groups.items():
        E0=lst[0][1]
        # подтвердить изогению точно
        conf=all(E0.is_isogenous(E) for _,E in lst[1:])
        lo,hi=[int(t) for t in pari(E0).ellrank()[:2]]
        out.append(dict(m=len(lst),confirmed=bool(conf),rank=[lo,hi],subsets=[list(map(int,S)) for S,_ in lst]))
    out.sort(key=lambda d:-d['m'])
    return out
if __name__=='__main__':
    sl=sys.argv[1]; res=run(sl)
    json.dump(res,open('mult_'+sl.replace('/','_')+'.json','w'))
    print(sl,'классов изогении:',len(res))
    for d in res[:12]: print(' m =',d['m'],'подтверждено' if d['confirmed'] else 'НЕ подтверждено','ранг',d['rank'],'наборы',d['subsets'][:4])
    print('Демьяненко применим (ранг < m):',[ (d['m'],d['rank']) for d in res if d['rank'][0]==d['rank'][1] and d['rank'][1] < d['m']])
