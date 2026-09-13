import json
from fractions import Fraction as F
from collections import Counter,defaultdict
from sympy import factorint
rows=json.load(open('grid_rows.json')); G={(x['m'],x['n']):x for x in rows}
def om(x): return sum(1 for p in factorint(x) if p!=2)
def n1(x): return sum(1 for p in factorint(x) if p%4==1)
# (a) таблица w2 mod 8
w2={}
for x in rows:
    m,n=x['m'],x['n']; e=om(m*n)+n1(n*n-m*m); w2.setdefault((m%8,n%8),set()).add(-x['w']*(-1)**e)
print('w2 по (m mod 8, n mod 8):')
for a in range(8):
    print(a, ' '.join(('%+d'%list(w2[(a,b)])[0]) if (a,b) in w2 else ' .' for b in range(8)))
# (b) ранг 0 / Sel тривиален против числа нечётных простых делителей mn(n^2-m^2)
def stat(name,keyf):
    d=defaultdict(lambda:[0,0,0,0])
    for x in rows:
        k=keyf(x); d[k][0]+=1; d[k][1]+=x['r2']==0; d[k][2]+=x['C']==2; d[k][3]+=x['w']==1
    print(f'\n{name}: ключ: всего, доля rank0, доля Sel2 трив, доля w=+1')
    for k in sorted(d,key=str):
        t,a,b,c=d[k]
        if t>=20: print(f'  {str(k):22s} {t:6d} {a/t:.3f} {b/t:.3f} {c/t:.3f}')
stat('ω_odd(m n (n²-m²))',lambda x:om(x['m']*x['n']*(x['n']**2-x['m']**2)))
stat('ω_odd(mn), ω_odd(n²-m²)',lambda x:(om(x['m']*x['n']),om(x['n']**2-x['m']**2)))
stat('(m,n) mod 4',lambda x:(x['m']%4,x['n']%4))
stat('(m,n) mod 8, только w=+1',lambda x:(x['m']%8,x['n']%8) if x['w']==1 else 'w=-1')
stat('(m,n) mod 3',lambda x:(x['m']%3,x['n']%3))
# при w=+1 и малом числе простых: доля Sel трив
stat('w=+1: ω_odd total',lambda x:om(x['m']*x['n']*(x['n']**2-x['m']**2)) if x['w']==1 else 'w=-1')
# (c) раскладка 711 наклонов по 12 парам
def pairs(r,s):
    P={'{r,s}':(r,s),'{r,s-r}':(r,s-r),'{s,s-r}':(s,s-r),'{s,r+s}':(s,r+s),'{r,r+s}':(r,r+s),
       'Q{s-r,r+s}':(s-r,r+s),'Q{s,|2r-s|}':(s,abs(2*r-s)),'Q{s,2r+s}':(s,2*r+s),'Q{r,2s-r}':(r,2*s-r),'Q{r,2s+r}':(r,2*s+r),'Q{s,2r}':(s,2*r),'Q{r,2s}':(r,2*s)}
    out={}
    for nm,(a,b) in P.items():
        from math import gcd
        g=gcd(a,b); a//=g; b//=g
        if a==b or a==0: out[nm]=None; continue
        out[nm]=G.get((min(a,b),max(a,b)))
    return out
slopes=sorted({F(r,s) for s in range(2,49) for r in range(1,s)})
closed_by=Counter(); nclosed=Counter(); only=Counter(); cubic_none=[]; allnone=[]
for k in slopes:
    r,s=k.numerator,k.denominator; P=pairs(r,s)
    cl=[nm for nm,x in P.items() if x and x['r2']==0]
    for nm in cl: closed_by[nm]+=1
    nclosed[len(cl)]+=1
    if len(cl)==1: only[cl[0]]+=1
    if not any(x and x['r2']==0 for nm,x in P.items() if not nm.startswith('Q')): cubic_none.append(str(k))
    if not cl: allnone.append(str(k))
print('\n(c) сколько наклонов (из 711) закрывает каждая пара:'); 
for nm,c in closed_by.most_common(): print(f'  {nm:14s} {c:4d}   единственная: {only[nm]}')
print('распределение числа закрывающих пар на наклон:',sorted(nclosed.items()))
print('без кубической пары ранга 0:',len(cubic_none)); print('без какой-либо пары ранга 0:',allnone)
