import json,sys
from collections import Counter,defaultdict
from sympy import factorint
G=json.load(open('grid_Emn_400.json'))
def omega_odd(x): return sum(1 for p in factorint(x) if p!=2)
def n1mod4(x): return sum(1 for p in factorint(x) if p%4==1)
rows=[]
for key,v in G.items():
    m,n=map(int,key.split(',')); rk=v['rk']
    if rk is None: continue
    r1,r2,s=rk; C=r2+2+s
    rows.append(dict(m=m,n=n,r1=r1,r2=r2,s=s,C=C,w=v['w'],N=v['N']))
print('пар',len(rows),'ellrank без ошибки; ранг определён точно у',sum(1 for x in rows if x['r1']==x['r2']))
print('распределение (r1,r2,s):',Counter((x['r1'],x['r2'],x['s']) for x in rows).most_common(12))
print('доля ранга 0 (r2=0):',sum(1 for x in rows if x['r2']==0)/len(rows))
print('доля Sel2 тривиален (C=2):',sum(1 for x in rows if x['C']==2)/len(rows))
print('доля w=+1:',sum(1 for x in rows if x['w']==1)/len(rows))
# формула знака: w =? -(w2) * (-1)^{omega_odd(mn)} * (-1)^{#p≡1(4) | m^2-n^2}
bad=Counter(); tot=Counter()
for x in rows:
    m,n=x['m'],x['n']
    e=omega_odd(m*n)+n1mod4((n*n-m*m))
    pred=-(-1)**e   # без w2
    cls=(m%16,n%16)
    tot[cls]+=1
    if pred!=x['w']: bad[cls]+=1
# w2 должен зависеть только от (m,n) mod 2^j: проверим, что внутри класса mod 16 несовпадение либо всегда, либо никогда
mixed=[(c,bad[c],tot[c]) for c in tot if 0<bad[c]<tot[c]]
print('классов mod16 с неоднозначным w2:',len(mixed), mixed[:10])
w2={c:(-1 if bad[c]==tot[c] else 1) for c in tot if bad[c] in (0,tot[c])}
# сжать до mod 8 / mod 4
for mod in (2,4,8):
    d=defaultdict(set)
    for (a,b),v in w2.items(): d[(a%mod,b%mod)].add(v)
    print(f'w2 определяется mod {mod}:', all(len(v)==1 for v in d.values()), dict(sorted((k,tuple(v)) for k,v in d.items())) if mod<=4 else '')
json.dump(rows,open('grid_rows.json','w'))
