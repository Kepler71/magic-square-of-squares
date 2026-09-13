import json,sys,time
from isogeny_descent import descent
N=int(sys.argv[1]); rows=json.load(open('grid_rows.json'))
out={}; bad=[]; t0=time.time(); cnt=0
for x in rows:
    m,n=x['m'],x['n']
    if n>N: continue
    D=descent(m,n); cnt+=1
    ru=D['rank_upper']; r1,r2,s=x['r1'],x['r2'],x['s']
    out[f'{m},{n}']=D
    if D['pure_rank0'] and not (r2==0 and s==0): bad.append(('pure but PARI',m,n,x))
    if ru<r2+s: bad.append(('upper<PARI Sel excess',m,n,ru,x))
    if ru<r1: bad.append(('upper<rank lower',m,n,ru,x))
print('пар',cnt,'время',round(time.time()-t0),'противоречий',len(bad)); print(bad[:10])
from collections import Counter
print('распределение (мой верх - (r2+s)):',Counter(out[f"{x['m']},{x['n']}"]['rank_upper']-(x['r2']+x['s']) for x in rows if x['n']<=N))
json.dump(out,open(f'descent_{N}.json','w'))
