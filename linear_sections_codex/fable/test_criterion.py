import json,random,time
from collections import Counter
from criterion import criterion
from isogeny_descent import descent
from fractions import Fraction as F
from math import gcd
rows=json.load(open('grid_rows.json')); G={(x['m'],x['n']):x for x in rows}
# (i) быстрый критерий против PARI на всей сетке
t0=time.time(); bad=[]; cnt=Counter(); crit={}
for x in rows:
    m,n=x['m'],x['n']; c=criterion(m,n); crit[(m,n)]=c
    ru=c['rank_upper']; sel_triv=(x['r2']==0 and x['s']==0)
    if c['pure_rank0'] and not sel_triv: bad.append(('pure, but PARI Sel nontrivial',m,n,x))
    if ru<x['r2']+x['s']: bad.append(('upper < PARI',m,n,ru,x))
    if ru<x['r1']: bad.append(('upper < r1',m,n,ru,x))
    cnt[('pure' if c['pure_rank0'] else 'notpure','SelTriv' if sel_triv else 'SelNontriv','rank0' if x['r2']==0 else 'rank>0')]+=1
print('(i) критерий vs PARI, пар',len(rows),'время',round(time.time()-t0),'противоречий',len(bad)); print(bad[:5])
for k,v in sorted(cnt.items()): print('   ',k,v)
print('   распределение (верх критерия − (r2+s)):',Counter(crit[(x['m'],x['n'])]['rank_upper']-(x['r2']+x['s']) for x in rows))
# (ii) быстрый критерий против полного решателя на независимой выборке 130<=n<=400
random.seed(7); sample=random.sample([(m,n) for (m,n) in G if n>=130],250); t0=time.time(); mism=[]
for m,n in sample:
    c=criterion(m,n); d=descent(m,n)
    if set(c['S_phi'])!=set(d['S_phi']) or set(c['S_phihat'])!=set(d['S_phihat']): mism.append((m,n,c,d))
print('(ii) критерий vs решатель Гензеля на',len(sample),'парах n>=130: расхождений',len(mism),'время',round(time.time()-t0)); print(mism[:3])
# (iii) наклоны со знаменателем <=200: 12 пар
def red(a,b):
    g=gcd(a,b); a//=g; b//=g; return (min(a,b),max(a,b))
def pairs(r,s):
    return {'{r,s}':(r,s),'{r,s-r}':(r,s-r),'{s,s-r}':(s,s-r),'{s,r+s}':(s,r+s),'{r,r+s}':(r,r+s),
       'Q{s-r,r+s}':(s-r,r+s),'Q{s,|2r-s|}':(s,abs(2*r-s)),'Q{s,2r+s}':(s,2*r+s),'Q{r,2s-r}':(r,2*s-r),'Q{r,2s+r}':(r,2*s+r),'Q{s,2r}':(s,2*r),'Q{r,2s}':(r,2*s)}
res={}
for lo,hi in ((2,48),(49,100),(101,200)):
    slopes=[F(r,s) for s in range(lo,hi+1) for r in range(1,s) if gcd(r,s)==1]
    c=Counter(); per_pair=Counter(); notclosed=[]
    for k in slopes:
        r,s=k.numerator,k.denominator; P=pairs(r,s); closed_th=False; closed_pari=False; unknown=False
        for nm,(a,b) in P.items():
            a,b=red(a,b)
            if a==b or a==0: continue
            cc=criterion(a,b)
            if cc['pure_rank0']: closed_th=True; per_pair[nm]+=1
            x=G.get((a,b))
            if x is None: unknown=True
            elif x['r2']==0: closed_pari=True
        c['theorem closes']+=closed_th
        c['PARI rank0 (known pairs)']+=closed_pari
        c['theorem closes but PARI not (impossible)']+=closed_th and (not closed_pari) and (not unknown)
        c['PARI closes, theorem not']+=closed_pari and not closed_th
        if not closed_th: notclosed.append(str(k))
    print(f'(iii) знаменатели {lo}..{hi}: наклонов {len(slopes)}',dict(c))
    print('     теорема закрывает по парам:',per_pair.most_common())
    res[f'{lo}-{hi}']={'n':len(slopes),'counts':dict(c),'not_closed_by_theorem':notclosed}
json.dump(res,open('theorem_coverage_200.json','w'),indent=1)
