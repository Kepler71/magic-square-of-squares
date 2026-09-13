import json
from fractions import Fraction as F
from math import gcd
from collections import Counter
rows=json.load(open('grid_rows.json')); G={(x['m'],x['n']):x for x in rows}
def red(a,b):
    g=gcd(a,b); a//=g; b//=g
    return (min(a,b),max(a,b))
def pairs(r,s):
    return {'{r,s}':(r,s),'{r,s-r}':(r,s-r),'{s,s-r}':(s,s-r),'{s,r+s}':(s,r+s),'{r,r+s}':(r,r+s),
       'Q{s-r,r+s}':(s-r,r+s),'Q{s,|2r-s|}':(s,abs(2*r-s)),'Q{s,2r+s}':(s,2*r+s),'Q{r,2s-r}':(r,2*s-r),'Q{r,2s+r}':(r,2*s+r),'Q{s,2r}':(s,2*r),'Q{r,2s}':(r,2*s)}
def status(k):
    r,s=k.numerator,k.denominator
    st={}
    for nm,(a,b) in pairs(r,s).items():
        a,b=red(a,b)
        if a==b or a==0: st[nm]='deg'; continue
        x=G.get((a,b))
        if x is None: st[nm]='?'; continue
        st[nm]='0' if x['r2']==0 else ('+' if x['r1']>=1 else 'u')   # 0: ранг 0 доказан; +: ранг>=1 доказан; u: не определён
    return st
for lo,hi in ((2,48),(49,100),(101,133),(134,200)):
    slopes=[F(r,s) for s in range(lo,hi+1) for r in range(1,s) if gcd(r,s)==1]
    c=Counter(); hopeless=[]; undet=[]; cubic_hopeless=0
    for k in slopes:
        st=status(k); v=list(st.values())
        cub=[st[n] for n in st if not n.startswith('Q')]
        if '0' in v: c['closed']+=1
        elif '?' in v: c['unknown(grid)']+=1; 
        elif 'u' in v: c['undetermined']+=1; undet.append(str(k))
        else: c['all rank>=1']+=1; hopeless.append(str(k))
        if '0' not in cub and '?' not in cub and 'u' not in cub: cubic_hopeless+=1
        if '0' in cub: c['closed by cubic']+=1
    print(f'знаменатели {lo}..{hi}: наклонов {len(slopes)}',dict(c),'кубические 5 пар все ранга>=1:',cubic_hopeless)
    print('   все 12 кривых ранга >=1 (эллиптический метод бессилен):',hopeless[:40], '...' if len(hopeless)>40 else '')
    print('   не определено (r1=0<r2 у некоторых, ранга 0 нет):',undet[:20])
# особые: 1/2,1/3,2/3,1/4,3/4 и квартические 56
print('\nОсобые наклоны:')
for ks in ['1/3','2/3','1/4','3/4','1/2']:
    print(ks,status(F(ks)))
cod=json.load(open('codex_table.json'))
qk=[k for k,w in cod['witness'].items() if w['kind']=='Q']
print('\nКвартические у Codex:',len(qk))
cnt=Counter()
for ks in qk:
    st=status(F(ks)); cnt[tuple(sorted(n for n,v in st.items() if v=='0'))]+=1
print('какие пары ранга 0 у квартических наклонов (комбинации):'); 
for comb,n in cnt.most_common(15): print('  ',n,comb)
print('число наклонов ≤48, у которых ровно один из 12 ранга 0:',sum(1 for s in range(2,49) for r in range(1,s) if gcd(r,s)==1 and list(status(F(r,s)).values()).count('0')==1))
