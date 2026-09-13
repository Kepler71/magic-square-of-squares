import json
from fractions import Fraction as F
from collections import Counter,defaultdict
T=json.load(open('full_table_2_48.json')); rows=T['rows']; CL=T['classes']
cod=json.load(open('codex_table.json'))
ks=sorted(rows,key=lambda s:F(s))
print('slopes',len(ks))
# A. сила классов
print('\n== Классы троек: число k с rank 0 (ellrank r_hi=0), с Sel2-избытком 0, с p0 не кручением, rk=None ==')
def stats(kind):
    out=[]
    for cls in CL[kind]:
        key=kind+'|'+','.join(cls)
        r0=s0=p0n=err=lo1=0
        for k in ks:
            e=rows[k][key]
            if e['rk'] is None: err+=1; continue
            if e['rk'][1]==0: r0+=1
            if e['rk'][2]==0: s0+=1
            if e['rk'][0]>=1: lo1+=1
            if not e['p0tors']: p0n+=1
        out.append((r0,s0,lo1,p0n,err,key))
    out.sort(reverse=True)
    for r0,s0,lo1,p0n,err,key in out: print(f'{key:32s} rank0={r0:4d} sel0={s0:4d} rank>=1={lo1:4d} p0nontors={p0n:4d} err={err}')
stats('C'); print(); stats('Q')
# B. свидетели Codex по классам
print('\n== Свидетели Codex по классам ==')
wc=Counter((w['kind'],w['cls']) for w in cod['witness'].values())
for (kind,cls),n in wc.most_common(): print(kind,cls,n)
# для каждого k: сколько классов дают rank 0
print('\n== Распределение числа классов с rank 0 на один k ==')
cnt=Counter(); cntQ=Counter(); none=[]
for k in ks:
    c=sum(1 for cls in CL['C'] if rows[k]['C|'+','.join(cls)]['rk'] and rows[k]['C|'+','.join(cls)]['rk'][1]==0)
    q=sum(1 for cls in CL['Q'] if rows[k]['Q|'+','.join(cls)]['rk'] and rows[k]['Q|'+','.join(cls)]['rk'][1]==0)
    cnt[c]+=1; cntQ[q]+=1
    if c==0: none.append((k,q))
print('cubic:',sorted(cnt.items())); print('quartic:',sorted(cntQ.items()))
print('k без кубического rank0:',len(none),none)
