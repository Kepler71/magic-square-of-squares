# Положительный контроль (скептик Claude/Opus, 2026-09-26): тот же вывод при ОСЛАБЛЕННОЙ гипотезе.
# Если убрать часть линий, классификатор должен находить пары PASS вне цели (иначе нули в основном
# прогоне ничего бы не значили). Заодно видно, какие линии реально нужны локально при p=2,3.
import json
CELLS=[(i,j) for i in (-1,0,1) for j in (-1,0,1)]
LINES=[[(i,-1),(i,0),(i,1)] for i in (-1,0,1)]+[[(-1,j),(0,j),(1,j)] for j in (-1,0,1)]
LINES+=[[(-1,-1),(0,0),(1,1)],[(-1,1),(0,0),(1,-1)]]
NAMES=['row-1','row0','row+1','col-1','col0','col+1','diag','anti']
def vint(n,p):
    v=0
    while n%p==0: n//=p; v+=1
    return v,n
def classify(p,k,m,B,C,lines):
    mod=p**k; need=3 if p==2 else 1; um=8 if p==2 else 3
    info={}
    for (i,j) in CELLS:
        if (i,j)==(0,0): info[(i,j)]=(False,m,1,10**6); continue
        N=(pow(p,m,mod)+i*B+j*C)%mod
        if N==0: info[(i,j)]=(True,0,1,0); continue
        v,u=vint(N,p); info[(i,j)]=(False,v,u%um,k-v)
    excl=False; allsq=True
    for ln in lines:
        if any(info[c][0] for c in ln): allsq=False; continue
        V=m+sum(info[c][1] for c in ln); U=1
        for c in ln: U=(U*info[c][2])%um
        known=min(info[c][3] for c in ln)>=need
        if V%2 or (known and U!=1): excl=True
        if not(V%2==0 and known and U==1): allsq=False
    return 'EXCL' if excl else ('PASS' if allsq else 'UND')
subsets={'all8':list(range(8))}
for t in range(8): subsets['drop_'+NAMES[t]]=[x for x in range(8) if x!=t]
subsets['central4']=[1,4,6,7]; subsets['rows_cols6']=[0,1,2,3,4,5]
res={}
for p,k in ((2,6),(3,4)):
    T=8 if p==2 else 3
    for name,sub in subsets.items():
        lines=[LINES[t] for t in sub]
        pn=un=0; ex=None
        for m in range(0,7):
            for B in range(p**k):
                for C in range(p**k):
                    if m>0 and B%p==0 and C%p==0: continue
                    st=classify(p,k,m,B,C,lines)
                    tgt=(m==0 and B%T==0 and C%T==0)
                    if st=='PASS' and not tgt:
                        pn+=1
                        if ex is None: ex=(m,B,C)
                    if st=='UND' and not tgt: un+=1
        res[f'p{p}_{name}']=dict(PASS_not_target=pn,UND_not_target=un,example_m_B_C=ex)
        print(f'p={p} k={k} {name:12s}: PASS вне цели={pn:6d}  UND вне цели={un:6d}  пример (m,B,C)={ex}',flush=True)
json.dump(res,open('p23_subsets.json','w'),indent=1)
