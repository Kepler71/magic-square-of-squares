# Полный перебор вычетов для шагов 3-4 теоремы об остаточных классах (скептик Claude/Opus, 2026-09-26).
# Утверждение (в сильной локальной форме, для ЛЮБЫХ b,c in Q_p, без структуры b=rz,c=sz):
#   если 8 произведений по линиям арифметической сетки 1+ib+jc - квадраты в Q_p (p=2 или 3),
#   то v_2(b),v_2(c) >= 3  (p=2),   v_3(b),v_3(c) >= 1  (p=3).
# Из этой формы следует нужная: одно из r,s - единица в Z_p, значит v_p(z)=min(v(b),v(c)).
# Параметризация: m = порядок полюса = -min(v(b),v(c),0); b=B/p^m, c=C/p^m, B,C in Z_p;
# при m>=1 хотя бы одно из B,C - единица. B,C известны по модулю p^k (перебор ВСЕХ вычетов).
# Числитель клетки N=p^m+iB+jC известен mod p^k; класс линии = класс(p^m * prod N).
# Статус линии: NSQ (не квадрат при ЛЮБОМ подъёме), SQ (квадрат при любом подъёме), UNK.
# Статус пары: EXCL (есть NSQ) / PASS (все 8 SQ) / UND.
# Опровержение = пара PASS (или UND, доведённая до PASS глубже), у которой b или c не делится на 8 (3).
import numpy as np, sys, json, time
p=int(sys.argv[1]); KMAX=int(sys.argv[2]); MMAX=int(sys.argv[3])
TGT=8 if p==2 else 3
CELLS=[(i,j) for i in (-1,0,1) for j in (-1,0,1)]
LINES=[[(i,-1),(i,0),(i,1)] for i in (-1,0,1)]+[[(-1,j),(0,j),(1,j)] for j in (-1,0,1)]
LINES+=[[(-1,-1),(0,0),(1,1)],[(-1,1),(0,0),(1,-1)]]

def unit_is_square(umod):   # umod: остаток единицы по модулю 8 (p=2) или 3 (p=3)
    return (umod==1)

def analyse(k,m,Bblk,C):
    mod=p**k
    B=Bblk[:,None]; Cc=C[None,:]
    shape=(len(Bblk),len(C))
    info={}
    for (i,j) in CELLS:
        if (i,j)==(0,0): continue
        N=(pow(p,m,mod)+i*B+j*Cc)%mod
        N=np.broadcast_to(N,shape).astype(np.int64)
        zero=(N==0)
        v=np.zeros(shape,dtype=np.int64)
        for t in range(1,k):
            v+=((N%(p**t))==0)
        pw=np.array([p**t for t in range(k+1)],dtype=np.int64)
        u=np.where(zero,1,N//pw[v])
        prec=k-v
        um=u%(8 if p==2 else 3)
        info[(i,j)]=(zero,v,um,prec)
    EXCL=np.zeros(shape,bool); ALLSQ=np.ones(shape,bool)
    for ln in LINES:
        anyzero=np.zeros(shape,bool); V=np.full(shape,m,dtype=np.int64)
        U=np.ones(shape,dtype=np.int64); P=np.full(shape,10**6,dtype=np.int64)
        for c in ln:
            if c==(0,0):
                V+=m; continue
            zero,v,um,prec=info[c]
            anyzero|=zero; V+=v; U=(U*um)%(8 if p==2 else 3); P=np.minimum(P,prec)
        known=(P>=3) if p==2 else (P>=1)
        nsq=(~anyzero)&(((V%2)==1)|(known&~unit_is_square(U)))
        sq=(~anyzero)&((V%2)==0)&known&unit_is_square(U)
        EXCL|=nsq; ALLSQ&=sq
    PASS=ALLSQ&~EXCL
    # индивидуальные клетки (для контроля вывода "все девять клеток - квадраты")
    cellsq=np.ones(shape,bool)
    for c,(zero,v,um,prec) in info.items():
        known=(prec>=3) if p==2 else (prec>=1)
        # класс клетки = класс(N/p^m) = класс(p^m N): оценка v+m
        cellsq&=(~zero)&(((v+m)%2)==0)&known&unit_is_square(um)
    return EXCL,PASS,cellsq

out=[]
t0=time.time()
for k in range(1,KMAX+1):
    mod=p**k
    C=np.arange(mod,dtype=np.int64)
    for m in range(0,MMAX+1):
        cnt=dict(total=0,EXCL=0,PASS=0,UND=0,PASS_not_target=0,UND_not_target=0,
                 EXCL_target=0,PASS_cell_not_sq=0)
        unresolved=[]
        blk=max(1,min(mod,2_000_000//mod))
        for b0 in range(0,mod,blk):
            Bblk=np.arange(b0,min(mod,b0+blk),dtype=np.int64)
            EX,PA,CS=analyse(k,m,Bblk,C)
            Bg=np.broadcast_to(Bblk[:,None],EX.shape); Cg=np.broadcast_to(C[None,:],EX.shape)
            valid=np.ones(EX.shape,bool) if m==0 else ((Bg%p)!=0)|((Cg%p)!=0)
            UN=valid&~EX&~PA
            if m==0 and k>=(3 if p==2 else 1):
                tgt=((Bg%TGT)==0)&((Cg%TGT)==0)
            else:
                tgt=np.zeros(EX.shape,bool)   # при полюсе цель недостижима по определению
            decidable=(m>0) or k>=(3 if p==2 else 1)
            cnt['total']+=int(valid.sum()); cnt['EXCL']+=int((valid&EX).sum())
            cnt['PASS']+=int((valid&PA).sum()); cnt['UND']+=int(UN.sum())
            if decidable:
                cnt['PASS_not_target']+=int((valid&PA&~tgt).sum())
                cnt['UND_not_target']+=int((UN&~tgt).sum())
                cnt['EXCL_target']+=int((valid&EX&tgt).sum())
                cnt['PASS_cell_not_sq']+=int((valid&PA&~CS).sum())
                if k==KMAX:
                    idx=np.argwhere(UN&~tgt)
                    for a,bb in idx[:20000]:
                        unresolved.append((int(Bg[a,bb]),int(Cg[a,bb])))
        rec=dict(p=p,k=k,m=m,**cnt,decidable=decidable)
        if k==KMAX: rec['unresolved_sample']=unresolved[:50]; rec['unresolved_n']=len(unresolved)
        out.append(rec)
        print(f"p={p} k={k} m={m}: {cnt} decidable={decidable} [{time.time()-t0:.1f}s]",flush=True)
        if k==KMAX and unresolved:
            with open(f'unresolved_p{p}_k{k}_m{m}.json','w') as f: json.dump(unresolved,f)
bad=sum(r['PASS_not_target'] for r in out)
print('ИТОГ p=%d: PASS вне цели (опровержение) = %d; EXCL внутри цели (ошибка кода) = %d; '
      'PASS с неквадратной клеткой = %d'%(p,bad,sum(r['EXCL_target'] for r in out),
      sum(r['PASS_cell_not_sq'] for r in out)))
with open(f'p23_full_p{p}.json','w') as f: json.dump(out,f,indent=1)
