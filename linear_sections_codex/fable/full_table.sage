# Полная таблица: для каждого k=r/s (0<k<1, s<=MAXDEN) и каждого класса троек/четвёрок клеток
# ellrank минимальной модели: [r_lo, r_hi, s_excess] (s_excess = dim Sel_2 - dim E[2] по PARI),
# плюс: является ли образ точки p=0 точкой кручения (иначе ранг >= 1 автоматически).
import json,sys,time
from itertools import combinations
MINDEN=int(sys.argv[1]); MAXDEN=int(sys.argv[2]); OUT=sys.argv[3]
KINDS=sys.argv[4] if len(sys.argv)>4 else 'CQ'
R=PolynomialRing(QQ,'p'); p=R.gen()
ORDER=['1','k','1+k','1-k','-1','-k','-(1+k)','-(1-k)']
NEG={'1':'-1','k':'-k','1+k':'-(1+k)','1-k':'-(1-k)'}; NEG.update({v:u for u,v in NEG.items()})
def canon(syms):
    a=tuple(sorted(syms,key=ORDER.index)); b=tuple(sorted((NEG[s] for s in syms),key=ORDER.index))
    return min(a,b,key=lambda t:[ORDER.index(x) for x in t])
CLASSES={'C':sorted({canon(t) for t in combinations(ORDER,3)},key=lambda t:[ORDER.index(x) for x in t]),
         'Q':sorted({canon(t) for t in combinations(ORDER,4)},key=lambda t:[ORDER.index(x) for x in t])}
assert len(CLASSES['C'])==28 and len(CLASSES['Q'])==38
def val(sym,k): return {'1':QQ(1),'k':k,'1+k':1+k,'1-k':1-k,'-1':QQ(-1),'-k':-k,'-(1+k)':-1-k,'-(1-k)':k-1}[sym]
def curve_and_p0(cls,k):
    cs=[val(s,k) for s in cls]
    if len(cs)==3:
        f=prod(1+c*p for c in cs); A=f[3]
        E=EllipticCurve([0,f[2],0,A*f[1],A*A]); P0=E(0,A)   # p=0,Y=1 -> X=0,V=A
    else:
        a=cs[0]; root=-1/a; x=p
        g=a*prod((1-c/a)*x+c for c in cs[1:]); A=g[3]
        E=EllipticCurve([0,g[2],0,A*g[1],A*A*g[0]])
        # p=0 -> x=1/(0-root)=-1/root=a ; w=Y/(p-root)^2 = 1/root^2 = a^2 ; X=A x, V=A w
        P0=E(A*a,A*a*a)
    return E,P0
cache={}; rows={}; t0=time.time()
slopes=sorted(set(QQ(a)/b for b in range(MINDEN,MAXDEN+1) for a in range(1,b) if gcd(a,b)==1))
for n,k in enumerate(slopes):
    row={}
    for kind in KINDS:
        for cls in CLASSES[kind]:
            try: E,P0=curve_and_p0(cls,k)
            except ArithmeticError:
                row[kind+'|'+','.join(cls)]={'key':None,'rk':None,'p0tors':None,'j':None}; continue
            M=E.minimal_model(); key=','.join(map(str,M.ainvs()))
            if key not in cache:
                try:
                    alarm(20); rk=M.pari_curve().ellrank(); cancel_alarm()
                    cache[key]=[int(rk[0]),int(rk[1]),int(rk[2])]
                except Exception as e:
                    cancel_alarm(); cache[key]=None
            row[kind+'|'+','.join(cls)]={'key':key,'rk':cache[key],'p0tors':bool(P0.has_finite_order()),'j':str(M.j_invariant())}
    rows[str(k)]=row
    if n%25==0:
        print(n,len(slopes),k,round(time.time()-t0),flush=True)
        json.dump({'rows':rows,'classes':CLASSES},open(OUT,'w'))
json.dump({'rows':rows,'classes':CLASSES},open(OUT,'w'))
print('done',len(rows),'curves',len(cache),round(time.time()-t0))
