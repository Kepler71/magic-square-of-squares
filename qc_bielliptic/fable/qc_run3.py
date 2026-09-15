# Fable 15.09.2026, после рецензии Codex. QC с точным Ω^nine (omega_exact), всеми наборами тегов и сертификатом Штрассмана.
# Использование: python3 qc_run3.py a b c d p n [genfile.json]
import sys, time, json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
import qc_load; from qc_load import *
qc_load.DEBUG_RHO = False
from models import fpoly
from omega_exact import omega_exact
a,b,c,d,p,n = map(int, sys.argv[1:7]); p=ZZ(p); n=ZZ(n)
f=fpoly(a,b,c); gens=None
if len(sys.argv)>7:
    G=json.load(open(sys.argv[7]))[f'{a},{b},{c}']
    gens=[[QQ(t) for t in sage_eval(g[0].replace(':',',').strip('()'))][:2] for g in (G['E1'],G['E2'])]
print(f'C: w^2 = {f};  p = {p}, n = {n}, cells = {[a,b,c,d]}',flush=True)
t0=time.time(); Om,tags,qs,per,per9,certs=omega_exact(f,p,n,[a,b,c,d]); print('Ω^nine за',round(time.time()-t0,1),'с',flush=True)
tag_inf=tuple(QQ(certs[q]['tailinf']['c']) for q in qs)
inf_idx=[l for l,tg in enumerate(tags) if any(tuple(QQ(c) for c in t_)==tag_inf for t_ in tg)]
qc_load.INF_OMEGA_INDEX = inf_idx[0] if len(inf_idx)==1 else None
print('тег точки ∞:',tag_inf,'— индекс ω:',qc_load.INF_OMEGA_INDEX,flush=True)
uncert_classes={int(q):c.get('uncertified_classes',[]) for q,c in certs.items() if c.get('uncertified_classes')}
print('несертифицированные классы перечисления:',uncert_classes,flush=True)
qc_load.STRASSMANN.clear(); t0=time.time()
rat,other = quadratic_chabauty_bielliptic(f,p,n,Omega=Om,gens=gens,omega_info=True)
dt=time.time()-t0; print('QC за',round(dt,1),'с',flush=True)
K=Qp(p,n)
def sq(v): return v!=0 and v.valuation()%2==0 and (v/K(p)**v.valuation()).is_square()
res=[]; nrat=noth=surv=0
for l in range(len(Om)):
    for P in rat[l]:
        nrat+=1; res.append(dict(omega=l,tags=[[str(t) for t in tg] for tg in tags[l]],rational=True,x=str(P[0]) if P[2]!=0 else 'inf'))
    for P in other[l]:
        noth+=1; x=P[0]; ok=all(sq((x+s*k)/x) for k in (a,b,c,d) for s in (1,-1)); surv+=ok
        res.append(dict(omega=l,tags=[[str(t) for t in tg] for tg in tags[l]],rational=False,x=str(x),xprec=int(x.precision_absolute()),nine_local=ok))
S=qc_load.STRASSMANN; ncert=sum(1 for t in S if t[5]); bad=[t for t in S if not t[5]]
print('рациональных (по Ω):',nrat,' нераспознанных:',noth,' выживших локально при p:',surv)
print('рациональные:',sorted(set(r['x'] for r in res if r['rational'])))
tagset=set(); [tagset.update(tuple(tg) for tg in r['tags']) for r in res if not r['rational'] and r['nine_local']]
print('выжившие теги (все наборы):',sorted(tagset))
print(f'Сертификат корней: (диск,ω)-пар {len(S)}, сертифицировано {ncert}, несертифицировано {len(bad)}:', bad[:10])
print('ИТОГ: корней ρ−ω в дисках нет, кроме', 'двух точек на бесконечности (z=0, вырожденный квадрат)' if nrat==2 and noth==0 else f'{nrat} рациональных и {noth} нераспознанных')
json.dump(dict(abc=[a,b,c],d=d,p=int(p),n=int(n),qs=[int(q) for q in qs],sec=float(dt),n_omega=len(Om),
               per={int(q):[str(t) for t in sorted(v)] for q,v in per.items()},per9={int(q):[str(t) for t in sorted(v)] for q,v in per9.items()},
               certs={int(q):{k:(v if not isinstance(v,dict) else {kk:str(vv) for kk,vv in v.items()}) for k,v in cc.items()} for q,cc in certs.items()},
               strassmann=dict(pairs=len(S),certified=ncert,uncertified=[list(t) for t in bad]),points=res),open(f'qc3_{a}_{b}_{c}_p{p}.json','w'),indent=1,default=str)
