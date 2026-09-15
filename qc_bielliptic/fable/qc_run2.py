# Fable 15.09.2026. Как qc_run.py, но с omega_info=True: точки по элементам Ω (с тегами (c_q)_q), для решета по двум простым.
# Использование: python3 qc_run2.py a b c d p n [genfile.json]   (d — четвёртое значение наклона, для проверки клеток 1±dz)
import sys, time, json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
from models import fpoly, E12
from omega_direct import omega_direct
a,b,c,d,p,n = map(int, sys.argv[1:7]); p=ZZ(p); n=ZZ(n)
f=fpoly(a,b,c)
gens=None
if len(sys.argv)>7:
    G=json.load(open(sys.argv[7]))[f'{a},{b},{c}']
    gens=[[QQ(t) for t in sage_eval(g[0].replace(':',',').strip('()'))][:2] for g in (G['E1'],G['E2'])]
print(f'C: w^2 = {f};  p = {p}, n = {n}, d = {d}',flush=True)
t0=time.time(); Om,per=omega_direct(f,p,n,cells=[a,b,c,d]); tags=omega_direct.last_tags; qs=omega_direct.last_qs
print('V_q^nine:',{int(q):sorted(v) for q,v in omega_direct.last_per9.items()},flush=True)
print('Ω_direct за',round(time.time()-t0,1),'с; |Ω| =',len(Om),'; совпадений ω:',sum(len(t)>1 for t in tags),flush=True)
t0=time.time()
rat,other = quadratic_chabauty_bielliptic(f,p,n,Omega=Om,gens=gens,omega_info=True)
dt=time.time()-t0; print('QC за',round(dt,1),'с',flush=True)
K=Qp(p,n)
def sq(v):  # v ∈ Q_p — ненулевой квадрат?
    if v==0: return False
    return v.valuation()%2==0 and (v/K(p)**v.valuation()).is_square()
res=[]; nrat=0; noth=0; surv=0
for l in range(len(Om)):
    for P in rat[l]:
        nrat+=1; res.append(dict(omega=l,tag=[str(t) for t in tags[l][0]],rational=True,x=str(P[0]) if P[2]!=0 else 'inf'))
    for P in other[l]:
        noth+=1; x=P[0]
        cells=[sq((x+s*k)/x) for k in (a,b,c,d) for s in (1,-1)]   # 1 ± kz = (x ± k)/x
        ok=all(cells); surv+=ok
        res.append(dict(omega=l,tag=[str(t) for t in tags[l][0]],rational=False,x=str(x),xprec=int(x.precision_absolute()),nine_local=ok))
print('рациональных (по Ω):',nrat,' нераспознанных:',noth,' из них с восемью клетками — квадратами в Q_p:',surv)
print('рациональные:',sorted(set(r['x'] for r in res if r['rational'])))
print('выжившие теги:',sorted(set(tuple(r['tag']) for r in res if not r['rational'] and r['nine_local'])))
json.dump(dict(abc=[a,b,c],d=d,p=int(p),n=int(n),qs=[int(q) for q in qs],sec=float(dt),n_omega=len(Om),points=res),open(f'qc2_{a}_{b}_{c}_p{p}.json','w'),indent=1)
