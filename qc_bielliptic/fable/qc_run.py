# Fable 15.09.2026. Квадратичный Шаботи (код Bianchi–Padurariu, патч под Sage 10.9) для биэллиптической кривой
#   C: w^2 = (x^2-a^2)(x^2-b^2)(x^2-c^2),  x = 1/z  (z — параметр девяти клеток 1 ± kz).
# Использование: python3 qc_run.py a b c p n [genfile.json]
# Известные точки: x = ∞ (z=0, две точки), (±a,0),(±b,0),(±c,0). Любая другая рациональная точка = кандидат на шесть квадратов.
import sys, time, json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
from models import fpoly, E12
a,b,c,p,n = map(int, sys.argv[1:6]); p=ZZ(p); n=ZZ(n)
f=fpoly(a,b,c); E1,E2=E12(f)
gens=None
if len(sys.argv)>6:
    G=json.load(open(sys.argv[6]))[f'{a},{b},{c}']
    gens=[sage_eval(G['E1'][0].replace(':',',').strip('()')), sage_eval(G['E2'][0].replace(':',',').strip('()'))]
    gens=[[QQ(t) for t in g][:2] for g in gens]
    print('образующие (минимальные модели):',gens,flush=True)
print(f'C: w^2 = {f};  p = {p}, n = {n}',flush=True)
print('E1min',E1.minimal_model().ainvs(),'E2min',E2.minimal_model().ainvs(),flush=True)
from omega_direct import omega_direct
t0=time.time()
Om,per=omega_direct(f,p,n)
print('Ω_direct построено за',round(time.time()-t0,1),'с; |Ω| =',len(Om),flush=True)
rat,other = quadratic_chabauty_bielliptic(f,p,n,Omega=Om,gens=gens)
dt=time.time()-t0
known=set([(oo,0)]+[(QQ(s*v),QQ(0)) for v in (a,b,c) for s in (1,-1)])
ratz=[]
for P in rat:
    if P[2]==0: ratz.append('z=0'); continue
    x=QQ(P[0]); ratz.append(f'z={1/x}' if x!=0 else 'z=inf')
new=[P for P in rat if P[2]!=0 and (QQ(P[0]),QQ(0)) not in known]
print('время',round(dt,1),'с')
print('рациональные точки, найденные QC:',rat)
print('в терминах z:',sorted(set(ratz)))
print('НОВЫЕ рациональные точки (не из известных восьми):',new)
print('нераспознанные p-адические точки:',len(other))
for P in other: print('  ',P)
json.dump(dict(abc=[a,b,c],p=int(p),n=int(n),sec=float(round(dt,1)),rational=[str(P) for P in rat],z=sorted(set(ratz)),new=[str(P) for P in new],
               other=[str(P) for P in other],n_other=len(other)),open(f'qc_{a}_{b}_{c}_p{p}.json','w'),indent=1)
