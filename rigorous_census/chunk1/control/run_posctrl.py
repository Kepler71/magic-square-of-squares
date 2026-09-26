# Claude 26.09: положительный контроль конвейера dem_rigorous.run_T на Y-точках (z*, −z* обе дают точки E_T ранга 1).
# Для каждой пары вычисляются n, n′ (P_T(±z*) = nP0+t), класс A/B; затем полный run_T; z* обязано быть в Z_A ∪ Z_B,
# в случае B (n′ = ±n) — среди корней многочленов случая B; в случае A — |n| ≤ M0.
from sage.all import *
import json, sys, time
sys.path.insert(0,'/home/kep/magicKube/rigorous_census/chunk1')
import dem_rigorous as D
Y=json.load(open('Ypoints_14.json'))
def index_of(Emin,P0,tors,Q):
    if Q.order()!=oo: return 0
    n=int(round(sqrt(float(Q.height())/float(P0.height()))))
    for m in (n,-n):
        for t in tors:
            if m*P0+t==Q: return m
    raise ValueError('не найден индекс')
out=[]; nA=nB=0; seen=set()
for y in Y:
    T=y['T']; z=QQ(y['z'])
    if prod(1+c*z for c in T)==0 or prod(1-c*z for c in T)==0: continue
    key=(tuple(T),abs(z))
    if key in seen: continue
    seen.add(key)
    E1,L,c0=D.model(T); Emin=E1.minimal_model(); iso=E1.isomorphism_to(Emin)
    M,_,_=D.mobius(E1,Emin,iso,L,c0)
    er=pari(Emin).ellrank(); P=[Emin([QQ(c) for c in p]) for p in er[3]]; P=[p for p in P if p.order()==oo][0]
    P0=Emin.saturation([P])[0][0]; tors=Emin.torsion_points()
    Qp=Emin.lift_x(D.xmin_of_z(M,z)); Qm=Emin.lift_x(D.xmin_of_z(M,-z))
    n=index_of(Emin,P0,tors,Qp); n2=index_of(Emin,P0,tors,Qm)
    cls='B' if n*n==n2*n2 else 'A'
    if cls=='A' and (nA>=6 or min(abs(n),abs(n2))==0 and nA>=3): continue
    if cls=='B' and nB>=6: continue
    S=sorted(set(T+[-c for c in T]))
    res=D.run_T(S,T,do_sallows_z=z)
    zB=set(QQ(v) for o in res['caseB'] for v in o.get('z',[]))
    inB=(z in zB) or (-z in zB)
    ok=res['control_z_found'] and (inB if cls=='B' else max(abs(n),abs(n2))<=res['M0'])
    if cls=='A': nA+=1
    else: nB+=1
    rec=dict(T=T,z=str(z),n=n,n_prime=n2,cls=cls,M0=res['M0'],B_bezout=res['B_bezout'],found=res['control_z_found'],in_caseB_roots=inB,PASS=bool(ok))
    print(rec,flush=True); out.append(rec)
    if nA>=6 and nB>=6: break
print('итог: A',nA,'B',nB,'все прошли:',all(r['PASS'] for r in out),flush=True)
json.dump(out,open('posctrl_results.json','w'),indent=1)
