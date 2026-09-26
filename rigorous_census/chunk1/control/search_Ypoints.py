# Claude 26.09: поиск положительного контроля уровня Y = C_T ×_{P¹} C_{−T}: T (3 клетки, |c| ≤ NMAX), E_T ранга 1,
# точки Q = nP0 + t (|n| ≤ 8), z = z(Q) ≠ 0, Q не кручение, ∏_T(1 − λz) — квадрат (тогда и P_T(−z) рациональна).
from sage.all import *
import itertools, json, time, sys
sys.path.insert(0,'/home/kep/magicKube/rigorous_census/chunk1')
from dem_rigorous import model, mobius, z_of_xmin, xcoord
NMAX=int(sys.argv[1]); TLIM=float(sys.argv[2])
vals=[c for c in range(-NMAX,NMAX+1) if c!=0]
t0=time.time(); found=[]; nrank1=0; ntot=0
for T in itertools.combinations(vals,3):
    if time.time()-t0>TLIM: break
    T=list(T)
    if sorted(-x for x in T)==T or gcd(T)!=1: continue
    try: E1,L,c0=model(T)
    except Exception: continue
    ntot+=1
    Emin=E1.minimal_model()
    er=pari(Emin).ellrank()
    if not (int(er[0])==int(er[1])==1): continue
    nrank1+=1
    iso=E1.isomorphism_to(Emin)
    M,_,_=mobius(E1,Emin,iso,L,c0)
    P=[Emin([QQ(c) for c in p]) for p in er[3]]; P=[p for p in P if p.order()==oo][0]
    tors=Emin.torsion_points()
    Q=Emin(0)
    for n in range(0,9):
        for t in tors:
            R=Q+t
            if n==0: continue
            z=z_of_xmin(M,xcoord(R))
            if z is None or z==0: continue
            pr=prod(1-c*z for c in T)
            if pr>=0 and pr.is_square():
                found.append(dict(T=[int(x) for x in T],z=str(z),n=n,t=str(t)))
                print('Y-точка:',found[-1],flush=True)
        Q=Q+P
print('наборов',ntot,'ранг 1:',nrank1,'Y-точек (нетривиальных):',len(found),'время',round(time.time()-t0),flush=True)
json.dump(found,open(f'Ypoints_{NMAX}.json','w'))
