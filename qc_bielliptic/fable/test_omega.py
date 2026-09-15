# Проверка Ω_direct: (1) направление to_min на рациональной точке; (2) Ω_direct ⊆ Ω_code попримно (BP22-пример, BD18, контроль 9,22,31);
# (3) QC с Ω_direct даёт те же точки, что с Ω_code (BP22-пример, p=5).
import sys,time; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from qc_load import *
from omega_direct import *
from models import fpoly
def code_perprime(f,p,n):
    a0=f[0]; a6=f[6]; a4=f[4]; a2=f[2]
    E1=EllipticCurve([0,a4,0,a2*a6,a0*a6**2]); E2=EllipticCurve([0,a2,0,a0*a4,a0**2*a6]); E1min=E1.minimal_model(); E2min=E2.minimal_model()
    K=Qp(p,n); bad_primes_1,W1=local_heights_at_bad_primes_new(E1min,E1,K); bad_primes_2,W2=local_heights_at_bad_primes_new(E2min,E2,K)
    bad_primes_3=f7(ZZ(a0).prime_factors()+ZZ(a6).prime_factors()); bps=f7(bad_primes_1+bad_primes_2+bad_primes_3); out={}
    for q in bps:
        if q in bad_primes_1: WE1q=W1[bad_primes_1.index(q)]
        else: WE1q=[] if ((q==2 and E1.Np(2)==1) or (q==3 and E1.Np(3)==1) or (q==2 and E1.has_split_multiplicative_reduction(2) and E1.kodaira_symbol(2)==KodairaSymbol(5))) else [0]
        if q in bad_primes_2: WE2q=W2[bad_primes_2.index(q)]
        else: WE2q=[] if ((q==2 and E2.Np(2)==1) or (q==3 and E2.Np(3)==1) or (q==2 and E2.has_split_multiplicative_reduction(2) and E2.kodaira_symbol(2)==KodairaSymbol(5))) else [0]
        evenq=[2*i for i in range(((-(a6).valuation(q))/2).ceil(), ((a0).valuation(q)/2).floor()+1)]
        Wq=f7([-w1-a0.valuation(q)*log(K(q))+O(p**n) for w1 in WE1q]+[-w[0]+w[1]-w[2]*log(K(q))+O(p**n) for w in itertools.product(WE1q,WE2q,evenq)]+[w2+a6.valuation(q)*log(K(q))+O(p**n) for w2 in WE2q])
        out[q]=Wq
    return out
R=PolynomialRing(QQ,'x'); x=R.gen()
# (1) направление изоморфизма
f=fpoly(9,22,31); E1=EllipticCurve([0,f[4],0,f[2]*f[6],f[0]*f[6]**2]); E1min=E1.minimal_model(); i1=E1.isomorphism_to(E1min)
P=E1min.gens(proof=False)[0]; Pn=i1.inverse()(P) if hasattr(i1,'inverse') else E1min.isomorphism_to(E1)(P)
print('to_min совпадает с iso:', to_min(i1,Pn[0],Pn[1])==(P[0],P[1]))
for name,f,p in [('BP22 99856',x**6+22*x**4-19*x**2+4,5),('BD18',x**6-2*x**4-x**2+1,5),('BP22 24025', -x**6+11*x**4-19*x**2+25,7)]:
    print('==',name,'p =',p,flush=True); t0=time.time()
    Om,per=omega_direct(f,ZZ(p),20)
    cp=code_perprime(f,ZZ(p),20); K=Qp(p,20)
    for q in sorted(per):
        mine=[c*K(q).log() for c in per[q]]; code=cp.get(q,[])
        sub=all(any(m==w for w in code) for m in mine)
        print(f'  q={q}: |V_q|={len(mine)} ⊆ код({len(code)}): {sub}')
    print('  время',round(time.time()-t0,1),'с',flush=True)
# (3) QC с Ω_direct на BP22-примере
f=x**6+22*x**4-19*x**2+4; Om,per=omega_direct(f,ZZ(5),20,verbose=False)
r1,o1=quadratic_chabauty_bielliptic(f,ZZ(5),ZZ(20)); r2,o2=quadratic_chabauty_bielliptic(f,ZZ(5),ZZ(20),Omega=Om)
print('QC(BP22, p=5): Ω_code:',len(r1),len(o1),' Ω_direct:',len(r2),len(o2),' rat равны:',sorted(r1)==sorted(r2),' other равны:',sorted(map(str,o1))==sorted(map(str,o2)))
