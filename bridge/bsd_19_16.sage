m=19;n=16;s=QQ(m^2+n^2)/2;b=s*m^2*n^2
R.<X>=QQ[]; f=(X+b)*(X+s*m^4)*(X+s*n^4)
E=EllipticCurve([0,f[2],0,f[1],f[0]]).minimal_model()
print("Emin:",E.ainvs())
N=E.conductor(); print("N =",N, "=", N.factor())
T=E.torsion_order(); print("torsion order:",T)
tam=[(p,E.tamagawa_number(p)) for p in N.prime_factors()]
print("tamagawa:",tam)
tp=prod(c for _,c in tam); print("prod c_p =",tp)
om=E.period_lattice().omega(); print("omega (real period, full) =",om)
print("real components:",E.real_components())
import time
t0=time.time()
L=E.lseries().dokchitser(prec=40)
print("L(1)  =", L(1))
print("L'(1) =", L.derivative(1,1))
print("L''(1)=", L.derivative(1,2))
print("time", time.time()-t0)
Lp=L.derivative(1,1)
Reg = Lp*T^2/(om*tp)
print("predicted Reg (=h(P) if rank1, Sha=1) =", Reg)
for sh in [1,4,9,16,25]:
    print("  if |Sha|=%d: h(P)="%sh, Reg/sh)
