import time
t0=time.time()
m=19;n=16
s=QQ(m^2+n^2)/2; b=s*m^2*n^2
R.<X>=QQ[]
f=(X+b)*(X+s*m^4)*(X+s*n^4)
E=EllipticCurve([0,f[2],0,f[1],f[0]])
Em=E.minimal_model()
print("Em",Em)
N=Em.conductor(); print("N =",N, factor(N))
print("tors order", Em.torsion_order())
print("tamagawa product", Em.tamagawa_product(), [ (p,Em.tamagawa_number(p)) for p in N.prime_factors()])
Om = Em.period_lattice().omega()   # real period (includes factor for #components)
print("omega (real period incl. components) =", Om)
print("t",time.time()-t0)
# L'(1) via pari lfun
pe = pari(Em)
L = pari('L = lfuninit(ellinit(%s),[0,3,1]); 0'%(list(Em.a_invariants()),))
v = pari('lfun(L,1,1)')
print("L'(1) =", v)
Lp = RR(v)
print("t",time.time()-t0)
tors = Em.torsion_order()
cprod = Em.tamagawa_product()
Reg_pred = Lp * tors^2 / (RR(Om) * cprod)
print("Predicted Reg (assuming Sha=1) =", Reg_pred)
for sh in [1,4,9,16,25]:
    print("  if Sha =",sh," =>  hhat(G) =", Reg_pred/sh)
print("t",time.time()-t0)
