import time
t0=time.time()
m=19;n=16
s=QQ(m^2+n^2)/2; b=s*m^2*n^2
R.<X>=QQ[]
f=(X+b)*(X+s*m^4)*(X+s*n^4)
E=EllipticCurve([0,f[2],0,f[1],f[0]])
Em=E.minimal_model()
print("Em",Em)
print("cond Em", factor(Em.conductor()))
# the untwisted small curve E' : w^2=(y+m^2n^2)(y+m^4)(y+n^4)
g=(X+m^2*n^2)*(X+m^4)*(X+n^4)
Ep=EllipticCurve([0,g[2],0,g[1],g[0]])
Epm=Ep.minimal_model()
print("E' =",Ep)
print("E' min:",Epm, "cond:",factor(Epm.conductor()))
print("check twist: E' twisted by 1234 isomorphic to E ?", Ep.quadratic_twist(1234).is_isomorphic(Em))
print("rank E' :", Ep.rank())
print("gens E' :", Ep.gens())
print("t",time.time()-t0)
print("--- 2-isogenous curves of Em ---")
for phi in Em.isogenies_prime_degree(2):
    C=phi.codomain().minimal_model()
    print("codomain min:",C, " cond:",factor(C.conductor()), " disc:", factor(C.discriminant()))
print("t",time.time()-t0)
