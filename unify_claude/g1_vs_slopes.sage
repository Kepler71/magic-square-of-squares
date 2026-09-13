# Гипотеза: фактор B семейства G1(m,n), Y^2 = X(X-(m^2-n^2)^2)(X-(m^2+n^2)^2),
# 2-изогенен кривой E_{m^2,n^2}: y^2 = x(x+m^4)(x+n^4) из разреза по наклонам.
import random
def EB(m,n):  return EllipticCurve(QQ,[0,-((m^2-n^2)^2+(m^2+n^2)^2),0,(m^2-n^2)^2*(m^2+n^2)^2,0])
def Emn(a,b): return EllipticCurve(QQ,[0,a^2+b^2,0,a^2*b^2,0])      # x(x+a^2)(x+b^2)
pairs=[(m,n) for m in range(2,21) for n in range(1,m) if gcd(m,n)==1]
bad=[]; rk_mismatch=[]
for m,n in pairs:
    B=EB(m,n); E=Emn(m^2,n^2)
    if not B.is_isogenous(E): bad.append((m,n)); continue
    phi=E.isogeny(E(0,0)); iso=phi.codomain().is_isomorphic(B)
    if not iso: bad.append((m,n,'изогенны, но не через (0,0)'))
print("пар G1:", len(pairs), "; не изогенны / не та изогения:", bad[:5], len(bad))
# ранги на выборке — изогения сохраняет ранг, но сверим с данными
random.seed(int(5)); samp=random.sample(pairs,int(25))
for m,n in samp:
    rB=pari(EB(m,n).minimal_model().a_invariants()).ellinit().ellrank(); rE=pari(Emn(m^2,n^2).minimal_model().a_invariants()).ellinit().ellrank()
    if (int(rB[0]),int(rB[1]))!=(int(rE[0]),int(rE[1])): rk_mismatch.append((m,n,rB[:2],rE[:2]))
print("выборка 25 пар: расхождений интервалов ранга:", rk_mismatch)
