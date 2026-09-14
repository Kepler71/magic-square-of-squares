# неопределённые [0,2] у кривых Фибоначчи C_k: L(E,1)≠0 ⇒ ранг 0; L(E,1)=0 при w=+1 ⇒ аналитический ранг ≥ 2
from sage.all import *
F=[1,1]
while len(F)<70: F.append(F[-1]+F[-2])
X=PolynomialRing(QQ,'z').gen()
for k in [8,30,41,43,45,47]:
    f=(1-F[k]*X)*(1-F[k+1]*X)*(1-F[k+2]*X)
    L=f.leading_coefficient(); h=(f(X/L)*L**2).monic(); co=h.list()
    E=EllipticCurve([0,co[2],0,co[1],co[0]]).minimal_model()
    try: Lr=E.lseries().L_ratio()
    except Exception as e: Lr='ошибка: '+str(e)[:60]
    print(k,'ellrank',E.pari_curve().ellrank()[:2],'L/Ω =',Lr,flush=True)
