# Fable, 14.09.2026. Теорема B (неизотривиальность): ни один эллиптический фактор из 12 классов и ни одна
# кривая рода 2 из пяти/шести клеток не повторяется на бесконечно многих наклонах.
# Плюс: тест простоты якобианов пятиклеточных кривых по многочленам Фробениуса (наблюдение) и PARI-ранг пар Бремнера–Саллоуса.
from sage.all import *
import json, sys
R = PolynomialRing(QQ,'k'); k = R.gen(); F = R.fraction_field()

def j_Emn(m,n):
    # E_{m,n}: y^2 = x(x+m^2)(x+n^2) — лежандрова форма с lambda = n^2/m^2
    lam = F(n)**2/F(m)**2
    return 256*(lam**2-lam+1)**3/(lam**2*(lam-1)**2)

r, s = k, F(1)            # наклон k = r/s: положим s = 1, r = k
pairs = {'(r,s)':(r,s),'(r,s-r)':(r,s-r),'(s,s-r)':(s,s-r),'(s,r+s)':(s,r+s),'(r,r+s)':(r,r+s),
         '(s-r,r+s)':(s-r,r+s),'(s,2r-s)':(s,2*r-s),'(s,2r+s)':(s,2*r+s),'(r,2s-r)':(r,2*s-r),
         '(r,2s+r)':(r,2*s+r),'(s,2r)':(s,2*r),'(r,2s)':(r,2*s)}
print('--- j-инварианты 12 классов как функции наклона k ---')
for name,(m,n) in pairs.items():
    j = j_Emn(m,n)
    num = j.numerator(); den = j.denominator()
    print(f'{name:10s} j непостоянен: {j.derivative()!=0}, deg num/den = {num.degree()}/{den.degree()}')
    assert j.derivative()!=0

print('--- инварианты Игузы пятиклеточной кривой y^2=(1+cz)(1-a^2z^2)(1-b^2z^2) ---')
P = PolynomialRing(QQ,'c,a,b'); c,a,b = P.gens(); Pz = PolynomialRing(P.fraction_field(),'z'); z = Pz.gen()
f5 = (1+c*z)*(1-a**2*z**2)*(1-b**2*z**2)
C5 = HyperellipticCurve(f5)
I2,I4,I6,I10 = C5.igusa_clebsch_invariants()
abs1 = I2**5/I10
print('I2^5/I10 =', factor(abs1.numerator()), '/', factor(abs1.denominator()))
# ограничение на наклон: (c;a,b) = (r; s, s+r) и др. — подставляем s=1, r=k
Pk = PolynomialRing(QQ,'kk'); kk = Pk.gen()
for (cc,aa,bb) in [(kk,1,1+kk),(1,kk,1+kk),(1+kk,kk,1-kk),(1-kk,kk,1+kk),(kk,1-kk,1+kk),(1,1-kk,1+kk)]:
    val = abs1.subs({c:cc,a:aa,b:bb}) if False else None
    g = Pk.fraction_field()(abs1.numerator().subs(c=cc,a=aa,b=bb))/Pk.fraction_field()(abs1.denominator().subs(c=cc,a=aa,b=bb))
    print(f'  (c;a,b)=({cc};{aa},{bb}): I2^5/I10 непостоянен по k: {g.derivative()!=0}')
    assert g.derivative()!=0
f6 = (1-c**2*z**2)*(1-a**2*z**2)*(1-b**2*z**2)
I2s,I4s,I6s,I10s = HyperellipticCurve(f6).igusa_clebsch_invariants()
g6 = I2s**5/I10s
for (cc,aa,bb) in [(kk,1,1+kk),(1,kk,1+kk),(1+kk,kk,1-kk)]:
    g = Pk.fraction_field()(g6.numerator().subs(c=cc,a=aa,b=bb))/Pk.fraction_field()(g6.denominator().subs(c=cc,a=aa,b=bb))
    print(f'  шесть клеток ({cc};{aa},{bb}): I2^5/I10 непостоянен по k: {g.derivative()!=0}')
    assert g.derivative()!=0

print('--- простота якобианов пятиклеточных кривых (многочлены Фробениуса, наблюдение) ---')
Rz = PolynomialRing(QQ,'x'); x = Rz.gen()
def frob_split_stats(f, pmax=300):
    C = HyperellipticCurve(f); D = ZZ(f.discriminant())
    nsplit=0; ntot=0; nirr=0
    for p in prime_range(7,pmax):
        if D%p==0: continue
        L = C.change_ring(GF(p)).frobenius_polynomial()
        fac = L.factor(); ntot+=1
        if len(fac)==1 and fac[0][1]==1: nirr+=1
        else: nsplit+=1
    return ntot,nsplit,nirr
tests = {'79/110 (79;31,189)':(79,31,189),'79/110 (79;110,189)':(79,110,189),'11/142 (131;142,153)':(131,142,153),
         '48/163 (115;48,163)':(115,48,163),'67/168 (235;67,101)':(235,67,101),'случайная (5;7,11)':(5,7,11)}
for name,(cc,aa,bb) in tests.items():
    f = (1+cc*x)*(1-aa**2*x**2)*(1-bb**2*x**2)
    ntot,nsplit,nirr = frob_split_stats(f)
    print(f'  пять клеток {name}: простых {ntot}, L_p неприводим над Q у {nirr}, приводим у {nsplit}')
for name,(cc,aa,bb) in list(tests.items())[:2]:
    f = (1-cc**2*x**2)*(1-aa**2*x**2)*(1-bb**2*x**2)
    ntot,nsplit,nirr = frob_split_stats(f)
    print(f'  контроль, шесть клеток {name}: простых {ntot}, неприводим у {nirr}, приводим у {nsplit} (ожидается: всегда приводим — якобиан расщеплён)')

print('--- PARI ellrank для пар Бремнера–Саллоуса (ожидается ранг >= 1) ---')
for (m,n) in [(247,578),(247,825),(578,825),(825,1072)]:
    E = EllipticCurve([0,m**2+n**2,0,m**2*n**2,0])
    try:
        rk = E.rank(only_use_mwrank=False)
        print(f'  ({m},{n}): rank = {rk}')
    except Exception as e:
        lo,hi = E.rank_bounds() if hasattr(E,'rank_bounds') else (None,None)
        print(f'  ({m},{n}): ранг не определён точно: {e}; границы {lo},{hi}')
