# Модель для кода Bianchi: x = 1/z, w = y/z^3:  w^2 = (x^2-a^2)(x^2-b^2)(x^2-c^2) = x^6 + a4 x^4 + a2 x^2 + a0 (моническая).
# E1 = y^2 = (x-a^2)(x-b^2)(x-c^2) = E^-;  E2 = y^2 = x^3 + a2 x^2 + a0 a4 x + a0^2 ~ E^+ (проверяется изоморфизмом ниже).
from sage.all import *
def fpoly(a,b,c):
    R=PolynomialRing(QQ,'x'); x=R.gen()
    return (x**2-a**2)*(x**2-b**2)*(x**2-c**2)
def E12(f):
    a6,a4,a2,a0=f[6],f[4],f[2],f[0]
    return EllipticCurve([0,a4,0,a2*a6,a0*a6**2]), EllipticCurve([0,a2,0,a0*a4,a0**2*a6])
def eplus(a,b,c):
    R=PolynomialRing(QQ,'u'); u=R.gen()
    pol=(1-a**2*u)*(1-b**2*u)*(1-c**2*u); lc=pol.leading_coefficient()
    g=(pol(u/lc)*lc**2).monic(); co=g.list()
    return EllipticCurve([0,co[2],0,co[1],co[0]])
if __name__=='__main__':
    import sys
    a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f)
    print('f =',f)
    print('E1 = E^- ?', E1.is_isomorphic(EllipticCurve([0,-(a**2+b**2+c**2),0,a**2*b**2+a**2*c**2+b**2*c**2,-a**2*b**2*c**2])))
    print('E2 ~ E^+ ?', E2.is_isomorphic(eplus(a,b,c)))
    for E,name in ((E1,'E1'),(E2,'E2')):
        Em=E.minimal_model(); r=pari(Em).ellrank(); print(name,'ellrank',r[:2],'точки',r[2], 'кручение',Em.torsion_order())
    good=[p for p in primes(5,80) if E1.has_good_reduction(p) and E2.has_good_reduction(p) and E1.is_ordinary(p) and E2.is_ordinary(p)]
    print('хорошие ординарные p:',good)
    print('плохие простые C:',factor(f.discriminant()))
    for E,name in ((E1,'E1'),(E2,'E2')):
        Em=E.minimal_model(); print(name,'min disc',factor(Em.discriminant()), 'тамагава',[(q,Em.tamagawa_number(q),str(Em.kodaira_symbol(q))) for q in ZZ(Em.discriminant()).prime_factors()])
