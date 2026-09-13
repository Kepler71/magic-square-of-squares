# Типы редукции и числа Тамагавы E и E' при p=2,3 по сетке; символьные тождества.
import json, glob
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
R.<a,b> = QQ[]
alpha = a^2+a*b+b^2; beta = a^2*b^2*(a+b)^2
print("4alpha^3-27beta == ((a-b)(a+2b)(2a+b))^2 :", 4*alpha^3-27*beta == ((a-b)*(a+2*b)*(2*a+b))^2)
Eg = EllipticCurve(R.fraction_field(), [0, alpha^2, 0, 2*alpha*beta, beta^2])
print("disc factor:", factor(Eg.discriminant()))
from collections import Counter, defaultdict
tab = defaultdict(Counter)
def v(n,p):
    return ZZ(n).valuation(p) if n != 0 else 99
pairs = [(aa,bb) for bb in range(2, 61) for aa in range(1, bb) if gcd(aa,bb) == 1]
for (aa,bb) in pairs:
    E, Ep, phi, al, be = curve_data(aa,bb)
    for p in (2,3):
        F = dict(a=aa,b=bb,apb=aa+bb,amb=aa-bb,ap2b=aa+2*bb,a2pb=2*aa+bb)
        divs = tuple(sorted((k,v(F[k],p)) for k in F if v(F[k],p)>0))
        ld = E.local_data(p); ldp = Ep.local_data(p)
        key = (p, divs)
        tab[key][(str(ld.kodaira_symbol()), ld.tamagawa_number(), str(ldp.kodaira_symbol()), ldp.tamagawa_number(), ld.conductor_valuation())] += 1
for key in sorted(tab, key=str):
    print(key, dict(tab[key]))
