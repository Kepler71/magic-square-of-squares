load('/home/kep/magicKube/corners/corner_search.sage')   # QuarticCurve (Коннелл), наш код
R = PolynomialRing(QQ, 'T'); T = R.gen()
m, n = 7, 1
F0 = m^2 + n^2*T^2; F4 = (m^2+n^2)*(T^2+1)/2; F8 = n^2 + m^2*T^2
a = QQ(m^2)/n^2 + 1 + QQ(n^2)/m^2
def rk(C):
    r = pari(C.minimal_model().ainvs()).ellinit().ellrank(); return ZZ(r[0]), ZZ(r[1])
E = EllipticCurve(QQ, [0, a, 0, a, 1]).minimal_model()
qs = {'A = Jac(F0·F4)': F0*F4, 'B = Jac(F0·F8)': F0*F8, 'C48 = Jac(F4·F8)': F4*F8}
curves = {}
for name, q in qs.items():
    qc = QuarticCurve(R(q), QQ(0))          # рациональная точка t = 0
    Em = qc.E.minimal_model(); curves[name] = Em
    lo, hi = rk(Em)
    print(f"{name}: {Em.ainvs()}  ранг [{lo},{hi}]  кручение {Em.torsion_order()}  проводник {Em.conductor().factor()}")
lo, hi = rk(E)
print(f"E (фактор J_H):   {E.ainvs()}  ранг [{lo},{hi}]  кручение {E.torsion_order()}")
A = curves['A = Jac(F0·F4)']; B = curves['B = Jac(F0·F8)']; C48 = curves['C48 = Jac(F4·F8)']
print("\nA ≅ C48 (проверка изоморфизма из записки):", A.is_isomorphic(C48))
ra, rb, re = rk(A)[0], rk(B)[0], rk(E)[0]
print(f"r(J_C) = 2·{ra} + {rb} + 2·{re} = {2*ra+rb+2*re}   (порог квадратичного Чабо для G1: r ≤ 10)")
