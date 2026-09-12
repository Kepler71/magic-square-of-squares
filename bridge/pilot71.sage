# Пилот Astra: G1, клетки 0,4,8, (m,n) = (7,1). Проверка тождеств и ранги факторов.
R.<t> = QQ[]
m, n = 7, 1
F0 = m^2 + n^2*t^2; F4 = (m^2+n^2)*(t^2+1)/2; F8 = n^2 + m^2*t^2
print("F0 =", F0, " F4 =", F4, " F8 =", F8)
# 1) тождества t²F0(1/t) = F8, t²F4(1/t) = F4
Ft = R.fraction_field()
print("t²F0(1/t) == F8:", Ft(t^2*F0.subs(t=1/t)) == Ft(F8))
print("t²F4(1/t) == F4:", Ft(t^2*F4.subs(t=1/t)) == Ft(F4))
# 2) H: Y² = F0F4F8 и нормализация Y = 35 y
P = F0*F4*F8; a = QQ(m)^2 + 1 + QQ(n)^2/m^2 * (n^2/m^2)^0 * 1
a = QQ(m^2)/n^2 + 1 + QQ(n^2)/m^2; b = QQ((m^2+n^2)*m^2*n^2)/2
print("a =", a, " b =", b, " b — квадрат:", b.is_square())
print("P == b*(t^6 + a t^4 + a t^2 + 1):", P == b*(t^6 + a*t^4 + a*t^2 + 1))
# 3) эллиптический фактор: X = t², E: v² = X³ + aX² + aX + 1
S.<X> = QQ[]
E = EllipticCurve(QQ, [0, a, 0, a, 1]).minimal_model()
Eastra = EllipticCurve(QQ, [0, 1+49+2401, 0, 1*49+1*2401+49*2401, 1*49*2401]).minimal_model()  # (x+1)(x+49)(x+2401)
print("E (из шестерички) =", E.ainvs(), " E (модель Astra) =", Eastra.ainvs(), " изоморфны:", E.is_isomorphic(Eastra))
def rk(C):
    r = pari(C.minimal_model().ainvs()).ellinit().ellrank()
    return ZZ(r[0]), ZZ(r[1])
# 4) A = Jac(v² = F0F4), B = Jac(v² = F0F8) — кривые рода 1 с точкой t=0
from sage.schemes.elliptic_curves.jacobian import Jacobian
A = Jacobian(R(F0*F4)).minimal_model()
B = Jacobian(R(F0*F8)).minimal_model()
C48 = Jacobian(R(F4*F8)).minimal_model()
print("\nA = Jac(F0F4):", A.ainvs(), " ранг", rk(A))
print("B = Jac(F0F8):", B.ainvs(), " ранг", rk(B))
print("C48 = Jac(F4F8):", C48.ainvs(), " ранг", rk(C48), " A ≅ C48:", A.is_isomorphic(C48))
print("E:", E.ainvs(), " ранг", rk(E), " кручение", E.torsion_order())
ra, rb, re = rk(A)[0], rk(B)[0], rk(E)[0]
print(f"\nr(J_C) = 2r(A) + r(B) + 2r(E) = 2·{ra} + {rb} + 2·{re} = {2*ra + rb + 2*re}")
print("порог квадратичного Чабо для G1 (rho >= 7, g = 5): r <= 10")
