# -*- coding: utf-8 -*-
# Шаг 0: построение E, минимальной модели, I,J, и проверка всей цепочки соответствий
# для пары G1 (m,n) = (11,4).  Ранг НЕ используется.
import sys, json
load('/home/kep/magicKube/bridge/ctp_cert_11_4_claude/ctp_quartic_snapshot.sage')

m, n = 11, 4
s = QQ(m^2 + n^2)/2
b = s * m^2 * n^2
print("m,n =", m, n, " s =", s, " b =", b)
assert s == QQ(137)/2 and b == 132616

e = [-b, -s*m^4, -s*n^4]
print("корни исходной E:", e)

# исходная E: V^2 = (X+b)(X+s m^4)(X+s n^4) = prod (X - e_i)
R.<X> = QQ[]
fE = prod(X - ei for ei in e)
E0 = EllipticCurve([0, fE.coefficients(sparse=False)[2], 0, fE.coefficients(sparse=False)[1], fE.coefficients(sparse=False)[0]])
print("E0 =", E0)

# целочисленная модель x = 4X + c: корни 4 e_i + c
# выбираем c так, чтобы получилась модель y^2 = x^3 + A x + B (сумма корней = 0)
c = -sum(4*ei for ei in e)/3
print("сдвиг c =", c, " (должен быть целым)")
assert c in ZZ
rootsM = [4*ei + c for ei in e]
print("корни модели M (в порядке e1,e2,e3):", rootsM)
fM = prod(X - r for r in rootsM)
print("fM =", fM)
cofs = fM.coefficients(sparse=False)
assert cofs[2] == 0
M = EllipticCurve([0, 0, 0, cofs[1], cofs[0]])
print("M =", M)
print("M минимальна?", M == M.minimal_model(), " минимальная модель:", M.minimal_model())

# масштаб 4 — квадрат, значит класс delta сохраняется
print("масштаб x - e_i^M = 4 (X - e_i);  4 — квадрат:", QQ(4).is_square())

I, J = IJ_of_curve(M)
print("I =", I)
print("J =", J)
assert I == 222926259931200, I
assert J == -6316841727504282240000, J

phis = [phi_of_root(M, r) for r in rootsM]
print("phi_i = -12 e_i^M :", phis)
Rx.<T> = QQ[]
cub = T^3 - 3*I*T + J
for p in phis:
    assert cub(p) == 0
print("все phi_i — корни T^3-3IT+J: да")

F = FisherCTP(QQ, I, J, verbose=True)
print("компоненты L:", [(cc['deg'], cc['r']) for cc in F.E.comps])
order = _comp_order(F, phis)
print("order (индекс корня phi для каждой компоненты L):", order)

# delta = (1, s, s) в порядке (e1,e2,e3);  s = 137/2 ~ 274 mod квадратов
delta_trip = [QQ(1), QQ(274), QQ(274)]
assert (QQ(274)/s).is_square(), "274 и s=137/2 в одном классе"
delta = F.E.from_comps([delta_trip[i] for i in order])
print("delta в L:", delta)
print("компоненты delta:", [F.E.comp(delta, i) for i in range(3)])

save((m, n, s, b, e, c, rootsM, I, J, phis, order, delta_trip), '/home/kep/magicKube/bridge/ctp_cert_11_4_claude/setup.sobj')
print("OK setup")
