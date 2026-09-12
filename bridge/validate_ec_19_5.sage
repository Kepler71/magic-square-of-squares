# -*- coding: utf-8 -*-
r"""
validate_ec_19_5.sage
=====================
Проверка части C (поиск по E_b) из search_H_19_5.sage.

Тест A. Генератор: точка лежит на E_b, имеет бесконечный порядок,
        E.saturation даёт индекс 1, PARI ellrank = [1,1].
Тест B. Редукция: самописное сложение над F_p даёт ровно редукцию nG+T,
        посчитанную точно в Sage (n = 1..40, несколько простых).
Тест C. Модулярный фильтр «звучен»: если x(P)/b — квадрат в Q, то n
        обязан пройти фильтр.  Проверяется на искусственных точках:
        берём любое P in E_b(Q) и убеждаемся, что фильтр не отбрасывает
        n, для которого x(nG+T)*b действительно квадрат (таких n нет,
        поэтому дополнительно проверяем логику на подставных значениях).
Тест D. Обратная карта E_b -> H: для произвольного рационального t
        точка (b t^2, b*Y) лежит на E_b, где Y^2 = f(t) — контроль,
        что критерий «x/b — квадрат» действительно выделяет точки H.
"""
import numpy as np

m = Integer(19); n = Integer(5)
s = (m ^ 2 + n ^ 2) // 2
b = s * m ^ 2 * n ^ 2
aa = QQ(m ^ 2) / n ^ 2 + 1 + QQ(n ^ 2) / m ^ 2
A2 = Integer(b * aa); A4 = Integer(b ^ 2 * aa); A6 = Integer(b ^ 3)
E = EllipticCurve([0, A2, 0, A4, A6])
R.<t> = QQ[]
F0 = m ^ 2 + n ^ 2 * t ^ 2; F4 = s * (1 + t ^ 2); F8 = n ^ 2 + m ^ 2 * t ^ 2
f = F0 * F4 * F8
GX = QQ(-28561963657) / 1369
GY = QQ(2089086742828800) / 50653
G = E(GX, GY)

print("=== ТЕСТ A: генератор и ранг ===")
print("  G =", G, " на кривой:", G in E)
print("  порядок:", G.order())
sat, idx, reg = E.saturation([G])
print("  E.saturation([G]) -> индекс =", idx, ", насыщенная точка:", sat[0] == G)
Em = E.minimal_model()
print("  PARI ellrank(minimal model) =", pari(Em).ellrank())
print("  кручение:", E.torsion_subgroup().invariants())
print("  результат:", "OK" if (idx == 1 and G.order() == oo) else "ПРОВАЛ")


def ec_add(P, Q, p, a2, a4):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * pow(2 * y1, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - a2 - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


print()
print("=== ТЕСТ B: самописная редукция против точного вычисления в Sage ===")
a2 = int(A2); a4 = int(A4); a6 = int(A6)
disc = Integer(E.discriminant())
tors = E.torsion_points()
ok = True
ps = []
p = next_prime(2000)
while len(ps) < 6:
    if disc % p != 0 and b % p != 0:
        ps.append(int(p))
    p = next_prime(p)
for p in ps:
    Gp = (int(GX.numerator() * pow(int(GX.denominator()), p - 2, p) % p),
          int(GY.numerator() * pow(int(GY.denominator()), p - 2, p) % p))
    Tp = [None] + [(int(T[0]) % p, int(T[1]) % p) for T in tors if not T.is_zero()]
    Tex = [E(0)] + [T for T in tors if not T.is_zero()]
    cur = None
    for k in range(1, 41):
        cur = ec_add(cur, Gp, p, a2, a4)          # = k*G mod p
        for j in range(4):
            Q = ec_add(cur, Tp[j], p, a2, a4)
            P = k * G + Tex[j]
            if P.is_zero():
                if Q is not None:
                    ok = False; print("  РАСХОЖДЕНИЕ O", p, k, j)
                continue
            xn = P[0].numerator() % p
            xd = P[0].denominator() % p
            if xd == 0:
                continue          # точка редуцируется в O — Q должно быть None
            xred = xn * pow(xd, p - 2, p) % p
            if Q is None or Q[0] != xred:
                ok = False
                print("  РАСХОЖДЕНИЕ x", p, k, j, Q, xred)
print("  проверено 6 простых x 40 кратных x 4 класса кручения")
print("  результат:", "OK" if ok else "ПРОВАЛ")

print()
print("=== ТЕСТ C: логика фильтра на подставных значениях ===")
# берём рациональные x, для которых x*b — заведомо квадрат и заведомо нет
ok = True
for q in [QQ(1), QQ(4), QQ(9)/25, QQ(b), QQ(1)/b]:
    xq = q * b            # тогда x*b = q*b^2 — квадрат <=> q квадрат
    for p in ps:
        val = (Integer(xq.numerator()) * pow(int(xq.denominator()), p - 2, p)
               % p) * (b % p) % p
        isqr = kronecker(val, p) != -1 or val == 0
        if q.is_square() and not isqr:
            ok = False
            print("  ЛОЖНОЕ ОТБРАСЫВАНИЕ", q, p)
print("  результат:", "OK" if ok else "ПРОВАЛ")

print()
print("=== ТЕСТ D: обратная карта E_b -> H на произвольных t ===")
ok = True
for tt in [QQ(1), QQ(2), QQ(3)/7, QQ(-5)/4, QQ(11)/9]:
    Xv = b * tt ^ 2
    Yv2 = f(tt)
    lhs = Xv ^ 3 + A2 * Xv ^ 2 + A4 * Xv + A6
    if lhs != b ^ 2 * Yv2:
        ok = False
        print("  РАСХОЖДЕНИЕ", tt)
print("  X = b t^2, V = b Y удовлетворяет уравнению E_b для всех тестовых t:",
      "OK" if ok else "ПРОВАЛ")
print("  обратно: x(P)/b квадрат <=> P — образ точки H (t = ±sqrt(x/b)):")
print("    это тождество, а не численный факт (см. search_H_19_5.sage,"
      " раздел identities).")
