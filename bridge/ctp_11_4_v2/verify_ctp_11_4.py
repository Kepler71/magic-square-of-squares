# -*- coding: utf-8 -*-
"""
НЕЗАВИСИМАЯ ПРОВЕРКА СЕРТИФИКАТА СПАРИВАНИЯ КАССЕЛСА-ТЕЙТА.  Пара G1 (m,n) = (11,4).

Цель.  Показать, что квадратный класс delta = (1, 274, 274) = (1, s, s) в (Q*/Q*^2)^3,
необходимый для любой конечной рациональной точки кривой C_(11,4) рода 5,
НЕ лежит в образе E(Q)/2E(Q) -> Sel^2(E/Q).
Инструмент: явное ненулевое спаривание Касселса-Тейта <[g1],[g2]>_CT = -1.
Никакая верхняя граница ранга E(Q) НЕ используется.

Формула: T. Fisher, "On binary quartics and the Cassels-Tate pairing",
arXiv:2208.14977, Research in Number Theory 8 (2022), Art. 74, Theorem 3.1
(+ Remark 3.3 о конечности набора мест, + Remark 3.2(v) об антисимметрии).

Скрипт использует ТОЛЬКО стандартную библиотеку Python (fractions, math, json,
random для самотестов).  Sage, PARI и код проекта не импортируются.
Вся арифметика точная (Fraction / int).  p-адических приближений нет.

Запуск:   python3 verify_ctp_11_4.py [путь/к/certificate.json]
Выход:    independent_result.json рядом с сертификатом; код возврата 0 = всё сошлось.
"""

import json
import random
import sys
from fractions import Fraction as Fr
from math import gcd, isqrt
from pathlib import Path

# ----------------------------------------------------------------- журнал проверок

CHECKS = []


def ok(cond, name, extra=""):
    CHECKS.append((bool(cond), name, extra))
    if not cond:
        print("[ПРОВАЛ] " + name + ("  " + extra if extra else ""))
        raise SystemExit("проверка не прошла: " + name)
    print("[ok] " + name + ("  " + extra if extra else ""))


# ----------------------------------------------------------------- базовая арифметика

def QQ(x):
    """Точное рациональное из строки/int/Fraction."""
    if isinstance(x, Fr):
        return x
    return Fr(str(x))


def is_square_Q(a):
    a = QQ(a)
    if a < 0:
        return False
    if a == 0:
        return True
    return isqrt(a.numerator) ** 2 == a.numerator and isqrt(a.denominator) ** 2 == a.denominator


def factor_int(n):
    """Полная факторизация |n| пробным делением. Числа сертификата малы."""
    n = abs(int(n))
    out = {}
    if n == 0:
        return out
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def primes_of(a):
    """Множество простых в числителе и знаменателе рационального a (a != 0)."""
    a = QQ(a)
    if a == 0:
        return set()
    return set(factor_int(a.numerator)) | set(factor_int(a.denominator))


def val_unit(a, p):
    """a = p^k * u, u — p-единица (рациональная).  Возвращает (k, u)."""
    a = QQ(a)
    assert a != 0
    num, den = a.numerator, a.denominator
    k = 0
    while num % p == 0:
        num //= p
        k += 1
    while den % p == 0:
        den //= p
        k -= 1
    return k, Fr(num, den)


def res_mod(u, mod):
    """Вычет рациональной p-единицы u по модулю mod (mod взаимно прост со знаменателем)."""
    u = QQ(u)
    return (u.numerator * pow(u.denominator, -1, mod)) % mod


def legendre(u, p):
    """Символ Лежандра (u|p) для нечётного p и p-единицы u:  +1 или -1."""
    r = pow(res_mod(u, p), (p - 1) // 2, p)
    return 1 if r == 1 else -1


def is_local_square(a, place):
    """a — ненулевой квадрат в Q_v ?  place: 'real' или строка/int простого p."""
    a = QQ(a)
    if a == 0:
        return False
    if place == "real":
        return a > 0
    p = int(place)
    k, u = val_unit(a, p)
    if k % 2:
        return False
    if p == 2:
        return res_mod(u, 8) == 1
    return legendre(u, p) == 1


def hilbert(a, b, place):
    """Символ Гильберта (a,b)_v в {+1,-1}.  Собственная реализация; ниже самотест."""
    a, b = QQ(a), QQ(b)
    assert a != 0 and b != 0
    if place == "real":
        return -1 if (a < 0 and b < 0) else 1
    p = int(place)
    al, u = val_unit(a, p)
    be, w = val_unit(b, p)
    if p == 2:
        uu, ww = res_mod(u, 8), res_mod(w, 8)
        eps_u, eps_w = (uu - 1) // 2, (ww - 1) // 2
        om_u, om_w = (uu * uu - 1) // 8, (ww * ww - 1) // 8
        e = eps_u * eps_w + al * om_w + be * om_u
        return -1 if e % 2 else 1
    s = 1
    if (al * be) % 2 and ((p - 1) // 2) % 2:
        s = -s
    if be % 2 and legendre(u, p) == -1:
        s = -s
    if al % 2 and legendre(w, p) == -1:
        s = -s
    return s


def selftest_hilbert(trials=400, seed=20260912):
    """Самотест символа Гильберта.

    (1) Закон взаимности Гильберта: произведение (a,b)_v по ВСЕМ местам = +1.
        Это ловит любую ошибку знака в отдельной формуле.
    (2) Прямое определение: (a,b)_v = +1  <=>  a x^2 + b y^2 = z^2 имеет
        нетривиальное решение в Q_v.  Проверяется перебором для вещественного места
        и малых p (полный перебор по вычетам + подъём по Гензелю через is_local_square).
    (3) Известные значения: (-1,-1)_2 = -1, (2,5)_5 = -1, (-1,-1)_real = -1.
    """
    rnd = random.Random(seed)
    bad = 0
    for _ in range(trials):
        a = Fr(rnd.choice([1, -1]) * rnd.randint(1, 400), rnd.randint(1, 60))
        b = Fr(rnd.choice([1, -1]) * rnd.randint(1, 400), rnd.randint(1, 60))
        if a == 0 or b == 0:
            continue
        pls = sorted(primes_of(a) | primes_of(b) | {2})
        prod = hilbert(a, b, "real")
        for p in pls:
            prod *= hilbert(a, b, p)
        if prod != 1:
            bad += 1
    ok(bad == 0, "самотест символа Гильберта: взаимность на %d случайных парах" % trials)
    ok(hilbert(-1, -1, 2) == -1 and hilbert(2, 5, 5) == -1 and hilbert(-1, -1, "real") == -1,
       "самотест символа Гильберта: контрольные значения")
    # прямое определение через разрешимость a x^2 + b y^2 = z^2
    def soluble(a, b, p):
        if p == "real":
            return not (a < 0 and b < 0)
        p = int(p)
        M = p ** (2 * (2 if p == 2 else 1) + 3)
        for x in range(M):
            for y in range(M):
                if x == 0 and y == 0:
                    continue
                v = a * x * x + b * y * y
                if v != 0 and is_local_square(v, p):
                    return True
        return False
    small = 0
    for a in [-3, -2, -1, 1, 2, 3, 5, 6, 7, 10, -5, 14, 15, -7]:
        for b in [-3, -2, -1, 1, 2, 3, 5, 6, 7, 10, -5, 14, 15, -7]:
            for p in ["real", 2, 3, 5, 7]:
                if hilbert(a, b, p) != (1 if soluble(Fr(a), Fr(b), p) else -1):
                    small += 1
    ok(small == 0, "самотест символа Гильберта: сверка с разрешимостью a x^2 + b y^2 = z^2")


# ----------------------------------------------------------------- бинарные квартики

def quartic_I(g):
    a, b, c, d, e = g
    return 12 * a * e - 3 * b * d + c * c


def quartic_J(g):
    a, b, c, d, e = g
    return 72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3


def quartic_hessian(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c,
            4 * (b * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e),
            3 * d * d - 8 * c * e]


def quartic_eval(g, x, z):
    a, b, c, d, e = g
    return a * x ** 4 + b * x ** 3 * z + c * x * x * z * z + d * x * z ** 3 + e * z ** 4


def z_invariant(g, phi):
    """z(g) = (4 a phi + 3 b^2 - 8 a c)/3  в компоненте L, отвечающей корню phi."""
    a, b, c = g[0], g[1], g[2]
    return (4 * a * phi + 3 * b * b - 8 * a * c) / 3


def H_form(g, phi, I):
    """Квадратичная форма H (Fisher, формула (4)) в компоненте phi.

    G = (4 phi g + h)/3,  H = (1/12) d^2 G/dx^2 + (2/9)(I - phi^2) z^2,
    т.е. H = [G0, G1/2, G2/6 + (2/9)(I - phi^2)] при базисе (x^2, xz, z^2).
    Тождество G(1,0) * G = H^2 проверяется вызывающей стороной.
    """
    h = quartic_hessian(g)
    G = [(4 * phi * gi + hi) / 3 for gi, hi in zip(g, h)]
    H = [G[0], G[1] / 2, G[2] / 6 + Fr(2, 9) * (I - phi * phi)]
    return G, H


def sq_of_quadratic(H):
    """Коэффициенты H^2 как квартики (x^4, x^3z, x^2z^2, xz^3, z^4)."""
    out = [Fr(0)] * 5
    for i in range(3):
        for j in range(3):
            out[i + j] += H[i] * H[j]
    return out


def eval_quadratic(Q3, x, z):
    return Q3[0] * x * x + Q3[1] * x * z + Q3[2] * z * z


# ----------------------------------------------------------------- сертификат

def main():
    cert_path = Path(sys.argv[1]) if len(sys.argv) > 1 else \
        Path(__file__).resolve().parent / "certificate.json"
    R = json.loads(cert_path.read_text())
    print("сертификат:", cert_path)
    print("пара (m,n) =", R["m"], R["n"])
    print()

    print("--- 0. самотесты арифметики -------------------------------------------")
    selftest_hilbert()
    print()

    # ---------------------------------------------------------------- 1. исходная C
    print("--- 1. кривая C рода 5, необходимый квадратный класс -------------------")
    m, n = int(R["m"]), int(R["n"])
    ok((m, n) == (11, 4), "пара (m,n) = (11,4)")
    ok(gcd(m, n) == 1, "gcd(m,n) = 1")
    s = QQ(R["s"])
    ok(s == Fr(m * m + n * n, 2), "s = (m^2+n^2)/2 = %s" % s)
    b = QQ(R["b"])
    ok(b == s * m * m * n * n, "b = s m^2 n^2 = %s" % b)
    ok(not is_square_Q(s), "s не квадрат в Q  =>  на бесконечности C точек нет "
                           "((u4/t)^2 = s невозможно)")

    # Тождества, задающие необходимый класс, — проверяются СИМВОЛЬНО через
    # многочлены от t^2 (точные рациональные коэффициенты при каждой степени).
    # X = b t^2;  F0 = m^2 + n^2 t^2, F4 = s(1+t^2), F8 = n^2 + m^2 t^2.
    e1, e2, e3 = -b, -s * m ** 4, -s * n ** 4
    ok([QQ(v) for v in R["orig_roots"]] == [e1, e2, e3],
       "корни E: (e1,e2,e3) = (-b, -s m^4, -s n^4)")
    # X - e1 = b t^2 + b            = (mn)^2 * F4
    # X - e2 = b t^2 + s m^4        = s m^2 * F0
    # X - e3 = b t^2 + s n^4        = s n^2 * F8
    for (rt, coef, F_const, F_t2, label) in [
            (e1, QQ((m * n) ** 2), s, s, "X-e1 = (mn)^2 * F4,  F4 = s(1+t^2)"),
            (e2, s * m * m, QQ(m * m), QQ(n * n), "X-e2 = s m^2 * F0,  F0 = m^2+n^2 t^2"),
            (e3, s * n * n, QQ(n * n), QQ(m * m), "X-e3 = s n^2 * F8,  F8 = n^2+m^2 t^2")]:
        # слева: (b) t^2 + (-rt);  справа: (coef*F_t2) t^2 + (coef*F_const)
        ok(b == coef * F_t2 and -rt == coef * F_const, "тождество " + label)
    ok(is_square_Q(QQ(1) / 1) and is_square_Q(s / QQ(274)) and is_square_Q(s / QQ(274)),
       "квадратные классы (1, s, s) = (1, 274, 274)")
    target = [QQ(t) for t in R["target_class"]]
    ok(target == [Fr(1), Fr(274), Fr(274)], "целевой класс сертификата = (1, 274, 274)")
    print("  => любая КОНЕЧНАЯ рациональная точка C даёт точку E(Q) класса (1,274,274);")
    print("     точек C над бесконечностью нет, т.к. s не квадрат.")
    print()

    # ---------------------------------------------------------------- 2. модель
    print("--- 2. минимальная модель и (I,J) -------------------------------------")
    A, B = int(R["M_A"]), int(R["M_B"])
    mr = [QQ(v) for v in R["minimal_roots"]]
    ok(len(set(mr)) == 3, "три различных корня минимальной модели")
    ok(sum(mr) == 0, "сумма корней = 0 (модель y^2 = x^3 + A x + B)")
    ok(all(e ** 3 + A * e + B == 0 for e in mr), "корни удовлетворяют x^3 + A x + B = 0")
    u2 = QQ(R["iso_xscale"])
    shift = QQ(R["iso_shift"])
    ok(all(u2 * X + shift == e for X, e in zip([e1, e2, e3], mr)),
       "изоморфизм x_M = %s * X + %s переводит корни в корни" % (u2, shift))
    ok(is_square_Q(u2), "масштаб %s — квадрат  =>  квадратные классы (X - e_i) СОХРАНЯЮТСЯ" % u2)
    # (u,r,s,t)-изоморфизм: x_M = u^-2 (X - r), y_M = u^-3 V; (u^-2)^3 = (u^-3)^2
    ok(u2 ** 3 == (QQ(R["iso_yscale"])) ** 2,
       "согласованность масштабов: (x-масштаб)^3 = (y-масштаб)^2")
    I, J = QQ(R["I"]), QQ(R["J"])
    ok(I == -48 * A, "I = c4(M) = -48 A")
    ok(J == -1728 * B, "J = 2 c6(M) = -1728 B")
    disc = Fr(16, 27) * (4 * I ** 3 - J ** 2)
    ok(disc != 0, "Delta = 16(4I^3 - J^2)/27 != 0")
    ok(disc == QQ(R["discriminant"]), "Delta совпадает с сертификатом")
    phis = [QQ(v) for v in R["phi_roots"]]
    ok(phis == [-12 * e for e in mr],
       "phi_i = -12 e_i  (x-координаты 2-кручения E_{I,J} равны -3 phi, X = 36 x)")
    ok(all(ph ** 3 - 3 * I * ph + J == 0 for ph in phis), "phi_i — корни X^3 - 3 I X + J")
    ok(len(set(phis)) == 3, "L = Q[phi]/(X^3-3IX+J) = Q x Q x Q, три различные компоненты")
    print()

    # ---------------------------------------------------------------- 3. квартики
    print("--- 3. три квартики, их инварианты и классы ---------------------------")
    gs = [[QQ(c) for c in g] for g in R["quartics"]]
    for k, g in enumerate(gs):
        ok(quartic_I(g) == I and quartic_J(g) == J,
           "квартика g%d имеет инварианты (I,J)" % (k + 1))
    zs = [[z_invariant(g, ph) for ph in phis] for g in gs]
    ok(zs == [[QQ(v) for v in row] for row in R["z_values"]], "z(g_i) совпадают с сертификатом")
    ok(all(all(c != 0 for c in row) for row in zs), "все z(g_i) — единицы в L (не делители нуля)")
    ok(all(is_square_Q(z / d) for z, d in zip(zs[0], target)),
       "КЛЮЧЕВОЕ: класс z(g1) в L*/(L*)^2 равен (1, 274, 274) — это и есть целевой класс")
    ms = [QQ(v) for v in R["m_values"]]
    ok(all(mv * mv == zs[0][i] * zs[1][i] * zs[2][i] for i, mv in enumerate(ms)),
       "z(g1) z(g2) z(g3) = m^2 в L  =>  [g1]+[g2]+[g3] = 0 в Sel^2")
    print()

    # ---------------------------------------------------------------- 4. ELS
    print("--- 4. всюду локальная разрешимость всех трёх квартик ------------------")
    for k, (g, loc) in enumerate(zip(gs, R["local_solubility"])):
        need = {2, 3} | primes_of(disc)
        for c in g:
            need |= primes_of(QQ(c.denominator))
        have = {v["place"] for v in loc}
        ok({str(p) for p in need} | {"real"} <= have,
           "g%d: свидетели есть во всех местах, не покрытых хорошей редукцией" % (k + 1),
           "нужно %s" % sorted(str(p) for p in need))
        for v in loc:
            x, z = QQ(v["x"]), QQ(v["z"])
            ok(not (x == 0 and z == 0), "g%d @ %s: (x:z) != (0:0)" % (k + 1, v["place"]))
            val = quartic_eval(g, x, z)
            ok(val == QQ(v["value"]) and val != 0,
               "g%d @ %s: g(x,z) = %s" % (k + 1, v["place"], val))
            ok(is_local_square(val, v["place"]),
               "g%d @ %s: g(x,z) — квадрат в Q_v" % (k + 1, v["place"]))
        # остальные места: p >= 5, p не делит Delta, квартика p-целая =>
        # y^2 = g(x,z) — гладкая кривая рода 1 над F_p; #C(F_p) >= p+1-2 sqrt(p) > 0;
        # гладкая точка поднимается по Гензелю.  content(g) — p-единица, т.к. p не делит Delta.
    print("  прочие p: p>=5, p не делит Delta, g p-целая  =>  хорошая редукция рода 1,")
    print("  оценка Хассе p+1-2sqrt(p) > 0 даёт точку, Гензель поднимает её в Q_p.")
    print()

    # ---------------------------------------------------------------- 5. спаривание
    print("--- 5. спаривание Касселса-Тейта (Fisher, Theorem 3.1) -----------------")
    products = {}
    for run in R["pairing_runs"]:
        tag = run["tag"]
        iH = int(run["H_from"])       # индекс квартики, по которой строится H  (0-based)
        ia = int(run["a_from"])       # индекс квартики, чей старший коэффициент — это a
        gH = gs[iH]
        a = QQ(run["a"])
        ok(a != 0, "%s: a != 0" % tag)
        ok(is_square_Q(a / gs[ia][0]) or a == gs[ia][0],
           "%s: a = g%d(1,0) с точностью до квадрата (%s ~ %s)" % (tag, ia + 1, a, gs[ia][0]))
        # G(1,0) * G = H^2 в каждой компоненте L
        Hs = []
        for ph in phis:
            G, H = H_form(gH, ph, I)
            ok(sq_of_quadratic(H) == [G[0] * v for v in G],
               "%s: тождество Гессиана G(1,0) G = H^2 (компонента phi = %s)" % (tag, ph))
            Hs.append(H)
        # fac = z(g_j) z(g_k) / m = m / z(g_H)   (т.к. z1 z2 z3 = m^2)
        fac = [ms[i] / zs[iH][i] for i in range(3)]
        for i in range(3):
            j, k = [t for t in range(3) if t != iH][0], [t for t in range(3) if t != iH][1]
            ok(fac[i] * zs[iH][i] == ms[i] and fac[i] * ms[i] == zs[j][i] * zs[k][i],
               "%s: fac = z(g_j)z(g_k)/m в компоненте %d" % (tag, i))
        # gamma = коэффициент при phi^2 элемента fac * H  (интерполяция Лагранжа в L = Q^3)
        gamma_raw = []
        for j in range(3):
            tot = Fr(0)
            for i in range(3):
                den = Fr(1)
                for kk in range(3):
                    if kk != i:
                        den *= (phis[i] - phis[kk])
                tot += fac[i] * Hs[i][j] / den
            gamma_raw.append(tot)
        ok(gamma_raw == [QQ(v) for v in run["gamma_raw"]],
           "%s: gamma_1 по формуле Фишера совпала с сертификатом" % tag)
        ok(any(c != 0 for c in gamma_raw), "%s: gamma_1 != 0" % tag)
        lam = QQ(run["gamma_scale"])
        ok(lam != 0, "%s: масштаб gamma не равен нулю" % tag)
        gamma = [c * lam for c in gamma_raw]
        ok(gamma == [QQ(v) for v in run["gamma"]], "%s: нормированная gamma" % tag)
        # набор мест: дополнение к достаточным условиям Remark 3.3
        need = {2, 3, 5, 7} | primes_of(disc) | primes_of(a)
        for c in list(gH) + list(gamma):
            if c != 0:
                need |= primes_of(QQ(c.denominator))
        nums = [abs(c.numerator) for c in gamma if c != 0]
        cont = nums[0]
        for t in nums[1:]:
            cont = gcd(cont, t)
        need |= primes_of(QQ(cont))
        places = [v["place"] for v in run["symbols"]]
        ok(len(places) == len(set(places)), "%s: места не повторяются" % tag)
        ok({str(p) for p in need} | {"real"} <= set(places),
           "%s: набор мест ПОЛОН по Remark 3.3" % tag,
           "необходимо: %s" % sorted(str(p) for p in need))
        prod = 1
        minus = []
        for v in run["symbols"]:
            x, z = QQ(v["x"]), QQ(v["z"])
            ok(not (x == 0 and z == 0), "%s @ %s: (x:z) != (0:0)" % (tag, v["place"]))
            gv = quartic_eval(gH, x, z)
            ok(gv == QQ(v["g_value"]) and gv != 0, "%s @ %s: g(x_v,z_v) = %s" % (tag, v["place"], gv))
            ok(is_local_square(gv, v["place"]),
               "%s @ %s: g(x_v,z_v) — квадрат в Q_v" % (tag, v["place"]))
            gam_v = eval_quadratic(gamma, x, z)
            ok(gam_v == QQ(v["gamma_value"]) and gam_v != 0,
               "%s @ %s: gamma(x_v,z_v) = %s != 0" % (tag, v["place"], gam_v))
            h = hilbert(a, gam_v, v["place"])
            ok(h == int(v["hilbert"]), "%s @ %s: (a, gamma)_v = %+d" % (tag, v["place"], h))
            prod *= h
            if h == -1:
                minus.append(v["place"])
        ok(prod == int(run["product"]), "%s: произведение символов = %+d" % (tag, prod))
        print("    %s: мест %d, отрицательный вклад в %s" % (tag, len(places), minus or "нигде"))
        products[tag] = prod
    print()

    # ---------------------------------------------------------------- 6. вывод
    print("--- 6. вывод -----------------------------------------------------------")
    vals = set(products.values())
    ok(len(vals) == 1, "все прогоны дали одно значение спаривания: %s" % products)
    val = vals.pop()
    ok(val == -1, "спаривание <[g1],[g2]>_CT = -1 (нетривиально)")
    fwd = [t for t in products if t.startswith("forward")]
    rev = [t for t in products if t.startswith("reverse")]
    ok(len(fwd) >= 1 and len(rev) >= 1,
       "есть и прямой, и ОБРАТНЫЙ порядок аргументов (формула Фишера по аргументам несимметрична)")
    alt = [t for t in products if "alt" in t]
    ok(len(alt) >= 1, "есть прогон с ДРУГИМИ локальными свидетелями")

    print()
    print("Спаривание Касселса-Тейта обращается в нуль на образе E(Q)/2E(Q) в Sel^2(E/Q)")
    print("(оно индуцировано спариванием на Ш, куда образ E(Q)/2E(Q) идёт в нуль).")
    print("Поэтому класс [g1], чей квадратный класс равен (1, 274, 274),")
    print("НЕ лежит в образе E(Q)/2E(Q).  Значит у E(Q) нет точки с (X-e1, X-e2, X-e3)")
    print("в классе (1, s, s), значит у C_(11,4) нет конечных рациональных точек,")
    print("а на бесконечности их нет, т.к. s не квадрат.   C_(11,4)(Q) = пусто.")
    print()
    print("НИ ОДНА проверка выше не использовала верхнюю границу ранга E(Q).")

    res = {
        "pair": "(11,4)",
        "verifier": "чистый Python (fractions/math/json), без Sage, PARI и кода проекта",
        "checks_total": len(CHECKS),
        "checks_failed": sum(1 for c in CHECKS if not c[0]),
        "target_square_class": ["1", "274", "274"],
        "z_g1_class_matches_target": True,
        "three_quartics_ELS": True,
        "hilbert_products": {k: v for k, v in products.items()},
        "pairing_value": -1,
        "rank_bound_used": False,
        "conclusion": "delta=(1,274,274) не в образе E(Q)/2E(Q); C_(11,4)(Q) пусто",
    }
    outp = cert_path.parent / "independent_result.json"
    outp.write_text(json.dumps(res, indent=2, ensure_ascii=False) + "\n")
    print("записано:", outp)
    return 0


if __name__ == "__main__":
    sys.exit(main())
