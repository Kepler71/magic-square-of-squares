# -*- coding: utf-8 -*-
# Локальная разрешимость кривой рода 2
#     H : Y^2 = F0(t)*F4(t)*F8(t),   (m,n) = (16,5)
# во ВСЕХ местах поля Q.
#
# Модель (целочисленная, эквивалентная над Q):
#     F0 = m^2 + n^2 t^2, F4 = s(1+t^2), F8 = n^2 + m^2 t^2, s = (m^2+n^2)/2
#     Y^2 = F0 F4 F8   <=>   W^2 = f(t) := 4 F0 F4 F8,  W = 2Y.
# Гладкая проективная модель: y^2 = F(x,z) = сумма f_i x^i z^(6-i) в P(1,3,1),
# точки P^1(Q_p) параметризуют (x:z); класс квадрата F(x,z) корректно определён,
# так как F(λx,λz) = λ^6 F(x,z).
#
# Запуск:  sage /home/kep/magicKube/bridge/local_H_16_5.sage

import sys

m, n = 16, 5
s = QQ(m^2 + n^2)/2

Rq.<t> = QQ[]
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
P  = F0*F4*F8

Rz.<T> = ZZ[]
f = Rz(4*P)                       # W^2 = f(T), W = 2Y

print("="*78)
print("H для (m,n) = (%d,%d):  Y^2 = F0*F4*F8" % (m, n))
print("  s  =", s)
print("  F0 =", F0)
print("  F4 =", F4)
print("  F8 =", F8)
print("  P  = F0*F4*F8 =", P)
print("  целая модель  W^2 = f(T) = 4*P,  W = 2Y")
print("  f  =", f)
print("  f  =", factor(f))
print("="*78)

# ---------- 0. Сверка модели с формулами задания ----------
print("\n[0] СВЕРКА МОДЕЛИ (подстановка конкретных t)")
ok_model = True
for tv in [QQ(0), QQ(1), QQ(2), QQ(-3), QQ(1)/2, QQ(-7)/3, QQ(5)/4]:
    lhs = 4*(F0(tv)*F4(tv)*F8(tv))
    rhs = f(tv)
    good = (lhs == rhs)
    ok_model = ok_model and good
    print("   t=%-8s 4*F0F4F8 = %-22s f(t) = %-22s %s"
          % (tv, lhs, rhs, "ok" if good else "РАСХОЖДЕНИЕ"))
# сверка с формой Codex: Y^2 = b (t^6 + a t^4 + a t^2 + 1)
b = s*m^2*n^2
a = QQ(m^2)/n^2 + 1 + QQ(n^2)/m^2
check_ba = (P == b*(t^6 + a*t^4 + a*t^2 + 1))
print("   P == b*(t^6+a t^4+a t^2+1) с b=%s, a=%s : %s" % (b, a, check_ba))
print("   класс квадрата b =", QQ(b).squarefree_part(), "(Codex: 562)")
print("   ВЫВОД: модель совпадает с формулами задания:", ok_model and check_ba)

# ---------- 1. Вещественное место ----------
print("\n[1] ВЕЩЕСТВЕННОЕ МЕСТО R")
print("   F0 = %s: коэффициенты %s > 0 => F0(t) > 0 для всех t в R" % (F0, [m^2, n^2]))
print("   F4 = %s: s = %s > 0 => F4(t) > 0 для всех t в R" % (F4, s))
print("   F8 = %s: коэффициенты %s > 0 => F8(t) > 0 для всех t в R" % (F8, [n^2, m^2]))
print("   => f(t) = 4*F0*F4*F8 > 0 для ВСЕХ t в R (нет вещественных корней:",
      len(f.change_ring(QQ).roots(RR)) == 0, ")")
print("   Явная точка: t=1, f(1) =", f(1), "= 4*281^3;  sqrt =", sqrt(RR(f(1))))
print("   Также старший коэффициент f_6 =", f.leading_coefficient(), "> 0 =>")
print("   и точки на бесконечности над R существуют.")
print("   СТАТУС: доказано  H(R) != пусто.")

# ---------- 2. Инфраструктура для Q_p ----------
disc = f.discriminant()
bad = [p for p, e in factor(disc)] + [2]
bad = sorted(set([ZZ(p) for p in bad if p > 0]))
print("\n[2] ДИСКРИМИНАНТ И ПЛОХИЕ ПРОСТЫЕ")
print("   disc(f) =", factor(disc))
print("   плохие простые (делители disc, плюс 2):", bad)
print("   старший коэффициент f_6 =", factor(f.leading_coefficient()))
print("   f палиндромичен (f_i = f_{6-i}):", all(f[i] == f[6-i] for i in range(7)))
print("   => F(1,z) = f(z), т.е. карта в бесконечности даёт ТОТ ЖЕ многочлен;")
print("      вопрос сводится к: существует ли t в Z_p с f(t) в (Q_p)^2 (включая 0).")


def unit_square_class_ok(hr, e, p):
    """p^e * hr, где hr -- p-адическая единица (целое, p не делит hr).
       True <=> это квадрат в Q_p^*."""
    if e % 2 == 1:
        return False              # нечётная валюация -> не квадрат
    if p == 2:
        return (hr % 8) == 1      # единица Z_2 -- квадрат <=> == 1 mod 8
    return kronecker(hr, p) == 1  # единица Z_p (p нечётно) -- квадрат <=> QR mod p


def prim_part(g, p):
    """g в ZZ[u], g != 0 -> (k, h) с g = p^k * h, h примитивен по p."""
    k = min(ZZ(c).valuation(p) for c in g.coefficients() if c != 0)
    if k == 0:
        return ZZ(0), g
    return ZZ(k), g.map_coefficients(lambda c: ZZ(c) // p**k)


def disk_solvable(h, e, p, depth, maxdepth, trace):
    """Решает: существует ли u в Z_p с p^e * h(u) квадратом в Q_p (0 разрешён)?
       h примитивен по p. Возвращает True / False / None (не решено на maxdepth).

       Корректность:
       * p нечётно, шаг p: если p не делит h(r), то для всех u = r mod p имеем
         h(u) = h(r) mod p, значит h(u) -- единица с тем же классом квадрата
         (класс единицы Z_p^* определяется вычетом mod p). Решение точное.
       * p = 2, шаг 8: если h(r) нечётно, то для u = r mod 8 имеем h(u) = h(r) mod 8
         (так как 8 | u-r и h целочислен), класс единицы Z_2^* определяется mod 8.
       * если p | h(r): рекурсия на поддиск u = r + step*u'.
    """
    step = 8 if p == 2 else p
    undecided = False
    for r in range(step):
        hr = ZZ(h(r))
        if hr % p != 0:
            if unit_square_class_ok(hr, e, p):
                trace.append((depth, r, "HENSEL: p^%d*h(%d) -- квадрат-единица" % (e, r)))
                return True
            continue                       # весь диск строго исключён
        # p | h(r): подставляем u -> r + step*u'
        u = h.parent().gen()
        sub = h(r + step*u)
        if sub.is_zero():
            trace.append((depth, r, "h тождественно 0 на диске -> точка с W=0"))
            return True
        k, h2 = prim_part(sub, p)
        e2 = (e + k) % 2
        if depth >= maxdepth:
            trace.append((depth, r, "ДОСТИГНУТ maxdepth -- НЕ РЕШЕНО"))
            undecided = True
            continue
        res = disk_solvable(h2, e2, p, depth + 1, maxdepth, trace)
        if res is True:
            return True
        if res is None:
            undecided = True
    return None if undecided else False


def local_solvable_p(f, p, maxdepth=40, verbose=True):
    """ПОЛНОЕ решение вопроса H(Q_p) != пусто для W^2 = f(T).
       Покрытие P^1(Q_p) = {t в Z_p} U {t = 1/z, z в p*Z_p}.
       (Здесь f палиндромичен, так что вторая карта даёт тот же многочлен,
        но проверяем её явно -- чтобы процедура не зависела от этого факта.)
       Процедура провалидирована против независимого перебора в
       /home/kep/magicKube/bridge/validate_decider.sage и validate_decider2.sage
       (600 случайных секстик, 0 расхождений, 43 подтверждённых ответа 'точек нет')."""
    # (a) корни f в Q_p -> точка Вейерштрасса
    K = Qp(p, 60)
    rts = f.change_ring(K).roots()
    if len(rts) > 0:
        if verbose:
            print("   p=%-4d корни f в Q_p: %d штук -> точка Вейерштрасса (W=0)" % (p, len(rts)))
        return True, "Weierstrass"
    # (b) карта t в Z_p
    k0, h0 = prim_part(f, p)
    trace = []
    res1 = disk_solvable(h0, k0 % 2, p, 0, maxdepth, trace)
    if res1 is True:
        return True, trace
    # (c) карта в бесконечности: t = 1/z, z = p*w, w в Z_p
    fr = f.reverse()
    w = fr.parent().gen()
    k1, h1 = prim_part(fr(p*w), p)
    res2 = disk_solvable(h1, k1 % 2, p, 0, maxdepth, trace)
    if res2 is True:
        return True, trace
    if res1 is None or res2 is None:
        return None, trace
    return False, trace


# ---------- 3. Точки на бесконечности ----------
print("\n[3] ТОЧКИ НА БЕСКОНЕЧНОСТИ")
lc = f.leading_coefficient()
print("   f_6 =", lc, "=", factor(lc))
print("   Точки на бесконечности гладкой модели: y^2 = f_6 при z=0.")
print("   f_6 -- квадрат в Q ?", ZZ(lc).is_square(), " => над Q точек на бесконечности НЕТ.")
print("   класс квадрата f_6 =", ZZ(lc).squarefree_part())
print("   Над Q_p они есть <=> f_6 -- квадрат в Q_p:")
inf_ok = []
for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 281]:
    v = ZZ(lc).valuation(p)
    u = ZZ(lc) // p**v
    isq = (v % 2 == 0) and (((u % 8) == 1) if p == 2 else (kronecker(u, p) == 1))
    inf_ok.append((p, isq))
    print("      p=%-4d v_p(f_6)=%d  единица=%s  квадрат в Q_p: %s" % (p, v, u, isq))

# ---------- 4. Проверка всех мест ----------
print("\n[4] ПРОВЕРКА H(Q_p) ДЛЯ ВСЕХ ПЛОХИХ ПРОСТЫХ И МАЛЫХ ПРОСТЫХ")
to_check = sorted(set(bad + [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]))
results = {}
for p in to_check:
    res, info = local_solvable_p(f, p)
    results[p] = res
    tag = "ЕСТЬ ТОЧКА" if res is True else ("ТОЧЕК НЕТ" if res is False else "НЕ РЕШЕНО")
    if info != "Weierstrass":
        print("   p=%-4d %s   (рекурсия: %d решающих узлов)" % (p, tag, len(info)))
        for d, r, msg in info[:3]:
            print("        depth=%d r=%d : %s" % (d, r, msg))
    else:
        print("   p=%-4d %s   (корень f в Q_p)" % (p, tag))

# ---------- 5. Независимая проверка: явные точки над Q_p ----------
print("\n[5] НЕЗАВИСИМАЯ ПРОВЕРКА: явный поиск t с f(t) -- квадрат в Q_p")


def is_square_Qp(x, p):
    if x == 0:
        return True
    x = QQ(x)
    v = x.valuation(p)
    if v % 2 != 0:
        return False
    u = x / QQ(p)**v
    num, den = u.numerator(), u.denominator()
    uu = num * inverse_mod(den, p**6) % (p**6)
    if p == 2:
        return (uu % 8) == 1
    return kronecker(uu, p) == 1


for p in to_check:
    found = None
    # карта t в Z_p
    for num in range(0, 400):
        if is_square_Qp(f(num), p):
            found = ("t=%d" % num, f(num))
            break
    if found is None:
        for num in range(-400, 0):
            if is_square_Qp(f(num), p):
                found = ("t=%d" % num, f(num))
                break
    if found is None:
        # рациональные t = u/v с малыми u,v (класс квадрата F(u,v) = v^6 f(u/v))
        for v in range(2, 60):
            for u in range(-60, 61):
                if gcd(u, v) != 1:
                    continue
                val = ZZ(v)**6 * f(QQ(u)/v)
                if is_square_Qp(val, p):
                    found = ("t=%d/%d" % (u, v), val)
                    break
            if found:
                break
    print("   p=%-4d %s" % (p, ("точка: %s, F=%s" % found) if found else
                            "явная точка в переборе НЕ найдена"))

# ---------- 5b. Третий независимый путь: перебор по Z/p^k + sqrt в Qp ----------
print("\n[5b] ТРЕТИЙ ПУТЬ: перебор t по Z/p^k, тест квадрата через Qp(p).is_square()")
print("     (внимание: content(f) = %s, поэтому v_p(f(t)) >= 1 при p=2,281;" % f.content())
print("      тест НЕ требует v_p = 0, а проверяет квадрат в Q_p при любой чётной валюации)")
for p, k in [(2, 12), (3, 8), (5, 6), (7, 6), (11, 5), (281, 2)]:
    K = Qp(p, 40)
    M = p**k
    hits = []
    for t0 in range(M):
        if K(f(t0)).is_square():
            hits.append(t0)
            if len(hits) >= 6:
                break
    dens_hits = []
    for t0 in range(min(M, p**min(k, 4))):
        if K(f(t0)).is_square():
            dens_hits.append(t0)
    print("     p=%-4d k=%-2d первые t с f(t) в (Q_p)^2: %s   (в [0,p^%d): %d/%d)"
          % (p, k, hits, min(k, 4), len(dens_hits), p**min(k, 4)))

# ---------- 6. Хорошие простые p >= 17: граница Вейля ----------
print("\n[6] ХОРОШИЕ ПРОСТЫЕ p, НЕ ВОШЕДШИЕ В СПИСОК")
print("   Плохие простые: %s. Для p вне этого списка редукция -- гладкая кривая" % bad)
print("   рода 2 над F_p, и #H(F_p) >= p+1-4*sqrt(p) (граница Хассе-Вейля).")
for p in [13, 17, 19, 23]:
    print("      p=%-4d p+1-4sqrt(p) = %.3f" % (p, p + 1 - 4*sqrt(RR(p))))
print("   => для всех хороших p >= 17 граница положительна, точки над F_p есть,")
print("      и они гладкие, значит поднимаются в Z_p по лемме Гензеля.")
print("   Единственное хорошее простое < 17 здесь -- p=13, оно проверено явно выше.")
# прямая проверка числа точек над F_p для хороших p
print("   Прямой подсчёт #H(F_p) для хороших p (контроль):")
Hq = HyperellipticCurve(f.change_ring(QQ))
for p in [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]:
    if p in bad:
        continue
    Hp = HyperellipticCurve(f.change_ring(GF(p)))
    print("      p=%-4d #H(F_p) = %d   (нижняя граница Вейля %.2f)"
          % (p, Hp.count_points(1)[0], p + 1 - 4*sqrt(RR(p))))

# ---------- 7. ИТОГ ----------
print("\n" + "="*78)
print("ИТОГ")
allp = sorted(results.keys())
bad_places = [p for p in allp if results[p] is False]
undec = [p for p in allp if results[p] is None]
print("   R            : точка есть (доказано)")
print("   проверено p  :", allp)
print("   H(Q_p)=пусто :", bad_places if bad_places else "НЕТ ТАКИХ")
print("   не решено    :", undec if undec else "нет")
print("   ПРЕПЯТСТВИЕ НАЙДЕНО:", bool(bad_places) )
print("="*78)
