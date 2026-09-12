# -*- coding: utf-8 -*-
# =====================================================================
#  ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ H : Y^2 = F0(t)*F4(t)*F8(t)  для (m,n) = (15,1)
# =====================================================================
#  s  = (m^2+n^2)/2 ;  F0 = m^2 + n^2 t^2 ;  F4 = s(1+t^2) ;  F8 = n^2 + m^2 t^2
#  (m,n)=(15,1):  s=113,  F0 = t^2+225,  F4 = 113(t^2+1),  F8 = 225t^2+1
#
#  Что делается:
#   0) сверка модели подстановкой конкретных t;
#   1) вещественное место;
#   2) дискриминант -> плохие простые;
#   3) ПОЛНЫЙ рекурсивный разбор Z_p (обе карты: t in Z_p и w=1/t in Z_p)
#      -- решающая процедура, не поиск;
#   4) независимая перепроверка: (а) явные сертификаты Гензеля,
#      (б) грубый перебор по модулю p^k;
#   5) точки на бесконечности;
#   6) хорошие простые p >= 17 -- оценка Хассе-Вейля.
#
#  Запуск:  sage /home/kep/magicKube/bridge/local_H_15_1.sage
# =====================================================================

R = PolynomialRing(QQ, 't'); t = R.gen()
RZ = PolynomialRing(ZZ, 'T'); T = RZ.gen()

m, n = 15, 1
s  = QQ(m^2 + n^2)/2
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
f  = R(F0*F4*F8)

print("="*72)
print("H : Y^2 = f(t),  (m,n) = (%d,%d),  s = %s" % (m, n, s))
print("  F0 =", F0, "   F4 =", F4, "   F8 =", F8)
print("  f  =", f)
assert f.denominator() == 1, "f должен быть целочисленным"
fZ = RZ(f)
lead = fZ.leading_coefficient()
const = fZ.constant_coefficient()
print("  deg f = %d,  старший коэф. = %s = %s,  свободный член = %s"
      % (fZ.degree(), lead, factor(lead), const))
print("  f бесквадратен:", fZ.is_squarefree())
assert fZ.degree() == 6 and fZ.is_squarefree()

# ---------------------------------------------------------------- 0) сверка
print("\n--- 0) СВЕРКА МОДЕЛИ (подстановка конкретных t) " + "-"*20)
ok_model = True
for tv in [QQ(0), QQ(1), QQ(2), QQ(-3), QQ(5)/7, QQ(-11)/4]:
    lhs = f(tv)
    rhs = (m^2 + n^2*tv^2) * (s*(1+tv^2)) * (n^2 + m^2*tv^2)
    same = (lhs == rhs)
    ok_model = ok_model and same
    print("   t=%-7s  F0=%-12s F4=%-14s F8=%-12s  f(t)=%s   [%s]"
          % (tv, m^2+n^2*tv^2, s*(1+tv^2), n^2+m^2*tv^2, lhs, "OK" if same else "РАСХОЖДЕНИЕ"))
print("   модель совпадает с формулами задания:", ok_model)
assert ok_model

# ---------------------------------------------------------- 1) R (вещ. место)
print("\n--- 1) ВЕЩЕСТВЕННОЕ МЕСТО " + "-"*44)
# каждый множитель строго положителен при всех вещественных t:
for name, poly in [("F0", R(F0)), ("F4/s", R(1+t^2)), ("F8", R(F8))]:
    d = poly.discriminant()
    print("   %-5s = %-16s  disc = %-8s  (<0 => нет вещ. корней),  значение в 0: %s"
          % (name, poly, d, poly(0)))
    assert d < 0 and poly(0) > 0
print("   s = %s > 0.  Значит f(t) > 0 для ВСЕХ t из R." % s)
print("   f(0) = %s > 0 ; точка H(R): (t,Y) = (0, %s...)" % (f(0), RR(sqrt(f(0)))))
real_ok = True
print("   ВЫВОД: H(R) != пусто  [доказано]")

# --------------------------------------------------- 2) дискриминант, плохие p
print("\n--- 2) ДИСКРИМИНАНТ И ПЛОХИЕ ПРОСТЫЕ " + "-"*33)
disc_f = fZ.discriminant()
print("   disc(f) =", disc_f)
print("   disc(f) =", factor(disc_f))
bad_from_disc = set(p for p, _ in factor(disc_f))
bad = sorted(bad_from_disc | set([2]) | set(p for p, _ in factor(lead)))
print("   простые | disc(f), плюс 2 и делители старшего коэф.:", bad)
# для рода 2 (y^2 = sextic) хорошая редукция вне {2} u {p | disc}
print("   => плохая редукция возможна только при p из", bad)

# ------------------------------------------------------------------------
#  РЕШАЮЩАЯ ПРОЦЕДУРА: существует ли x in Z_p с f(x) in (Q_p^*)^2 u {0} ?
#  Рекурсия по дискам a + p^k Z_p.
#  h(T) = f(a + p^k T) in Z[T];  m = min val коэффициентов;  g = h/p^m.
#   * m нечётно  -> нужно v(g(T)) нечётно  => g(T) = 0 mod p необходимо;
#                   рекурсия по корням g mod p; если корней нет -> диск пуст.
#   * m чётно    -> нужно g(T) квадрат.
#       p нечётное: если есть t0 с g(t0) != 0 mod p и (g(t0)|p) = +1 -> ЕСТЬ ТОЧКА
#                   (Гензель: g(T) = единица == g(t0) mod p на всём подддиске);
#                   t0 с (g(t0)|p) = -1 -> подддиск пуст;  g(t0)=0 -> рекурсия.
#       p = 2:     смотрим t0 mod 8; g(t0) = 1 mod 8 -> ЕСТЬ ТОЧКА;
#                   g(t0) нечётно и != 1 mod 8 -> подддиск пуст;
#                   g(t0) чётно -> рекурсия (шаг k+3).
#  Завершаемость: f бесквадратен => рекурсия конечна. Глубина контролируется.
# ------------------------------------------------------------------------
MAXDEPTH = 60

class DepthExceeded(Exception):
    pass

def _min_val(h, p):
    return min(ZZ(c).valuation(p) for c in h.coefficients() if c != 0)

def soluble_chart(fpoly, p, maxdepth=MAXDEPTH, verbose=False):
    """Есть ли x in Z_p с fpoly(x) квадратом в Q_p (включая 0)?
       Возвращает (answer:bool, witness or None, stats dict)."""
    stats = {'nodes': 0, 'maxdepth': 0, 'capped': False}
    hZ = RZ(fpoly)

    def rec(a, k, depth):
        stats['nodes'] += 1
        stats['maxdepth'] = max(stats['maxdepth'], depth)
        if depth > maxdepth:
            stats['capped'] = True
            raise DepthExceeded()
        h = hZ(a + p**k * T)
        if h == 0:
            return (True, (a, k, 'f тождественно 0'))
        mval = _min_val(h, p)
        g = RZ([ZZ(c)//p**mval for c in h.list()])
        if mval % 2 == 1:
            # нужно v(g) нечётно => g(t0) = 0 mod p
            for t0 in range(p):
                if g(t0) % p == 0:
                    r = rec(a + p**k*t0, k+1, depth+1)
                    if r[0]:
                        return r
            return (False, None)
        else:
            if p != 2:
                roots = []
                for t0 in range(p):
                    v = ZZ(g(t0)) % p
                    if v == 0:
                        roots.append(t0)
                    elif kronecker(v, p) == 1:
                        return (True, (a + p**k*t0, 'Hensel: v(f)=%d чёт., unit — КВ mod %d' % (mval, p)))
                for t0 in roots:
                    r = rec(a + p**k*t0, k+1, depth+1)
                    if r[0]:
                        return r
                return (False, None)
            else:
                evens = []
                for t0 in range(8):
                    v = ZZ(g(t0)) % 8
                    if v == 1:
                        return (True, (a + p**k*t0, 'Hensel: v(f)=%d чёт., unit = 1 mod 8' % mval))
                    if v % 2 == 0:
                        evens.append(t0)
                    # v нечётно и != 1 mod 8 -> подддиск пуст
                for t0 in evens:
                    r = rec(a + p**k*t0, k+3, depth+1)
                    if r[0]:
                        return r
                return (False, None)

    try:
        ans, wit = rec(0, 0, 0)
    except DepthExceeded:
        return (None, None, stats)     # None = НЕ ОПРЕДЕЛЕНО
    return (ans, wit, stats)

def reverse_poly(fpoly, d=6):
    """f*(w) = w^d f(1/w) — карта в бесконечности."""
    c = RZ(fpoly).list()
    c = c + [0]*(d+1-len(c))
    return RZ(list(reversed(c)))

fstar = reverse_poly(fZ, 6)
print("\n   карта бесконечности: f*(w) = w^6 f(1/w) =", fstar)
print("   f* бесквадратен:", fstar.is_squarefree(), " f*(0) =", fstar(0), "(= старший коэф. f)")
assert fstar.is_squarefree()

# --------------------------- самотест решающей процедуры на известных примерах
print("\n--- ПРОВЕРКА САМОЙ ПРОЦЕДУРЫ на контрольных примерах " + "-"*17)
def brute_force_modpk(poly, p, k):
    """Есть ли (x,y) mod p^k с y^2 = poly(x)?  Отсутствие => нет точек в этой карте."""
    M = p**k
    sq = set((y*y) % M for y in range(M))
    for x in range(M):
        if ZZ(poly(x)) % M in sq:
            return x
    return None

selftests = [
    # (полином, p, ожидаемый ответ, комментарий)
    (RZ(3),                3, False, "y^2=3 над Q_3: v(3)=1 нечётна"),
    (RZ(3),                2, False, "y^2=3 над Q_2: 3 != 1 mod 8"),
    (RZ(2),                7, True,  "y^2=2 над Q_7: 2 — КВ mod 7"),
    (RZ(T^6+1),            5, True,  "y^2=x^6+1: точка (0,1)"),
    (RZ(5*(T^2+1)),        5, True,  "y^2=5(x^2+1) над Q_5: f(2)=25 — КВАДРАТ (мой ручной прогноз был ошибочен)"),
    (RZ(3*(T^6+1)),        3, False, "СЕКСТИКА БЕЗ Q_3-ТОЧЕК: x^6+1 != 0 mod 3 всегда => v=1 нечётна"),
    (RZ(3*(T^6+1)),        2, False, "СЕКСТИКА БЕЗ Q_2-ТОЧЕК: x чёт -> 3 mod 8; x нечёт -> v=1"),
    (RZ(7*(T^6+T^3+1)),    7, False, "y^2=7*g, g != 0 mod 7 при всех x => v нечётна"),
]
selftest_ok = True
for poly, pp, expect, cmt in selftests:
    a, w, st = soluble_chart(poly, pp)
    good = (a == expect)
    selftest_ok = selftest_ok and good
    # независимый контроль: грубый перебор mod p^k
    kk = 6 if pp <= 3 else 4
    bf = brute_force_modpk(poly, pp, kk)
    bf_consistent = (bf is not None) if a else True   # a=True => перебор обязан найти
    if (not a) and (bf is not None):
        # рекурсия говорит "точек нет", а перебор mod p^k нашёл — это НЕ противоречие само по
        # себе только если перебор ещё не достиг нужной глубины; отмечаем явно
        bf_note = "перебор mod %d^%d ещё находит x=%s (нужна большая глубина)" % (pp, kk, bf)
    elif (not a):
        bf_note = "перебор mod %d^%d НЕ находит решений — подтверждает" % (pp, kk)
    else:
        bf_note = "перебор mod %d^%d находит x=%s — подтверждает" % (pp, kk, bf)
    print("   p=%-4d %-22s -> %-6s [ожид %-6s %s]  (узлов %d, глуб %d, обрез %s)"
          % (pp, str(poly)[:22], str(a), str(expect), "OK" if good else "!!! РАСХОЖДЕНИЕ",
             st['nodes'], st['maxdepth'], st['capped']))
    print("        %s ; %s" % (cmt, bf_note))
print("   ВСЕ САМОТЕСТЫ ПРОЙДЕНЫ:", selftest_ok)

# ---------------------------------------------- 3) перебор всех мест p
print("\n--- 3) ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ ПО ПРОСТЫМ " + "-"*29)
primes_to_check = sorted(set(bad) | set(primes(2, 20)))
print("   проверяем p in", primes_to_check)
print("   (плохие простые + все p < 20; для хороших p >= 17 ниже — Хассе-Вейль)\n")

results = {}
for p in primes_to_check:
    a1, w1, st1 = soluble_chart(fZ, p)
    a2, w2, st2 = soluble_chart(fstar, p)
    if a1 is None or a2 is None:
        verdict = None
    else:
        verdict = bool(a1 or a2)
    results[p] = (verdict, a1, w1, st1, a2, w2, st2)
    print("   p = %-4d  карта t in Z_p : %-5s  (узлов %3d, глуб %2d, обрез %s)"
          % (p, str(a1), st1['nodes'], st1['maxdepth'], st1['capped']))
    print("            карта w in Z_p : %-5s  (узлов %3d, глуб %2d, обрез %s)"
          % (str(a2), st2['nodes'], st2['maxdepth'], st2['capped']))
    if w1: print("            свидетель (t-карта):", w1)
    elif w2: print("            свидетель (w-карта):", w2)
    print("            => H(Q_%d) != пусто : %s" % (p, verdict))

# --------------------------------------- 4) НЕЗАВИСИМЫЕ ЯВНЫЕ СЕРТИФИКАТЫ
print("\n--- 4) НЕЗАВИСИМЫЕ ЯВНЫЕ СЕРТИФИКАТЫ ГЕНЗЕЛЯ " + "-"*26)
print("   ищем целое t0 с f(t0) = p^(2a)*u, u in Z_p^*, u — квадрат в Z_p^*")
print("   (p нечёт.: (u|p)=+1 ; p=2: u = 1 mod 8). Это даёт Y in Q_p явно.\n")

def hensel_cert(p, trange=200):
    for t0 in range(0, trange):
        for sgn in ([1] if t0 == 0 else [1, -1]):
            x = sgn*t0
            val = ZZ(fZ(x))
            if val == 0:
                return (x, 0, 0, "f(t0)=0, точка (t0,0)")
            a = val.valuation(p)
            if a % 2: continue
            u = val // p**a
            if p != 2:
                if kronecker(u, p) == 1:
                    return (x, a, u, "(u|%d)=+1" % p)
            else:
                if u % 8 == 1:
                    return (x, a, u, "u = 1 mod 8")
    return None

cert_all_ok = True
cert_table = []
for p in primes_to_check:
    c = hensel_cert(p)
    if c is None:
        cert_all_ok = False
        print("   p = %-4d  явного целого t0 с |t0|<200 не найдено (см. решающую процедуру выше)" % p)
        cert_table.append((p, None))
    else:
        x, a, u, why = c
        print("   p = %-4d  t0 = %-4d  f(t0) = %-22s = %d^%d * %-18s  [%s]"
              % (p, x, fZ(x), p, a, u, why))
        cert_table.append((p, (x, a, u)))

# 4b) ТРЕТЬЯ, НЕЗАВИСИМАЯ ПРОВЕРКА: явное извлечение корня в Q_p средствами Sage
print("\n   независимая проверка через Qp(p, prec): Y = sqrt(f(t0)) в Q_p, затем Y^2 - f(t0) = 0")
for p, c in cert_table:
    if c is None: continue
    x, a, u = c
    K = Qp(p, 60)
    val = K(fZ(x))
    issq = val.is_square()
    Y = val.sqrt() if issq else None
    resid = (Y^2 - val) if Y is not None else None
    print("     p=%-4d t0=%-5d f(t0) is_square в Q_p: %-5s  Y = %s ; Y^2-f(t0) = %s"
          % (p, x, issq, (str(Y)[:34] + "...") if Y is not None else "-",
             resid if resid is None else resid.add_bigoh(50)))
    assert issq, "сертификат для p=%d не подтвердился Sage-ом!" % p

# 4c) структурное замечание: при p = 1 mod 4 у f есть КОРЕНЬ в Z_p
print("\n   структурно: 225 = 15^2, поэтому корни f — это ±15i, ±i, ±i/15.")
print("   Все три условия равносильны '-1 — квадрат в Q_p', т.е. p = 1 mod 4.")
for p in primes_to_check:
    if p == 2: continue
    has_root = (p % 4 == 1)
    if has_root:
        K = Qp(p, 40); r = K(-1).sqrt()
        chk = K(fZ(0))  # просто метка
        print("     p=%-4d = 1 mod 4 -> i in Z_p, f(i)=0 -> точка (i,0) на H(Q_%d)  [v_p(i)=%d]"
              % (p, p, r.valuation()))
    else:
        print("     p=%-4d = 3 mod 4 -> корней f в Q_p нет, нужен сертификат Гензеля выше" % p)

# грубая перепроверка перебором mod p^k (независимо от рекурсии)
print("\n   перепроверка грубым перебором mod p^k (наличие решения y^2 = f(x)):")
for p in primes_to_check:
    k = 6 if p <= 3 else (4 if p <= 7 else (3 if p <= 13 else 2))
    if p == 113: k = 2
    xf = brute_force_modpk(fZ, p, k)
    xg = brute_force_modpk(fstar, p, k)
    print("     p=%-4d k=%d : решение в t-карте x=%s ; в w-карте w=%s" % (p, k, xf, xg))

# ------------------------------------------------- 5) точки на бесконечности
print("\n--- 5) ТОЧКИ НА БЕСКОНЕЧНОСТИ " + "-"*40)
print("   На гладкой модели Y^2 = f(t) (deg 6) точки над t=inf существуют над полем K")
print("   тогда и только тогда, когда старший коэффициент — квадрат в K.")
print("   старший коэф. lc = %s = %s" % (lead, factor(lead)))
print("   бесквадратная часть lc = %s ; lc — квадрат в Q: %s"
      % (ZZ(lead).squarefree_part(), ZZ(lead).is_square()))
for p in primes_to_check:
    v = ZZ(lead).valuation(p)
    u = ZZ(lead) // p**v
    if v % 2:
        issq = False
    elif p != 2:
        issq = (kronecker(u, p) == 1)
    else:
        issq = (u % 8 == 1)
    print("     p=%-4d : v_p(lc)=%d, unit=%-8s -> lc квадрат в Q_p: %s   (точки на бескон. над Q_%d: %s)"
          % (p, v, u, issq, p, issq))
print("     R      : lc = %s > 0 -> квадрат в R: True" % lead)
print("   NB: f(0) = %s = lc — свободный член и старший коэф. совпадают" % const)

# ------------------------------------- 6) хорошие простые: Хассе-Вейль
print("\n--- 6) ХОРОШИЕ ПРОСТЫЕ p >= 17 " + "-"*40)
print("   Плохие простые: %s. Все они < 17, кроме p = 113." % bad)
print("   Для хорошей редукции (род 2): #H(F_p) >= p + 1 - 4*sqrt(p).")
for p in [11, 13, 17, 19, 23]:
    print("     p=%-4d : p+1-4*sqrt(p) = %.3f  %s" % (p, p+1-4*sqrt(float(p)),
          "> 0 => гладкая точка в F_p => Гензель" if p+1-4*sqrt(float(p)) > 0 else "<= 0 => нужна прямая проверка"))
print("   => все хорошие p >= 17 автоматически разрешимы; p = 113 проверено выше явно.")

# избыточная, но дешёвая подстраховка: явный сертификат Гензеля для ВСЕХ p <= 500
print("\n   подстраховка: явный целый t0 с f(t0) in (Q_p*)^2 для всех p <= 500")
missing = []
for p in primes(2, 500):
    c = hensel_cert(p, trange=60)
    if c is None:
        missing.append(p)
print("     простые p <= 500 без явного t0 (|t0|<60):", missing if missing else "нет — сертификат найден для каждого p")

# =====================================================================
#  7) ПОЛНЫЙ РАЗБОР СЛУЧАЕВ ДЛЯ ВСЕХ ПРОСТЫХ  (это и есть ДОКАЗАТЕЛЬСТВО;
#     секции 3-4 выше — вычислительная проверка, но перебор p<=500 сам по
#     себе доказательством для ВСЕХ p не является).
# =====================================================================
print("\n--- 7) ПОЛНЫЙ РАЗБОР СЛУЧАЕВ ДЛЯ ВСЕХ ПРОСТЫХ " + "-"*26)
print("   f = 113*(t^2+1)*(t^2+225)*(225t^2+1);  disc(f) = %s" % factor(disc_f))
assert sorted(bad_from_disc | {2}) == [2, 3, 5, 7, 113], "неожиданный набор плохих простых"
print("   плохие простые РОВНО: {2,3,5,7,113};  f палиндромичен: %s" % (fstar == fZ))
assert fstar == fZ            # => инволюция t -> 1/t, карта бескон. = карта нуля

print("\n   A) p = 1 mod 4.  Тогда i = sqrt(-1) in Z_p и f(i) = 0 (множитель t^2+1).")
print("      f бесквадратен => (i,0) — гладкая точка => (i,0) in H(Q_p).  [доказано]")
A_ok = True
for p in [q for q in primes(3, 400) if q % 4 == 1]:
    K = Qp(p, 40); i = K(-1).sqrt()
    A_ok = A_ok and (f.change_ring(K)(i).valuation() > 30)
print("      численная сверка f(i)=0 для всех p = 1 mod 4, p < 400:", A_ok)
assert A_ok
print("      => покрыты 5, 13, 17, 29, 37, ..., и плохое p = 113 (113 = 1 mod 4).")

print("\n   B) p = 2.  f(0) = %s = 15^2 * 113,  113 mod 8 = %d  => 113 in (Z_2^*)^2."
      % (fZ(0), 113 % 8))
B_ok = Qp(2, 60)(fZ(0)).is_square()
print("      Qp(2)(f(0)).is_square() =", B_ok, " => точка (0, 15*sqrt(113)) in H(Q_2).  [доказано]")
assert B_ok

print("\n   C) p in {3,7,11}  — это ВСЕ простые p = 3 mod 4 с p < 17.")
C_ok = True
for p, t0 in [(3, 1), (7, 0), (11, 0)]:
    v = ZZ(fZ(t0)); a = v.valuation(p); u = v // p**a
    sq = Qp(p, 60)(v).is_square()
    C_ok = C_ok and sq
    print("      p=%-3d t0=%d  f(t0)=%-10s=%-20s v_p=%d  unit mod p=%-3d (unit|p)=%+d  квадрат в Q_p: %s"
          % (p, t0, v, str(factor(v)), a, u % p, kronecker(u, p), sq))
assert C_ok
print("      [доказано: Гензель для y^2 = единица, единица — КВ mod p, p нечётно]")

print("\n   D) p = 3 mod 4 и p >= 19.  Плохие простые = {2,3,5,7,113}; из них")
print("      сравнимы с 3 mod 4 только 3 и 7, оба < 19 => все такие p ХОРОШЕЙ редукции.")
print("      Вейль (род 2): #H(F_p) >= p + 1 - 4*sqrt(p) > 0 при p >= 17.")
print("      Гладкая модель над Z_p + гладкая F_p-точка => Гензель => H(Q_p) != пусто. [доказано]")
D_ok = True
worst = None
for p in [q for q in primes(19, 300) if q % 4 == 3]:
    fp = fZ.change_ring(GF(p))
    good = fp.is_squarefree() and fp.degree() == 6
    N = HyperellipticCurve(fp).count_points(1)[0]
    bnd = p + 1 - 4*sqrt(float(p))
    D_ok = D_ok and good and (N > 0) and (float(N) >= bnd)
    if worst is None or N - bnd < worst[1]:
        worst = (p, N - bnd)
    if p < 50:
        print("      p=%-5d хорошая редукция: %-5s  #H(F_p) = %-5d >= %.2f" % (p, good, N, bnd))
print("      ... все p = 3 mod 4 в [19,300) проверены явным счётом точек:", D_ok)
print("      минимальный запас над границей Вейля на этом отрезке:", worst)
print("      p=17: p+1-4*sqrt(p) = %.4f > 0 (граница положительна начиная с 17)"
      % (17 + 1 - 4*sqrt(17.0)))
assert D_ok

print("\n   ПОКРЫТИЕ A u B u C u D = все простые числа:")
uncovered = [p for p in primes(2, 2000)
             if not (p % 4 == 1 or p == 2 or p in (3, 7, 11) or (p % 4 == 3 and p >= 19))]
print("      простые p < 2000 вне всех четырёх случаев:", uncovered if uncovered else "НЕТ")
assert not uncovered
print("      (аргумент чисто арифметический: любое нечётное p = 1 или 3 mod 4;")
print("       при 3 mod 4 либо p in {3,7,11}, либо p >= 19.)")
case_split_ok = A_ok and B_ok and C_ok and D_ok and not uncovered
print("\n   ИТОГ СЕКЦИИ 7: H(Q_p) != пусто для ВСЕХ простых p :", case_split_ok)

# ------------------------------------------------------------ ИТОГ
print("\n" + "="*72)
print("ИТОГ")
undecided = [p for p in primes_to_check if results[p][0] is None]
noloc     = [p for p in primes_to_check if results[p][0] is False]
print("  места с ДОКАЗАННЫМ отсутствием точек : %s" % (noloc if noloc else "нет"))
print("  места, где процедура НЕ завершилась   : %s" % (undecided if undecided else "нет"))
print("  вещественное место                    : разрешимо (f > 0 всюду)")
print("  разбор случаев секции 7 (ВСЕ простые)   : %s" % case_split_ok)
if not noloc and not undecided and case_split_ok:
    print("  ВЫВОД [доказано]: H локально разрешима во ВСЕХ местах Q (R и все p).")
    print("  => локального препятствия НЕТ; отсутствие рациональных точек (если оно есть)")
    print("     должно объясняться Брауэром-Манином или просто пустотой H(Q).")
else:
    print("  ВЫВОД: см. список выше.")
print("="*72)
