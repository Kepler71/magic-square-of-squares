# local_H_15_8.sage
#
# Локальная разрешимость кривой рода 2
#     H : Y^2 = F0(t)*F4(t)*F8(t)
# для пары G1 (m,n) = (15,8).
#
#     s  = (m^2+n^2)/2
#     F0 = m^2 + n^2 t^2
#     F4 = s (1 + t^2)
#     F8 = n^2 + m^2 t^2
#
# Проверяются:
#   (1) вещественное место;
#   (2) все простые плохой редукции (через дискриминант секстики) + контрольные хорошие простые;
#   (3) точки на бесконечности (u = 1/t, u -> 0).
#
# Метод для Q_p: полный рекурсивный разбор P^1(Q_p) с точным критерием
# «квадрат в Q_p» и явной оценкой точности (лемма Гензеля в форме
# «значение постоянно по модулю достаточной степени p на диске»).
# Процедура даёт ТРИ ответа: True (точка предъявлена), False (доказано, что
# на диске точек нет), 'UNKNOWN' (исчерпана глубина рекурсии).
#
# Запуск:  sage /home/kep/magicKube/bridge/local_H_15_8.sage

import sys

print("="*78)
print("ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ H : Y^2 = F0*F4*F8,  (m,n) = (15,8)")
print("="*78)

R = PolynomialRing(QQ, 't')
t = R.gen()

m, n = 15, 8
assert gcd(m, n) == 1
s = QQ(m**2 + n**2) / 2
F0 = m**2 + n**2 * t**2
F4 = s * (1 + t**2)
F8 = n**2 + m**2 * t**2
f = F0 * F4 * F8

print("\n[0] Модель и сверка с формулами Codex")
print("  m, n      =", m, n)
print("  s         =", s)
print("  F0        =", F0)
print("  F4        =", F4)
print("  F8        =", F8)
print("  f = F0F4F8=", f)

# сверка с записью H : Y^2 = b (t^6 + a t^4 + a t^2 + 1)
b = s * m**2 * n**2
a = QQ(m**2)/n**2 + 1 + QQ(n**2)/m**2
assert f == b * (t**6 + a*t**4 + a*t**2 + 1), "модель не совпала с b(t^6+at^4+at^2+1)"
print("  b = s m^2 n^2 =", b, " = ", factor(b))
print("  a             =", a)
print("  сверка f == b(t^6+a t^4+a t^2+1): OK")

# точечная сверка подстановкой конкретных t (контроль против опечаток)
for tt in [QQ(0), QQ(1), QQ(-1), QQ(2), QQ(3)/5, QQ(-7)/4]:
    lhs = f(tt)
    rhs = (m**2 + n**2*tt**2) * s * (1+tt**2) * (n**2 + m**2*tt**2)
    assert lhs == rhs
print("  поточечная сверка f(t) в t = 0, 1, -1, 2, 3/5, -7/4: OK")

# --- целая модель, Q-изоморфная H -------------------------------------------
# f = (289/2) * P(t),  P = (225+64t^2)(1+t^2)(64+225t^2)
# Y^2 = f  <=>  (2Y/17)^2 = 2*P(t).   Множитель f/g = (17/2)^2 -- квадрат в Q^*,
# поэтому g и f задают Q-изоморфные кривые (Y |-> Y*(17/2)).
P = (m**2 + n**2*t**2) * (1 + t**2) * (n**2 + m**2*t**2)
g = 2 * P
g = R(g)
ratio = f / g
print("\n[0b] Целая модель g(t) = 2*P(t):")
print("  g =", g)
print("  f/g =", ratio, " квадрат в Q^*? ", ratio.is_square(), " sqrt =", sqrt(ratio))
assert ratio.is_square()
assert all(c in ZZ for c in g.coefficients())
gZ = PolynomialRing(ZZ, 't')(g)
print("  g ∈ Z[t]:", gZ)
print("  cont(g) =", gcd(gZ.coefficients()))
print("  g squarefree:", gZ.is_squarefree(), " => род 2 (deg 6, без кратных корней)")

# самовзаимность: t^6 g(1/t) = g(t)
g_rev = R(gZ.reverse())
print("  g реципрокен (t^6 g(1/t) = g(t)):", g_rev == g)

# ============================================================================
# (1) ВЕЩЕСТВЕННОЕ МЕСТО
# ============================================================================
print("\n" + "="*78)
print("[1] ВЕЩЕСТВЕННОЕ МЕСТО")
print("="*78)
print("  F0 = 225 + 64 t^2 > 0 для всех t ∈ R  (сумма положительных)")
print("  F4 = (289/2)(1 + t^2) > 0")
print("  F8 = 64 + 225 t^2 > 0")
print("  => f(t) > 0 для всех t ∈ R; H(R) ≠ ∅.")
print("  Формальная проверка: число вещественных корней f =",
      len(f.roots(RR)), " (должно быть 0)")
print("  min f на [-50,50] по сетке:",
      min([f(QQ(k)/20) for k in range(-1000, 1001)]))
print("  старший коэффициент g:", gZ.leading_coefficient(),
      "> 0 => две вещественные точки на бесконечности (гладкая модель).")
print("  ВЫВОД (доказано): H(R) ≠ ∅.")

# ============================================================================
# Дискриминант и плохие простые
# ============================================================================
print("\n" + "="*78)
print("[2] ДИСКРИМИНАНТ И ПЛОХИЕ ПРОСТЫЕ")
print("="*78)
disc = gZ.discriminant()
print("  disc(g) =", disc)
print("  factor  =", factor(disc))
bad_from_disc = sorted([p for p, e in factor(disc)])
# для модели y^2 = g(x) к плохим всегда добавляем 2 и делители старшего коэфф.
lead_primes = sorted([p for p, e in factor(gZ.leading_coefficient())])
bad = sorted(set(bad_from_disc) | set(lead_primes) | {2})
print("  простые | disc      :", bad_from_disc)
print("  простые | старш.коэф:", lead_primes)
print("  ИТОГО кандидаты в плохие простые:", bad)

# ============================================================================
# Аппарат: квадраты в Q_p и рекурсивный разбор
# ============================================================================

def is_square_qp(x, p):
    """x ∈ Q^*, решить, квадрат ли x в Q_p^*. Точный ответ."""
    x = QQ(x)
    if x == 0:
        return True
    v = x.valuation(p)
    u = x / QQ(p)**v          # p-адическая единица (рациональная)
    if v % 2 != 0:
        return False
    num = u.numerator()
    den = u.denominator()
    if p == 2:
        # u ≡ num * den^{-1} (mod 8); den нечётен
        inv = inverse_mod(den % 8, 8)
        return (num * inv) % 8 == 1
    else:
        inv = inverse_mod(den % p, p)
        return kronecker((num * inv) % p, p) == 1


def _needed_prec(p):
    """Сколько разрядов единичной части нужно, чтобы решить вопрос о квадрате."""
    return 3 if p == 2 else 1


def disc_solvable(poly, p, a, k, depth, maxdepth, witness, budget=None):
    """
    Решает: существует ли x ∈ a + p^k Z_p с poly(x) ∈ (Q_p^*)^2 ∪ {0}?
    poly ∈ Z[t], a ∈ Z, k ≥ 0.
    Возвращает True / False / 'UNKNOWN'. При True кладёт свидетеля в witness.

    Обоснование вердикта False на диске a + p^k Z_p:
      для x = a + p^k s имеем poly(x) = poly(a) + Σ_{j≥1} c_j p^{jk} s^j,
      где c_j = poly^{(j)}(a)/j! ∈ Z. Пусть v0 = v_p(poly(a)),
      mprec = min_j (v_p(c_j) + jk). Если mprec > v0 и mprec - v0 ≥ e
      (e = 1 при p нечётном, e = 3 при p = 2), то на всём диске
      v_p(poly(x)) = v0, а единичная часть poly(x)/p^{v0} сравнима с
      единичной частью poly(a)/p^{v0} по модулю p^{mprec-v0}, т.е. с точностью,
      достаточной для критерия квадрата. Значит квадратичность poly(x)
      одинакова на всём диске и определяется значением в a.
    """
    if budget is not None:
        budget[0] -= 1
        if budget[0] <= 0:
            return 'UNKNOWN'
    val = poly(a)
    if val == 0:
        witness.append(('weierstrass', a, 0))
        return True
    if is_square_qp(val, p):
        witness.append(('point', a, val))
        return True

    # попытка доказать, что на всём диске значений-квадратов нет
    v0 = val.valuation(p)
    # poly(a + p^k s) = val + sum_{j>=1} c_j p^{jk} s^j
    ders = []
    dpoly = poly
    fact = 1
    for j in range(1, poly.degree() + 1):
        dpoly = dpoly.derivative()
        fact *= j
        c_j = ZZ(dpoly(a) / fact)         # биномиальный коэффициент Тейлора, целый
        if c_j != 0:
            ders.append(c_j.valuation(p) + j * k)
    mprec = min(ders) if ders else Infinity

    if mprec - v0 >= _needed_prec(p) and mprec > v0:
        # на всём диске v_p(poly(x)) = v0, единичная часть ≡ const mod p^{mprec-v0}
        # => квадратичность определена значением в a => нигде не квадрат
        return False

    if depth >= maxdepth:
        return 'UNKNOWN'

    res_unknown = False
    for i in range(p):
        r = disc_solvable(poly, p, a + i * p**k, k + 1, depth + 1, maxdepth,
                          witness, budget)
        if r is True:
            return True
        if r == 'UNKNOWN':
            res_unknown = True
    return 'UNKNOWN' if res_unknown else False


def locally_solvable(polyZ, p, maxdepth=40, verbose=False, budget_size=400000):
    """
    Полная проверка H: y^2 = polyZ(x), deg polyZ = 6, над Q_p, включая бесконечность.
    Разбиение P^1(Q_p) = {x ∈ Z_p} ∪ {x = 1/u, u ∈ p Z_p}.
    u = 0 отвечает двум точкам на бесконечности (существуют <=> старший коэф. квадрат).
    """
    assert polyZ.degree() == 6, "процедура написана для секстики (deg = 6)"
    w1 = []
    b1 = [budget_size]
    r1 = disc_solvable(polyZ, p, 0, 0, 0, maxdepth, w1, b1)
    # x = 1/u, y = w/u^3  =>  w^2 = u^6 polyZ(1/u) = rev(u);  нужно u ∈ p Z_p
    rev = PolynomialRing(ZZ, 't')(polyZ.reverse(degree=6))
    w2 = []
    b2 = [budget_size]
    r2 = disc_solvable(rev, p, 0, 1, 0, maxdepth, w2, b2)
    if r1 is True:
        return True, ('x ∈ Z_p', w1[-1])
    if r2 is True:
        return True, ('x = 1/u, u ∈ pZ_p', w2[-1])
    if r1 == 'UNKNOWN' or r2 == 'UNKNOWN':
        return 'UNKNOWN', (r1, r2)
    return False, (r1, r2)


def brute_force_point(polyZ, p, jmax=4, kmax=None):
    """
    НЕЗАВИСИМЫЙ прямой поиск точки: перебор x = a/p^j, a mod p^k, j = 0..jmax,
    плюс отдельно ветка бесконечности. Проверка «квадрат в Q_p» — точная.
    Возвращает найденное x (или 'inf'), либо None.
    """
    if kmax is None:
        kmax = {2: 14, 3: 9, 5: 7, 7: 6, 17: 4, 23: 4}.get(p, 4)
    rev = PolynomialRing(ZZ, 't')(polyZ.reverse(degree=6))
    if is_square_qp(rev(0), p):
        return 'inf'
    for j in range(0, jmax + 1):
        for a in range(0, p**kmax):
            x = QQ(a) / p**j
            val = polyZ(x)
            if val == 0 or is_square_qp(val, p):
                return x
    return None


# --- самотестирование аппарата ---------------------------------------------
print("\n" + "="*78)
print("[2a] САМОТЕСТИРОВАНИЕ АППАРАТА")
print("="*78)
ZT = PolynomialRing(ZZ, 't')
u = ZT.gen()

print("\n  (i) is_square_qp против Sage Qp(p).is_square():")
mism = 0
set_random_seed(1)
for p in [2, 3, 5, 7, 17, 23]:
    K = Qp(p, 60)
    for _ in range(400):
        x = QQ(ZZ.random_element(-10**6, 10**6)) / ZZ.random_element(1, 10**4)
        if x == 0:
            continue
        if is_square_qp(x, p) != K(x).is_square():
            mism += 1
            print("     РАСХОЖДЕНИЕ p=%s x=%s" % (p, x))
print("     расхождений:", mism)
assert mism == 0, "is_square_qp неверен"

print("\n  (ii) Кривые с ответом, разобранным вручную:")
tests = [
    (ZT(-1 - u**6), 3, True,
     "y^2=-1-x^6 /Q_3: x=1 -> -2, а -2 квадрат в Q_3"),
    (ZT(3*(u**6 + 1)), 3, False,
     "y^2=3(x^6+1) /Q_3: x^6+1 ≡ 1 или 2 (mod 3), v(3·u)=1 нечётна"),
    (ZT(3*(u**6 + 1)), 2, False,
     "y^2=3(x^6+1) /Q_2: x чётн -> 3 mod 8; x нечётн -> v=1 нечётна"),
    (ZT(7*(u**6 + 1)), 7, False,
     "y^2=7(x^6+1) /Q_7: x^6 ∈ {0,1} mod 7, v(7·u)=1 нечётна"),
    (ZT(u**6 + 1), 5, True,  "y^2=x^6+1 /Q_5: x=0 -> 1"),
    (ZT(2*u**6 + 2), 2, True, "y^2=2x^6+2 /Q_2: x=1 -> 4"),
    (ZT(u**6 + 2), 2, True,
     "y^2=x^6+2 /Q_2: аффинных нет, но старш.коэф 1 — квадрат => точки на бесконечности"),
    (ZT(5*u**6 + 5), 5, True, "y^2=5(x^6+1) /Q_5: x=2 -> 5·65=325=5^2·13, 13 кв. mod 5"),
]
for poly, pp, expect, note in tests:
    r, wit = locally_solvable(poly, pp, maxdepth=25)
    ok = "OK" if r == expect else "!!! РАСХОЖДЕНИЕ !!!"
    print("     p=%-3s %-66s -> %-7s %s" % (pp, note[:66], r, ok))
    assert r == expect, note

print("\n  (iii) Проверка НАДЁЖНОСТИ вердикта False (главное направление):")
print("       на случайных секстиках: если прямой перебор нашёл точку,")
print("       процедура обязана вернуть True. Обратное (True без свидетеля) — тоже ошибка.")
set_random_seed(20260912)
bad_cases = 0
n_tested = 0
n_false = 0
for p in [2, 3, 5, 7]:
    for trial in range(60):
        cs = [ZZ.random_element(-30, 31) for _ in range(6)] + [ZZ.random_element(1, 31)]
        poly = ZT(list(cs))
        if poly.degree() != 6 or not poly.is_squarefree():
            continue
        n_tested += 1
        r, wit = locally_solvable(poly, p, maxdepth=30)
        bf = brute_force_point(poly, p, jmax=2,
                               kmax={2: 11, 3: 7, 5: 5, 7: 4}[p])
        if r is False:
            n_false += 1
            if bf is not None:
                bad_cases += 1
                print("     !!! p=%s poly=%s : процедура False, перебор нашёл x=%s" % (p, poly, bf))
        if r is True and bf is None:
            bad_cases += 1
            print("     ??? p=%s poly=%s : процедура True, перебор не нашёл (проверь свидетеля %s)" % (p, poly, wit))
print("     протестировано секстик: %d, из них вердикт False: %d, противоречий: %d"
      % (n_tested, n_false, bad_cases))
assert bad_cases == 0, "процедура несогласована с прямым перебором"

# ============================================================================
# (2) ПРОВЕРКА ВО ВСЕХ ПЛОХИХ ПРОСТЫХ + КОНТРОЛЬНЫЕ ХОРОШИЕ
# ============================================================================
print("\n" + "="*78)
print("[3] ПРОВЕРКА H(Q_p) ДЛЯ ПЛОХИХ И КОНТРОЛЬНЫХ ПРОСТЫХ")
print("="*78)

check_primes = sorted(set(bad) | set(primes(2, 60)))
results = {}
for p in check_primes:
    r, wit = locally_solvable(gZ, p, maxdepth=40)
    results[p] = r
    tag = "ПЛОХОЕ " if p in bad else "хорошее"
    print("  p = %-4d %s  H(Q_p) ≠ ∅ ? %-8s  свидетель: %s" % (p, tag, r, wit))
    sys.stdout.flush()

print("\n  Сводка по плохим простым:")
for p in bad:
    print("    p =", p, "->", results[p])

# ============================================================================
# (3) ТОЧКИ НА БЕСКОНЕЧНОСТИ
# ============================================================================
print("\n" + "="*78)
print("[4] ТОЧКИ НА БЕСКОНЕЧНОСТИ")
print("="*78)
lc = gZ.leading_coefficient()
print("  старший коэффициент g:", lc, "=", factor(lc))
print("  lc — квадрат в Q?", ZZ(lc).is_square())
sqfree = ZZ(lc).squarefree_part()
print("  бесквадратная часть lc:", sqfree)
print("  => на гладкой модели две точки над Q_v <=> %s — квадрат в Q_v" % sqfree)
print("  R:   %s > 0 => ДА (две вещественные точки на бесконечности)" % sqfree)
for p in check_primes:
    print("    p = %-4d : %s квадрат в Q_p ? %s" % (p, sqfree, is_square_qp(sqfree, p)))

# ============================================================================
# (4) КОНТРОЛЬ: ГЛОБАЛЬНАЯ ТОЧКА
# ============================================================================
print("\n" + "="*78)
print("[5] КОНТРОЛЬ: ГЛОБАЛЬНАЯ РАЦИОНАЛЬНАЯ ТОЧКА НА H")
print("="*78)
for tt in [QQ(1), QQ(-1)]:
    val = f(tt)
    print("  t = %-3s : F0=%s, F4=%s, F8=%s" % (tt, F0(tt), F4(tt), F8(tt)))
    print("            f(t) = %s = %s ; квадрат? %s ; Y = ±%s"
          % (val, factor(val), val.is_square(), sqrt(val)))
    assert val.is_square()
    # на целой модели g
    vg = g(tt)
    print("            g(t) = %s ; квадрат? %s ; Y_g = ±%s" % (vg, vg.is_square(), sqrt(vg)))

print("\n  => H(Q) ⊇ {(1,±4913), (-1,±4913)} ≠ ∅.")
print("     Следовательно H(Q_v) ≠ ∅ для КАЖДОГО места v поля Q (доказано,")
print("     глобальная точка вкладывается в любое пополнение).")
print("     Замечание: эта точка вырожденная — F0=F4=F8=289, т.е. u0=u4=u8=17,")
print("     все три клетки равны, магического квадрата из РАЗЛИЧНЫХ квадратов не даёт.")

# независимая проверка через Sage HyperellipticCurve
print("\n  Независимая проверка моделью Sage HyperellipticCurve:")
try:
    Hc = HyperellipticCurve(gZ)
    print("    H:", Hc)
    print("    genus =", Hc.genus())
    pt = Hc(QQ(1), sqrt(g(QQ(1))))
    print("    точка на кривой:", pt, " лежит? OK")
    for p in [2, 3, 5, 7, 17]:
        try:
            Cp = Hc.change_ring(GF(p))
            print("    #H(F_%d) (аффинная модель, если гладкая) = %s" % (p, Cp.count_points(1)))
        except Exception as e:
            print("    p=%d: редукция особая/недоступна (%s)" % (p, type(e).__name__))
except Exception as e:
    print("    HyperellipticCurve недоступна:", e)

print("\n" + "="*78)
print("ИТОГ")
print("="*78)
allsolv = all(results[p] is True for p in results)
print("  Все проверенные простые дали H(Q_p) ≠ ∅ :", allsolv)
print("  Вещественное место: H(R) ≠ ∅ (доказано).")
print("  Глобальная точка (1, 4913) => локальная разрешимость всюду (доказано).")
print("  ПРЕПЯТСТВИЯ НЕТ.")
