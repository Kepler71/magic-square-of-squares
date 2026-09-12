# -*- coding: utf-8 -*-
# СТАТУС: независимая проверка чужого заявления
# ДЛЯ: Claude / Codex / пользователь
# ИТОГ: прямой перебор всех 18 точек P^1(F_17) для полной системы девяти клеток G1 при (m,n)=(13,8)
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЕНО: клетки выведены символически ЗАНОВО из условий магического квадрата, не скопированы
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: проверяете локальное исключение (13,8) mod 17
#
# Запуск:  sage /home/kep/magicKube/bridge/local_13_8_p17_прямая.sage
#
# План:
#   §1  вывести девять клеток G1 из условий магического квадрата (символически, без копирования)
#   §2  проверить магичность и автоматическую квадратность c1,c3,c5,c7 тождественно
#   §3  прямой перебор ВСЕХ p+1 точек P^1(F_p) для p=17, (m,n)=(13,8)
#   §4  разбор вырождений: нули клеток и совпадения клеток — опирается ли вывод на них
#   §5  контроли: (15,8) t=1 должно проходить; (7,1) t=0 должно проходить на трёх клетках
#   §6  17-адическое подтверждение по Гензелю (избыточно, но независимо от §3)
#   §7  сканирование других p и всех (m:n) mod 17

import itertools
from fractions import Fraction as Fr

def hdr(s):
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)

# ============================================================================
# §1. ВЫВОД ДЕВЯТИ КЛЕТОК ИЗ УСЛОВИЙ МАГИЧЕСКОГО КВАДРАТА
# ============================================================================
hdr("§1. Вывод девяти клеток G1 заново, из условий магического квадрата")

R = PolynomialRing(QQ, ['m', 'n', 't'])
m, n, t = R.gens()
Fld = R.fraction_field()

# Схема G1 (BREMNER16 §2.2): полные противоположные пары — рёберные {1,7} и {3,5}.
# Параметризация Сегре квадрики P^2+Q^2 = R^2+S^2 = 2c при q=1, p=t:
P = m*t + n
Q = m - n*t
Rr = m*t - n
S = m + n*t
# автоматические квадраты
c1 = P**2
c7 = Q**2
c3 = Rr**2
c5 = S**2

# центр: 2*c4 = c1 + c7
c4 = Fld(c1 + c7) / 2

# Оставшиеся четыре клетки восстанавливаем ИЗ МАГИЧНОСТИ, не подставляя готовых формул.
# Неизвестные c0, c2, c6, c8. Уравнения:
#   строка 0:      c0 + c1 + c2 = 3 c4
#   строка 2:      c6 + c7 + c8 = 3 c4
#   столбец 0:     c0 + c3 + c6 = 3 c4
#   столбец 2:     c2 + c5 + c8 = 3 c4
#   диагональ:     c0 + c4 + c8 = 3 c4
#   антидиагональ: c2 + c4 + c6 = 3 c4
M = Matrix(Fld, [
    [1, 1, 0, 0],
    [0, 0, 1, 1],
    [1, 0, 1, 0],
    [0, 1, 0, 1],
    [1, 0, 0, 1],
    [0, 1, 1, 0],
])
rhs = vector(Fld, [3*c4 - c1, 3*c4 - c7, 3*c4 - c3, 3*c4 - c5, 2*c4, 2*c4])
sol = M.solve_right(rhs)
c0, c2, c6, c8 = sol
print("Решение линейной системы магичности единственно (ранг матрицы = %d из 4):" % M.rank())
print("  c0 =", R(c0))
print("  c2 =", Fld(c2))
print("  c6 =", Fld(c6))
print("  c8 =", R(c8))
print("  c4 =", Fld(c4))

s_sym = Fld(m**2 + n**2) / 2
print()
print("Сверка с обозначениями проекта (s = (m^2+n^2)/2):")
checks = [
    ("c0 = F0 = m^2 + n^2 t^2", c0 - (m**2 + n**2*t**2)),
    ("c4 = F4 = s (1 + t^2)",   c4 - s_sym*(1 + t**2)),
    ("c8 = F8 = n^2 + m^2 t^2", c8 - (n**2 + m**2*t**2)),
    ("c2 = L  = F4 - 2 m n t",  c2 - (s_sym*(1+t**2) - 2*m*n*t)),
    ("c6 = U  = F4 + 2 m n t",  c6 - (s_sym*(1+t**2) + 2*m*n*t)),
]
for name, diff in checks:
    print("  %-28s разность = %s   %s" % (name, Fld(diff), "OK" if Fld(diff) == 0 else "!!! РАСХОЖДЕНИЕ"))

cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]

# ============================================================================
# §2. ТОЖДЕСТВЕННЫЕ ПРОВЕРКИ
# ============================================================================
hdr("§2. Тождественные проверки (символически, над Q(m,n,t))")

lines = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
ok = True
for L in lines:
    d = Fld(sum(cells[i] for i in L) - 3*c4)
    if d != 0:
        ok = False
        print("  !!! линия", L, "не сходится:", d)
print("  все 8 линий = 3*центр:", "ДА" if ok else "НЕТ")

auto = {1: P, 3: Rr, 5: S, 7: Q}
for i, root in auto.items():
    d = Fld(cells[i] - root**2)
    print("  c%d = (%s)^2 тождественно: %s" % (i, root, "ДА" if d == 0 else "НЕТ"))
print("  => неавтоматических квадратностей ровно пять: c0, c2, c4, c6, c8")

print()
print("  §2a. Устойчивость к переобозначению: в параметризации Сегре есть свобода выбора,")
print("       какая из P,Q,R,S садится в клетку 3, а какая в клетку 5. Проверяю ВТОРОЙ вариант")
print("       (c3 = S^2, c5 = R^2) — он мог бы дать ДРУГИЕ пять форм.")
c1b, c7b, c3b, c5b = P**2, Q**2, S**2, Rr**2
c4b = Fld(c1b + c7b) / 2
rhs_b = vector(Fld, [3*c4b - c1b, 3*c4b - c7b, 3*c4b - c3b, 3*c4b - c5b, 2*c4b, 2*c4b])
c0b, c2b, c6b, c8b = M.solve_right(rhs_b)
five_a = sorted([Fld(c0), Fld(c2), Fld(c4), Fld(c6), Fld(c8)], key=str)
five_b = sorted([Fld(c0b), Fld(c2b), Fld(c4b), Fld(c6b), Fld(c8b)], key=str)
print("       вариант 1, пять форм:", [str(x) for x in five_a])
print("       вариант 2, пять форм:", [str(x) for x in five_b])
print("       МУЛЬТИМНОЖЕСТВА СОВПАДАЮТ:", five_a == five_b,
      "=> пять квадратностей не зависят от этого выбора (переставлены только позиции)")
print("       Аналогично n -> -n и t -> -t лишь меняют L <-> U:",
      Fld(c2).subs({n: -n}) == Fld(c6) and Fld(c2).subs({t: -t}) == Fld(c6))

# ============================================================================
# §3. ПРЯМОЙ ПЕРЕБОР P^1(F_p)
# ============================================================================

def homog_cells(mm, nn, a, b, K):
    """девять клеток, ОДНОРОДИЗОВАННЫЕ умножением на b^2 (t = a/b);
       умножение на квадрат b^2 не меняет квадратности.
       Работает в любом кольце K, где 2 обратима."""
    two = K(2)
    s = K(mm*mm + nn*nn) / two
    a = K(a); b = K(b); mm = K(mm); nn = K(nn)
    F0 = mm**2*b**2 + nn**2*a**2
    F4 = s*(a**2 + b**2)
    F8 = nn**2*b**2 + mm**2*a**2
    L  = F4 - 2*mm*nn*a*b
    U  = F4 + 2*mm*nn*a*b
    C1 = (mm*a + nn*b)**2
    C3 = (mm*a - nn*b)**2
    C5 = (mm*b + nn*a)**2
    C7 = (mm*b - nn*a)**2
    return [F0, C1, L, C3, F4, C5, U, C7, F8]

def squares_mod(p):
    return set((x*x) % p for x in range(p))       # содержит 0

def brute_force_Fp(mm, nn, p, verbose=True):
    """Перебор всех p+1 точек P^1(F_p).
       Возвращает список прошедших точек и подробную таблицу."""
    K = GF(p)
    QR = squares_mod(p)
    pts = [(x, 1) for x in range(p)] + [(1, 0)]
    survivors = []
    table = []
    for (a, b) in pts:
        cs = homog_cells(mm, nn, a, b, K)
        vals = [int(x) for x in cs]
        five = [vals[i] for i in (0, 4, 8, 2, 6)]      # F0, F4, F8, L, U
        flags = [v in QR for v in five]
        # автоматические клетки — квадраты тождественно, проверим и это
        auto_ok = all(vals[i] in QR for i in (1, 3, 5, 7))
        table.append(((a, b), vals, flags, auto_ok))
        if all(flags):
            survivors.append(((a, b), vals))
    return survivors, table, QR

hdr("§3. (m,n) = (13,8), p = 17 — прямой перебор ВСЕХ 18 точек P^1(F_17)")

MM, NN, PP = 13, 8, 17
print("  m^2 = %d, n^2 = %d, 2mn = %d, s = (m^2+n^2)/2 = %d  (всё mod %d)"
      % (MM*MM % PP, NN*NN % PP, 2*MM*NN % PP, int(GF(PP)(MM*MM+NN*NN)/2), PP))
print("  QR_%d (с нулём) = %s" % (PP, sorted(squares_mod(PP))))
print()

surv, tab, QR17 = brute_force_Fp(MM, NN, PP)

print("  t      | F0  F4  F8  |  L   U   | квадраты (F0,F4,F8,L,U)      | итог")
print("  " + "-"*88)
for (ab, vals, flags, auto_ok) in tab:
    a, b = ab
    tname = str(a) if b == 1 else "inf"
    F0, C1, L, C3, F4, C5, U, C7, F8 = vals
    mark = "ПРОХОДИТ" if all(flags) else "отказ"
    # какая именно квадратность ломается
    names = ["F0", "F4", "F8", "L", "U"]
    bad = [names[i] for i, f in enumerate(flags) if not f]
    badvals = [[F0, F4, F8, L, U][i] for i, f in enumerate(flags) if not f]
    why = "" if all(flags) else "  (не квадрат: " + ", ".join("%s=%d" % (x, y) for x, y in zip(bad, badvals)) + ")"
    print("  %-6s | %-3d %-3d %-3d | %-3d %-3d | %-28s | %s%s"
          % (tname, F0, F4, F8, L, U, str([int(f) for f in flags]), mark, why))

print()
print("  проверено точек P^1(F_17): %d (ожидалось %d)" % (len(tab), PP + 1))
print("  прошедших все пять квадратностей: %d" % len(surv))
if surv:
    print("  !!! ЛОКАЛЬНОГО ПРЕПЯТСТВИЯ НЕТ, прошли:", surv)
else:
    print("  => НИ ОДНА точка P^1(F_17) не даёт всех пяти квадратов")
    print("  => [доказано] полная система девяти клеток G1 при (m,n)=(13,8) не имеет решений в Q_17")

# ============================================================================
# §4. ОПИРАЕТСЯ ЛИ ВЫВОД НА ПОЛОЖИТЕЛЬНОСТЬ / РАЗЛИЧНОСТЬ / ЗАПРЕТ НУЛЕЙ
# ============================================================================
hdr("§4. Вырождения: опирается ли отказ на нули клеток или на совпадения клеток")

print("  Критерий в §3 — САМЫЙ СЛАБЫЙ из возможных: значение годится, если оно квадрат в F_17,")
print("  ВКЛЮЧАЯ нуль. Нули разрешены явно (ненулевой рациональный квадрат может обнулиться mod p).")
print("  Положительность и попарная различность НЕ использовались вовсе.")
print()
print("  Где встречаются нули среди пяти клеток и мешает ли это:")
zero_rows = 0
for (ab, vals, flags, auto_ok) in tab:
    a, b = ab
    tname = str(a) if b == 1 else "inf"
    F0, C1, L, C3, F4, C5, U, C7, F8 = vals
    five = [("F0", F0), ("F4", F4), ("F8", F8), ("L", L), ("U", U)]
    zs = [nm for nm, v in five if v == 0]
    if zs:
        zero_rows += 1
        print("    t=%-4s нули в %s — приняты как квадраты (не причина отказа)" % (tname, zs))
print("    строк с нулями:", zero_rows)
print()
print("  Совпадения клеток mod 17 (для сведения; НА ВЫВОД НЕ ВЛИЯЮТ):")
for (ab, vals, flags, auto_ok) in tab[:0]:
    pass
coll = 0
for (ab, vals, flags, auto_ok) in tab:
    a, b = ab
    tname = str(a) if b == 1 else "inf"
    eq = [(i, j) for i in range(9) for j in range(i+1, 9) if vals[i] == vals[j]]
    if eq:
        coll += 1
print("    точек, где есть хоть одно совпадение клеток mod 17: %d из %d" % (coll, len(tab)))
print("    (совпадение mod p НИЧЕГО не значит для различности над Q — здесь оно и не использовалось)")

print()
print("  Контрольный перебор с ЗАПРЕТОМ нулей (более сильное требование) — должен дать не больше:")
surv_nz = []
for (ab, vals, flags, auto_ok) in tab:
    F0, C1, L, C3, F4, C5, U, C7, F8 = vals
    five = [F0, F4, F8, L, U]
    if all(v in QR17 and v != 0 for v in five):
        surv_nz.append(ab)
print("    прошло при запрете нулей: %d; при разрешённых нулях: %d" % (len(surv_nz), len(surv)))
print("    => вывод НЕ усилен запретом нулей: он и так пуст при самом слабом требовании")

print()
print("  §4b. Положительность в этом семействе вообще НЕ является ограничением (доказательство):")
print("       F0 = m^2 b^2 + n^2 a^2 >= 0, F8 = n^2 b^2 + m^2 a^2 >= 0, F4 = s(a^2+b^2) >= 0;")
print("       L = s(a^2+b^2) - 2mn ab,  s - mn = (m-n)^2/2 >= 0, поэтому")
print("       L >= mn(a^2+b^2) - 2mn ab = mn(a-b)^2 >= 0, и то же для U с (a+b)^2.")
R4 = PolynomialRing(QQ, ['M', 'N', 'A', 'B'])
MM4, NN4, AA, BB = R4.gens()
F4h = (MM4**2 + NN4**2)*(AA**2 + BB**2)/2
Lh = F4h - 2*MM4*NN4*AA*BB
Uh = F4h + 2*MM4*NN4*AA*BB
certL = MM4*NN4*(AA - BB)**2 + (MM4 - NN4)**2/2*(AA**2 + BB**2)
certU = MM4*NN4*(AA + BB)**2 + (MM4 - NN4)**2/2*(AA**2 + BB**2)
for nm, expr, cert, sgn in [("L", Lh, certL, "-"), ("U", Uh, certU, "+")]:
    d = expr - cert
    print("       тождество %s = mn(a%sb)^2 + ((m-n)^2/2)(a^2+b^2):  разность = %s  %s"
          % (nm, sgn, d, "OK" if d == 0 else "!!! РАСХОЖДЕНИЕ"))
print("       => при m,n > 0 все пять форм неотрицательны при ЛЮБОМ вещественном t;")
print("          требование положительности пусто, и опираться на него отказ не может.")

# ============================================================================
# §5. КОНТРОЛИ — процедура должна ПРОПУСКАТЬ известные точки
# ============================================================================
hdr("§5. Контроли: тот же код на случаях, где точка ТОЧНО есть")

print("  A. (15,8), t=1 — известная вырожденная точка (девять квадратов с повторами).")
K = QQ
cs = homog_cells(15, 8, 1, 1, QQ)
names9 = ["c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]
print("     клетки над Q при t=1:", dict(zip(names9, [QQ(x) for x in cs])))
allsq = all(QQ(x) >= 0 and QQ(x).is_square() for x in cs)
print("     все девять — рациональные квадраты:", allsq)
print("     суммы линий =", set(sum(QQ(cs[i]) for i in L) for L in lines), " 3*центр =", 3*QQ(cs[4]))
print("     различны ли девять клеток:", len(set(QQ(x) for x in cs)) == 9, "(нет — точка вырожденная)")
print()
print("     этот код на (15,8) НЕ должен давать локального препятствия ни при каком p:")
bad = []
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101]:
    sv, _, _ = brute_force_Fp(15, 8, p)
    if not sv:
        bad.append(p)
print("     p, где перебор пуст (должно быть пусто):", bad if bad else "нет — контроль пройден")
print()
print("  B. (7,1), t=0 — контрпример Codex к смешению C7 и полной девятки.")
cs = homog_cells(7, 1, 0, 1, QQ)
print("     клетки над Q при t=0:", dict(zip(names9, [QQ(x) for x in cs])))
three = [QQ(cs[0]), QQ(cs[4]), QQ(cs[8])]
print("     F0,F4,F8 =", three, "— квадраты:", [x.is_square() for x in three])
five = [QQ(cs[0]), QQ(cs[4]), QQ(cs[8]), QQ(cs[2]), QQ(cs[6])]
print("     L,U =", [QQ(cs[2]), QQ(cs[6])], "— квадраты:", [QQ(cs[2]).is_square(), QQ(cs[6]).is_square()])
print("     => C7 выполнено, полная девятка — нет. Различение Codex подтверждено.")

# ============================================================================
# §6. 17-АДИЧЕСКОЕ ПОДТВЕРЖДЕНИЕ (независимо от §3)
# ============================================================================
hdr("§6. Независимое подтверждение: перебор по дискам в Z_17 (без языка «редукция mod p»)")

def vp_int(x, p):
    x = int(x); p = int(p)
    if x == 0:
        return None
    v = 0
    while x % p == 0:
        x //= p; v += 1
    return v

def is_sq_Qp_int(x, p):
    """x — целое, x != 0: квадрат ли в Q_p (p нечётное)"""
    x = int(x); p = int(p)
    v = vp_int(x, p)
    if v % 2:
        return False
    y = x // p**v
    return pow(y % p, (p-1)//2, p) == 1

def cells4(mm, nn, a, b):
    """ЧЕТЫРЁХКРАТНЫЕ однородные клетки — ровно целые числа, без Fraction.
       cell = G/4, а 4 — квадрат, поэтому класс квадратов у G и у cell один и тот же.
         4*F0 = 4(m^2 b^2 + n^2 a^2)
         4*F4 = 2(m^2+n^2)(a^2+b^2)
         4*F8 = 4(n^2 b^2 + m^2 a^2)
         4*L  = 4*F4 - 8 m n a b
         4*U  = 4*F4 + 8 m n a b
    """
    mm = int(mm); nn = int(nn); a = int(a); b = int(b)
    G4 = 2*(mm*mm + nn*nn)*(a*a + b*b)
    G0 = 4*(mm*mm*b*b + nn*nn*a*a)
    G8 = 4*(nn*nn*b*b + mm*mm*a*a)
    GL = G4 - 8*mm*nn*a*b
    GU = G4 + 8*mm*nn*a*b
    return [G0, G4, G8, GL, GU]

def disk_search(mm, nn, p, KMAX=6):
    """Ветвление по дискам t = r mod p^k в Z_p плюс окрестность бесконечности.
       Ветвь закрывается ТОЛЬКО когда значение ОПРЕДЕЛЕНО на всём диске
       (v_p(G(r)) < k, потому что возмущение r на p^k Z_p меняет G не более чем на p^k)
       и при этом не квадрат. Иначе спускаемся глубже.
       Возвращает (есть_точка / None=НЕ ЗАВЕРШЕНО, список неразрешённых ветвей)."""
    p = int(p)
    undecided = []
    def go(chart, r, k):
        # chart=0: t = r + p^k Z_p ; chart=1: t = 1/u, u = r + p^k Z_p
        a, b = (r, 1) if chart == 0 else (1, r)
        vals = cells4(mm, nn, a, b)
        need_deeper = False
        for x in vals:
            v = vp_int(x, p)              # None, если x == 0 ровно
            if v is not None and v < k:
                if not is_sq_Qp_int(x, p):
                    return False
            else:
                need_deeper = True
        if not need_deeper:
            return True
        if k >= KMAX:
            undecided.append((chart, r, k))
            return None
        res = False
        for z in range(p):
            out = go(chart, r + z*p**k, k+1)
            if out is True:
                return True
            if out is None:
                res = None
        return res
    o1 = go(0, 0, 0)
    if o1 is True:
        return True, undecided
    o2 = go(1, 0, 1)      # u in p Z_p
    if o2 is True:
        return True, undecided
    if o1 is None or o2 is None:
        return None, undecided
    return False, undecided

res, und = disk_search(13, 8, 17)
print("  (13,8) в Q_17, поиск по дискам глубины 6:", res, " неразрешённых ветвей:", len(und))
print("  (False = точки нет; None = НЕ ЗАВЕРШЕНО, это НЕ доказательство)")
res2, und2 = disk_search(15, 8, 17)
print("  контроль (15,8) в Q_17:", res2, " неразрешённых ветвей:", len(und2), " (ожидается True)")
res3, und3 = disk_search(13, 8, 2, KMAX=10)
print("  контроль: (13,8) в Q_2 (там препятствия быть не должно):", res3, "неразрешённых:", len(und3))

# ---------------------------------------------------------------------------
# §6a. Сверка промежуточного шага Codex: сколько конечных t проходят ПЕРВЫЕ ТРИ
# ---------------------------------------------------------------------------
print()
print("  §6a. Сверка промежуточного шага Codex «остаются ровно четыре конечных параметра»:")
three_ok = []
for (ab, vals, flags, auto_ok) in tab:
    a, b = ab
    if b == 0:
        continue
    if flags[0] and flags[1] and flags[2]:      # F0, F4, F8
        F0, C1, L, C3, F4, C5, U, C7, F8 = vals
        three_ok.append((a, (F0, F4, F8), (L, U)))
print("    t, проходящие C7 (F0,F4,F8 квадраты) mod 17:", [x[0] for x in three_ok])
for a, f, lu in three_ok:
    print("      t=%-3d (F0,F4,F8)=%s  (L,U)=%s   L квадрат: %s, U квадрат: %s"
          % (a, f, lu, lu[0] in QR17, lu[1] in QR17))
print("    таблица Codex: t in {2,8,9,15}; (L,U) = (5,4),(1,14),(14,1),(4,5) —",
      "СОВПАЛА" if [x[0] for x in three_ok] == [2, 8, 9, 15]
      and [x[2] for x in three_ok] == [(5, 4), (1, 14), (14, 1), (4, 5)] else "РАСХОЖДЕНИЕ")

# ---------------------------------------------------------------------------
# §6b. Независимый рациональный контроль малой высоты (НЕ доказательство)
# ---------------------------------------------------------------------------
print()
print("  §6b. Прямой рациональный перебор t = a/b, |a|,|b| <= 200, gcd=1 — контроль, НЕ доказательство:")
found = []
H = 200
for b in range(1, H+1):
    for a in range(-H, H+1):
        if gcd(a, b) != 1:
            continue
        G = cells4(13, 8, a, b)
        if all(Integer(x) >= 0 and Integer(x).is_square() for x in G):
            found.append((a, b))
print("     найдено рациональных t с пятью квадратами:", found if found else "нет")
print("     (пустой перебор сам по себе НИЧЕГО не доказывает; доказывает §3)")

# ============================================================================
# §7. ШИРЕ: другие p для (13,8) и все отношения (m:n) mod 17
# ============================================================================
hdr("§7. Шире: какие ещё p ловят (13,8), и таблица запрещённых отношений mod 17")

killers = []
for p in [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]:
    sv, _, _ = brute_force_Fp(13, 8, p)
    if not sv:
        killers.append(p)
print("  простые p <= 113, при которых система пуста над F_p для (13,8):", killers)

print()
print("  Все отношения (m:n) в P^1(F_17), при которых полная система пуста над F_17:")
K = GF(17)
forbidden = []
for lam in range(17):
    sv, _, _ = brute_force_Fp(lam, 1, 17)
    if not sv:
        forbidden.append(lam)
sv0, _, _ = brute_force_Fp(1, 0, 17)      # n = 0 mod 17
print("    lambda = m/n in", forbidden)
print("    n = 0 mod 17 (точка (1:0) отношений): система пуста?", not sv0)
print("    13/8 mod 17 =", int(K(13)/K(8)), "-> в списке:", int(K(13)/K(8)) in forbidden)

hdr("КОНЕЦ")
