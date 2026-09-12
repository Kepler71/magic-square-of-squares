# -*- coding: utf-8 -*-
# СТАТУС: скептическая перепроверка чужого утверждения (попытка опровержения)
# ДЛЯ: Codex (автор утверждения), пользователь, Grok
# ИТОГ: см. вывод скрипта; вывод скрипта = данные, статусы расставлены в отчёте
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЯЕТ: PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md, §«(13,8), p=17»
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: сомневаетесь в локальном исключении (13,8) mod 17
#
# Запуск:  sage /home/kep/magicKube/bridge/local_13_8_p17_скептик.sage
#
# ВНИМАНИЕ. Рядом лежит local_13_8_p17_скептик.py другого скептика (22:02) — он НЕ читался
# при написании этого файла, чтобы проверки остались независимыми.
#
# Что делает скрипт (в порядке нарастания скепсиса):
#   A. выводит девять клеток G1 С НУЛЯ: из линейных аксиом магического квадрата + полной
#      рациональной параметризации квадрики P^2+Q^2=R^2+S^2, а НЕ из формул Codex;
#   B. сверяет их с пятью формами Codex (F0,F4,F8,L,U) как многочлены;
#   C. проверяет, что параметризация ПОЛНАЯ (обратное восстановление m,n,t по клеткам)
#      и что масштабной свободы, портящей редукцию, нет;
#   D. перебирает ВСЕ 18 точек P^1(F_17) при разных соглашениях о нуле;
#   E. проверяет отсутствие точки уже над Q_17 (Гензель), а не только над F_17;
#   F. положительные контроли: (7,1) t=0 и (15,8) t=1 — машина ОБЯЗАНА их не исключать;
#   G. проверяет корректность выбора p (делимости знаменателей) и перебирает другие p;
#   H. независимый рациональный поиск t=a/b (контроль, НЕ доказательство).

import itertools
from fractions import Fraction as Fr

OUT = []
def say(s=""):
    OUT.append(str(s))
    print(s)

say("=" * 100)
say("СКЕПТИК: локальное исключение (m,n) = (13,8) по модулю 17")
say("=" * 100)

# ---------------------------------------------------------------------------
# A. НЕЗАВИСИМЫЙ ВЫВОД ДЕВЯТИ КЛЕТОК G1
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("A. Вывод девяти клеток G1 с нуля (не из формул Codex)")
say("-" * 100)

R9 = PolynomialRing(QQ, ['m', 'n', 't', 'A', 'B', 'C'])
m, n, t, A, B, C = R9.gens()

# A1. Общий магический квадрат 3x3 над Q: 8 линейных условий на 9 клеток.
# Решаем систему честно, а не берём готовую параметризацию Бремнера.
V = PolynomialRing(QQ, ['x%d' % i for i in range(9)])
xs = V.gens()
rows = [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
cols = [(0, 3, 6), (1, 4, 7), (2, 5, 8)]
dias = [(0, 4, 8), (2, 4, 6)]
lines = rows + cols + dias
Mrows = []
for L in lines[1:]:
    v = [0] * 9
    for i in L:
        v[i] += 1
    for i in lines[0]:
        v[i] -= 1
    Mrows.append(v)
Mat = Matrix(QQ, Mrows)
ker = Mat.right_kernel().basis()
say("  Размерность пространства магических квадратов 3x3 над Q: %d (ожидается 3)" % len(ker))
say("  Базис ядра (как векторы клеток 0..8):")
for v in ker:
    say("    %s" % (list(v),))

# A2. Проверяем следствие: c_i + c_{8-i} = 2*c_4 для всех i — на базисе ядра.
ok_pairs = all(all(v[i] + v[8 - i] == 2 * v[4] for i in range(9)) for v in ker)
say("  Тождество c_i + c_{8-i} = 2*c_4 выполняется на всём пространстве: %s" % ok_pairs)

# A3. Конфигурация G1: квадратны пары клеток {1,7} и {3,5}.
#     c1 = P^2, c7 = Q^2, c3 = R^2, c5 = S^2, и P^2+Q^2 = R^2+S^2 = 2c.
#     Полная рациональная параметризация этой квадрики в P^3 (она P^1 x P^1):
P = m * t + n * 1
Q = m * 1 - n * t
Rr = m * t - n * 1
S = m * 1 + n * t
say()
say("  Квадрика P^2+Q^2 = R^2+S^2, параметризация (P,Q,R,S) = (mt+n, m-nt, mt-n, m+nt):")
say("    P^2+Q^2 - (R^2+S^2) = %s   (должно быть 0)" % (P**2 + Q**2 - Rr**2 - S**2))
say("    P^2+Q^2 = %s" % (P**2 + Q**2))

# A4. Восстанавливаем ВСЕ девять клеток из четырёх известных, решая линейную систему.
# Неизвестные: c0..c8. Известно: c1=P^2, c3=R^2, c5=S^2, c7=Q^2 и магичность.
# Решаем символически: c4 = (P^2+Q^2)/2, далее по линиям.
c4 = (P**2 + Q**2) / 2
c1 = P**2
c7 = Q**2
c3 = Rr**2
c5 = S**2
# строка 1: c3+c4+c5 = 3c4  -> проверка совместности
consist = c3 + c4 + c5 - 3 * c4
say("  Совместность (c3+c4+c5 = 3c4): %s  (должно быть 0)" % consist)
# столбец 1: c1+c4+c7 = 3c4 -> проверка
say("  Совместность (c1+c4+c7 = 3c4): %s  (должно быть 0)" % (c1 + c4 + c7 - 3 * c4))
# теперь c0,c2,c6,c8 из линий:
# строка 0: c0+c1+c2 = 3c4 ; столбец 0: c0+c3+c6 = 3c4 ; диагональ: c0+c4+c8 = 3c4
# антидиагональ: c2+c4+c6 = 3c4 ; c8 = 2c4-c0 ; c6 = 2c4-c2
# Из диагонали: c8 = 2c4 - c0 (уже). Нужно ещё одно уравнение, связывающее c0 и c2:
# столбец 0: c0 + c3 + c6 = 3c4, c6 = 2c4 - c2  =>  c0 + c3 + 2c4 - c2 = 3c4 => c0 - c2 = c4 - c3
# строка 0: c0 + c1 + c2 = 3c4  =>  c0 + c2 = 3c4 - c1
c0 = ((c4 - c3) + (3 * c4 - c1)) / 2
c2 = ((3 * c4 - c1) - (c4 - c3)) / 2
c8 = 2 * c4 - c0
c6 = 2 * c4 - c2
cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]

say()
say("  Девять клеток как многочлены от m, n, t (мой независимый вывод):")
names = ["c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"]
for nm, cc in zip(names, cells):
    say("    %s = %s" % (nm, cc))

# A5. Контроль: получившийся квадрат действительно магический.
bad = []
for L in lines:
    ssum = sum(cells[i] for i in L)
    if ssum != 3 * c4:
        bad.append((L, ssum))
say()
say("  Все 8 линий дают 3*c4: %s   (нарушения: %s)" % (not bad, bad))

# A6. Контроль: четыре клетки — тождественные квадраты линейных форм.
auto = {1: P, 3: Rr, 5: S, 7: Q}
say("  Автоматические квадраты: " + ", ".join(
    "c%d = (%s)^2 : %s" % (i, f, cells[i] - f**2 == 0) for i, f in auto.items()))

# ---------------------------------------------------------------------------
# B. СВЕРКА С ФОРМАМИ CODEX
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("B. Сверка моих клеток с пятью формами Codex (F0, F4, F8, L, U)")
say("-" * 100)

Rab = PolynomialRing(QQ, ['m', 'n', 'a', 'b'])
mm, nn, aa, bb = Rab.gens()
s_sym = (mm**2 + nn**2) / 2
F0_cx = mm**2 * bb**2 + nn**2 * aa**2
F4_cx = s_sym * (aa**2 + bb**2)
F8_cx = nn**2 * bb**2 + mm**2 * aa**2
L_cx = F4_cx - 2 * mm * nn * aa * bb
U_cx = F4_cx + 2 * mm * nn * aa * bb

# гомогенизируем мои клетки: t = a/b, умножаем на b^2
def homog(poly):
    """подставить t -> a/b и умножить на b^2"""
    num = poly.numerator() if hasattr(poly, 'numerator') else poly
    F = Rab.fraction_field()
    sub = poly.subs({m: Rab.fraction_field()(mm), n: Rab.fraction_field()(nn),
                     t: Rab.fraction_field()(aa) / Rab.fraction_field()(bb)})
    return Rab.fraction_field()(sub) * Rab.fraction_field()(bb)**2

pairs_check = [("c0 <-> F0", c0, F0_cx), ("c4 <-> F4", c4, F4_cx), ("c8 <-> F8", c8, F8_cx),
               ("c2 <-> L", c2, L_cx), ("c6 <-> U", c6, U_cx)]
allmatch = True
for nm, mine, theirs in pairs_check:
    d = homog(mine) - Rab.fraction_field()(theirs)
    ok = (d == 0)
    allmatch = allmatch and ok
    say("  %-12s  b^2 * (моя клетка) - (форма Codex) = %s   -> %s" % (nm, d, "СОВПАЛО" if ok else "РАЗОШЛОСЬ"))
say()
say("  ИТОГ B: все пять неавтоматических клеток совпали: %s" % allmatch)
say("  (обозначения Codex L = c2, U = c6 — «красные» клетки; при t -> -t они меняются местами)")

# ---------------------------------------------------------------------------
# C. ПОЛНОТА ПАРАМЕТРИЗАЦИИ И ОТСУТСТВИЕ ЛИШНЕЙ МАСШТАБНОЙ СВОБОДЫ
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("C. Полна ли параметризация? Нет ли множителя, который портит редукцию?")
say("-" * 100)
say("  Обратное восстановление: n = (P-R)/2, m = (Q+S)/2, t = (P+R)/(Q+S).")
say("    n - (P-R)/2 = %s" % (n - (P - Rr) / 2))
say("    m - (Q+S)/2 = %s" % (m - (Q + S) / 2))
say("    t*(Q+S) - (P+R) = %s" % (t * (Q + S) - (P + Rr)))
say("  => любой рациональный магический квадрат конфигурации G1 получается при НЕКОТОРЫХ")
say("     рациональных m,n,t точно, без дополнительного множителя.")
say("  Масштабирование (m,n) -> λ(m,n) умножает КАЖДУЮ клетку на λ^2 (проверка):")
lam = R9.gens()[3]  # A
scaled = [cc.subs({m: lam * m, n: lam * n}) - lam**2 * cc for cc in cells]
say("    c_i(λm, λn, t) - λ^2 c_i(m,n,t) = %s  (все нули => класс квадратов не меняется)" % scaled)
say("  Значит переход к целым взаимно простым представителям (m:n) законен,")
say("  и проверка «клетка — квадрат» корректно определена на проективном классе.")

# ---------------------------------------------------------------------------
# D. ПЕРЕБОР ВСЕХ ТОЧЕК P^1(F_p)
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("D. Полный перебор P^1(F_17) для (m,n) = (13,8) — попытка найти пропущенное t")
say("-" * 100)

def forms_mod(mv, nv, p):
    """пять однородных форм как функции (a,b) -> элемент F_p; s = (m^2+n^2)/2 обращается в F_p"""
    Fp = GF(p)
    M, N = Fp(mv), Fp(nv)
    S = (M**2 + N**2) / Fp(2)
    def F0(a, b): return M**2 * b**2 + N**2 * a**2
    def F4(a, b): return S * (a**2 + b**2)
    def F8(a, b): return N**2 * b**2 + M**2 * a**2
    def Lf(a, b): return F4(a, b) - 2 * M * N * a * b
    def Uf(a, b): return F4(a, b) + 2 * M * N * a * b
    return [F0, F4, F8, Lf, Uf], Fp

def qr_set(p):
    return set(GF(p)(x)**2 for x in range(p))

def proj_points(p):
    Fp = GF(p)
    pts = [(Fp(x), Fp(1)) for x in range(p)] + [(Fp(1), Fp(0))]
    return pts

def enumerate_all(mv, nv, p, zero_ok=True, verbose=False):
    fs, Fp = forms_mod(mv, nv, p)
    QR = qr_set(p)
    QRnz = QR - {Fp(0)}
    good = []
    survive3 = []
    rows_out = []
    for (a, b) in proj_points(p):
        vals = [f(a, b) for f in fs]
        def issq(x):
            return (x in QR) if zero_ok else (x in QRnz)
        ok3 = all(issq(vals[i]) for i in (0, 1, 2))
        ok5 = all(issq(v) for v in vals)
        label = ("%s" % a) if b != 0 else "oo"
        rows_out.append((label, [int(v) for v in vals], ok3, ok5))
        if ok3:
            survive3.append(label)
        if ok5:
            good.append(label)
    return good, survive3, rows_out

mv, nv, p = 13, 8, 17
say("  Предварительно: s = (m^2+n^2)/2 = %s;  m^2+n^2 = %d" % (QQ(13**2 + 8**2) / QQ(2), 13**2 + 8**2))
say("  Редукция mod 17: s = %s, m^2 = %s, n^2 = %s, 2mn = %s"
    % (GF(17)(233) / GF(17)(2), GF(17)(169), GF(17)(64), GF(17)(208)))
say("  QR(17) (с нулём) = %s" % sorted(int(x) for x in qr_set(17)))
say()
say("  Все 18 проективных точек, формы (F0,F4,F8,L,U) mod 17:")
say("    %-4s %-28s %-8s %-8s" % ("t", "(F0,F4,F8,L,U)", "три?", "пять?"))
good, surv3, rows_out = enumerate_all(mv, nv, p, zero_ok=True)
for label, vals, ok3, ok5 in rows_out:
    mark = "  <-- прошла три" if ok3 else ""
    say("    %-4s %-28s %-8s %-8s%s" % (label, tuple(vals), ok3, ok5, mark))
say()
say("  Прошли первые три квадратности (F0,F4,F8): %s   (всего %d)" % (surv3, len(surv3)))
say("  Прошли все пять: %s" % (good if good else "НЕТ НИ ОДНОЙ"))

# Сверка с таблицей Codex
codex_tab = {"2": (0, 13, 9, 5, 4), "8": (15, 16, 0, 1, 14),
             "9": (15, 16, 0, 14, 1), "15": (0, 13, 9, 4, 5)}
say()
say("  Сверка с таблицей Codex (t : F0,F4,F8,L,U):")
mine_tab = {lab: tuple(v) for lab, v, ok3, ok5 in rows_out if ok3}
for k in sorted(codex_tab, key=lambda z: int(z)):
    say("    t=%-3s Codex %s ; мой расчёт %s ; совпадение: %s"
        % (k, codex_tab[k], mine_tab.get(k), mine_tab.get(k) == codex_tab[k]))
say("  Множество t, прошедших три квадратности, совпало с Codex: %s"
    % (set(mine_tab) == set(codex_tab)))
say("  На бесконечности: (F0,F4,F8,L,U) = %s" % (mine_tab.get("oo") or
    [tuple(v) for lab, v, _, _ in rows_out if lab == "oo"][0],))

# D2. Соглашения о нуле и про девять клеток
say()
say("  D2. Устойчивость вывода к соглашениям:")
for zo in (True, False):
    g, s3, _ = enumerate_all(mv, nv, p, zero_ok=zo)
    say("     ноль считается квадратом = %-5s : прошли три = %2d, прошли пять = %d"
        % (zo, len(s3), len(g)))
say("     (Codex берёт zero_ok=True — это КОНСЕРВАТИВНЫЙ выбор: он пропускает БОЛЬШЕ точек,")
say("      значит препятствие при таком соглашении сильнее, а не слабее.)")

# девять клеток: c1,c3,c5,c7 — тождественные квадраты, добавление их условий ничего не меняет
def enumerate_nine(mv, nv, p):
    Fp = GF(p)
    M, N = Fp(mv), Fp(nv)
    S = (M**2 + N**2) / Fp(2)
    QR = qr_set(p)
    cnt = 0
    for (a, b) in proj_points(p):
        c = [M**2 * b**2 + N**2 * a**2,
             (M * a + N * b)**2,
             S * (a**2 + b**2) - 2 * M * N * a * b,
             (M * a - N * b)**2,
             S * (a**2 + b**2),
             (M * b + N * a)**2,
             S * (a**2 + b**2) + 2 * M * N * a * b,
             (M * b - N * a)**2,
             N**2 * b**2 + M**2 * a**2]
        if all(x in QR for x in c):
            cnt += 1
    return cnt
say("     точек P^1(F_17), где ВСЕ ДЕВЯТЬ клеток — квадраты (вкл. 0): %d" % enumerate_nine(13, 8, 17))
say("     (ожидание: столько же, сколько для пяти — клетки 1,3,5,7 тождественно квадраты)")

# D3. Квадратичный твист: а если вся девятка умножена на неквадратную константу?
say()
say("  D3. Проверка на «а вдруг весь квадрат умножен на неквадратную константу?»")
say("      По §C такой свободы нет (m,n восстанавливаются точно), но проверим и этот случай:")
Fp = GF(17)
nonres = [x for x in range(1, 17) if GF(17)(x) not in qr_set(17)]
for cst in [1] + nonres[:4]:
    fs, _ = forms_mod(13, 8, 17)
    QR = qr_set(17)
    cnt = 0
    for (a, b) in proj_points(17):
        if all(GF(17)(cst) * f(a, b) in QR for f in fs):
            cnt += 1
    say("      множитель %2d (квадрат: %s): точек P^1(F_17) с пятью квадратами = %d"
        % (cst, GF(17)(cst) in qr_set(17), cnt))

# ---------------------------------------------------------------------------
# E. ОТ F_17 К Q_17 (Гензель) — независимая проверка
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("E. Есть ли точка над Q_17 (а не только над F_17)? Независимый рекурсивный подъём")
say("-" * 100)

def vp(x, p):
    if x == 0:
        return None
    num, den = x.numerator(), x.denominator()
    v = 0
    while num % p == 0:
        num //= p; v += 1
    while den % p == 0:
        den //= p; v -= 1
    return v

def is_sq_Qp(x, p):
    v = vp(x, p)
    if v is None:
        return True
    if v % 2:
        return False
    y = x / QQ(p)**v
    u = (y.numerator() * inverse_mod(y.denominator(), p)) % p
    return kronecker(u, p) == 1

def cells5_Q(mv, nv, a, b):
    S = QQ(mv**2 + nv**2) / 2
    return [QQ(mv**2 * b**2 + nv**2 * a**2),
            S * (a**2 + b**2),
            QQ(nv**2 * b**2 + mv**2 * a**2),
            S * (a**2 + b**2) - 2 * mv * nv * a * b,
            S * (a**2 + b**2) + 2 * mv * nv * a * b]

def solve_disk(mv, nv, p, r, k, KMAX, budget):
    """t in r + p^k Z_p; True/False/None"""
    if budget[0] <= 0:
        return None
    budget[0] -= 1
    vals = cells5_Q(mv, nv, QQ(r), QQ(1))
    undec = False
    for x in vals:
        v = vp(x, p)
        if v is not None and v < k:
            if not is_sq_Qp(x, p):
                return False
        else:
            undec = True
    if not undec:
        return True
    if k >= KMAX:
        return None
    res = False
    for z in range(p):
        out = solve_disk(mv, nv, p, r + z * p**k, k + 1, KMAX, budget)
        if out is True:
            return True
        if out is None:
            res = None
    return res

def solve_inf(mv, nv, p, KMAX, budget):
    """окрестность бесконечности: t = 1/u, u in pZ_p; формы (a,b)=(1,u)"""
    def rec(r, k):
        if budget[0] <= 0:
            return None
        budget[0] -= 1
        vals = cells5_Q(mv, nv, QQ(1), QQ(r))
        undec = False
        for x in vals:
            v = vp(x, p)
            if v is not None and v < k:
                if not is_sq_Qp(x, p):
                    return False
            else:
                undec = True
        if not undec:
            return True
        if k >= KMAX:
            return None
        res = False
        for z in range(p):
            out = rec(r + z * p**k, k + 1)
            if out is True:
                return True
            if out is None:
                res = None
        return res
    return rec(0, 1)

for (mv2, nv2, pp) in [(13, 8, 17), (7, 1, 17), (15, 8, 17), (5, 16, 41), (16, 5, 41)]:
    b1 = [200000]
    r1 = solve_disk(mv2, nv2, pp, 0, 0, 8, b1)
    b2 = [200000]
    r2 = solve_inf(mv2, nv2, pp, 8, b2)
    verdict = "ЕСТЬ точка" if (r1 is True or r2 is True) else (
        "НЕТ точки (доказано перебором дисков)" if (r1 is False and r2 is False)
        else "НЕ РЕШЕНО (бюджет/глубина)")
    say("  (m,n)=(%2d,%2d), p=%2d : Z_p -> %-5s, окрестность oo -> %-5s  ==> %s"
        % (mv2, nv2, pp, r1, r2, verdict))

# ---------------------------------------------------------------------------
# F. ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("F. Положительные контроли: машина обязана НЕ исключать пары с известной точкой")
say("-" * 100)

def cells9_Q(mv, nv, tv):
    tv = QQ(tv)
    S = QQ(mv**2 + nv**2) / 2
    return [mv**2 + nv**2 * tv**2,
            (mv * tv + nv)**2,
            S * (1 + tv**2) - 2 * mv * nv * tv,
            (mv * tv - nv)**2,
            S * (1 + tv**2),
            (mv + nv * tv)**2,
            S * (1 + tv**2) + 2 * mv * nv * tv,
            (mv - nv * tv)**2,
            nv**2 + mv**2 * tv**2]

def is_rational_square(x):
    x = QQ(x)
    return x >= 0 and QQ(x).is_square()

for (mv2, nv2, tv) in [(7, 1, 0), (15, 8, 1), (15, 8, -1), (13, 8, 0), (13, 8, 1)]:
    c = cells9_Q(mv2, nv2, tv)
    sq = [is_rational_square(x) for x in c]
    magic = all(sum(c[i] for i in L) == 3 * c[4] for L in lines)
    distinct = len(set(c)) == 9
    positive = all(x > 0 for x in c)
    say("  (m,n,t)=(%2d,%2d,%3s): клетки %s" % (mv2, nv2, tv, [str(x) for x in c]))
    say("      магический: %s; квадратов из 9: %d; все девять квадраты: %s; различны: %s; положительны: %s"
        % (magic, sum(sq), all(sq), distinct, positive))

say()
say("  Контроль локальной процедуры на этих парах (mod 17):")
for (mv2, nv2) in [(7, 1), (15, 8), (13, 8), (16, 5), (5, 16), (8, 13)]:
    g, s3, _ = enumerate_all(mv2, nv2, 17, zero_ok=True)
    say("     (m,n)=(%2d,%2d): точек P^1(F_17) с пятью квадратами = %2d  %s"
        % (mv2, nv2, len(g), "ИСКЛЮЧЕНА mod 17" if not g else ""))
say("  (важно: (7,1) и (15,8) имеют глобальную точку => их исключить нельзя; проверяем, что не исключены)")
say("  (важно: (8,13) — та же проективная пара, что (13,8), с переставленными m,n)")

# ---------------------------------------------------------------------------
# G. ВЫБОР ПРОСТОГО, ЗНАМЕНАТЕЛИ, ДРУГИЕ p
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("G. Корректность выбора p = 17 и перебор других простых")
say("-" * 100)
say("  Знаменатели, которые могли бы испортить редукцию: только 2 (в s = (m^2+n^2)/2).")
say("  17 | 2 ? %s;  17 | m=13 ? %s;  17 | n=8 ? %s;  17 | m^2+n^2=233 ? %s"
    % (2 % 17 == 0, 13 % 17 == 0, 8 % 17 == 0, 233 % 17 == 0))
say("  233 простое: %s" % is_prime(233))
say("  => 17 нечётно и не делит ни одного знаменателя: редукция форм корректна.")
say()
say("  Какие ещё простые дают препятствие для (13,8) (пять клеток, ноль разрешён):")
killers = []
for pp in prime_range(3, 200):
    if pp == 2:
        continue
    g, s3, _ = enumerate_all(13, 8, pp, zero_ok=True)
    if not g:
        killers.append(pp)
say("     %s" % killers)
say("  (если бы 17 было единственным и «случайным», это был бы повод для подозрения;")
say("   наличие/отсутствие других свидетелей — просто данные)")

say()
say("  Проверка опубликованного списка запрещённых λ = m/n mod 17:")
forb = []
for lam in range(17):
    g, _, _ = enumerate_all(lam, 1, 17, zero_ok=True)
    if not g:
        forb.append(lam)
g_inf, _, _ = enumerate_all(1, 0, 17, zero_ok=True)
say("     мой список λ (n=1): %s" % forb)
say("     Codex:              [2, 3, 6, 8, 9, 11, 14, 15]")
say("     совпадение: %s" % (forb == [2, 3, 6, 8, 9, 11, 14, 15]))
say("     случай n=0 mod 17 (λ = oo) исключён? %s (Codex: не входит в запрещённые)" % (not g_inf))
say("     (13,8): λ = 13/8 mod 17 = %s  -> в запрещённом списке: %s"
    % (GF(17)(13) / GF(17)(8), int(GF(17)(13) / GF(17)(8)) in forb))

# ---------------------------------------------------------------------------
# H. НЕЗАВИСИМЫЙ РАЦИОНАЛЬНЫЙ ПОИСК (контроль, НЕ доказательство)
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("H. Рациональный поиск t = a/b для (13,8) — КОНТРОЛЬ, не доказательство")
say("-" * 100)
best = (0, None)
N = 200
cnt_pts = 0
for b0 in range(1, N + 1):
    for a0 in range(-N, N + 1):
        if gcd(a0, b0) != 1:
            continue
        cnt_pts += 1
        tv = QQ(a0) / QQ(b0)
        c = cells9_Q(13, 8, tv)
        k = sum(1 for x in [c[0], c[2], c[4], c[6], c[8]] if is_rational_square(x))
        if k > best[0]:
            best = (k, tv)
say("  Перебрано %d несократимых t = a/b, |a|,|b| <= %d" % (cnt_pts, N))
say("  Максимум квадратов среди пяти неавтоматических клеток: %d при t = %s" % (best[0], best[1]))
say("  СТАТУС: отсутствие находки здесь НИЧЕГО не доказывает; это лишь согласуется с локальным запретом.")
say("  (для (13,8) локальный запрет mod 17 уже запрещает 3 квадрата из пяти в любой точке —")
say("   сверим: максимум, который допускает F_17 для (13,8):)")
fs, Fp = forms_mod(13, 8, 17)
QR = qr_set(17)
mx = 0
for (a, b) in proj_points(17):
    k = sum(1 for f in fs if f(a, b) in QR)
    mx = max(mx, k)
say("  Максимум квадратных форм из пяти по всем 18 точкам P^1(F_17): %d" % mx)

# ---------------------------------------------------------------------------
# I. САМОЕ ГЛАВНОЕ: ОПРЕДЕЛЁННОСТЬ ОТКАЗА В КАЖДОЙ ИЗ 18 ТОЧЕК
# ---------------------------------------------------------------------------
say()
say("-" * 100)
say("I. В каждой ли точке отказ ОПРЕДЕЛЁН, то есть засвидетельствован НЕНУЛЕВЫМ невычетом?")
say("-" * 100)
say("  Если в какой-то точке единственные «плохие» значения были бы нулями, вывод бы не прошёл:")
say("  ноль mod p ничего не говорит о квадратности в Q_p. Проверяем явно.")
fs, Fp17 = forms_mod(13, 8, 17)
QR17 = qr_set(17)
fnames = ["F0", "F4", "F8", "L", "U"]
undetermined = []
for (a, b) in proj_points(17):
    vals = [f(a, b) for f in fs]
    witness = [(fnames[i], int(vals[i])) for i in range(5)
               if vals[i] != Fp17(0) and vals[i] not in QR17]
    zeros = [fnames[i] for i in range(5) if vals[i] == Fp17(0)]
    label = ("t=%s" % a) if b != 0 else "t=oo"
    if not witness:
        undetermined.append(label)
    say("    %-6s нулевые клетки: %-10s ; ненулевые невычеты-свидетели: %s"
        % (label, zeros if zeros else "нет", witness if witness else "НЕТ — ДЫРА!"))
say()
say("  Точек без определённого свидетеля: %s" % (undetermined if undetermined else "НЕТ — препятствие определено везде"))
say("  => вывод НЕ зависит от того, считать ли ноль квадратом: в каждой из 18 точек")
say("     есть клетка с ненулевым значением, которое не является квадратом в F_17,")
say("     а значит (по 17-целостности форм) не является квадратом и в Q_17.")

say()
say("  I2. Насколько сильно это сверх нужного: сколько клеток из пяти минимум ломается")
say("      в каждой точке (чем больше, тем труднее испортить вывод одной опечаткой):")
cnts = []
for (a, b) in proj_points(17):
    vals = [f(a, b) for f in fs]
    k = sum(1 for v in vals if v != Fp17(0) and v not in QR17)
    cnts.append(k)
say("      распределение числа ненулевых невычетов по 18 точкам: %s"
    % {k: cnts.count(k) for k in sorted(set(cnts))})
say("      минимум по точкам: %d (если бы был 0 — препятствия не было бы)" % min(cnts))

say()
say("  I3. Устойчивость к опечатке в одной клетке: если ВЫЧЕРКНУТЬ любую одну из пяти форм,")
say("      сохраняется ли пустота? (проверка того, что вывод не держится на одной клетке)")
for drop in range(5):
    idx = [i for i in range(5) if i != drop]
    cnt = 0
    for (a, b) in proj_points(17):
        if all(fs[i](a, b) in QR17 for i in idx):
            cnt += 1
    say("      без %-2s : точек с четырьмя квадратами = %d" % (fnames[drop], cnt))
say("      (ожидание: без L или без U появляются точки — значит «красные» клетки")
say("       действительно несут препятствие, а F0,F4,F8 сами по себе его не дают)")

# ---------------------------------------------------------------------------
say()
say("=" * 100)
say("КОНЕЦ. Итоговые статусы расставляются в отчёте, а не в этом выводе.")
say("=" * 100)

with open("/home/kep/magicKube/bridge/local_13_8_p17_скептик.log", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
