# -*- coding: utf-8 -*-
# check_kummer_11_4_модель.sage
#
# НЕЗАВИСИМАЯ ПОПЫТКА СЛОМАТЬ заявление Codex: C_{11,4}(Q) = пусто.
# Угол атаки: МОДЕЛЬ И ТОЖДЕСТВА.  Всё строится с нуля, чужие числа не копируются;
# числа Codex используются ТОЛЬКО в финальном блоке сверки (раздел J).
#
#  A. Символьные тождества в Q(m,n)[t] — общие m,n, без подстановки чисел.
#  B. Специализация (11,4): неособость, различие корней, полное 2-кручение.
#  C. Собственная реализация delta (полный 2-спуск) + контроль гомоморфности.
#  D. КОНТРОЛЬНЫЙ ЭКСПЕРИМЕНТ: пара (15,8), у которой на C ЕСТЬ явная точка t=1.
#     Это эмпирически фиксирует, какой корень обязан быть e1 (класс 1).
#  E. Точки на E, образ delta, требуемый класс.
#  F. Ранг: несколько независимых источников верхней границы.
#  G. Бесконечность — в проективных координатах, все 8 ветвей.
#  H. Перестановки корней: где вывод ломается и почему.
#  I. Локальный образ delta_v и вопрос «лежит ли (1,s,s) в группе Селмера».
#  J. Сверка с числами Codex.

import sys, itertools, time

def hdr(t):
    print("\n" + "=" * 78); print(t); print("=" * 78); sys.stdout.flush()

VERD = {}
def rec(key, ok, note=""):
    VERD[key] = (bool(ok), note)
    print("   [%s] %s%s" % ("OK  " if ok else "ПРОВАЛ", key, ("  -- " + note) if note else ""))
    sys.stdout.flush()

def sqcls(q):
    """Представитель квадратного класса ненулевого q из Q*: бесквадратное целое."""
    q = QQ(q)
    if q == 0:
        raise ValueError("нулевой квадратный класс")
    return ZZ(q.numerator() * q.denominator()).squarefree_part()

# ============================================================ A
hdr("A. СИМВОЛЬНЫЕ ТОЖДЕСТВА В Q(m,n)[t] — БЕЗ ПОДСТАНОВКИ ЧИСЕЛ")

Rmn = PolynomialRing(QQ, ['m', 'n'])
Fmn = Rmn.fraction_field()
m, n = Fmn.gens()
Rt.<t> = PolynomialRing(Fmn)

s_sym = (m**2 + n**2) / 2
F0 = m**2 + n**2 * t**2
F4 = s_sym * (1 + t**2)
F8 = n**2 + m**2 * t**2
b_sym = s_sym * m**2 * n**2

e1 = -b_sym
e2 = -s_sym * m**4
e3 = -s_sym * n**4
X = b_sym * t**2

print("   s   = %s" % s_sym)
print("   b   = %s" % b_sym)
print("   e1  = %s" % e1)
print("   e2  = %s" % e2)
print("   e3  = %s" % e3)

id1 = (X - e1) - (m * n)**2 * F4
id2 = (X - e2) - s_sym * m**2 * F0
id3 = (X - e3) - s_sym * n**2 * F8
rec("A1  X-e1 == (mn)^2 * F4  тождественно в Q(m,n)[t]", id1 == 0, "остаток %s" % id1)
rec("A2  X-e2 == s*m^2 * F0  тождественно", id2 == 0, "остаток %s" % id2)
rec("A3  X-e3 == s*n^2 * F8  тождественно", id3 == 0, "остаток %s" % id3)

# уравнение E выполняется: V^2 = (X-e1)(X-e2)(X-e3), V = b*u0*u4*u8
lhs = (b_sym)**2 * F0 * F4 * F8
rhs = (X - e1) * (X - e2) * (X - e3)
rec("A4  (b u0u4u8)^2 == (X-e1)(X-e2)(X-e3)", lhs == rhs, "разность %s" % (lhs - rhs))

# КРИТИЧЕСКОЕ: какой множитель у каждого корня, ЯВНО.
# X-e1 = (mn)^2 * u4^2  -> квадрат.       класс 1
# X-e2 = s * m^2 * u0^2 -> s * квадрат.   класс [s]
# X-e3 = s * n^2 * u8^2 -> s * квадрат.   класс [s]
cof1 = (X - e1) / F4          # должно быть (mn)^2, константа по t
cof2 = (X - e2) / F0
cof3 = (X - e3) / F8
ok_c = (cof1 == (m * n)**2) and (cof2 == s_sym * m**2) and (cof3 == s_sym * n**2)
rec("A5  коэффициенты (X-ei)/F_j: (mn)^2, s m^2, s n^2", ok_c,
    "cof1=%s cof2=%s cof3=%s" % (cof1, cof2, cof3))
print("   ВЫВОД A: класс 1 принадлежит ИМЕННО корню e1=-b=-s m^2 n^2, потому что")
print("            X-e1 = b(t^2+1) = (mn)^2 * s(1+t^2) = (mn)^2 u4^2.")
print("            e2 и e3 оба дают класс [s]; их взаимный порядок НЕ важен.")

# Проверим, что нули F_i не могут быть рациональными (т.е. u_i != 0 на C(Q))
disc_checks = []
for nm, F in [("F0", F0), ("F4", F4), ("F8", F8)]:
    # F как квадратичный по t: дискриминант
    c = F.coefficients(sparse=False)
    disc = c[1]**2 - 4 * c[0] * c[2] if len(c) == 3 else None
    disc_checks.append((nm, disc))
    print("   %s: старший=%s, свободный=%s, дискриминант по t = %s" % (nm, c[2], c[0], disc))
print("   Все дискриминанты = -4*(положительное) < 0 при m,n>0 => корни мнимые =>")
print("   при вещественном t все F_i > 0, в частности u_i != 0. (см. B6)")

# ============================================================ B
hdr("B. СПЕЦИАЛИЗАЦИЯ (m,n) = (11,4)")

M, N = 11, 4
rec("B0  gcd(m,n)=1", gcd(M, N) == 1, "gcd=%s" % gcd(M, N))

s = QQ(M**2 + N**2) / 2
b = s * M**2 * N**2
E1 = -b
E2 = -s * M**4
E3 = -s * N**4
print("   s = %s      (= %s)" % (s, s.factor() if s.denominator() == 1 else "137/2"))
print("   b = %s" % b)
print("   e1 = %s" % E1)
print("   e2 = %s" % E2)
print("   e3 = %s" % E3)
rec("B1  b целое", b.denominator() == 1, "b=%s" % b)
rec("B2  корни попарно различны", len(set([E1, E2, E3])) == 3)

# кривая: V^2 = (X-e1)(X-e2)(X-e3) = X^3 + a2 X^2 + a4 X + a6
Rx.<Xv> = PolynomialRing(QQ)
cub = (Xv - E1) * (Xv - E2) * (Xv - E3)
a2 = cub[2]; a4 = cub[1]; a6 = cub[0]
E = EllipticCurve([0, a2, 0, a4, a6])
print("   E: %s" % E)
print("   дискриминант E = %s" % E.discriminant())
rec("B3  E неособа (disc != 0)", E.discriminant() != 0)
rec("B4  кубика E совпадает с (X-e1)(X-e2)(X-e3)",
    Xv**3 + a2 * Xv**2 + a4 * Xv + a6 == cub)

T = E.torsion_subgroup()
tors_pts = sorted([P for P in E.torsion_points()], key=lambda P: str(P))
two_tors = [P for P in E.torsion_points() if P.order() == 2]
roots_from_curve = sorted([P[0] for P in two_tors])
rec("B5  E[2] полностью рационально (3 точки порядка 2)", len(two_tors) == 3,
    "структура кручения %s" % T.invariants().__str__())
rec("B5b корни 2-кручения = {e1,e2,e3}", set(roots_from_curve) == set([E1, E2, E3]),
    "из кривой %s" % roots_from_curve)
print("   E(Q)_tors = %s, порядок %s" % (T.invariants(), T.order()))

# положительность F_i при вещественном t
F0n = M**2 + N**2 * Xv**2
F4n = s * (1 + Xv**2)
F8n = N**2 + M**2 * Xv**2
allpos = all(len(f.roots(RR)) == 0 and f(0) > 0 for f in [F0n, F4n, F8n])
rec("B6  F0,F4,F8 не имеют вещественных корней и положительны", allpos)

req_class = (1, sqcls(s), sqcls(s))
print("   ТРЕБУЕМЫЙ КЛАСС (мой расчёт): %s   [ s = %s, sqfree(s) = %s ]"
      % (req_class.__str__(), s, sqcls(s)))

# ============================================================ C
hdr("C. СОБСТВЕННАЯ РЕАЛИЗАЦИЯ delta И КОНТРОЛЬ ЕЁ СВОЙСТВ")

ROOTS = [E1, E2, E3]

def delta(P, roots=ROOTS):
    """Полный 2-спуск: delta(P) = [X-e1, X-e2, X-e3] в (Q*/Q*^2)^3,
       с заменой нулевой координаты для 2-кручения и (1,1,1) для O."""
    if P.is_zero():
        return (1, 1, 1)
    x = P[0]
    out = []
    for i in range(3):
        v = x - roots[i]
        if v == 0:
            j, k = [z for z in range(3) if z != i]
            v = (roots[i] - roots[j]) * (roots[i] - roots[k])
        out.append(sqcls(v))
    return tuple(out)

# C1: произведение координат всегда квадрат
def prod_is_square(d):
    return ZZ(d[0] * d[1] * d[2]).squarefree_part() == 1

# C2: гомоморфность — численная проверка на многих точках
def mulcls(d, dd):
    return tuple(sqcls(QQ(d[i]) * dd[i]) for i in range(3))

pts_test = []
for P in E.torsion_points():
    pts_test.append(P)
srch = E.point_search(11)
for P in srch:
    pts_test.append(P); pts_test.append(-P)
# ещё точек: кратные и суммы
extra = []
for P in pts_test[:]:
    for Q in pts_test[:]:
        R = P + Q
        if R not in extra:
            extra.append(R)
pts_test = list(set(pts_test + extra))
print("   точек E(Q) для теста гомоморфности: %s" % len(pts_test))

homok = True
bad = None
cnt = 0
for P in pts_test:
    for Q in pts_test:
        if mulcls(delta(P), delta(Q)) != delta(P + Q):
            homok = False; bad = (P, Q); break
        cnt += 1
    if not homok:
        break
rec("C1  delta — гомоморфизм (проверено %d пар)" % cnt, homok, "контрпример %s" % (bad,))
rec("C2  произведение координат delta — всегда квадрат",
    all(prod_is_square(delta(P)) for P in pts_test))
rec("C3  delta(O) = (1,1,1)", delta(E(0)) == (1, 1, 1))
# ядро: delta(2P) = (1,1,1)
rec("C4  delta(2P) = (1,1,1) для всех тестовых P",
    all(delta(2 * P) == (1, 1, 1) for P in pts_test))

# ============================================================ D
hdr("D. КОНТРОЛЬНЫЙ ЭКСПЕРИМЕНТ: ПАРА (15,8), У КОТОРОЙ НА C ЕСТЬ ЯВНАЯ ТОЧКА")
print("   Если порядок корней перепутан, требуемый класс окажется перестановкой")
print("   и вывод перевернётся.  Нужен случай, где точка C ИЗВЕСТНА, а [s] != 1.")

def build(mm, nn):
    ss = QQ(mm**2 + nn**2) / 2
    bb = ss * mm**2 * nn**2
    rr = [-bb, -ss * mm**4, -ss * nn**4]
    cc = (Xv - rr[0]) * (Xv - rr[1]) * (Xv - rr[2])
    EE = EllipticCurve([0, cc[2], 0, cc[1], cc[0]])
    return ss, bb, rr, EE

def C_point(mm, nn, tt):
    """Проверяет, лежит ли t=tt на C_{mm,nn}; возвращает (u0,u4,u8) или None."""
    ss = QQ(mm**2 + nn**2) / 2
    vals = [mm**2 + nn**2 * tt**2, ss * (1 + tt**2), nn**2 + mm**2 * tt**2]
    us = []
    for v in vals:
        v = QQ(v)
        if not v.is_square():
            return None
        us.append(v.sqrt())
    return tuple(us)

ctrl = []
for (mm, nn) in [(15, 8), (20, 21), (7, 24), (119, 120), (44, 117)]:
    if gcd(mm, nn) != 1:
        continue
    u = C_point(mm, nn, QQ(1))
    if u is None:
        continue
    ss, bb, rr, EE = build(mm, nn)
    Xp = bb * QQ(1)**2
    Vp = bb * u[0] * u[1] * u[2]
    on = (Vp**2 == (Xp - rr[0]) * (Xp - rr[1]) * (Xp - rr[2]))
    P = EE(Xp, Vp)
    d = delta(P, rr)
    want = (1, sqcls(ss), sqcls(ss))
    ctrl.append((mm, nn, ss, d, want, on, d == want))
    print("   (m,n)=(%d,%d)  t=1  u=(%s,%s,%s)  s=%s  [s]=%s" % (mm, nn, u[0], u[1], u[2], ss, sqcls(ss)))
    print("        точка на E: %s ;  лежит на E: %s" % (P, on))
    print("        delta(P) = %s      требуемый (1,[s],[s]) = %s   -> %s"
          % (d.__str__(), want.__str__(), "СОВПАЛО" if d == want else "НЕ СОВПАЛО"))

rec("D1  контрольные пары дают delta = (1,[s],[s]) точно в этом порядке",
    len(ctrl) > 0 and all(c[6] for c in ctrl),
    "проверено пар: %d" % len(ctrl))
rec("D2  у контрольных пар [s] != 1, значит тест РАЗЛИЧАЕТ позицию 1 и позиции 2,3",
    all(sqcls(c[2]) != 1 for c in ctrl))

# дополнительный контроль: широкий поиск точек C при t=p/q
hdr("D'. ПРЯМОЙ ПОИСК РАЦИОНАЛЬНЫХ ТОЧЕК НА C_{11,4} (попытка опровергнуть в лоб)")
found = []
HT = 400
t0 = time.time()
for q in range(1, HT + 1):
    for p in range(-HT, HT + 1):
        if gcd(abs(p), q) != 1:
            continue
        tt = QQ(p) / q
        # быстрый предварительный фильтр: F0 квадрат
        num0 = M**2 * q**2 + N**2 * p**2
        if not ZZ(num0).is_square():
            continue
        num8 = N**2 * q**2 + M**2 * p**2
        if not ZZ(num8).is_square():
            continue
        # F4 = s(1+t^2) = (137/2)(p^2+q^2)/q^2 -> квадрат <=> 2*137*(p^2+q^2) квадрат
        if not ZZ(2 * 137 * (p**2 + q**2)).is_square():
            continue
        found.append(tt)
print("   перебор |p|,q <= %d, время %.1f с" % (HT, time.time() - t0))
rec("D3  прямой поиск точек C_{11,4} ничего не нашёл (это НЕ доказательство)",
    len(found) == 0, "найдено: %s" % found)

# ============================================================ E
hdr("E. ТОЧКИ E(Q), ОБРАЗ delta, ТРЕБУЕМЫЙ КЛАСС")

t0 = time.time()
gens = None
Emin = E.minimal_model()
iso = Emin.isomorphism_to(E)
print("   минимальная модель Emin: %s" % Emin)
try:
    gmin = Emin.gens()
    gens = [iso(g) for g in gmin]
except Exception as ex:
    print("   Emin.gens() не отработал: %s" % ex)
    try:
        gens = E.gens()
    except Exception as ex2:
        print("   E.gens() не отработал: %s" % ex2)
print("   образующие (перенесённые на E) = %s   (%.1f с)" % (gens, time.time() - t0))

# собственный поиск точек, независимо от gens()
pool = set()
for P in E.torsion_points():
    pool.add(P)
t0 = time.time()
for h in [12, 16, 20]:
    try:
        for P in Emin.point_search(h):
            pool.add(iso(P)); pool.add(-iso(P))
    except Exception as ex:
        print("   point_search(%s) сбой: %s" % (h, ex))
print("   углублённый поиск точек до высоты 20 на Emin: %.1f с, точек %s" % (time.time() - t0, len(pool)))
if gens:
    for g in gens:
        pool.add(g)
        for k in range(-4, 5):
            pool.add(k * g)
# замыкаем по сложению один раз
pool2 = set(pool)
for P in list(pool):
    for Q in list(pool):
        pool2.add(P + Q)
pool = pool2
print("   всего точек в пуле: %s" % len(pool))

img = set()
for P in pool:
    img.add(delta(P))
img = sorted(img)
print("   различных классов delta от найденных точек: %s" % len(img))
for d in img:
    print("      %s" % (d.__str__()))

rec("E1  найдено ровно 8 различных классов", len(img) == 8, "найдено %s" % len(img))
rec("E2  требуемый класс %s ОТСУТСТВУЕТ среди найденных" % (req_class.__str__()),
    req_class not in img)

# проверим, что образ — подгруппа
def is_subgroup(S):
    S = set(S)
    for a in S:
        for c in S:
            if mulcls(a, c) not in S:
                return False
    return True
rec("E3  найденный образ замкнут по умножению (подгруппа)", is_subgroup(img))

# перестановки требуемого класса — присутствуют ли?
perms_req = set()
for sg in itertools.permutations(range(3)):
    perms_req.add(tuple(req_class[i] for i in sg))
print("   перестановки требуемого класса: %s" % sorted(perms_req).__str__())
present_perm = [d for d in img if d in perms_req]
print("   ИЗ НИХ ПРИСУТСТВУЮТ В ОБРАЗЕ: %s" % present_perm.__str__())
rec("E4  ОПАСНОСТЬ: некоторая перестановка требуемого класса В ОБРАЗЕ ЕСТЬ",
    len(present_perm) > 0,
    "это значит, что путаница в порядке корней ПЕРЕВЕРНУЛА БЫ ВЫВОД")

# ============================================================ F
hdr("F. РАНГ E(Q): НЕСКОЛЬКО НЕЗАВИСИМЫХ ИСТОЧНИКОВ ВЕРХНЕЙ ГРАНИЦЫ")

print("   ЛОГИКА: |E(Q)/2E(Q)| = 2^(r+2) при полном 2-кручении.")
print("   8 различных классов доказывают ТОЛЬКО r >= 1.  Полнота образа требует r <= 1.")

res_rank = {}
try:
    lo, hi = E.rank_bounds()
    res_rank['eclib rank_bounds'] = (lo, hi)
    print("   eclib rank_bounds: [%s, %s]" % (lo, hi))
except Exception as ex:
    print("   rank_bounds не отработал: %s" % ex)

try:
    sr = E.selmer_rank()
    res_rank['selmer_rank'] = sr
    print("   2-Selmer rank (eclib) = %s  =>  r <= %s - 2 = %s" % (sr, sr, sr - 2))
except Exception as ex:
    print("   selmer_rank не отработал: %s" % ex)

try:
    r_pari = E.rank(algorithm='pari', only_use_mwrank=False)
    res_rank['pari rank'] = r_pari
    print("   PARI ellrank -> rank = %s" % r_pari)
except Exception as ex:
    print("   PARI rank не отработал: %s" % ex)

try:
    ar = E.analytic_rank()
    res_rank['analytic_rank'] = ar
    print("   аналитический ранг = %s" % ar)
    print("   знак функционального уравнения = %s" % E.root_number())
    print("   кондуктор = %s" % E.conductor().factor())
except Exception as ex:
    print("   analytic_rank не отработал: %s" % ex)

ok_rank_1_unconditional = False
note = ""
if 'eclib rank_bounds' in res_rank:
    lo, hi = res_rank['eclib rank_bounds']
    if hi == 1:
        ok_rank_1_unconditional = True
        note = "eclib даёт верхнюю границу 1"
    else:
        note = "eclib верхняя граница = %s > 1" % hi
rec("F1  БЕЗУСЛОВНАЯ верхняя граница ранга = 1 получена 2-спуском", ok_rank_1_unconditional, note)

print("\n   --- F2: сырой вывод PARI ellrank (что именно доказано) ---")
try:
    Emin = E.minimal_model()
    print("   минимальная модель: %s" % Emin)
    pe = pari(Emin).ellrank()
    print("   PARI ellrank(Emin) = %s" % pe)
    print("   формат PARI: [нижняя граница, верхняя граница, s, точки]")
    lo_p = ZZ(pe[0]); hi_p = ZZ(pe[1])
    res_rank['pari bounds'] = (lo_p, hi_p)
    rec("F2  PARI ellrank доказывает rank = 1 (lo == hi == 1)", lo_p == 1 and hi_p == 1,
        "PARI: [%s, %s]" % (lo_p, hi_p))
except Exception as ex:
    print("   PARI ellrank не отработал: %s" % ex)
    rec("F2  PARI ellrank доказывает rank = 1", False, str(ex))

print("\n   --- F3: строгая дорожка Гросс-Загир + Колывагин ---")
print("   Знак функционального уравнения -1 => ord_{s=1} L(E,s) НЕЧЁТЕН (безусловно,")
print("   по модулярности и функциональному уравнению).  Достаточно показать L'(E,1) != 0.")
try:
    Emin = E.minimal_model()
    Lser = Emin.lseries()
    val, err = Lser.deriv_at1(4000)
    print("   L'(E,1) ~ %s   с оценкой ошибки %s" % (val, err))
    rig = (abs(val) > 2 * err) and (err > 0)
    rec("F3  L'(E,1) != 0 с запасом по оценке ошибки", rig,
        "|L'| = %s, error = %s, отношение %s" % (val, err, (abs(val) / err) if err > 0 else "inf"))
    if rig:
        print("   => ord = 1 => (Гросс-Загир + Колывагин) rank E(Q) = 1 и Sha конечна.")
        print("   ЭТО БЕЗУСЛОВНАЯ ТЕОРЕМА, если численная оценка ошибки корректна.")
except Exception as ex:
    print("   deriv_at1 не отработал: %s" % ex)
    rec("F3  L'(E,1) != 0 с запасом", False, str(ex))

print("\n   --- F4: что даёт 2-спуск сам по себе ---")
if 'selmer_rank' in res_rank:
    sr = res_rank['selmer_rank']
    print("   dim_F2 S^2(E/Q) = %s = r + 2 + dim Sha[2].  При r=1 => dim Sha[2] = %s." % (sr, sr - 3))
    print("   Индекс образа E(Q)/2E(Q) в Селмере при r=1 равен 2^%s = %s." % (sr - 3, 2**(sr - 3)))
    print("   ЗНАЧИТ: чистый 2-спуск НЕ закрывает (11,4); нужен внешний аргумент про ранг.")

# ============================================================ G
hdr("G. БЕСКОНЕЧНОСТЬ: ПРОЕКТИВНО, ВСЕ ВЕТВИ")

print("   C аффинно: u0^2=F0, u4^2=F4, u8^2=F8.  Гладкая проективная модель добавляет")
print("   только точки над t=inf (аффинная часть при конечном t уже гладкая: см. G1).")

# G1: гладкость аффинной части
Rq = PolynomialRing(QQ, ['tt', 'x0', 'x4', 'x8'])
tt, x0, x4, x8 = Rq.gens()
g0 = x0**2 - (M**2 + N**2 * tt**2)
g4 = x4**2 - s * (1 + tt**2)
g8 = x8**2 - (N**2 + M**2 * tt**2)
Jac = matrix(Rq, [[g.derivative(v) for v in Rq.gens()] for g in [g0, g4, g8]])
# особая точка требует ранг < 3, т.е. все 3x3 миноры = 0 вместе с g_i.
I = Rq.ideal([g0, g4, g8] + list(Jac.minors(3)))
dimI = I.dimension()
print("   идеал (уравнения + все 3x3 миноры якобиана): размерность = %s" % dimI)
rec("G1  аффинная C гладкая над Q~ (особое множество пусто, dim = -1)", dimI == -1)

# G2: замена z = 1/t, U_i = u_i/t
Rz.<z> = PolynomialRing(QQ)
G0 = N**2 + M**2 * z**2      # (u0/t)^2
G4 = s * (1 + z**2)          # (u4/t)^2
G8 = M**2 + N**2 * z**2      # (u8/t)^2
print("   U0^2 = (u0/t)^2 = m^2/t^2 + n^2 = %s   -> при z=0: %s" % (G0, G0(0)))
print("   U4^2 = (u4/t)^2 = s(1/t^2 + 1)   = %s   -> при z=0: %s" % (G4, G4(0)))
print("   U8^2 = (u8/t)^2 = n^2/t^2 + m^2 = %s   -> при z=0: %s" % (G8, G8(0)))
# символьная проверка подстановки
u0sq = M**2 + N**2 * t**2  # в Rt, но с числами
chk0 = (M**2 + N**2 * (1 / z)**2) * z**2 - G0
chk4 = (s * (1 + (1 / z)**2)) * z**2 - G4
chk8 = (N**2 + M**2 * (1 / z)**2) * z**2 - G8
FFz = Rz.fraction_field()
rec("G2  подстановка t=1/z, U_i=u_i/t корректна",
    FFz(chk0) == 0 and FFz(chk4) == 0 and FFz(chk8) == 0)

print("   Над z=0 система распадается: U0^2=n^2=%s, U4^2=s=%s, U8^2=m^2=%s" % (N**2, s, M**2))
print("   U0 = +-%s (рационально), U8 = +-%s (рационально), U4 = +-sqrt(%s)" % (N, M, s))
rec("G3  s НЕ квадрат в Q", not QQ(s).is_square(), "s = %s" % s)
print("   => все 8 точек над t=inf определены над Q(sqrt(s)) = Q(sqrt(%s)) и разбиты" % sqcls(s))
print("      на 4 пары сопряжённых; рациональных среди них НЕТ.")

# G4: проверим, что над z=0 ровно 8 точек и они неразветвлены
print("   Число точек над z=0: три квадратики G0,G4,G8 при z=0 дают ненулевые значения")
print("   %s, %s, %s -> ни одна не ветвится, точек 2*2*2 = 8." % (G0(0), G4(0), G8(0)))
rec("G4  над t=inf ровно 8 неразветвлённых точек", G0(0) != 0 and G4(0) != 0 and G8(0) != 0)

# G5: род
print("   Ветвление: F0,F4,F8 имеют по 2 различных корня, все 6 различны?")
allroots = []
for f in [F0n, F4n, F8n]:
    allroots += [r for r, _ in f.roots(QQbar)]
print("   6 точек ветвления: %s" % [CC(r) for r in allroots])
distinct6 = len(set(allroots)) == 6
rec("G5  6 точек ветвления попарно различны", distinct6)
print("   Риман-Гурвиц: 2g-2 = 8*(-2) + 6*4 = 8 => g = 5.  (степень накрытия 8)")

# G6: независимость F0,F4,F8 по модулю квадратов в Q(t)* (иначе C приводима)
Rt2.<T2v> = PolynomialRing(QQ)
Fs = {1: M**2 + N**2 * T2v**2, 2: s * (1 + T2v**2), 4: N**2 + M**2 * T2v**2}
indep = True
for mask in range(1, 8):
    prod_f = Rt2(1)
    for k in [1, 2, 4]:
        if mask & k:
            prod_f *= Fs[k]
    if prod_f.is_square():
        indep = False
        print("   ПОДОЗРЕНИЕ: произведение маски %s — квадрат в Q(t)" % mask)
rec("G6  F0,F4,F8 независимы mod квадратов => C неприводима, степень 8", indep)

# ============================================================ H
hdr("H. ПЕРЕСТАНОВКИ КОРНЕЙ")

print("   delta зависит от упорядочения (e1,e2,e3).  При перестановке sg образ и")
print("   требуемый класс переставляются ОДНОЙ И ТОЙ ЖЕ sg -> вывод инвариантен,")
print("   ЕСЛИ и только если обе величины считаны в ОДНОМ порядке.")
inv_ok = True
for sg in itertools.permutations(range(3)):
    rr = [ROOTS[i] for i in sg]
    img_s = set(delta(P, rr) for P in pool)
    # требуемый класс в том же порядке
    base = {0: 1, 1: sqcls(s), 2: sqcls(s)}
    req_s = tuple(base[i] for i in sg)
    inside = req_s in img_s
    print("   sg=%s  корни=%s  требуемый=%s  в образе: %s"
          % (sg.__str__(), [str(x) for x in rr], req_s.__str__(), inside))
    if inside:
        inv_ok = False
rec("H1  ни при какой перестановке (согласованной!) требуемый класс не попадает в образ", inv_ok)

print("   А теперь НЕСОГЛАСОВАННЫЙ случай — ровно та ошибка, которую ищем:")
for sg in itertools.permutations(range(3)):
    wrong = tuple(req_class[i] for i in sg)
    if wrong in img:
        print("   *** если требуемый класс ошибочно записать как %s (перестановка sg=%s"
              % (wrong.__str__(), sg.__str__()))
        print("       при НЕпереставленных корнях), он В ОБРАЗЕ ЕСТЬ и вывод рушится.")

# ============================================================ I
hdr("I. ЛОКАЛЬНЫЙ ОБРАЗ delta_v И ГРУППА СЕЛМЕРА")

def loc_class(q, p):
    """Канонический представитель квадратного класса q в Q_p*/Q_p*^2."""
    q = QQ(q)
    a = q.valuation(p)
    u = q / p**a
    if p == 2:
        un = ZZ(u.numerator()); ud = ZZ(u.denominator())
        uu = (un * ud) % 8          # ud нечётно, ud^2=1 mod 8, so un*ud ~ un/ud
        return (a % 2, uu % 8)
    else:
        un = ZZ(u.numerator()) % p; ud = ZZ(u.denominator()) % p
        val = (un * pow(ud, p - 2, p)) % p
        return (a % 2, 1 if kronecker(val, p) == 1 else -1)

def loc_class_R(q):
    return 1 if QQ(q) > 0 else -1

def loc_delta(x, p):
    return tuple(loc_class(x - ROOTS[i], p) for i in range(3))

def loc_size(p):
    # |E(Q_p)/2E(Q_p)| = |E(Q_p)[2]| / |2|_p ; полное 2-кручение => 4 / |2|_p
    return 4 if p != 2 else 8

BADP = sorted(set([2] + [q for q, _ in ZZ(E.discriminant().numerator()).factor()]
                      + [q for q, _ in ZZ(2 * M * N * (M**2 - N**2) * (M**2 + N**2)).factor()]))
print("   плохие простые (мой расчёт): %s" % BADP)

sel_report = {}
set_random_seed(1)
for p in BADP:
    target = loc_size(p)
    seen = {}
    # образы кручения и O
    seen[(loc_class(1, p),) * 3 if False else tuple([loc_class(1, p)] * 3)] = "O"
    for i in range(3):
        x = ROOTS[i]
        cls = []
        for j in range(3):
            v = x - ROOTS[j]
            if v == 0:
                a, c = [z for z in range(3) if z != j]
                v = (ROOTS[j] - ROOTS[a]) * (ROOTS[j] - ROOTS[c])
            cls.append(loc_class(v, p))
        seen[tuple(cls)] = "T%d" % (i + 1)
    # случайный поиск X из Q_p с f(X) квадратом
    tries = 0
    while len(seen) < target and tries < 400000:
        tries += 1
        e = ZZ.random_element(-3, 4)
        a = ZZ.random_element(1, p**4 + 1)
        x = QQ(a) * p**e
        f = (x - ROOTS[0]) * (x - ROOTS[1]) * (x - ROOTS[2])
        if f == 0:
            continue
        # f — квадрат в Q_p ?
        a_ = f.valuation(p)
        if a_ % 2 != 0:
            continue
        u = f / p**a_
        if p == 2:
            un = ZZ(u.numerator()) * ZZ(u.denominator())
            if un % 8 != 1:
                continue
        else:
            un = ZZ(u.numerator()); ud = ZZ(u.denominator())
            if kronecker((un * ud) % p, p) != 1:
                continue
        seen[loc_delta(x, p)] = "x=%s" % x
    got = len(seen)
    req_loc = tuple(loc_class(QQ(req_class[i]), p) for i in range(3))
    inside = req_loc in seen
    sel_report[p] = (got, target, inside)
    print("   p=%-5s локальных классов найдено %s из ожидаемых %s;  требуемый класс локально %s"
          % (p, got, target, "ЕСТЬ" if inside else "НЕТ"))

# вещественное место
Rsize = 2
seenR = set()
seenR.add((1, 1, 1))
for i in range(3):
    cls = []
    for j in range(3):
        v = ROOTS[i] - ROOTS[j]
        if v == 0:
            a, c = [z for z in range(3) if z != j]
            v = (ROOTS[j] - ROOTS[a]) * (ROOTS[j] - ROOTS[c])
        cls.append(loc_class_R(v))
    seenR.add(tuple(cls))
for _ in range(20000):
    x = QQ(ZZ.random_element(-10**8, 10**8)) / ZZ.random_element(1, 1000)
    f = (x - ROOTS[0]) * (x - ROOTS[1]) * (x - ROOTS[2])
    if f > 0:
        seenR.add(tuple(loc_class_R(x - ROOTS[i]) for i in range(3)))
reqR = tuple(loc_class_R(QQ(req_class[i])) for i in range(3))
print("   место R: локальных классов %s (ожидалось %s); требуемый %s -> %s"
      % (len(seenR), Rsize, reqR.__str__(), "ЕСТЬ" if reqR in seenR else "НЕТ"))
sel_report['R'] = (len(seenR), Rsize, reqR in seenR)

all_complete = all(sel_report[k][0] == sel_report[k][1] for k in sel_report)
all_inside = all(sel_report[k][2] for k in sel_report)
rec("I1  локальные образы посчитаны ПОЛНОСТЬЮ во всех местах", all_complete,
    str({k: sel_report[k] for k in sel_report}))
rec("I2  требуемый класс ЛОКАЛЬНО РАЗРЕШИМ ВЕЗДЕ (=> лежит в группе Селмера)", all_inside)
if all_inside and all_complete:
    print("   СЛЕДСТВИЕ: локального препятствия НЕТ.  Значит исключение (11,4) держится")
    print("   ИСКЛЮЧИТЕЛЬНО на верхней границе ранга r <= 1.  Если r = 3, образ = Селмер")
    print("   и требуемый класс окажется В ОБРАЗЕ -> заявление ложно.")

# ============================================================ J
hdr("J. СВЕРКА С ЧИСЛАМИ CODEX (только здесь используются его величины)")

codex = {
    's': QQ(137) / 2, 'b': 132616,
    'e': (QQ(-132616), QQ(-2005817) / 2, QQ(-17536)),
    'G': (QQ(70664), QQ(138738600)),
    'req': (1, 274, 274),
    'img': [(1, 1, 1), (-1, 28770, -28770), (-28770, 137, -210), (28770, 210, 137),
            (105, 210, 2), (-105, 137, -14385), (-274, 28770, -105), (274, 1, 274)],
}
rec("J1  s совпадает", codex['s'] == s)
rec("J2  b совпадает", codex['b'] == b)
rec("J3  корни совпадают И В ТОМ ЖЕ ПОРЯДКЕ", codex['e'] == (E1, E2, E3))
rec("J4  требуемый класс совпадает", codex['req'] == req_class)
try:
    Gc = E(codex['G'][0], codex['G'][1])
    rec("J5  точка G Codex лежит на МОЕЙ кривой E", True, "G=%s, delta(G)=%s" % (Gc, delta(Gc).__str__()))
    rec("J6  G имеет бесконечный порядок", Gc.order() == Infinity)
except Exception as ex:
    rec("J5  точка G Codex лежит на моей E", False, str(ex))
rec("J7  список из 8 классов Codex совпадает с моим образом",
    sorted(codex['img']) == sorted(img), "мой образ: %s" % sorted(img).__str__())

# ============================================================ ИТОГ
hdr("ИТОГ")
nfail = [k for k in VERD if not VERD[k][0]]
for k in sorted(VERD):
    print("   %-6s %s%s" % ("OK" if VERD[k][0] else "ПРОВАЛ", k, ("  -- " + VERD[k][1]) if VERD[k][1] else ""))
print("\n   провалов: %s" % len(nfail))
print("   %s" % nfail)
