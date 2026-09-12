# -*- coding: utf-8 -*-
# ПРЯМАЯ ПРОВЕРКА заявленного ЛОКАЛЬНОГО исключения пары (m,n) = (16,5) по модулю 41.
#
# Источник проверяемого утверждения: PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md (Codex),
#   «полная система девяти клеток для (16,5) не имеет решений уже над F_41».
#
# Что делает скрипт (всё считается заново, формулы клеток НЕ берутся на веру):
#   §0  символический вывод девяти клеток из аксиом магического квадрата
#       (сумма линии = 3*центр; противоположные клетки дают 2*центр) + четыре автоматических квадрата;
#   §1  контроли процедуры: параметры с ИЗВЕСТНЫМ глобальным решением обязаны выживать при всех p;
#   §2  главный перебор: все 42 точки P^1(F_41), пять неавтоматических клеток;
#   §3  независимое второе кодирование: перебор всех 1680 пар (a,b) != (0,0) в F_41^2,
#       девять клеток, квадратность проверяется таблицей квадратов (без символа Лежандра);
#   §4  сверка с таблицей Codex построчно;
#   §5  вырождения: нули клеток, совпадения клеток, знаки — опирается ли вывод на них;
#   §6  подъём в Q_41 (Гензель по дискам) — независимо от редукционного аргумента;
#   §7  побочное: минимальное p, исключающее пару; список запрещённых λ = m/n при p = 41 и p = 17.
#
# Запуск: sage /home/kep/magicKube/bridge/local_16_5_p41_прямая.sage
# (файл написан так, что работает и как обычный python3 — Sage используется только в §0)

from fractions import Fraction as Fr
from math import gcd, isqrt
import itertools, sys

M, N, P = int(16), int(5), int(41)   # int(): Sage-препарсер иначе делает Integer, Fraction его не ест

def head(s):
    print()
    print("=" * 100)
    print(s)
    print("=" * 100)

# ----------------------------------------------------------------------------------------------
head("§0. СИМВОЛИЧЕСКИЙ ВЫВОД ДЕВЯТИ КЛЕТОК ИЗ АКСИОМ (формулы не берутся из чужих файлов)")
# ----------------------------------------------------------------------------------------------
# Аксиомы: c4 = C — центр; все восемь линий = 3C; каждая пара противоположных клеток = 2C.
# Отсюда c0 = C+α, c8 = C-α, c1 = C+β, c7 = C-β, c2 = C+γ, c6 = C-γ, c3 = C+δ, c5 = C-δ,
# и связи: α+β+γ = 0 (строка 0), δ = γ-α (столбец 0). Свободных параметров два.
#
# Семейство G1 задаётся четырьмя АВТОМАТИЧЕСКИМИ квадратами на рёбрах:
#   c1 = (am+bn)^2, c3 = (am-bn)^2, c5 = (an+bm)^2, c7 = (an-bm)^2,  t = a/b.
# Тогда C вынужден: c1+c7 = c3+c5 = (a^2+b^2)(m^2+n^2) = 2C.
try:
    R = PolynomialRing(QQ, ['a', 'b', 'm', 'n'])          # noqa: F821  (Sage)
    a, b, m, n = R.gens()
    HAVE_SAGE = True
except NameError:
    import sympy as sp
    a, b, m, n = sp.symbols('a b m n')
    HAVE_SAGE = False

def expand(e):
    # в Sage элементы кольца многочленов уже в нормальной форме; в sympy нужен expand
    return e if HAVE_SAGE else __import__('sympy').expand(e)

c1 = (a * m + b * n) ** 2
c3 = (a * m - b * n) ** 2
c5 = (a * n + b * m) ** 2
c7 = (a * n - b * m) ** 2

# центр из ассоциативности (проверяем, что обе пары рёбер дают ОДИН И ТОТ ЖЕ центр)
C_from_17 = (c1 + c7) / 2
C_from_35 = (c3 + c5) / 2
print("  c1+c7 = c3+c5 (обе рёберные пары дают один центр):", expand(C_from_17 - C_from_35) == 0)
C = expand(C_from_17)
print("  => c4 = C =", C, "   т.е. (a^2+b^2)(m^2+n^2)/2 = F4:",
      expand(C - (a * a + b * b) * (m * m + n * n) / 2) == 0)

beta  = expand(c1 - C)
delta = expand(c3 - C)
# из α+β+γ = 0 и δ = γ-α:  2γ = δ-β
gamma = expand((delta - beta) / 2)
alpha = expand(-beta - gamma)
c0 = expand(C + alpha); c8 = expand(C - alpha)
c2 = expand(C + gamma); c6 = expand(C - gamma)
c4 = C
cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]

print("  ВЫВЕДЕНО:")
print("    c0 =", c0)
print("    c2 =", c2, "   (= F4 - 2mn*ab:", expand(c2 - (C - 2 * m * n * a * b)) == 0, ")")
print("    c6 =", c6, "   (= F4 + 2mn*ab:", expand(c6 - (C + 2 * m * n * a * b)) == 0, ")")
print("    c8 =", c8)
print("  сверка с формулами проекта (EVEN_N_2ADIC §1, PARALLEL_RESULT_FULL_NINE_126):")
print("    c0 = n^2 a^2 + m^2 b^2 :", expand(c0 - (n * n * a * a + m * m * b * b)) == 0)
print("    c8 = m^2 a^2 + n^2 b^2 :", expand(c8 - (m * m * a * a + n * n * b * b)) == 0)

LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
ok_lines = all(expand(cells[i] + cells[j] + cells[k] - 3 * c4) == 0 for i, j, k in LINES)
ok_pairs = all(expand(cells[i] + cells[8 - i] - 2 * c4) == 0 for i in range(4))
print("  все восемь линий = 3*c4:", ok_lines, "; все четыре пары противоположных = 2*c4:", ok_pairs)

# какие клетки квадраты ТОЖДЕСТВЕННО (корректная проверка: полный квадрат в кольце многочленов)
def is_poly_square(e):
    if HAVE_SAGE:
        try:
            f = R(e * 4)          # умножаем на 4, чтобы уйти от знаменателя 2, квадратность не меняется
        except Exception:
            return False
        if f == 0:
            return True
        fac = f.factor()
        return f.is_square()
    else:
        import sympy as sp
        r = sp.sqrt(sp.expand(e))
        return r.is_polynomial(a, b, m, n) if hasattr(r, 'is_polynomial') else False

auto = [i for i in range(9) if is_poly_square(cells[i])]
print("  тождественные (автоматические) квадраты: клетки", auto, " -> ожидались рёбра [1,3,5,7]")
NONAUTO = [i for i in range(9) if i not in auto]
print("  НЕАВТОМАТИЧЕСКИЕ клетки (их и надо проверять):", NONAUTO, " -> ожидались [0,2,4,6,8]")

# ----------------------------------------------------------------------------------------------
# Числовые клетки: пять неавтоматических форм от (a,b) при данных (m,n). Точная рациональная арифметика.
# ----------------------------------------------------------------------------------------------
def cells_Q(mm, nn, aa, bb):
    """все девять клеток как Fraction при t = aa/bb (однородно, домножено на bb^2)"""
    mm, nn, aa, bb = int(mm), int(nn), int(aa), int(bb)   # Sage-препарсер даёт Integer, Fraction его не ест
    F4 = Fr((mm * mm + nn * nn) * (aa * aa + bb * bb), int(2))
    return [Fr(nn * nn * aa * aa + mm * mm * bb * bb),
            Fr((aa * mm + bb * nn) ** 2),
            F4 - int(2) * mm * nn * aa * bb,
            Fr((aa * mm - bb * nn) ** 2),
            F4,
            Fr((aa * nn + bb * mm) ** 2),
            F4 + int(2) * mm * nn * aa * bb,
            Fr((aa * nn - bb * mm) ** 2),
            Fr(mm * mm * aa * aa + nn * nn * bb * bb)]

def five_Q(mm, nn, aa, bb):
    c = cells_Q(mm, nn, aa, bb)
    return [c[i] for i in (0, 2, 4, 6, 8)]

def is_rat_square(x):
    """точный тест «x — квадрат рационального» через isqrt (никаких float: коэффициенты бывают огромны)"""
    if x < 0: return False
    nu, de = int(x.numerator), int(x.denominator)
    return isqrt(nu) ** 2 == nu and isqrt(de) ** 2 == de

# ----------------------------------------------------------------------------------------------
head("§1. КОНТРОЛИ ПРОЦЕДУРЫ (обязательны: без них пустой перебор ничего не значит)")
# ----------------------------------------------------------------------------------------------
def sq_or_zero_mod(x_frac, p):
    """редукция Fraction в F_p и ответ: квадрат-или-ноль? (знаменатель обязан быть обратим)"""
    num, den = x_frac.numerator, x_frac.denominator
    if den % p == 0:
        return None                       # не p-целое; для наших форм при p ∤ 2 не бывает
    v = (num % p) * pow(den % p, p - 2, p) % p
    if v == 0:
        return True
    return pow(v, (p - 1) // 2, p) == 1

CTRL_PRIMES = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113]
controls = [
    ("(7,1), t=0   [вырожденная глобальная точка, REVIEW §1]", 7, 1, 0, 1),
    ("(15,8), t=1  [вырожденная глобальная точка]",           15, 8, 1, 1),
    ("(15,8), t=-1 [вырожденная глобальная точка]",           15, 8, -1, 1),
    ("(7,1), t=inf [точка на бесконечности, (a:b)=(1:0)]",     7, 1, 1, 0),
]
ctrl_ok = True
for name, mm, nn, aa, bb in controls:
    vals = five_Q(mm, nn, aa, bb)
    allsq = all(is_rat_square(v) for v in vals)
    bad = [p for p in CTRL_PRIMES if not all(sq_or_zero_mod(v, p) for v in vals)]
    print(f"  {name}")
    print(f"     пять клеток над Q = {[str(v) for v in vals]}; все квадраты в Q: {allsq}")
    print(f"     редукция даёт квадрат-или-ноль при ВСЕХ p из списка: {not bad}"
          + (f"   ПРОВАЛ при p = {bad}" if bad else ""))
    if allsq and bad:
        ctrl_ok = False
print("  ИТОГ §1:", "контроли пройдены — процедура не режет настоящие решения" if ctrl_ok
      else "!!! КОНТРОЛЬ ПРОВАЛЕН — всё ниже недействительно")

# ----------------------------------------------------------------------------------------------
head(f"§2. ГЛАВНОЕ: (m,n) = ({M},{N}), p = {P}. ПОЛНЫЙ перебор всех {P+1} точек P^1(F_{P})")
# ----------------------------------------------------------------------------------------------
m2, n2, mn2 = M * M % P, N * N % P, 2 * M * N % P
inv2 = pow(2, P - 2, P)
s_mod = (m2 + n2) * inv2 % P
QRs = sorted({x * x % P for x in range(1, P)})
print(f"  константы mod {P}: m^2={m2}, n^2={n2}, 2mn={mn2}, m^2+n^2={(m2+n2)%P}, s=(m^2+n^2)/2={s_mod}")
print(f"  ненулевые квадраты F_{P} ({len(QRs)} шт.): {QRs}")
print(f"  сверка списка QR со списком Codex: "
      f"{QRs == [1,2,4,5,8,9,10,16,18,20,21,23,25,31,32,33,36,37,39,40]}")

def five_mod(aa, bb):
    """пять неавтоматических клеток в F_p для проективной точки (aa:bb)"""
    F0 = (n2 * aa * aa + m2 * bb * bb) % P
    F4 = s_mod * (aa * aa + bb * bb) % P
    F8 = (m2 * aa * aa + n2 * bb * bb) % P
    L  = (F4 - mn2 * aa * bb) % P
    U  = (F4 + mn2 * aa * bb) % P
    return F0, F4, F8, L, U

SQ_OR_0 = set(QRs) | {0}
points = [(t, 1) for t in range(P)] + [(1, 0)]      # все p+1 точек P^1(F_p), включая бесконечность
assert len(points) == P + 1

surv3, surv5, table = [], [], []
for (aa, bb) in points:
    F0, F4, F8, L, U = five_mod(aa, bb)
    ok3 = all(v in SQ_OR_0 for v in (F0, F4, F8))
    ok5 = ok3 and all(v in SQ_OR_0 for v in (L, U))
    if ok3: surv3.append((aa, bb))
    if ok5: surv5.append((aa, bb))
    table.append(((aa, bb), (F0, F4, F8), (L, U), ok3, ok5))

def pname(aa, bb):
    return "inf" if bb == 0 else str(aa)

print(f"\n  (а) после ТРЁХ квадратностей (F0,F4,F8) выживает {len(surv3)} точек: "
      f"{[pname(*pt) for pt in surv3]}")
print(f"      Codex: «остаются ровно четыре конечных параметра 4,10,31,37»  -> совпадение: "
      f"{sorted(pname(*pt) for pt in surv3) == sorted(['10','31','37','4'])}")

print("\n  (б) подробная таблица выживших троек и что с ними делают красные L,U:")
print("          t |   F0   F4   F8 |    L    U | L квадрат? U квадрат? | вердикт")
for pt, (F0, F4, F8), (L, U), ok3, ok5 in table:
    if not ok3: continue
    print(f"      {pname(*pt):>5} | {F0:4} {F4:4} {F8:4} | {L:4} {U:4} | "
          f"{str(L in SQ_OR_0):>9} {str(U in SQ_OR_0):>10} | {'ПРОХОД' if ok5 else 'отказ'}")

print(f"\n  (в) точка на бесконечности (1:0): (F0,F4,F8,L,U) = {five_mod(1,0)}; "
      f"F4 = s = {s_mod}, квадрат? {s_mod in SQ_OR_0}")
print(f"\n  ИТОГ §2: точек P^1(F_{P}) со всеми ПЯТЬЮ квадратностями: {len(surv5)} "
      f"(из проверенных {len(points)})")

# ----------------------------------------------------------------------------------------------
head("§3. НЕЗАВИСИМОЕ ВТОРОЕ КОДИРОВАНИЕ: все 1680 пар (a,b) != (0,0), ДЕВЯТЬ клеток, "
     "квадратность — по таблице квадратов")
# ----------------------------------------------------------------------------------------------
# Здесь не используются ни символ Лежандра, ни имена F0/F4/F8/L/U: клетки берутся из §0-вывода
# и вычисляются заново; квадратность — «существует ли x in F_p с x^2 = v».
sq_table = {}
for x in range(P):
    sq_table.setdefault(x * x % P, []).append(x)

def cells_mod_from_scratch(aa, bb):
    F4 = (M * M + N * N) % P * inv2 % P * ((aa * aa + bb * bb) % P) % P
    return [(N * N * aa * aa + M * M * bb * bb) % P,
            (aa * M + bb * N) ** 2 % P,
            (F4 - 2 * M * N * aa * bb) % P,
            (aa * M - bb * N) ** 2 % P,
            F4,
            (aa * N + bb * M) ** 2 % P,
            (F4 + 2 * M * N * aa * bb) % P,
            (aa * N - bb * M) ** 2 % P,
            (M * M * aa * aa + N * N * bb * bb) % P]

good_pairs = []
auto_always = True
for aa in range(P):
    for bb in range(P):
        if aa == 0 and bb == 0: continue
        cs = cells_mod_from_scratch(aa, bb)
        for i in (1, 3, 5, 7):
            if cs[i] not in sq_table: auto_always = False
        if all(cs[i] in sq_table for i in range(9)):
            good_pairs.append((aa, bb))
print(f"  проверено пар (a,b): {P*P-1}")
print(f"  рёберные клетки c1,c3,c5,c7 всегда квадраты (контроль вывода §0): {auto_always}")
print(f"  пар со ВСЕМИ ДЕВЯТЬЮ клетками-квадратами(-или-нулями): {len(good_pairs)}")
if good_pairs:
    print(f"  примеры: {good_pairs[:10]}")
    proj = sorted({(aa * pow(bb, P - 2, P) % P if bb else 'inf') for aa, bb in good_pairs})
    print(f"  соответствующие проективные t: {proj}")
print(f"  согласие §2 и §3: {(len(good_pairs) == 0) == (len(surv5) == 0)}")

# ----------------------------------------------------------------------------------------------
head("§4. ПОСТРОЧНАЯ СВЕРКА С ТАБЛИЦЕЙ CODEX (PARALLEL_RESULT_FULL_NINE_126, раздел «(16,5), p=41»)")
# ----------------------------------------------------------------------------------------------
codex_rows = {4: ((0,31,21), (6,15)), 10: ((9,25,0), (24,26)),
              31: ((9,25,0), (26,24)), 37: ((0,31,21), (15,6))}
allmatch = True
for t in sorted(codex_rows):
    F0, F4, F8, L, U = five_mod(t, 1)
    ok = (F0, F4, F8) == codex_rows[t][0] and (L, U) == codex_rows[t][1]
    allmatch &= ok
    print(f"  t={t:>2}: мой расчёт (F0,F4,F8)=({F0},{F4},{F8}), (L,U)=({L},{U}); "
          f"Codex {codex_rows[t][0]}, {codex_rows[t][1]} -> {'совпадает' if ok else 'РАСХОЖДЕНИЕ'}")
print(f"  Codex: «на бесконечности F4 = s = 38 не квадрат» -> s={s_mod}, не квадрат: {s_mod not in SQ_OR_0}")
print(f"  Codex: «оба не квадраты» в каждой из четырёх строк -> проверка: "
      f"{all((L not in SQ_OR_0) and (U not in SQ_OR_0) for t in codex_rows for L, U in [five_mod(t,1)[3:]])}")
print(f"  ИТОГ §4: таблица Codex воспроизведена полностью: {allmatch}")

# ----------------------------------------------------------------------------------------------
head("§5. ВЫРОЖДЕНИЯ: опирается ли вывод на положительность, различность или запрет нулей")
# ----------------------------------------------------------------------------------------------
print("  (а) нули среди клеток у четырёх выживших троек:")
for t in sorted(codex_rows):
    cs = cells_mod_from_scratch(t, 1)
    zeros = [i for i in range(9) if cs[i] == 0]
    dupl = [(i, j) for i in range(9) for j in range(i + 1, 9) if cs[i] == cs[j]]
    print(f"      t={t:>2}: девять клеток mod {P} = {cs}; нулевые клетки: {zeros or 'нет'}; "
          f"совпадающие пары: {len(dupl)}")
print("  (б) нули РАЗРЕШЕНЫ в тесте (ненулевой рациональный квадрат может редуцироваться в 0).")
print("      Это ослабляет тест; значит пустота при разрешённых нулях — более сильный результат.")
strict = [pt for pt in points if all(v in set(QRs) for v in five_mod(*pt))]   # нули ЗАПРЕЩЕНЫ
print(f"      контроль в обратную сторону (нули запрещены, тест строже): выживших точек {len(strict)}")
print("  (в) положительность и попарная различность в §2/§3 НЕ используются вовсе:")
print("      тест — только «значение клетки есть квадрат или ноль в F_p». Ни одного сравнения >0,")
print("      ни одной проверки c_i != c_j в решающем переборе нет (см. §2, §3).")
print("  (г) проверка: не спасают ли пару вырожденные точки над Q (их тест обязан пропускать).")
for t in [Fr(int(0)), Fr(int(1)), Fr(int(-1)), Fr(M, N), Fr(N, M), Fr(int(3), int(2))]:
    vals = five_Q(M, N, t.numerator, t.denominator)
    print(f"      t={str(t):>6}: пять клеток = {[str(v) for v in vals]}; все квадраты в Q: "
          f"{all(is_rat_square(v) for v in vals)}")

# ----------------------------------------------------------------------------------------------
head("§6. НЕЗАВИСИМЫЙ ПОДЪЁМ В Q_41 (рекурсия по дискам + Гензель, без ссылки на §2)")
# ----------------------------------------------------------------------------------------------
def vp(x, p):
    if x == 0: return None
    nu, de = x.numerator, x.denominator
    v = 0
    while nu % p == 0: nu //= p; v += 1
    while de % p == 0: de //= p; v -= 1
    return v

def is_sq_Qp(x, p):
    v = vp(x, p)
    if v is None: return True
    if v % 2: return False
    y = x / Fr(int(p)) ** int(v)
    u = (y.numerator % p) * pow(y.denominator % p, p - 2, p) % p
    return pow(u, (p - 1) // 2, p) == 1

def solve_disk(mm, nn, p, r, k, KMAX, at_inf, budget):
    if budget[0] <= 0: return None
    budget[0] -= 1
    vals = five_Q(mm, nn, 1, r) if at_inf else five_Q(mm, nn, r, 1)
    undec = False
    for x in vals:
        v = vp(x, p)
        if v is not None and v < k:
            if not is_sq_Qp(x, p): return False
        else:
            undec = True
    if not undec: return True
    if k >= KMAX: return None
    res = False
    for z in range(p):
        out = solve_disk(mm, nn, p, r + z * p ** k, k + 1, KMAX, at_inf, budget)
        if out is True: return True
        if out is None: res = None
    return res

bud = [200000]
o_fin = solve_disk(M, N, P, 0, 0, 5, False, bud)
o_inf = solve_disk(M, N, P, 0, 1, 5, True, bud)
print(f"  ветка t in Z_{P}:           {o_fin}   (False = точек нет)")
print(f"  ветка t = 1/u, u in {P}Z_{P}: {o_inf}")
print(f"  ИТОГ §6: Q_{P}-точка пяти квадратностей существует: "
      f"{o_fin is True or o_inf is True}; доказано её отсутствие: {o_fin is False and o_inf is False}")

# ----------------------------------------------------------------------------------------------
head("§7. ПОБОЧНОЕ: минимальное исключающее p; запрещённые классы λ = m/n")
# ----------------------------------------------------------------------------------------------
def empty_at(mm, nn, p, k=5):
    m2_, n2_, mn2_ = mm * mm % p, nn * nn % p, 2 * mm * nn % p
    i2 = pow(2, p - 2, p); s_ = (m2_ + n2_) * i2 % p
    S = {x * x % p for x in range(p)}
    for (aa, bb) in [(t, 1) for t in range(p)] + [(1, 0)]:
        F0 = (n2_ * aa * aa + m2_ * bb * bb) % p
        F4 = s_ * (aa * aa + bb * bb) % p
        F8 = (m2_ * aa * aa + n2_ * bb * bb) % p
        vals = (F0, F4, F8) if k == 3 else (F0, F4, F8, (F4 - mn2_ * aa * bb) % p, (F4 + mn2_ * aa * bb) % p)
        if all(v in S for v in vals): return False
    return True

PR = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
k5 = [p for p in PR if empty_at(M, N, p, 5)]
k3 = [p for p in PR if empty_at(M, N, p, 3)]
print(f"  (16,5): простые p <= 97, где ПЯТЬ клеток локально невозможны: {k5 or 'нет'}")
print(f"  (16,5): простые p <= 97, где ТРИ клетки (кривая C7) локально невозможны: {k3 or 'нет'}")
print(f"  => минимальное исключающее p для девяти клеток: {k5[0] if k5 else 'нет'}")
print(f"  => C7 локально НЕ исключается ни одним из этих p: {not k3}  "
      f"(то есть пустота C7(Q) этим методом не доказана — и не заявляется)")

for p in (41, 17):
    bad_lam = []
    for lam in range(p):
        if empty_at(lam, 1, p, 5): bad_lam.append(lam)
    inf_empty = empty_at(1, 0, p, 5)
    print(f"  p={p}: запрещённые λ = m/n (m=λ, n=1): {bad_lam}; класс n=0 mod p исключён: {inf_empty}")
codex_41 = [5,8,13,19,22,28,33,36]; codex_17 = [2,3,6,8,9,11,14,15]
lam41 = [l for l in range(41) if empty_at(l, 1, 41, 5)]
lam17 = [l for l in range(17) if empty_at(l, 1, 17, 5)]
print(f"  сверка с Codex: p=41 {lam41 == codex_41}; p=17 {lam17 == codex_17}")
lam_of_pair = M * pow(N, P - 2, P) % P
print(f"  λ = m/n = {M}/{N} mod {P} = {lam_of_pair}; попадает в запрещённый список: {lam_of_pair in lam41}")

# ----------------------------------------------------------------------------------------------
head("ФИНАЛЬНЫЙ ВЕРДИКТ")
# ----------------------------------------------------------------------------------------------
print(f"  контроли §1 пройдены:                                  {ctrl_ok}")
print(f"  точек P^1(F_41) со всеми пятью квадратностями (§2):    {len(surv5)}")
print(f"  пар (a,b) с девятью квадратами-или-нулями (§3):        {len(good_pairs)}")
print(f"  таблица Codex воспроизведена (§4):                     {allmatch}")
print(f"  Q_41-точка отсутствует, доказано подъёмом (§6):        {o_fin is False and o_inf is False}")
verdict = ctrl_ok and len(surv5) == 0 and len(good_pairs) == 0
print(f"\n  ЛОКАЛЬНОЕ ПРЕПЯТСТВИЕ ДЛЯ (16,5) ПРИ p=41 ПОДТВЕРЖДЕНО: {verdict}")
print("  Область утверждения: нет точки (a:b) in P^1(Q_41), при которой все пять неавтоматических")
print("  клеток — квадраты в Q_41. Следовательно нет и рациональной. Это НЕ утверждение C7(Q)=пусто.")
