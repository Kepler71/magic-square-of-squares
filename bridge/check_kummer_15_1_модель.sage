# -*- coding: utf-8 -*-
# =====================================================================
#  АТАКА НА ЗАЯВЛЕНИЕ CODEX:  C_{15,1}(Q) = пусто  (семейство G1)
#  Угол атаки: МОДЕЛЬ И ТОЖДЕСТВА.
#  Автор проверки: Claude (независимая пересборка, числа Codex
#  используются ТОЛЬКО в финальном разделе сверки L).
#
#  Установка: НЕ подтвердить, а СЛОМАТЬ.
#
#  A. Общие символьные тождества в Q[m,n,t]
#  B. Специализация (15,1): E, дискриминант, корни, 2-кручение
#  C. Собственная реализация delta + калибровка на известных кривых
#  D. Ранг E(Q) несколькими независимыми инструментами
#  E. Восемь классов, различность, требуемый класс
#  F. ХРУПКОСТЬ ПОРЯДКА КОРНЕЙ: все 6 перестановок
#  G. Бесконечность: проективно, все 8 геометрических точек
#  H. Гладкость аффинной модели C (нормализация не добавляет Q-точек)
#  I. Положительные контроли: (4,3) и (15,8) — метод НЕ должен исключать
#  J. Прямой поиск точек на C и на торсоре D
#  K. Локальная разрешимость C (нетривиальность исключения)
#  L. Сверка с числами Codex
# =====================================================================

import sys
from sage.all import *

FAIL = []
UNVER = []

def hdr(s):
    print("\n" + "=" * 72)
    print(s)
    print("=" * 72)
    sys.stdout.flush()

def ck(name, cond, extra=""):
    ok = bool(cond)
    print(("  [OK]     " if ok else "  [ПРОВАЛ] ") + name + (("   " + extra) if extra else ""))
    sys.stdout.flush()
    if not ok:
        FAIL.append(name)
    return ok

def unver(name, why=""):
    print("  [НЕ ПРОВЕРЕНО] " + name + (("   " + why) if why else ""))
    UNVER.append(name + ((" | " + why) if why else ""))


# ---------------------------------------------------------------------
def sqclass(x):
    """Представитель квадратного класса рационального x != 0:
       бесквадратное целое со знаком."""
    x = QQ(x)
    if x == 0:
        raise ValueError("sqclass(0)")
    z = ZZ(x.numerator() * x.denominator())   # x * (denom)^2 -> тот же класс
    return ZZ(z.sign()) * ZZ(z.abs().squarefree_part())


def delta_of_point(P, e):
    """Полный 2-спусковый гомоморфизм для y^2=(x-e0)(x-e1)(x-e2).
       e — СПИСОК корней в выбранном порядке."""
    if P.is_zero():
        return (ZZ(1), ZZ(1), ZZ(1))
    X = P[0]; V = P[1]
    d = [QQ(X - e[i]) for i in range(3)]
    for i in range(3):
        if d[i] == 0:
            j, k = [q for q in range(3) if q != i]
            d[i] = QQ((e[i] - e[j]) * (e[i] - e[k]))
    return tuple(sqclass(v) for v in d)


def curve_from_roots(e):
    """E: y^2 = (x-e0)(x-e1)(x-e2)."""
    R = PolynomialRing(QQ, 'x'); x = R.gen()
    f = (x - e[0]) * (x - e[1]) * (x - e[2])
    c = f.coefficients(sparse=False)   # c0 + c1 x + c2 x^2 + x^3
    return EllipticCurve([0, c[2], 0, c[1], c[0]])


# =====================================================================
hdr("A. ОБЩИЕ СИМВОЛЬНЫЕ ТОЖДЕСТВА в Q[m,n,t]  (никакой специализации)")

R3 = PolynomialRing(QQ, ['m', 'n', 't'])
m, n, t = R3.gens()

s_sym = (m**2 + n**2) / 2          # полуцелое допустимо: коэффициенты в Q
F0 = m**2 + n**2 * t**2
F4 = s_sym * (1 + t**2)
F8 = n**2 + m**2 * t**2
b_sym = s_sym * m**2 * n**2

e1 = -b_sym
e2 = -s_sym * m**4
e3 = -s_sym * n**4
X = b_sym * t**2

ck("F0 + F8 = 2 F4  (диагональная арифметическая прогрессия)", F0 + F8 - 2 * F4 == 0)
ck("X - e1 == (m n)^2 * F4", (X - e1) - (m * n)**2 * F4 == 0)
ck("X - e2 == s * m^2 * F0", (X - e2) - s_sym * m**2 * F0 == 0)
ck("X - e3 == s * n^2 * F8", (X - e3) - s_sym * n**2 * F8 == 0)
ck("(X-e1)(X-e2)(X-e3) == b^2 * F0 F4 F8   (т.е. V = b u0u4u8 корректно)",
   (X - e1) * (X - e2) * (X - e3) - b_sym**2 * F0 * F4 * F8 == 0)

# множители при u_i^2 — квадратные классы координат delta
print("\n  Множители: X-e1 = (mn)^2 * u4^2  -> класс 1")
print("             X-e2 = s*m^2 * u0^2   -> класс s  (т.к. (m u0)^2 квадрат)")
print("             X-e3 = s*n^2 * u8^2   -> класс s  (т.к. (n u8)^2 квадрат)")
print("  => НЕОБХОДИМЫЙ класс delta = (1, s, s) при порядке (e1,e2,e3)=(-b, -s m^4, -s n^4)")

# различность корней в общем виде
ck("e1 != e2 в общем виде (разность = s m^2 (n^2-m^2), ненулевая при m!=n)",
   (e1 - e2) == s_sym * m**2 * (m**2 - n**2))
ck("e1 != e3 в общем виде", (e1 - e3) == s_sym * n**2 * (n**2 - m**2))
ck("e2 != e3 в общем виде", (e2 - e3) == s_sym * (n**4 - m**4))


# =====================================================================
hdr("B. СПЕЦИАЛИЗАЦИЯ (m,n) = (15,1)")

M, N = ZZ(15), ZZ(1)
ck("gcd(m,n) = 1", gcd(M, N) == 1)
S = QQ(M**2 + N**2) / 2
B = S * M**2 * N**2
print("  s = %s   b = %s" % (S, B))
ck("s = 113 (моим счётом)", S == 113)
ck("b = 25425 (моим счётом)", B == 25425)

E1 = -B
E2 = -S * M**4
E3 = -S * N**4
E_ROOTS = [QQ(E1), QQ(E2), QQ(E3)]
print("  e1 = %s   e2 = %s   e3 = %s" % (E1, E2, E3))

E = curve_from_roots(E_ROOTS)
print("  E: %s" % E)
print("  минимальная модель: %s" % (E.minimal_model().a_invariants(),))
print("  кондуктор: %s" % E.conductor())

ck("дискриминант E != 0", E.discriminant() != 0, "disc = %s" % E.discriminant().factor())
ck("три корня попарно различны", len(set(E_ROOTS)) == 3)
ck("кубика справа раскладывается над Q (=> полное 2-кручение)",
   all(E.is_on_curve(r, 0) for r in E_ROOTS))
T = E.torsion_subgroup()
print("  E(Q)_tors = %s" % (T.invariants(),))
ck("E(Q)[2] = (Z/2)^2  (ровно 4 точки порядка 1 или 2)",
   sum(1 for P in E.torsion_points() if P.order() in (1, 2)) == 4)

# требуемый класс
REQ = (sqclass(1), sqclass(S), sqclass(S))
print("\n  ТРЕБУЕМЫЙ класс (моим счётом) = %s" % (REQ,))
ck("произведение координат требуемого класса — квадрат (необходимое условие)",
   QQ(REQ[0] * REQ[1] * REQ[2]).is_square())

# прямая проверка тождеств на специализации через кольцо Q[t]
Rt = PolynomialRing(QQ, 't'); tt = Rt.gen()
f0 = M**2 + N**2 * tt**2
f4 = S * (1 + tt**2)
f8 = N**2 + M**2 * tt**2
Xt = B * tt**2
ck("(15,1): X-e1 = 225*F4 как многочлены", Xt - E1 - (M * N)**2 * f4 == 0)
ck("(15,1): X-e2 = 113*225*F0 как многочлены", Xt - E2 - S * M**2 * f0 == 0)
ck("(15,1): X-e3 = 113*1*F8 как многочлены", Xt - E3 - S * N**2 * f8 == 0)
ck("F0,F4,F8 не имеют вещественных корней (u_i != 0 всюду на C(R))",
   len(f0.roots(RR)) == 0 and len(f4.roots(RR)) == 0 and len(f8.roots(RR)) == 0)


# =====================================================================
hdr("C. КАЛИБРОВКА собственной реализации delta на известных кривых")

# rank 0, полное 2-кручение: y^2 = x(x-1)(x+1)
Ec0 = curve_from_roots([QQ(0), QQ(1), QQ(-1)])
im0 = set()
for P in Ec0.torsion_points():
    im0.add(delta_of_point(P, [QQ(0), QQ(1), QQ(-1)]))
ck("y^2=x(x-1)(x+1): rank=0 => |im|=4", Ec0.rank() == 0 and len(im0) == 4,
   "получено %d" % len(im0))

# rank 1 (конгруэнтное число 5): y^2 = x(x-25)(x+25)... используем x(x^2-25)
r5 = [QQ(0), QQ(5), QQ(-5)]
Ec1 = curve_from_roots(r5)
im1 = set()
gens1 = Ec1.gens()
base = list(Ec1.torsion_points()) + gens1
pts = set()
for a in range(-3, 4):
    for P in Ec1.torsion_points():
        pts.add(a * gens1[0] + P)
for P in pts:
    im1.add(delta_of_point(P, r5))
ck("y^2=x(x-5)(x+5): rank=1 => |im|=8", Ec1.rank() == 1 and len(im1) == 8,
   "получено %d" % len(im1))

# гомоморфность delta — численный контроль на E_{15,1} (после построения точек — ниже)


# =====================================================================
hdr("D. РАНГ E(Q) — несколько независимых инструментов")

Emin = E.minimal_model()
print("  минимальная модель: %s" % Emin)

try:
    rb = E.rank_bounds()
    print("  Sage E.rank_bounds() = %s" % (rb,))
except Exception as ex:
    rb = None
    print("  rank_bounds не сработал: %s" % ex)

try:
    rk = E.rank(only_use_mwrank=True, proof=True)
    ck("mwrank (proof=True): rank = 1", rk == 1, "rank = %s" % rk)
except Exception as ex:
    unver("mwrank proof=True", str(ex))
    rk = None

try:
    rk2 = E.rank(only_use_mwrank=False)
    ck("Sage rank() (pari/ellrank путь): rank = 1", rk2 == 1, "rank = %s" % rk2)
except Exception as ex:
    unver("Sage rank() второй путь", str(ex))

# прямой вызов PARI ellrank
try:
    pe = pari(Emin.a_invariants()).ellinit()
    er = pe.ellrank()
    print("  PARI ellrank -> %s" % er)
    lo = ZZ(er[0]); hi = ZZ(er[1])
    ck("PARI ellrank: верхняя граница ранга = 1", hi == 1, "[lo,hi]=[%s,%s]" % (lo, hi))
except Exception as ex:
    unver("PARI ellrank", str(ex))

# 2-Selmer
try:
    sr = E.selmer_rank()
    print("  E.selmer_rank() = %s" % sr)
    print("  (ожидание: при rank=1 и нетривиальной Sha[2] селмеров ранг > 1)")
except Exception as ex:
    unver("selmer_rank", str(ex))

# аналитический ранг — ЧИСЛЕННЫЙ кросс-контроль, не доказательство
try:
    ar = E.analytic_rank()
    print("  аналитический ранг (численно) = %s" % ar)
except Exception as ex:
    unver("analytic_rank", str(ex))

try:
    G_list = E.gens()
    print("  Sage E.gens() = %s" % (G_list,))
    ck("ровно один свободный генератор", len(G_list) == 1)
except Exception as ex:
    G_list = []
    unver("E.gens()", str(ex))

# независимый поиск точек: не найдётся ли ВТОРОЙ независимой?
try:
    found = E.point_search(12, rank_bound=2)
    Esat, idx, reg = E.saturation(found) if found else ([], 1, 0)
    print("  point_search(height 12): найдено %d точек, после насыщения ранг подгруппы = %d"
          % (len(found), len(Esat)))
    ck("поиск не даёт ранг >= 2", len(Esat) <= 1, "|базис| = %d" % len(Esat))
except Exception as ex:
    unver("point_search/saturation", str(ex))


# =====================================================================
hdr("E. ВОСЕМЬ КЛАССОВ И ТРЕБУЕМЫЙ КЛАСС")

T1 = E(E_ROOTS[0], 0)
T2 = E(E_ROOTS[1], 0)
T3 = E(E_ROOTS[2], 0)
ck("T1+T2+T3 = O (три 2-кручения)", (T1 + T2 + T3).is_zero())

if G_list:
    Gpt = G_list[0]
else:
    Gpt = None

ck("G — точка E и не кручение", Gpt is not None and Gpt.order() == oo)

rows = []
image = set()
for a in [0, 1]:
    for c in [0, 1]:
        for d in [0, 1]:
            P = a * Gpt + c * T1 + d * T2
            dd = delta_of_point(P, E_ROOTS)
            rows.append(((a, c, d), P, dd))
            image.add(dd)

print("\n  (a,c,d) |        delta (мой счёт)        |  X(P)")
for (acd, P, dd) in rows:
    xs = "O" if P.is_zero() else str(P[0])
    print("  %s |  %-28s | %s" % (str(acd), str(dd), xs))

ck("восемь классов ПОПАРНО РАЗЛИЧНЫ", len(image) == 8, "|image| = %d" % len(image))
ck("каждый класс имеет квадратное произведение координат",
   all(QQ(u * v * w).is_square() for (u, v, w) in image))

# гомоморфность delta на всех парах из 8 точек
homok = True
for i in range(8):
    for j in range(8):
        Pi = rows[i][1]; Pj = rows[j][1]
        di = rows[i][2]; dj = rows[j][2]
        Ps = Pi + Pj
        ds = delta_of_point(Ps, E_ROOTS)
        prod = tuple(sqclass(di[k] * dj[k]) for k in range(3))
        if prod != ds:
            homok = False
            print("    гомоморфность нарушена: %s + %s" % (rows[i][0], rows[j][0]))
ck("delta гомоморфен на всех 64 парах (численный контроль)", homok)

ck("требуемый класс %s ОТСУТСТВУЕТ в образе" % (REQ,), REQ not in image)

# как выглядят «почти совпадающие» классы
print("\n  Проверка соседних (переставленных) классов:")
from itertools import permutations
for perm in permutations([0, 1, 2]):
    cand = tuple(REQ[perm[k]] for k in range(3))
    mark = "ЕСТЬ В ОБРАЗЕ" if cand in image else "нет"
    print("    перестановка %s -> %-20s : %s" % (str(perm), str(cand), mark))


# =====================================================================
hdr("F. ХРУПКОСТЬ ПОРЯДКА КОРНЕЙ: все 6 перестановок (СОГЛАСОВАННО)")

print("  Если И таблицу, И требуемый класс переставить ОДНОЙ и той же")
print("  перестановкой, вывод обязан не меняться.\n")
allsame = True
for perm in permutations([0, 1, 2]):
    eperm = [E_ROOTS[perm[k]] for k in range(3)]
    Eperm = curve_from_roots(eperm)
    imp = set()
    T1p = Eperm(eperm[0], 0); T2p = Eperm(eperm[1], 0)
    Gp = Eperm(Gpt[0], Gpt[1])
    for a in [0, 1]:
        for c in [0, 1]:
            for d in [0, 1]:
                imp.add(delta_of_point(a * Gp + c * T1p + d * T2p, eperm))
    reqp = tuple(REQ[perm[k]] for k in range(3))
    ok = (len(imp) == 8) and (reqp not in imp)
    if not ok:
        allsame = False
    print("    perm %s : |image|=%d, требуемый %s -> %s"
          % (str(perm), len(imp), str(reqp), "ОТСУТСТВУЕТ" if reqp not in imp else "ПРИСУТСТВУЕТ"))
ck("вывод инвариантен относительно СОГЛАСОВАННОЙ перестановки корней", allsame)

print("\n  А теперь НЕСОГЛАСОВАННО: таблица в порядке (e1,e2,e3),")
print("  а требуемый класс ошибочно переставлен — точка обрушения:")
for perm in permutations([0, 1, 2]):
    cand = tuple(REQ[perm[k]] for k in range(3))
    if cand in image:
        print("    !!! ошибка порядка %s дала бы ЛОЖНЫЙ вывод «не исключено»: %s"
              % (str(perm), str(cand)))


# =====================================================================
hdr("G. БЕСКОНЕЧНОСТЬ — проективно, все 8 геометрических точек")

Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
# t = 1/z, U_i = u_i / t = z*u_i  =>  U_i^2 = z^2 F_i(1/z)
g0 = (M**2 * z**2 + N**2)
g4 = S * (z**2 + 1)
g8 = (N**2 * z**2 + M**2)
# честная проверка: подставляем t=1/z в поле рациональных функций и умножаем на z^2
Fz = FractionField(Rz)
zz = Fz.gen()
ck("U0^2 = z^2 * F0(1/z) = m^2 z^2 + n^2", Fz(z**2) * Fz(f0(1 / zz)) == Fz(g0))
ck("U4^2 = z^2 * F4(1/z) = s(z^2+1)", Fz(z**2) * Fz(f4(1 / zz)) == Fz(g4))
ck("U8^2 = z^2 * F8(1/z) = n^2 z^2 + m^2", Fz(z**2) * Fz(f8(1 / zz)) == Fz(g8))
print("  U0^2 = %s ;  U4^2 = %s ;  U8^2 = %s" % (g0, g4, g8))
print("  при z = 0:  U0^2 = %s,  U4^2 = %s,  U8^2 = %s" % (g0(0), g4(0), g8(0)))
ck("U0^2|_{z=0} = n^2 = 1 — рационально", QQ(g0(0)).is_square())
ck("U8^2|_{z=0} = m^2 = 225 — рационально", QQ(g8(0)).is_square())
ck("U4^2|_{z=0} = s = 113 — НЕ квадрат в Q", not QQ(g4(0)).is_square())
ck("все три значения при z=0 ненулевые => накрытие неразветвлено над t=oo",
   g0(0) != 0 and g4(0) != 0 and g8(0) != 0)
K113 = QuadraticField(113, 'a')
ck("все 8 геометрических точек над t=oo имеют поле определения Q(sqrt(113)) != Q",
   not QQ(113).is_square() and K113.degree() == 2)
print("  8 точек = (U0,U4,U8) = (±1, ±sqrt(113), ±15): каждая определена над Q(sqrt113),")
print("  Галуа переставляет знак U4 => ни одна не рациональна.")


# =====================================================================
hdr("H. ГЛАДКОСТЬ АФФИННОЙ МОДЕЛИ C (нормализация не добавляет Q-точек)")

res04 = f0.resultant(f4)
res08 = f0.resultant(f8)
res48 = f4.resultant(f8)
ck("F0,F4 не имеют общих корней (res != 0)", res04 != 0, "res=%s" % res04)
ck("F0,F8 не имеют общих корней", res08 != 0, "res=%s" % res08)
ck("F4,F8 не имеют общих корней", res48 != 0, "res=%s" % res48)
ck("каждая F_i без кратных корней (disc != 0)",
   f0.discriminant() != 0 and f4.discriminant() != 0 and f8.discriminant() != 0)
print("  => аффинная кривая {u_i^2 = F_i} гладка; нормализация — изоморфизм над ней;")
print("  => C(Q) = (аффинные Q-точки) U (Q-точки над t=oo).")
# род
print("  Риман–Гурвиц (степень 8, 6 точек ветвления, e=2 в каждой, 4 точки над каждой):")
print("    2g-2 = 8*(-2) + 6*4 = 8  =>  g = 5")
ck("род 5", (8 * (-2) + 6 * 4) == 8)


# =====================================================================
hdr("I. ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ: метод НЕ должен исключать (4,3) и (15,8)")

def full_test(mm, nn, hbound=12):
    mm = ZZ(mm); nn = ZZ(nn)
    ss = QQ(mm**2 + nn**2) / 2
    bb = ss * mm**2 * nn**2
    ee = [QQ(-bb), QQ(-ss * mm**4), QQ(-ss * nn**4)]
    Ex = curve_from_roots(ee)
    rq = (sqclass(1), sqclass(ss), sqclass(ss))
    rkx = Ex.rank()
    gx = Ex.gens()
    T1x = Ex(ee[0], 0); T2x = Ex(ee[1], 0)
    im = set()
    ranges = [range(-1, 2)] * len(gx)
    import itertools
    for coef in itertools.product(*ranges):
        P0 = Ex(0)
        for i, cc in enumerate(coef):
            P0 = P0 + cc * gx[i]
        for c in [0, 1]:
            for d in [0, 1]:
                im.add(delta_of_point(P0 + c * T1x + d * T2x, ee))
    return ss, rkx, len(im), 2**(rkx + 2), rq, (rq in im)

for (mm, nn) in [(4, 3), (15, 8)]:
    ss, rkx, nim, full, rq, present = full_test(mm, nn)
    print("  (%d,%d): s=%s rank=%d, классов найдено %d из %d, требуемый %s -> %s"
          % (mm, nn, ss, rkx, nim, full, str(rq), "ПРИСУТСТВУЕТ" if present else "отсутствует"))
    # у обоих m^2+n^2 — квадрат, значит t=1 даёт точку C(Q)
    sq = ZZ(mm**2 + nn**2)
    ck("(%d,%d): m^2+n^2 = %s — квадрат, т.е. t=1 даёт ТОЧКУ на C(Q)" % (mm, nn, sq), sq.is_square())
    ck("(%d,%d): метод НЕ исключает (требуемый класс присутствует)" % (mm, nn), present)

ck("(15,1): m^2+n^2 = 226 НЕ квадрат (t=1 не даёт точки)", not ZZ(226).is_square())


# =====================================================================
hdr("J. ПРЯМОЙ ПОИСК точек: на C_{15,1} и на торсоре D")

# C: нужно 225q^2+p^2, 113(p^2+q^2), q^2+225p^2 — все квадраты
NB = 3000
hits = []
# сначала условие 113 | p^2+q^2 и (p^2+q^2)/113 квадрат (самое жёсткое)
cnt = 0
for q in range(1, NB + 1):
    q2 = q * q
    for p in range(0, NB + 1):
        v = p * p + q2
        if v % 113 != 0:
            continue
        w = v // 113
        if not ZZ(w).is_square():
            continue
        cnt += 1
        if ZZ(225 * q2 + p * p).is_square() and ZZ(q2 + 225 * p * p).is_square():
            hits.append((p, q))
print("  перебор t=p/q, 0<=p<=%d, 1<=q<=%d: кандидатов по F4: %d" % (NB, NB, cnt))
ck("рациональных точек на C_{15,1} малой высоты НЕ найдено", len(hits) == 0,
   "найдено: %s" % str(hits[:5]))

# Торсор D для класса (1,113,113):
#   X+25425 = r^2 ; X+5720625 = 113*z2^2 ; X+113 = 113*z3^2
#   => r^2 + 5695200 = 113 z2^2 ,  r^2 - 25312 = 113 z3^2
# Ищем рациональные r = a/c
print("\n  Торсор D (мой вывод):  r^2 - 113 z2^2 = -5695200 ,  r^2 - 113 z3^2 = 25312")
# r^2 = X-e1, 113 z2^2 = X-e2, 113 z3^2 = X-e3  =>
#   r^2 - 113 z2^2 = e2-e1 ,  r^2 - 113 z3^2 = e3-e1
ck("D получен из требуемого класса (проверка констант)",
   (E2 - E1) == -5695200 and (E3 - E1) == 25312,
   "e2-e1 = %s, e3-e1 = %s" % (E2 - E1, E3 - E1))

DB = 2500
dhits = []
for c in range(1, DB + 1):
    c2 = c * c
    A1 = 5695200 * c2
    A2 = 25312 * c2
    for a in range(0, DB + 1):
        a2 = a * a
        v1 = a2 + A1
        if v1 % 113 != 0:
            continue
        if not ZZ(v1 // 113).is_square():
            continue
        v2 = a2 - A2
        if v2 < 0 or v2 % 113 != 0:
            continue
        if ZZ(v2 // 113).is_square():
            dhits.append((a, c))
ck("рациональных точек на торсоре D малой высоты НЕ найдено", len(dhits) == 0,
   "найдено: %s" % str(dhits[:5]))

# Поиск точек E с требуемым классом среди найденных поиском точек
try:
    srch = E.point_search(14, rank_bound=2)
    bad = []
    for P in srch:
        for k in range(-6, 7):
            Q0 = k * P
            if Q0.is_zero():
                continue
            for c in [0, 1]:
                for d in [0, 1]:
                    Q1 = Q0 + c * T1 + d * T2
                    if Q1.is_zero():
                        continue
                    if delta_of_point(Q1, E_ROOTS) == REQ:
                        bad.append(Q1)
    ck("среди точек поиска (высота 14) нет точки с требуемым классом", len(bad) == 0,
       "нашлось %d" % len(bad))
except Exception as ex:
    unver("расширенный поиск точек E", str(ex))


# =====================================================================
hdr("K. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ C_{15,1} (нетривиальность исключения)")

def is_sq_Qp(x, p, prec=40):
    x = QQ(x)
    if x == 0:
        return False
    v = x.valuation(p)
    if v % 2 != 0:
        return False
    u = x / QQ(p)**v
    num = ZZ(u.numerator()); den = ZZ(u.denominator())
    if p == 2:
        r = (num * inverse_mod(den, 8)) % 8
        return r == 1
    else:
        r = (num * inverse_mod(den, p)) % p
        return r != 0 and kronecker(r, p) == 1

good_p = []
bad_p = []
for p in prime_range(200):
    ok = False
    for tv in [QQ(x) for x in range(-40, 41)] + [QQ(1) / QQ(k) for k in range(2, 40)]:
        if (is_sq_Qp(f0(tv), p) and is_sq_Qp(f4(tv), p) and is_sq_Qp(f8(tv), p)):
            ok = True
            good_p.append((p, tv))
            break
    if not ok:
        bad_p.append(p)
print("  простых с найденной локальной точкой: %d" % len(good_p))
print("  простые БЕЗ найденной точки в этом диапазоне t: %s" % bad_p)
ck("C(R) != пусто (t=0 даёт положительные F_i)", f0(0) > 0 and f4(0) > 0 and f8(0) > 0)
if bad_p:
    unver("локальная разрешимость при p in %s" % bad_p, "перебор t ограничен; не доказано отсутствие")
else:
    print("  => для всех p < 200 локальная точка предъявлена явно (перебор t, не доказательство для всех p)")
unver("локальная разрешимость для ВСЕХ p (Хассе–Вейль для p>=101)",
      "аргумент Codex через #C(F_p)>=p+1-10sqrt(p) мной не пересчитан символьно; "
      "проверено только явными t для p<200")


# =====================================================================
hdr("L. СВЕРКА С ЧИСЛАМИ CODEX (только здесь появляются его числа)")

ck("Codex s=113 совпадает", S == 113)
ck("Codex b=25425 совпадает", B == 25425)
ck("Codex e=(-25425,-5720625,-113) совпадает", E_ROOTS == [QQ(-25425), QQ(-5720625), QQ(-113)])
ck("Codex E: V^2=(X+25425)(X+5720625)(X+113) совпадает с моей",
   curve_from_roots([QQ(-25425), QQ(-5720625), QQ(-113)]).a_invariants() == E.a_invariants())
codexG = None
try:
    codexG = E(-75825, 146764800)
    ck("Codex G=(-75825,146764800) лежит на E", True)
    ck("Codex G — неторсионная", codexG.order() == oo)
    print("  delta(Codex G) моим счётом = %s" % (delta_of_point(codexG, E_ROOTS),))
    # эквивалентна ли она моему генератору по модулю 2E и кручения?
    ck("delta(Codex G) входит в мой образ", delta_of_point(codexG, E_ROOTS) in image)
except Exception as ex:
    ck("Codex G лежит на E", False, str(ex))

codex_table = [(1, 1, 1), (-1582, 226, -7), (-1, 1582, -1582), (1582, 7, 226),
               (-14, 2, -7), (113, 113, 1), (14, 791, 226), (-113, 14, -1582)]
codex_set = set(tuple(ZZ(x) for x in row) for row in codex_table)
ck("множество 8 классов Codex СОВПАДАЕТ с моим", codex_set == image)
if codex_set != image:
    print("    у Codex лишние: %s" % (codex_set - image))
    print("    у меня лишние:  %s" % (image - codex_set))
ck("Codex: требуемый (1,113,113) отсутствует — подтверждаю своим счётом", REQ not in image)
ck("ВНИМАНИЕ: (113,113,1) — присутствует (образ транспозиции e1<->e3)",
   (ZZ(113), ZZ(113), ZZ(1)) in image)

cminA = [0, -1, 0, -42422006041, 3362853259097305]
ck("минимальная модель Codex совпадает с моей", list(E.minimal_model().a_invariants()) == cminA,
   "моя: %s" % (E.minimal_model().a_invariants(),))


# =====================================================================
hdr("ИТОГ")
print("  ПРОВАЛОВ: %d" % len(FAIL))
for x in FAIL:
    print("    - " + x)
print("  НЕ ПРОВЕРЕНО: %d" % len(UNVER))
for x in UNVER:
    print("    - " + x)
print("")
if not FAIL:
    print("  Сломать заявление НЕ УДАЛОСЬ: модель, тождества, порядок корней,")
    print("  полнота образа при rank<=1 и аргумент про бесконечность выдержали.")
else:
    print("  ЕСТЬ ПРОВАЛЫ — см. список выше.")

# =====================================================================
#  ИТОГИ ВСЕХ РАУНДОВ АТАКИ (Claude, 2026-09-12)
#
#  Сопутствующие скрипты и логи в этом же каталоге:
#    check_kummer_15_1_ранг_КЛОД2.sage     / .log   (раунд 2)
#    check_kummer_15_1_модель_раунд3.sage  / .log   (раунд 3)
#    check_kummer_15_1_раунд4.sage         / .log   (раунд 4, секция P2 НЕ ДОСЧИТАНА)
#    check_kummer_15_1_Lряд_КЛОД.sage      / .log   (раунд 5)
#
#  ЧТО ПОДТВЕРЖДЕНО СВОИМ СЧЁТОМ:
#   1. Тождества X-e1=(mn)^2 F4, X-e2=s m^2 F0, X-e3=s n^2 F8 верны
#      СИМВОЛЬНО в Q[m,n,t] при e=(-b,-s m^4,-s n^4), b=s m^2 n^2.
#      Отсюда необходимый класс (1,s,s); для (15,1) это (1,113,113).
#   2. E неособа, disc = 2^36*3^4*5^4*7^6*113^8, три корня различны,
#      E(Q)_tors = (Z/2)^2, полное рациональное 2-кручение.
#   3. Мои восемь классов совпали с таблицей Codex посимвольно;
#      (1,113,113) среди них НЕТ; delta гомоморфен на всех 64 парах.
#   4. Вывод ИНВАРИАНТЕН относительно согласованной перестановки корней.
#   5. Над t=oo: U4^2 = s = 113 — не квадрат; все 8 геометрических точек
#      определены над Q(sqrt113); рациональных точек нет.
#   6. Аффинная модель C гладка (резултанты и дискриминанты != 0),
#      значит нормализация не добавляет рациональных точек.
#   7. Положительные контроли (4,3) и (15,8): метод НЕ исключает — верно,
#      у обоих m^2+n^2 квадрат, значит t=1 даёт реальную точку C(Q).
#   8. Поиск: точек C нет при t=p/q до 20000; точек торсора D нет.
#
#  ГДЕ РЕАЛЬНО ДЕРЖИТСЯ ДОКАЗАТЕЛЬСТВО (важно):
#   * Первый 2-спуск даёт ТОЛЬКО rank <= 3, и это так для ВСЕХ ЧЕТЫРЁХ
#     кривых изогенного класса (см. раунд 3, секция N2). Граница rank<=1
#     держится на втором спуске eclib и на спаривании Касселса–Тейта PARI.
#   * НЕЗАВИСИМОЕ подтверждение (раунд 5): w(E) = -1, значит L(E,1)=0;
#     L'(E,1) = 8.24611772899398609...  при оценке ошибки 3.1e-30.
#     Аналитический ранг = 1 => по Гросс–Загиру и Колывагину rank E(Q)=1
#     и Sha(E/Q) конечна. Этот путь НЕ использует спуск вообще.
#
#  ХРУПКОСТЬ, КОТОРУЮ НАДО ЗНАТЬ:
#   * (113,113,1) — ЕСТЬ в образе. Это в точности образ требуемого
#     (1,113,113) под транспозицией e1<->e3. Одна ошибка в метке корня —
#     и вывод молча переворачивается на «не исключено».
#   * У (11,4) та же картина: требуемый (1,274,274), а в образе (274,1,274),
#     т.е. транспозиция e1<->e2. Проверять метки корней надо в каждой паре.
# =====================================================================
