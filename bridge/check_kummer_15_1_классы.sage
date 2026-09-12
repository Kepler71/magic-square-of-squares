# -*- coding: utf-8 -*-
# check_kummer_15_1_классы.sage  — РАУНД 5. АТАКА на заявление Codex о (m,n)=(15,1).
#
# Заявление Codex: C_{15,1}(Q) = пусто.
#   s=113, b=25425, E: V^2=(X+25425)(X+5720625)(X+113),
#   требуемый класс (1,113,113), rank E = 1, 8 различных классов, требуемого нет.
#
# Установка: НЕ подтверждать, а ломать. Все числа считаю сам.
#
# ГЛАВНАЯ ЛОВУШКА, которую проверяю первым делом: в таблице Codex присутствует
# класс (113,113,1) — ПЕРЕСТАНОВКА требуемого (1,113,113). Значит весь вывод
# держится ТОЛЬКО на правильном порядке корней e1,e2,e3. Любая путаница в порядке
# переворачивает ответ. Порядок фиксируется символьными тождествами (раздел 2).
#
# Разделы:
#   1. Параметры s,b,e_i — своим счётом.
#   2. Символьная проверка X-e_i = c_i * F_j и вывод необходимого класса.
#      Проверка ловушки-перестановки.
#   3. Своя неторсионная точка: (a) грубый перебор по x, (b) E.gens() — сверка.
#   4. Свой групповой закон (чистый Python, Fraction) + своя delta.
#   5. Таблица 2^(r+2) классов, различность, замкнутость, гомоморфность.
#   6. Независимые границы ранга (Sage / mwrank / PARI) + 2-Selmer.
#   7. Бесконечность на C.
#   8. Прямой поиск рациональных точек на C (попытка контрпримера).
#   9. Сверка с числами Codex — В САМОМ КОНЦЕ.
#
# Автор: Claude (аудит-атака, раунд 5), 2026-09-12.

import sys, time
from fractions import Fraction
from math import isqrt

T0 = time.time()
def H(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)

FAILS = []
def ck(cond, msg):
    if cond:
        print("  [OK]   " + msg)
    else:
        print("  [FAIL] " + msg)
        FAILS.append(msg)
    return bool(cond)

# ==========================================================================
H("1. ПАРАМЕТРЫ ДЛЯ (m,n)=(15,1) — СЧИТАЮ САМ")
# ==========================================================================
m = Integer(15); n = Integer(1)
ck(gcd(m, n) == 1, "gcd(m,n)=1")
ck((m*m + n*n) % 2 == 0, "m^2+n^2 чётно => s целое")
s  = (m*m + n*n) / 2
b  = s * m*m * n*n
e1 = -b
e2 = -s * m**4
e3 = -s * n**4
print("  s = (m^2+n^2)/2 =", s)
print("  b = s*m^2*n^2   =", b)
print("  (e1,e2,e3) = (-b, -s*m^4, -s*n^4) =", (e1, e2, e3))
ck(s == 113, "s = 113")
ck(b == 25425, "b = 25425")
ck((e1, e2, e3) == (-25425, -5720625, -113), "(e1,e2,e3) = (-25425,-5720625,-113)")
ck(len(set([e1, e2, e3])) == 3, "корни попарно различны")
ck(not QQ(s).is_square(), "s=113 не квадрат в Q")

# ==========================================================================
H("2. СИМВОЛЬНЫЕ ТОЖДЕСТВА И НЕОБХОДИМЫЙ КЛАСС (порядок корней ФИКСИРОВАН)")
# ==========================================================================
R.<t> = QQ[]
F0 = m**2 + n**2 * t**2
F4 = s * (1 + t**2)
F8 = n**2 + m**2 * t**2
X  = b * t**2
print("  F0 =", F0, "   F4 =", F4, "   F8 =", F8)
print("  X  =", X)

c1 = (m*n)**2          # X - e1 = (mn)^2 * F4
c2 = s * m**2          # X - e2 = s*m^2 * F0
c3 = s * n**2          # X - e3 = s*n^2 * F8
ck((X - e1) - c1*F4 == 0, "X - e1 = (mn)^2 * F4 = %s * F4   тождественно" % c1)
ck((X - e2) - c2*F0 == 0, "X - e2 = s*m^2 * F0 = %s * F0   тождественно" % c2)
ck((X - e3) - c3*F8 == 0, "X - e3 = s*n^2 * F8 = %s * F8   тождественно" % c3)

# произведение должно быть квадратом * (u0u4u8)^2 => V = b*u0u4u8
ck(c1*c2*c3 == b**2, "c1*c2*c3 = %s = b^2  => V = b*u0*u4*u8 корректно" % (c1*c2*c3))

def sqfree(x):
    """квадратсвободный представитель класса x в Q*/Q*^2, со знаком"""
    x = QQ(x)
    if x == 0:
        raise ValueError("ноль не имеет квадратного класса")
    num = x.numerator(); den = x.denominator()
    z = num * den          # тот же класс, целое
    sgn = 1 if z > 0 else -1
    z = abs(z)
    r = 1
    for p, k in factor(z):
        if k % 2 == 1:
            r *= p
    return Integer(sgn * r)

req = (sqfree(c1), sqfree(c2), sqfree(c3))
print("  классы множителей c_i:", (c1, c2, c3), "->", req)
ck(req == (1, 113, 113), "НЕОБХОДИМЫЙ класс delta = (1, 113, 113)")

print()
print("  --- ЛОВУШКА-ПЕРЕСТАНОВКА ---")
print("  У Codex в таблице стоит класс (113,113,1) — это ПЕРЕСТАНОВКА требуемого.")
print("  Значит вывод целиком зависит от того, какой корень назван e1.")
print("  Проверяю ещё раз, что e1 — это именно -b (корень, отвечающий u4/F4):")
ck((X - (-b)) == c1*F4, "e1 = -b  <->  F4 (средняя клетка), множитель 225 = квадрат")
ck((X - (-s*m**4)) == c2*F0, "e2 = -s*m^4  <->  F0, множитель 25425 = 113*15^2 -> класс 113")
ck((X - (-s*n**4)) == c3*F8, "e3 = -s*n^4  <->  F8, множитель 113 -> класс 113")
print("  => единственная координата с ТРИВИАЛЬНЫМ классом — ПЕРВАЯ (при e1=-b).")

# Явная проверка на конкретном t, что порядок не перепутан:
tv = QQ(7)/3
Xv = b*tv**2
print("  контроль при t = %s:" % tv)
print("    X-e1 =", Xv-e1, " класс", sqfree(Xv-e1), " (должен быть класс F4 =", sqfree(s*(1+tv**2)), ")")
print("    X-e2 =", Xv-e2, " класс", sqfree(Xv-e2), " (должен быть класс 113*F0 =", sqfree(113*(225+tv**2)), ")")
print("    X-e3 =", Xv-e3, " класс", sqfree(Xv-e3), " (должен быть класс 113*F8 =", sqfree(113*(1+225*tv**2)), ")")
ck(sqfree(Xv-e1) == sqfree(s*(1+tv**2)), "числовой контроль коорд.1")
ck(sqfree(Xv-e2) == sqfree(113*(225+tv**2)), "числовой контроль коорд.2")
ck(sqfree(Xv-e3) == sqfree(113*(1+225*tv**2)), "числовой контроль коорд.3")

# ==========================================================================
H("3. МОЯ СОБСТВЕННАЯ НЕТОРСИОННАЯ ТОЧКА")
# ==========================================================================
E = EllipticCurve([0, -(e1+e2+e3), 0, e1*e2+e1*e3+e2*e3, -e1*e2*e3])
# y^2 = (x-e1)(x-e2)(x-e3)
print("  E: y^2 = (x-e1)(x-e2)(x-e3) =", E)
Rx.<xx> = QQ[]
fx = (xx-e1)*(xx-e2)*(xx-e3)
print("  f(x) =", fx)
ck(fx == xx**3 + Rx(E.a2())*xx**2 + Rx(E.a4())*xx + Rx(E.a6()),
   "модель Sage совпадает с (x-e1)(x-e2)(x-e3)")
ck(E.a1() == 0 and E.a3() == 0, "a1=a3=0 (короткая форма)")
print("  disc =", E.discriminant().factor())
print("  conductor =", E.conductor().factor())
tors = E.torsion_subgroup()
print("  torsion =", tors.invariants(), " порядок", tors.order())
ck(E.torsion_order() == 4, "E(Q)_tors порядка 4 => E(Q)[2]=(Z/2)^2, |T/2T|=4")

E1 = e1; E2 = e2; E3 = e3

print()
print("  (a) СОБСТВЕННЫЙ ГРУБЫЙ ПЕРЕБОР по целым x (без чужих чисел)")
sys.stdout.flush()
# сито: f(x) должно быть квадратичным вычетом по каждому малому модулю
mods = [9, 5, 7, 11, 13, 16]
MM = 1
for M in mods:
    MM *= M
iE1, iE2, iE3 = int(E1), int(E2), int(E3)
allow = {}
for M in mods:
    sq = set((i*i) % M for i in range(M))
    allow[M] = set(r for r in range(M)
                   if ((r - iE1)*(r - iE2)*(r - iE3)) % M in sq)
sieve = bytearray(MM)
for r in range(MM):
    good = 1
    for M in mods:
        if (r % M) not in allow[M]:
            good = 0; break
    sieve[r] = good
dens = sum(sieve) / MM
print("  сито по модулю %d: пропускает %.4f долю классов" % (MM, dens))
sys.stdout.flush()

def brute(lo, hi):
    found = []
    for x in range(lo, hi):
        if not sieve[x % MM]:
            continue
        v = (x - iE1)*(x - iE2)*(x - iE3)
        if v < 0:
            continue
        r = isqrt(v)
        if r*r == v:
            found.append((x, r))
    return found

# точки существуют при x >= e2=-5720625 в двух интервалах: [e2,e3] U [e1?..] —
# упорядочим корни по величине
roots_sorted = sorted([int(E1), int(E2), int(E3)])
print("  корни по возрастанию:", roots_sorted)
print("  f(x)>=0 на [%d, %d] и на [%d, +inf)" % (roots_sorted[0], roots_sorted[1], roots_sorted[2]))
mine = []
mine += brute(roots_sorted[0], roots_sorted[1] + 1)
print("  перебор левого интервала [%d,%d] завершён, найдено %d" % (roots_sorted[0], roots_sorted[1], len(mine)))
sys.stdout.flush()
mine += brute(roots_sorted[2], roots_sorted[2] + 4000000)
print("  перебор [%d, %d) завершён" % (roots_sorted[2], roots_sorted[2]+4000000))
sys.stdout.flush()
print("  ВСЕ найденные целые точки (x,|y|):")
for (x, y) in mine:
    print("     x =", x, "  y =", y, "  на кривой:", E.is_on_curve(QQ(x), QQ(y)))

nontors_brute = []
for (x, y) in mine:
    if y == 0:
        continue
    P = E(QQ(x), QQ(y))
    if P.order() == Infinity:
        nontors_brute.append(P)
print("  неторсионных из перебора:", len(nontors_brute))
for P in nontors_brute:
    print("     ", P, "  высота", RR(P.height()))

print()
print("  (b) Sage E.gens() — независимый источник той же группы")
sys.stdout.flush()
gens = E.gens()
print("  E.gens() =", gens)
Gsage = [g for g in gens if g.order() == Infinity]
print("  неторсионные образующие:", Gsage)

# выбираю СВОЮ рабочую точку: предпочитаю найденную перебором
if nontors_brute:
    G = nontors_brute[0]
    src = "собственный перебор по x"
else:
    G = Gsage[0]
    src = "E.gens()"
print("  РАБОЧАЯ ТОЧКА G =", G, "   источник:", src)
ck(E.is_on_curve(G[0], G[1]), "G лежит на E (прямая подстановка)")
ck(G.order() == Infinity, "G неторсионна")

T1 = E(E1, 0); T2 = E(E2, 0); T3 = E(E3, 0)
ck(T1 + T2 == T3, "T1+T2=T3 (полное 2-кручение)")

# ==========================================================================
H("4. СВОЙ ГРУППОВОЙ ЗАКОН (чистый Python, Fraction) И СВОЯ delta")
# ==========================================================================
_py = r'''
from fractions import Fraction as Fr
from math import isqrt

# ВНИМАНИЕ: препарсер Sage превращает литерал 0 в Integer(0); Fraction(Integer(0))
# создаёт битый объект (numerator становится методом). Поэтому все константы
# Fraction строятся ЗДЕСЬ, в чистом Python-блоке, и только через int().
ZERO = Fr(0)
ONE  = Fr(1)
def frac(a, b=1):
    return Fr(int(a), int(b))
def fr_of(q):
    """Sage Rational/Integer -> Fraction"""
    return Fr(int(q.numerator()), int(q.denominator()))

def add(P, Q, e1, e2, e3):
    """сложение на y^2 = (x-e1)(x-e2)(x-e3); P,Q = None (O) или (x,y) из Fr"""
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    A = -(e1+e2+e3)                 # y^2 = x^3 + A x^2 + B x + C
    B = e1*e2 + e1*e3 + e2*e3
    C = -e1*e2*e3
    if x1 == x2:
        if y1 != y2 or y1 == 0:
            return None
        lam = (3*x1*x1 + 2*A*x1 + B) / (2*y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    x3 = lam*lam - A - x1 - x2
    y3 = lam*(x1 - x3) - y1
    return (x3, y3)

def neg(P):
    if P is None: return None
    return (P[0], -P[1])

def on_curve(P, e1, e2, e3):
    if P is None: return True
    x, y = P
    return y*y == (x-e1)*(x-e2)*(x-e3)

def sqcls_int(z):
    """z != 0 рационально (Fr) -> целое num*den (тот же класс)"""
    return z.numerator * z.denominator

def same_class(a, b):
    """a,b целые != 0: равны ли квадратные классы -> a*b точный квадрат"""
    p = a*b
    if p < 0: return False
    r = isqrt(p)
    return r*r == p

def delta_raw(P, e1, e2, e3):
    """возвращает тройку ЦЕЛЫХ представителей класса (не приведённых)"""
    if P is None:
        return (1, 1, 1)
    x, y = P
    es = (e1, e2, e3)
    out = []
    for i in range(3):
        d = x - es[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = Fr(es[i] - es[j]) * Fr(es[i] - es[k])
        out.append(sqcls_int(d))
    return tuple(out)
'''
exec(_py, globals())

Fr = globals()['Fr']; frac = globals()['frac']; fr_of = globals()['fr_of']
ZERO = globals()['ZERO']

pe1, pe2, pe3 = frac(e1), frac(e2), frac(e3)
pG  = (fr_of(QQ(G[0])), fr_of(QQ(G[1])))
pT1 = (pe1, ZERO)
pT2 = (pe2, ZERO)
# самоконтроль: объекты Fraction настоящие (numerator — число, не метод)
ck(isinstance(pe1.numerator, int) and isinstance(pG[0].numerator, int),
   "Fraction-объекты корректны (numerator — int, не метод Sage)")

ck(on_curve(pG, pe1, pe2, pe3), "мой групповой закон: G на кривой")
ck(on_curve(pT1, pe1, pe2, pe3), "мой групповой закон: T1 на кривой")
ck(on_curve(pT2, pe1, pe2, pe3), "мой групповой закон: T2 на кривой")
ck(add(pT1, pT1, pe1, pe2, pe3) is None, "2*T1 = O")
ck(add(pT2, pT2, pe1, pe2, pe3) is None, "2*T2 = O")
ck(add(pT1, pT2, pe1, pe2, pe3) == (pe3, ZERO), "T1+T2 = T3")

# сверка моего сложения с Sage на нескольких примерах
def to_sage(P):
    return E(0) if P is None else E(QQ(P[0].numerator)/QQ(P[0].denominator),
                                    QQ(P[1].numerator)/QQ(P[1].denominator))
chk_ok = True
cur = pG
for k in range(1, 7):
    if k > 1:
        cur = add(cur, pG, pe1, pe2, pe3)
    if not on_curve(cur, pe1, pe2, pe3):
        chk_ok = False; break
    if to_sage(cur) != k*G:
        chk_ok = False; break
ck(chk_ok, "мой групповой закон совпадает с Sage на k*G, k=1..6")
for Pp, name in [(add(pG, pT1, pe1, pe2, pe3), "G+T1"),
                 (add(pG, pT2, pe1, pe2, pe3), "G+T2"),
                 (add(add(pG, pT1, pe1, pe2, pe3), pT2, pe1, pe2, pe3), "G+T1+T2")]:
    ck(on_curve(Pp, pe1, pe2, pe3), "мой групповой закон: %s на кривой" % name)

def delta_pretty(P):
    raw = delta_raw(P, pe1, pe2, pe3)
    return tuple(sqfree(Integer(r)) for r in raw)

# ==========================================================================
H("5. ТАБЛИЦА ВСЕХ 2^(r+2) КЛАССОВ delta")
# ==========================================================================
pts = {}
for a in range(2):
    for c in range(2):
        for d in range(2):
            P = None
            if a: P = add(P, pG, pe1, pe2, pe3)
            if c: P = add(P, pT1, pe1, pe2, pe3)
            if d: P = add(P, pT2, pe1, pe2, pe3)
            pts[(a, c, d)] = P

rows = []
print("  %-9s | %-28s | %s" % ("(a,c,d)", "точка (x,y)", "delta = (d1,d2,d3)"))
print("  " + "-" * 74)
table = {}
for key in sorted(pts):
    P = pts[key]
    dl = delta_pretty(P)
    table[key] = dl
    if P is None:
        ps = "O"
    else:
        xs = str(P[0]); ys = str(P[1])
        ps = "(%s, %s)" % (xs if len(xs) < 14 else xs[:11]+"...",
                           ys if len(ys) < 14 else ys[:11]+"...")
    print("  %-9s | %-28s | %s" % (str(key), ps, str(dl)))
    rows.append((key, dl))

print()
vals = [table[k] for k in sorted(table)]
distinct = (len(set(vals)) == 8)
ck(distinct, "все 8 классов ПОПАРНО РАЗЛИЧНЫ (иначе образ НЕ полон)")
if not distinct:
    from collections import Counter
    cc = Counter(vals)
    for v, k in cc.items():
        if k > 1:
            print("     ПОВТОР:", v, "встречается", k, "раз")

# замкнутость по покоординатному умножению классов (образ = подгруппа)
def mul(u, v):
    return tuple(sqfree(Integer(u[i]) * Integer(v[i])) for i in range(3))
S = set(vals)
closed = all(mul(u, v) in S for u in S for v in S)
ck(closed, "множество из 8 классов ЗАМКНУТО по умножению (подгруппа (Z/2)^3)")

# гомоморфность delta на всех 64 парах в моей группе <G,T1,T2>
def addp(P, Q): return add(P, Q, pe1, pe2, pe3)
hom_bad = []
for k1 in sorted(pts):
    for k2 in sorted(pts):
        Psum = addp(pts[k1], pts[k2])
        lhs = delta_pretty(Psum)
        rhs = mul(table[k1], table[k2])
        if lhs != rhs:
            hom_bad.append((k1, k2, lhs, rhs))
ck(len(hom_bad) == 0, "гомоморфность delta проверена на всех 64 парах")
for z in hom_bad[:5]:
    print("     НАРУШЕНИЕ:", z)

# дополнительная гомоморфность на случайных точках вне <G,T1,T2>
import random
rnd_bad = 0; rnd_tot = 0
for _ in range(40):
    a1 = int(random.randint(-6, 6)); c1r = int(random.randint(0, 1)); d1r = int(random.randint(0, 1))
    a2 = int(random.randint(-6, 6)); c2r = int(random.randint(0, 1)); d2r = int(random.randint(0, 1))
    def build(a, c, d):
        P = None
        if a != 0:
            Q = pG
            for _i in range(int(abs(a)) - 1):
                Q = addp(Q, pG)
            if a < 0: Q = neg(Q)
            P = Q
        if c: P = addp(P, pT1)
        if d: P = addp(P, pT2)
        return P
    A = build(a1, c1r, d1r); B = build(a2, c2r, d2r)
    Ssum = addp(A, B)
    try:
        lhs = delta_pretty(Ssum); rhs = mul(delta_pretty(A), delta_pretty(B))
    except Exception:
        continue
    rnd_tot += 1
    if lhs != rhs: rnd_bad += 1
ck(rnd_bad == 0, "гомоморфность на %d случайных парах aG+cT1+dT2, |a|<=6" % rnd_tot)

# 2-делимость: delta(P)=(1,1,1) <=> P in 2E(Q) — контроль ядра на удвоениях
ker_ok = True
for k in sorted(pts):
    P = pts[k]
    P2 = addp(P, P)
    if delta_pretty(P2) != (1, 1, 1):
        ker_ok = False
ck(ker_ok, "delta(2P) = (1,1,1) для всех 8 P (ядро содержит 2E(Q))")

print()
print("  ТРЕБУЕМЫЙ КЛАСС:", req)
in_image = (req in S)
ck(not in_image, "требуемый класс (1,113,113) ОТСУТСТВУЕТ в образе")
if in_image:
    print("     !!! КЛАСС ЕСТЬ В ОБРАЗЕ — ЗАЯВЛЕНИЕ CODEX ЛОЖНО !!!")

# проверка: а перестановки требуемого класса — в образе?
from itertools import permutations
print()
print("  Перестановки требуемого класса и их наличие в образе (диагностика ловушки):")
for pm in sorted(set(permutations(req))):
    print("     %-18s -> в образе: %s" % (str(pm), pm in S))

# ==========================================================================
H("6. НЕЗАВИСИМЫЕ ГРАНИЦЫ РАНГА И 2-SELMER")
# ==========================================================================
sys.stdout.flush()
Emin = E.minimal_model()
print("  минимальная модель:", Emin.ainvs())
try:
    rb = E.rank_bounds()
    print("  Sage E.rank_bounds() =", rb)
    ck(rb[1] <= 1, "верхняя граница ранга <= 1 (Sage rank_bounds)")
except Exception as ex:
    print("  rank_bounds исключение:", ex)

try:
    r_sage = E.rank()
    print("  Sage E.rank() =", r_sage)
    ck(r_sage == 1, "Sage: rank = 1")
except Exception as ex:
    print("  E.rank() исключение:", ex)

try:
    sel = E.selmer_rank()
    print("  Sage E.selmer_rank() (2-Selmer, без 2-кручения) =", sel)
except Exception as ex:
    print("  selmer_rank исключение:", ex)

try:
    import subprocess
    gpcode = ("E=ellinit([%s,%s,%s,%s,%s]); r=ellrank(E); print(r);"
              % tuple(str(a) for a in Emin.ainvs()))
    out = subprocess.run(["gp", "-q"], input=gpcode, capture_output=True, text=True, timeout=900)
    print("  PARI ellrank:", out.stdout.strip()[:400], out.stderr.strip()[:200])
except Exception as ex:
    print("  PARI недоступен/исключение:", ex)

try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    mw = mwrank_EllipticCurve(list(Emin.ainvs()))
    mw.two_descent(verbose=False)
    print("  mwrank rank =", mw.rank(), " bounds =", mw.rank_bound(), " selmer =", mw.selmer_rank())
    ck(mw.rank_bound() <= 1, "eclib/mwrank: верхняя граница ранга <= 1")
except Exception as ex:
    print("  mwrank исключение:", ex)

print()
print("  ЛОГИКА ЗАКРЫТИЯ: |E(Q)/2E(Q)| = 2^r * |E(Q)[2]| = 2^(r+2).")
print("  E(Q)[2] = (Z/2)^2 (проверено), значит при r<=1 имеем |E(Q)/2E(Q)| <= 8.")
print("  Предъявлено 8 различных классов => образ delta ИСЧЕРПАН.")

# ==========================================================================
H("7. БЕСКОНЕЧНОСТЬ НА C")
# ==========================================================================
print("  z = 1/t, U_i = u_i/t:  U0^2 = n^2 + m^2 z^2, U4^2 = s(z^2+1), U8^2 = m^2 + n^2 z^2")
print("  при z=0:  U0^2 = n^2 = 1,  U4^2 = s = 113,  U8^2 = m^2 = 225")
ck(not QQ(113).is_square(), "U4^2 = 113 не имеет рационального решения => точек над бесконечностью нет")
print("  контроль: u0^2/t^2 = (m^2+n^2 t^2)/t^2 = n^2 + m^2 (1/t)^2  ->  при 1/t=0 даёт n^2 =", n**2)
print("            u4^2/t^2 = s(1+t^2)/t^2 = s(1 + (1/t)^2)          ->  при 1/t=0 даёт s   =", s)
print("            u8^2/t^2 = (n^2+m^2 t^2)/t^2 = m^2 + n^2(1/t)^2   ->  при 1/t=0 даёт m^2 =", m**2)
print("  ВНИМАНИЕ: при z=0 первая координата даёт n^2 (квадрат), третья m^2 (квадрат),")
print("            и только СРЕДНЯЯ даёт s=113 — она и блокирует бесконечность.")

# ==========================================================================
H("8. ПРЯМОЙ ПОИСК РАЦИОНАЛЬНЫХ ТОЧЕК НА C (попытка контрпримера)")
# ==========================================================================
sys.stdout.flush()
# t = p/q, gcd(p,q)=1. Нужны квадраты:
#   225 q^2 + p^2,  113 (q^2+p^2),  q^2 + 225 p^2   (все / q^2)
# сильный фильтр: p^2+q^2 = 113 * square
from math import gcd as _gcd
BOUND = 20000
found_C = []
def issq(v):
    if v < 0: return False
    r = isqrt(v)
    return r*r == v
# p^2 + q^2 == 0 mod 113 => p == +-c*q mod 113, где c^2 == -1 mod 113
croots = [c for c in range(113) if (c*c + 1) % 113 == 0]
print("  корни c^2=-1 mod 113:", croots)
tested = 0
for q in range(1, BOUND + 1):
    q2 = q*q
    resid = set((c*q) % 113 for c in croots)
    for r0 in resid:
        p = r0
        while p <= BOUND:
            if p == 0 and q != 1:
                p += 113; continue
            if _gcd(p, q) == 1:
                tested += 1
                w = p*p + q2
                if issq(113*w) and issq(225*q2 + p*p) and issq(q2 + 225*p*p):
                    found_C.append((p, q))
            p += 113
print("  перебор t=p/q, 0<=p<=%d, 1<=q<=%d, gcd=1, 113|p^2+q^2 : проверено %d, найдено %d"
      % (BOUND, BOUND, tested, len(found_C)))
for pq in found_C[:20]:
    print("     t =", pq)
ck(len(found_C) == 0, "прямой перебор контрпримера НЕ нашёл (это НЕ доказательство отсутствия)")

# ==========================================================================
H("9. СВЕРКА С ЧИСЛАМИ CODEX (в самом конце)")
# ==========================================================================
codex_table = {
    (0,0,0): (1, 1, 1),
    (0,0,1): (-1582, 226, -7),
    (0,1,0): (-1, 1582, -1582),
    (0,1,1): (1582, 7, 226),
    (1,0,0): (-14, 2, -7),
    (1,0,1): (113, 113, 1),
    (1,1,0): (14, 791, 226),
    (1,1,1): (-113, 14, -1582),
}
codex_set = set(codex_table.values())
print("  множество классов Codex == моё множество :", codex_set == S)
ck(codex_set == S, "МНОЖЕСТВА классов совпадают (обязаны, если оба образа полны)")
if codex_set != S:
    print("     только у Codex:", sorted(codex_set - S))
    print("     только у меня :", sorted(S - codex_set))
Gcodex = E(QQ(-75825), QQ(146764800))
print("  точка Codex G_c =", Gcodex, " на кривой:", E.is_on_curve(Gcodex[0], Gcodex[1]))
print("  моя точка G     =", G)
print("  G_c == G ?", Gcodex == G, "   G_c == -G ?", Gcodex == -G)
try:
    print("  разность G_c - G =", Gcodex - G, " порядок:", (Gcodex - G).order())
except Exception as ex:
    print("  разность: исключение", ex)
print("  высоты: h(G) =", RR(G.height()), "  h(G_c) =", RR(Gcodex.height()))

# ==========================================================================
H("ИТОГ")
# ==========================================================================
print("  время: %.1f c" % (time.time() - T0))
if FAILS:
    print("  ПРОВАЛЕННЫЕ ПРОВЕРКИ (%d):" % len(FAILS))
    for x in FAILS:
        print("    -", x)
else:
    print("  все проверки пройдены")
