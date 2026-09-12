# -*- coding: utf-8 -*-
# СТАТУС: расчёт
# ДЛЯ: Claude / Codex / Grok / пользователь
# ИТОГ: три «дешёвых» способа рёберного случая (локальные препятствия, фактор ранга 0,
#       необходимый класс Куммера) в УГЛОВОМ случае неприменимы В ПРИНЦИПЕ; причина одна
#       и посчитана: у каждого углового слоя есть рациональная точка (вырожденный латинский квадрат).
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЕНО: девять клеток G2 выведены заново из магичности; сверены с bremner16/engine2.sage;
#            все отрицательные результаты посчитаны перебором, а не выведены рассуждением
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: вы хотите знать, можно ли закрыть угловой случай дёшево
#
# Запуск:  sage /home/kep/magicKube/bridge/corner_cheap.sage
#          CORNER_FAST=1 sage ...   — сокращённый прогон рангов (25 пар вместо 127)
#
# План:
#   §1  вывод девяти клеток УГЛОВОГО семейства G2 из магичности (символически, заново)
#   §2  ГЛАВНОЕ ТОЖДЕСТВО: четыре рациональных t, при которых все девять клеток — квадраты
#   §3  СПОСОБ 1. Локальные препятствия: перебор P^1(F_p), p < 200, все 127 пар + контроли
#   §4  СПОСОБ 2. Фактор ранга 0: эллиптические факторы углового семейства, кручение, ранг
#   §5  СПОСОБ 3. Необходимый класс Куммера: угловой аналог (1,s,s)
#   §6  Угловое СЕЧЕНИЕ (аналог теоремы об исключении отношений) — то же препятствие
#   §7  Что осталось и чего делать не надо

import functools, itertools, os, sys, time
print = functools.partial(print, flush=True)
FAST = os.environ.get('CORNER_FAST', '') == '1'   # CORNER_FAST=1 sage corner_cheap.sage
T00 = time.time()

def hdr(s):
    print(); print("=" * 100); print(s); print("=" * 100)

def sub(s):
    print(); print("-" * 100); print(s); print("-" * 100)

PAIRS = [(m, n) for m in range(2, 21) for n in range(1, m) if gcd(m, n) == 1]

# ============================================================================
# §1. ВЫВОД ДЕВЯТИ КЛЕТОК УГЛОВОГО СЕМЕЙСТВА G2
# ============================================================================
hdr("§1. Угловое семейство G2: девять клеток выведены заново из условий магичности")

print("""Разметка   a b c / d e f / g h i  =  c0..c8.
Пары: УГЛОВЫЕ (c0,c8) и (c2,c6); РЁБЕРНЫЕ (c1,c7) и (c3,c5).
G1 (рёберное, вся арифметика проекта): полные пары {1,7} и {3,5}; свободны 0,2,4,6,8 («икс»).
G2 (УГЛОВОЕ, предмет этой записки): полные пары {0,8} и {2,6}; свободны 1,3,4,5,7 («крест»).
Параметризация Сегре квадрики P^2+Q^2=R^2+S^2=2c:  P=mt+n, Q=m-nt, R=mt-n, S=m+nt.""")

R3 = PolynomialRing(QQ, ['m', 'n', 't'])
mm, nn, tt = R3.gens()
Fld = R3.fraction_field()

P_ = mm*tt + nn
Q_ = mm - nn*tt
R_ = mm*tt - nn
S_ = mm + nn*tt
# УГЛОВОЕ размещение: c0=P^2, c8=Q^2, c2=R^2, c6=S^2
g0, g8, g2, g6 = P_**2, Q_**2, R_**2, S_**2
g4 = Fld(g0 + g8) / 2

# оставшиеся четыре клетки (рёберные) — ИЗ МАГИЧНОСТИ, без готовых формул
# неизвестные (c1, c3, c5, c7)
Mmat = Matrix(Fld, [[1, 0, 0, 0],   # строка 0:      c0+c1+c2 = 3c4
                    [0, 0, 0, 1],   # строка 2:      c6+c7+c8 = 3c4
                    [0, 1, 0, 0],   # столбец 0:     c0+c3+c6 = 3c4
                    [0, 0, 1, 0],   # столбец 2:     c2+c5+c8 = 3c4
                    [1, 0, 0, 1],   # столбец 1:     c1+c4+c7 = 3c4
                    [0, 1, 1, 0]])  # строка 1:      c3+c4+c5 = 3c4
rhs = vector(Fld, [3*g4 - g0 - g2, 3*g4 - g6 - g8, 3*g4 - g0 - g6, 3*g4 - g2 - g8, 2*g4, 2*g4])
assert Mmat.rank() == 4
g1, g3, g5, g7 = Mmat.solve_right(rhs)

cellsG2 = [g0, g1, g2, g3, g4, g5, g6, g7, g8]
LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
ok_lines = all(Fld(sum(cellsG2[i] for i in L) - 3*g4) == 0 for L in LINES)
print("\n  все 8 линий = 3*центр (тождественно над Q(m,n,t)):", ok_lines)
assert ok_lines

s_sym = Fld(mm**2 + nn**2) / 2
alpha = Fld(3*nn**2 - mm**2) / 2      # старший коэффициент c1
beta  = Fld(3*mm**2 - nn**2) / 2      # свободный член c1
checks = [
    ("c4 = s(1+t^2)",                g4 - s_sym*(1 + tt**2)),
    ("c3 = c4 - 4mn t",              g3 - (s_sym*(1+tt**2) - 4*mm*nn*tt)),
    ("c5 = c4 + 4mn t",              g5 - (s_sym*(1+tt**2) + 4*mm*nn*tt)),
    ("c1 = alpha t^2 + beta",        g1 - (alpha*tt**2 + beta)),
    ("c7 = beta t^2 + alpha",        g7 - (beta*tt**2 + alpha)),
]
for nm, d in checks:
    print("  %-26s разность = %s" % (nm, Fld(d)))
    assert Fld(d) == 0
print("  где s=(m^2+n^2)/2,  alpha=(3n^2-m^2)/2,  beta=(3m^2-n^2)/2,  alpha+beta = m^2+n^2 = 2s")
print("  автоматические квадраты — ЧЕТЫРЕ УГЛА c0,c2,c6,c8; свободны c1,c3,c4,c5,c7 (крест).")

sub("§1a. Сверка с независимым кодом проекта bremner16/engine2.sage (cells_G('G2',m,n))")
load('/home/kep/magicKube/bremner16/engine2.sage')
agree = True
for (m0, n0) in [(3,2),(5,2),(11,4),(13,8),(16,5),(19,5),(20,19)]:
    cc, auto, free4 = cells_G('G2', m0, n0)
    # сравниваем как функции: значения в 13 точках однозначно определяют многочлен степени <= 2
    tvals = [QQ(k)/7 for k in range(-6, 7)]
    same = all(all(QQ(Fld(cellsG2[i]).subs({mm: m0, nn: n0, tt: tv})) == QQ(cc[i](tv))
                   for tv in tvals) for i in range(9))
    agree = agree and same and sorted(auto) == [0, 2, 6, 8] and sorted(free4) == [1, 3, 5, 7]
    print("  (m,n)=(%2d,%2d): совпадение всех девяти клеток: %s" % (m0, n0, same))
print("  ИТОГ сверки:", agree)
assert agree

sub("§1b. Для контраста — рёберное семейство G1 (то, на чём стоит вся арифметика проекта)")
e1_, e7_, e3_, e5_ = P_**2, Q_**2, R_**2, S_**2       # рёберное размещение
e4_ = Fld(e1_ + e7_)/2
rhsE = vector(Fld, [3*e4_ - e1_, 3*e4_ - e7_, 3*e4_ - e3_, 3*e4_ - e5_, 2*e4_, 2*e4_])
MmatE = Matrix(Fld, [[1,1,0,0],[0,0,1,1],[1,0,1,0],[0,1,0,1],[1,0,0,1],[0,1,1,0]])
e0_, e2_, e6_, e8_ = MmatE.solve_right(rhsE)
cellsG1 = [e0_, e1_, e2_, e3_, e4_, e5_, e6_, e7_, e8_]
assert all(Fld(sum(cellsG1[i] for i in L) - 3*e4_) == 0 for L in LINES)
print("  c0 = F0 = m^2+n^2 t^2 :", Fld(e0_ - (mm**2 + nn**2*tt**2)) == 0)
print("  c4 = F4 = s(1+t^2)   :", Fld(e4_ - s_sym*(1+tt**2)) == 0)
print("  c8 = F8 = n^2+m^2 t^2:", Fld(e8_ - (nn**2 + mm**2*tt**2)) == 0)
print("  c2 = L  = F4-2mn t   :", Fld(e2_ - (s_sym*(1+tt**2) - 2*mm*nn*tt)) == 0)
print("  c6 = U  = F4+2mn t   :", Fld(e6_ - (s_sym*(1+tt**2) + 2*mm*nn*tt)) == 0)
print("  свободны c0,c2,c4,c6,c8 — совпадает с RESULTS_FOR_CLAUDE_2026-09-12_G1.md §2/§4.")

# ============================================================================
# §2. ГЛАВНОЕ ТОЖДЕСТВО
# ============================================================================
hdr("§2. ГЛАВНОЕ ТОЖДЕСТВО: у КАЖДОГО углового слоя есть рациональная точка с девятью квадратами")

def sqnum(x):
    x = QQ(x)
    return x >= 0 and x.is_square()

sub("§2a. Символически над Q(m,n): четыре значения t дают девять квадратов")
tstars = [(Fld(mm-nn)/Fld(mm+nn), "t = (m-n)/(m+n)"),
          (-Fld(mm-nn)/Fld(mm+nn), "t = -(m-n)/(m+n)"),
          (Fld(mm+nn)/Fld(mm-nn), "t = (m+n)/(m-n)"),
          (-Fld(mm+nn)/Fld(mm-nn), "t = -(m+n)/(m-n)")]
for tv, nm in tstars:
    vals = [Fld(c).subs({tt: tv}) for c in cellsG2]
    allsq = all((Fld(v).numerator()*Fld(v).denominator()).is_square() for v in vals)
    print("  %-18s : все девять клеток — квадраты в Q(m,n):  %s" % (nm, allsq))
    assert allsq

v0 = [Fld(c).subs({tt: Fld(mm-nn)/Fld(mm+nn)}) for c in cellsG2]
print("""
  Явно, при t=(m-n)/(m+n) и после умножения на квадрат (m+n)^2 девять клеток равны

        (m^2+n^2)^2     (m^2+2mn-n^2)^2  (m^2-2mn-n^2)^2
        (m^2-2mn-n^2)^2 (m^2+n^2)^2      (m^2+2mn-n^2)^2
        (m^2+2mn-n^2)^2 (m^2-2mn-n^2)^2  (m^2+n^2)^2

  — ЛАТИНСКИЙ вырожденный магический квадрат из девяти квадратов, но лишь ТРЁХ различных
  значений (сумма 3(m^2+n^2)^2 сходится, потому что (A-B)^2+(A+B)^2 = 2(m^2+n^2)^2
  при A=m^2-n^2, B=2mn).  Это классический тривиальный пример, не новый квадрат.""")

lat = [(mm**2+nn**2)**2, (mm**2+2*mm*nn-nn**2)**2, (mm**2-2*mm*nn-nn**2)**2]
den2 = (mm+nn)**2
print("  проверка явной формы (умножение на (m+n)^2):",
      all(Fld(v0[i]*den2 - lat[j]) == 0 for i, j in
          [(0,0),(1,1),(2,2),(3,2),(4,0),(5,1),(6,1),(7,2),(8,0)]))
X_ = mm**2 - nn**2; Y_ = 2*mm*nn
print("  тождество (A-B)^2+(A+B)^2 = 2(m^2+n^2)^2:",
      R3((X_-Y_)**2 + (X_+Y_)**2 - 2*(mm**2+nn**2)**2) == 0)

sub("§2b. То же численно для всех 127 пар 1<=n<m<=20, gcd=1 (перебор, не рассуждение)")
bad = []
ndistinct = set()
for (m0, n0) in PAIRS:
    t0 = QQ(m0-n0)/QQ(m0+n0)
    vals = [QQ(Fld(c).subs({mm: m0, nn: n0, tt: t0})) for c in cellsG2]
    if not all(sqnum(v) for v in vals):
        bad.append((m0, n0))
    ndistinct.add(len(set(vals)))
print("  пар, где НЕ все девять клеток квадраты при t=(m-n)/(m+n):", len(bad), bad)
assert not bad
print("  => 127 из 127.  Число различных значений в квадрате (по всем парам):", sorted(ndistinct))

sub("§2b'. Поиск НЕвырожденных точек в угловом семействе (t=p/q, |p|,|q| <= 60)")
import math
def g2_cells_int(m0, n0, a, b):
    S = m0*m0 + n0*n0; aa = a*a; bb = b*b; ab = a*b
    C4 = 2*S*(aa+bb); W = m0*m0*aa + n0*n0*bb
    return [(m0*a+n0*b)**2, 3*C4-8*W, (m0*a-n0*b)**2, C4-16*m0*n0*ab, C4,
            C4+16*m0*n0*ab, (m0*b+n0*a)**2, 8*W-C4, (m0*b-n0*a)**2]
def is_sq_int(x):
    if x < 0:
        return False
    r = math.isqrt(x)
    return r*r == x
HB = 60
tcand = [(p, q) for q in range(1, HB+1) for p in range(-HB, HB+1) if gcd(abs(p), q) == 1]
t0_ = time.time()
best_global = 0; found7 = []; hist = {}
for (m0, n0) in PAIRS:
    best = 0
    for (a, b) in tcand:
        cs = g2_cells_int(m0, n0, a, b)
        k = sum(1 for v in cs if is_sq_int(v))
        if k >= 7:
            nd = len(set(cs))
            if k > best:
                best = k
            if nd >= 7:
                found7.append((m0, n0, a, b, k, nd))
    hist[best] = hist.get(best, 0) + 1
    best_global = max(best_global, best)
print("  перебрано t-кандидатов на пару: %d; время %.0f с" % (len(tcand), time.time()-t0_))
print("  распределение «максимум квадратных клеток при 7 и более» по парам:", hist)
print("  найдено точек с >= 7 квадратными И >= 7 РАЗЛИЧНЫМИ значениями:", len(found7), found7[:5])
print("""  [проверено численно, ограниченный поиск] В этом ящике все точки с семью и более
  квадратными клетками — вырожденные (три различных значения). Отсутствие находки НЕ есть
  доказательство отсутствия: ящик мал (|p|,|q|<=60) и семейство лишь одно из трёх.""")

sub("§2c. У РЁБЕРНОГО семейства G1 такой точки почти никогда НЕТ — вот где ломается симметрия")
print("""  В G1 вырожденные точки с девятью квадратами тоже бывают, но лишь при УСЛОВИИ на (m,n):
    t = +-1 : c0=c4=c8=m^2+n^2, c2=(m-n)^2, c6=(m+n)^2 => нужен квадратный m^2+n^2;
    t = 0   : c0=m^2, c8=n^2, c2=c4=c6=s          => нужен квадратный s=(m^2+n^2)/2.
  В G2 никакого условия нет — латинская точка есть тождественно.""")
pyth = [(m0, n0) for (m0, n0) in PAIRS if ZZ(m0**2 + n0**2).is_square()]
halfsq = [(m0, n0) for (m0, n0) in PAIRS if QQ(m0**2 + n0**2)/2 in ZZ
          and ZZ(QQ(m0**2+n0**2)/2).is_square()]
haveG1 = []
for (m0, n0) in PAIRS:
    for t0 in [QQ(1), QQ(-1), QQ(0), QQ(m0-n0)/QQ(m0+n0), QQ(m0+n0)/QQ(m0-n0),
               -QQ(m0-n0)/QQ(m0+n0), -QQ(m0+n0)/QQ(m0-n0)]:
        vals = [QQ(Fld(c).subs({mm: m0, nn: n0, tt: t0})) for c in cellsG1]
        if all(sqnum(v) for v in vals):
            haveG1.append((m0, n0)); break
print("  m^2+n^2 — квадрат у пар:", pyth)
print("  s=(m^2+n^2)/2 — квадрат у пар:", halfsq)
print("  пар, у которых G1 даёт девять квадратов в одной из семи вырожденных точек:",
      len(haveG1), haveG1)
print("  объединение двух критериев совпадает со списком:",
      sorted(set(pyth) | set(halfsq)) == sorted(haveG1))
G1_HAVE_POINT = set(haveG1)
print("""  Итог сравнения [доказано]:
    G2 (углы):  девять квадратов есть при ЛЮБЫХ (m,n) — 127/127;
    G1 (рёбра): только при m^2+n^2=кв. или (m^2+n^2)/2=кв. — %d/127.
  Именно поэтому в G1 бывают локальные препятствия и факторы ранга 0, а в G2 их быть не может.""" % len(haveG1))

# ============================================================================
# §3. СПОСОБ 1 — ЛОКАЛЬНЫЕ ПРЕПЯТСТВИЯ
# ============================================================================
hdr("§3. СПОСОБ 1. Локальные препятствия: перебор P^1(F_p) для всех нечётных p < 200")

def g2_mod(m0, n0, a, b, p):
    """девять клеток G2 mod p, однородизованные умножением на квадрат;
       квадратность значения не зависит от умножения на квадрат."""
    S = (m0*m0 + n0*n0) % p
    aa = (a*a) % p; bb = (b*b) % p; ab = (a*b) % p
    C0 = pow((m0*a + n0*b) % p, 2, p)
    C2 = pow((m0*a - n0*b) % p, 2, p)
    C6 = pow((m0*b + n0*a) % p, 2, p)
    C8 = pow((m0*b - n0*a) % p, 2, p)
    C4 = (2*S*(aa + bb)) % p                                  # = 4*c4*b^2
    C3 = (C4 - 16*m0*n0*ab) % p
    C5 = (C4 + 16*m0*n0*ab) % p
    W  = (m0*m0*aa + n0*n0*bb) % p
    C1 = (3*C4 - 8*W) % p
    C7 = (8*W - C4) % p
    return [C0, C1, C2, C3, C4, C5, C6, C7, C8]

def g1_mod(m0, n0, a, b, p):
    S = (m0*m0 + n0*n0) % p
    aa = (a*a) % p; bb = (b*b) % p; ab = (a*b) % p
    F0 = (4*(m0*m0*bb + n0*n0*aa)) % p
    F4 = (2*S*(aa + bb)) % p
    F8 = (4*(n0*n0*bb + m0*m0*aa)) % p
    L  = (F4 - 8*m0*n0*ab) % p
    U  = (F4 + 8*m0*n0*ab) % p
    C1 = pow((m0*a + n0*b) % p, 2, p); C3 = pow((m0*a - n0*b) % p, 2, p)
    C5 = pow((m0*b + n0*a) % p, 2, p); C7 = pow((m0*b - n0*a) % p, 2, p)
    return [F0, C1, L, C3, F4, C5, U, C7, F8]

QRC = {}
def qrset(p):
    if p not in QRC:
        QRC[p] = set((x*x) % p for x in range(p))
    return QRC[p]

def survivors(cellfun, m0, n0, p, free_idx):
    QR = qrset(p); cnt = 0; wit = None
    for (a, b) in [(x, 1) for x in range(p)] + [(1, 0)]:
        cs = cellfun(m0, n0, a, b, p)
        if all(cs[i] in QR for i in free_idx):
            cnt += 1
            if wit is None:
                wit = (a, b)
    return cnt, wit

FREE_G2 = (1, 3, 4, 5, 7)
FREE_G1 = (0, 2, 4, 6, 8)
PRIMES200 = [int(p) for p in prime_range(3, 200)]

sub("§3a. КОНТРОЛЬ кода на рёберном случае: (13,8) mod 17 и (16,5) mod 41 должны давать 0 точек")
for (m0, n0, p) in [(13, 8, 17), (16, 5, 41)]:
    c, w = survivors(g1_mod, m0, n0, p, FREE_G1)
    print("  G1 (m,n)=(%2d,%2d), p=%2d: выживших точек P^1(F_p) = %d   %s"
          % (m0, n0, p, c, "ПУСТО — известное препятствие воспроизведено" if c == 0 else "!!! НЕ воспроизведено"))
    assert c == 0
c, _ = survivors(g1_mod, 15, 8, 17, FREE_G1)
print("  контроль-наоборот: G1 (15,8) p=17 выживших =", c, "(> 0, там есть рациональная точка t=1)")
assert c > 0

sub("§3b. Полный скан G1 (рёберное) — эталон: сколько пар ловится локально")
t0_ = time.time(); g1_killed = {}
for (m0, n0) in PAIRS:
    for p in PRIMES200:
        if m0 % p == 0 and n0 % p == 0:
            continue
        c, _ = survivors(g1_mod, m0, n0, p, FREE_G1)
        if c == 0:
            g1_killed[(m0, n0)] = p
            break
print("  G1: локально исключено пар из 127: %d   (%.0f c)" % (len(g1_killed), time.time()-t0_))
byp = {}
for k, p in g1_killed.items():
    byp.setdefault(p, []).append(k)
for p in sorted(byp):
    print("     первое препятствие при p=%3d : %2d пар" % (p, len(byp[p])))
rest = sorted(set(PAIRS) - set(g1_killed))
print("  остались локально-разрешимыми (нечётные p<200): %d пар" % len(rest))
print("   ", rest if len(rest) <= 40 else rest[:40])
print("  ТОЧНОЕ СОВПАДЕНИЕ со списком пар, у которых есть рациональная точка (§2c):",
      sorted(rest) == sorted(G1_HAVE_POINT), " —", sorted(G1_HAVE_POINT))
assert sorted(rest) == sorted(G1_HAVE_POINT)
print("""  [проверено перебором] В G1 локальное препятствие отсутствует РОВНО у тех пар, у которых
  есть рациональная точка. Это и подтверждает механизм: препятствие исчезает тогда и только
  тогда (в этом диапазоне), когда точка есть. В G2 точка есть всегда => препятствий нет нигде.""")

sub("§3b'. Семь пар, которые в проекте остались открытыми: их пустые p (девятиклеточная система)")
SEVEN = [(11,4),(13,8),(15,1),(15,8),(16,5),(19,5),(19,16)]
for (m0, n0) in SEVEN:
    empt = [p for p in PRIMES200 if (m0 % p or n0 % p) and survivors(g1_mod, m0, n0, p, FREE_G1)[0] == 0]
    print("  (%2d,%2d): пустые p < 200: %s" % (m0, n0, empt))
print("""  Контроль: для (13,8) проект получил [17,47,73,79,103] при p<=113 (local_13_8_p17_прямая.log §7)
  и для (16,5) — первое препятствие 41. Совпало.
  ПОБОЧНАЯ НАХОДКА (относится к РЁБЕРНОМУ случаю, не к угловому): те же дешёвые сравнения
  закрывают также (11,4), (15,1), (19,5), (19,16) — но ТОЛЬКО для девятиклеточной системы.
  Это слабее, чем доказанное в проекте C(Q)=пусто для ТРЁХклеточной кривой, и не заменяет его;
  зато для вопроса «есть ли в этом слое магический квадрат из девяти квадратов» достаточно.
  Единственная из семи, где девятиклеточная система локально разрешима всюду, — (15,8):
  у неё 15^2+8^2=17^2, то есть ровно та самая рациональная точка при t=+-1.""")

sub("§3c. ГЛАВНЫЙ СКАН G2 (угловое): те же p < 200, те же 127 пар")
t0_ = time.time(); g2_killed = {}; minsurv = {}
for (m0, n0) in PAIRS:
    best = None
    for p in PRIMES200:
        if m0 % p == 0 and n0 % p == 0:
            continue
        c, _ = survivors(g2_mod, m0, n0, p, FREE_G2)
        if best is None or c < best[0]:
            best = (c, p)
        if c == 0:
            g2_killed[(m0, n0)] = p
            break
    minsurv[(m0, n0)] = best
print("  G2: локально исключено пар из 127: %d  %s   (%.0f c)"
      % (len(g2_killed), g2_killed, time.time()-t0_))
mv = min(v[0] for v in minsurv.values())
print("  минимум по всем (m,n) и всем p<200 числа выживших точек P^1(F_p):", mv)
print("  где достигается (пара -> (выживших, p)):",
      [(k, v) for k, v in minsurv.items() if v[0] == mv][:10])
print("""
  [доказано перебором] Ни одна из 127 угловых пар не исключается локально ни при одном
  нечётном p < 200. И это не случайность: по §2 у каждой пары есть РАЦИОНАЛЬНАЯ точка
  t=(m-n)/(m+n); она редуцируется в точку над F_p при любом p хорошей редукции, поэтому
  выживших всегда >= 1. Локальное препятствие невозможно НИ ПРИ КАКОМ p — в принципе,
  а не «до 200».""")

sub("§3d. Прямая проверка: редукция латинской точки выживает при каждом p")
okred = True; checked = 0
for (m0, n0) in PAIRS:
    a0, b0 = m0 - n0, m0 + n0
    for p in PRIMES200:
        if a0 % p == 0 and b0 % p == 0:
            continue
        cs = g2_mod(m0, n0, a0, b0, p); QR = qrset(p)
        checked += 1
        if not all(v in QR for v in cs):
            okred = False
            print("   !!! (m,n)=(%d,%d) p=%d — редукция латинской точки не квадратна" % (m0, n0, p))
print("  проверено пар-простых: %d; редукция даёт девять квадратов в F_p всюду: %s" % (checked, okred))
assert okred

sub("§3e. p=2 и вещественное место")
def surv_Z2(cellfun_int, m0, n0, k, free_idx):
    Mod = 2**k; cnt = 0
    for a in range(Mod):
        for b in range(Mod):
            if a % 2 == 0 and b % 2 == 0:
                continue
            cs = cellfun_int(m0, n0, a, b)
            good = True
            for i in free_idx:
                x = ZZ(cs[i])
                if x == 0:
                    continue
                v = x.valuation(2)
                if v + 3 <= k:
                    if v % 2 or (ZZ(x >> v) % 8) != 1:
                        good = False; break
            if good:
                cnt += 1
    return cnt
def g2_int(m0, n0, a, b):
    S = m0*m0 + n0*n0; aa = a*a; bb = b*b; ab = a*b
    C4 = 2*S*(aa+bb); W = m0*m0*aa + n0*n0*bb
    return [(m0*a+n0*b)**2, 3*C4-8*W, (m0*a-n0*b)**2, C4-16*m0*n0*ab, C4,
            C4+16*m0*n0*ab, (m0*b+n0*a)**2, 8*W-C4, (m0*b-n0*a)**2]
def g1_int(m0, n0, a, b):
    S = m0*m0 + n0*n0; aa = a*a; bb = b*b; ab = a*b
    F4 = 2*S*(aa+bb)
    return [4*(m0*m0*bb+n0*n0*aa), (m0*a+n0*b)**2, F4-8*m0*n0*ab, (m0*a-n0*b)**2, F4,
            (m0*b+n0*a)**2, F4+8*m0*n0*ab, (m0*b-n0*a)**2, 4*(n0*n0*bb+m0*m0*aa)]
for (m0, n0) in [(13, 8), (16, 5), (11, 4), (5, 2), (3, 2)]:
    print("  (m,n)=(%2d,%2d): 2-адически выживших mod 2^6:  G2 = %5d,  G1 = %5d"
          % (m0, n0, surv_Z2(g2_int, m0, n0, 6, FREE_G2), surv_Z2(g1_int, m0, n0, 6, FREE_G1)))
print("\n  вещественное место. В G1 все пять свободных клеток положительно определены,")
print("  в G2 — НЕТ: c1,c7 имеют старший коэффициент (3n^2-m^2)/2 (может быть < 0),")
print("  а c3,c5 имеют вещественные корни при 4mn > m^2+n^2. Но область положительности")
print("  всегда непуста, потому что содержит латинскую точку:")
for (m0, n0) in [(5, 2), (13, 8), (16, 5), (20, 19)]:
    t0 = QQ(m0-n0)/QQ(m0+n0)
    vals = [QQ(Fld(c).subs({mm: m0, nn: n0, tt: t0})) for c in cellsG2]
    a_i = QQ(3*n0**2 - m0**2)/2
    print("    (%2d,%2d): alpha=(3n^2-m^2)/2 = %-8s все клетки при t* положительны: %s"
          % (m0, n0, a_i, all(v > 0 for v in vals)))

# ============================================================================
# §4. СПОСОБ 2 — ФАКТОР РАНГА 0
# ============================================================================
hdr("§4. СПОСОБ 2. Эллиптический фактор ранга 0 и его кручение")

print("""В рёберном случае берут квартику z^2 = F0*F8 (произведение свободной УГЛОВОЙ пары),
её якобиан E_B: Y^2 = X(X-(m^2-n^2)^2)(X-(m^2+n^2)^2), кручение ВСЕГДА Z/2 x Z/4,
и rank=0 => t in {0,+-1,oo} => в слое невырожденного квадрата нет. Так закрылись 62 из 127 пар.

Угловой аналог: свободные клетки G2 — крест, в нём две РЁБЕРНЫЕ пары (c3,c5) и (c1,c7),
значит две квартики:
   B_df: z^2 = 4*c3*c5 = S^2(1+t^2)^2 - 64 m^2 n^2 t^2,      S = m^2+n^2
   B_bh: z^2 = 4*c1*c7 = (A t^2 + B)(B t^2 + A),  A = 3n^2-m^2, B = 3m^2-n^2""")

RT = PolynomialRing(QQ, 't'); T = RT.gen()

def Eroots(r1, r2, r3):
    return EllipticCurve([0, -(r1+r2+r3), 0, r1*r2+r1*r3+r2*r3, -r1*r2*r3])

def jac_biquad(a, c, e):
    """якобиан квартики z^2 = a t^4 + c t^2 + e:  Y^2 = X(X^2 - 2cX + (c^2-4ae)).
       ЗНАК ВАЖЕН: с +2c получается квадратичный твист на -1 (проверено счётом точек в §4a')."""
    return EllipticCurve([0, -2*c, 0, c**2 - 4*a*e, 0])

def corner_E_df(m0, n0):
    """якобиан B_df, масштаб X -> X/4: корни 0, S^2-16m^2n^2, -16m^2n^2"""
    S = m0**2 + n0**2
    return Eroots(0, S**2 - 16*m0**2*n0**2, -16*m0**2*n0**2)

def corner_E_bh(m0, n0):
    """якобиан B_bh, масштаб X -> X/4: корни 0, S^2, 4(m^2-n^2)^2"""
    S = m0**2 + n0**2
    return Eroots(0, S**2, 4*(m0**2 - n0**2)**2)

def edge_E_B(m0, n0):
    """якобиан z^2=F0*F8: корни 0, (m^2-n^2)^2, (m^2+n^2)^2 — в точности E_B проекта"""
    return Eroots(0, (m0**2 - n0**2)**2, (m0**2 + n0**2)**2)

sub("§4a. Сверка явных моделей с якобианами квартик")
for (m0, n0) in [(3, 2), (5, 2), (13, 8), (16, 5), (19, 5)]:
    S = m0**2 + n0**2
    E1 = jac_biquad(S**2, 2*S**2 - 64*m0**2*n0**2, S**2)
    E2 = corner_E_df(m0, n0)
    A_, B_ = 3*n0**2 - m0**2, 3*m0**2 - n0**2
    E3 = jac_biquad(A_*B_, A_**2 + B_**2, A_*B_)
    E4 = corner_E_bh(m0, n0)
    E5 = jac_biquad(m0**2*n0**2, m0**4 + n0**4, m0**2*n0**2)
    E6 = edge_E_B(m0, n0)
    print("  (%2d,%2d): E_df==Jac(B_df): %-5s ; E_bh==Jac(B_bh): %-5s ; E_B==Jac(F0F8): %-5s"
          % (m0, n0, E1.is_isomorphic(E2), E3.is_isomorphic(E4), E5.is_isomorphic(E6)))
    assert E1.is_isomorphic(E2) and E3.is_isomorphic(E4) and E5.is_isomorphic(E6)
print("""  Последний столбец — КОНТРОЛЬ: та же формула на рёберной квартике z^2=F0*F8 обязана
  дать в точности E_B: Y^2=X(X-(m^2-n^2)^2)(X-(m^2+n^2)^2) из REPORT_G1_SIX_KUMMER. Даёт.
  корни E_df: 0, S^2-16m^2n^2 = (m-n)^2(m^2+4mn+n^2), -16m^2n^2
  корни E_bh: 0, S^2, 4(m^2-n^2)^2""")

sub("§4a'. Контроль знака: #C(F_p) должно равняться #Jac(F_p) (Ленг), а не #твиста")
RTq = PolynomialRing(QQ, 'x'); xq = RTq.gen()
def count_quartic_Fp(q, p):
    K = GF(p); qq = q.change_ring(K); n = 0
    for x in K:
        v = qq(x)
        if v == 0:
            n += 1
        elif v.is_square():
            n += 2
    lc = K(q.leading_coefficient())
    if lc != 0 and lc.is_square():
        n += 2
    return n
for (m0, n0) in [(2, 1), (5, 2), (13, 8)]:
    S = m0**2 + n0**2
    q = S**2*xq**4 + (2*S**2 - 64*m0**2*n0**2)*xq**2 + S**2
    Egood = corner_E_df(m0, n0)
    Ebad = EllipticCurve([0, 2*(2*S**2 - 64*m0**2*n0**2), 0,
                          (2*S**2-64*m0**2*n0**2)**2 - 4*S**4, 0])   # твист на -1
    row = []
    for p in [7, 11, 13, 19, 23, 29, 31, 37]:
        if ZZ(q.discriminant()) % p == 0 or ZZ(Egood.discriminant()) % p == 0:
            continue
        row.append((p, count_quartic_Fp(q, p),
                    Egood.change_ring(GF(p)).cardinality(),
                    Ebad.change_ring(GF(p)).cardinality()))
    ok = all(r[1] == r[2] for r in row) and not all(r[1] == r[3] for r in row)
    print("  (%2d,%2d): (p, #C, #E_df, #твист) = %s" % (m0, n0, row[:5]))
    print("           #C == #E_df всюду: %s ; знак якобиана выбран верно: %s"
          % (all(r[1] == r[2] for r in row), ok))
    assert all(r[1] == r[2] for r in row)

sub("§4b. Кручение: угловые факторы против рёберного эталона (все 127 пар)")
def dist(d):
    out = {}
    for v in d.values():
        out[v] = out.get(v, 0) + 1
    return out
tors_df = {}; tors_bh = {}; tors_edge = {}
for (m0, n0) in PAIRS:
    tors_df[(m0,n0)] = tuple(corner_E_df(m0, n0).torsion_subgroup().invariants())
    tors_bh[(m0,n0)] = tuple(corner_E_bh(m0, n0).torsion_subgroup().invariants())
    tors_edge[(m0,n0)] = tuple(edge_E_B(m0, n0).torsion_subgroup().invariants())
print("  кручение E_B  (РЁБЕРНОЕ, эталон)  :", dist(tors_edge))
print("  кручение E_df (угловое, пара c3,c5):", dist(tors_df))
print("  кручение E_bh (угловое, пара c1,c7):", dist(tors_bh))
print("""
  ЭТО ВАЖНОЕ ОТЛИЧИЕ, и это прямой ответ на вопрос задания «каково кручение».
  В рёберном случае Z/2 x Z/4 УНИВЕРСАЛЬНО: разности корней 0, (m^2-n^2)^2, (m^2+n^2)^2
  равны (m^2-n^2)^2, (m^2+n^2)^2 и (2mn)^2 — все три квадраты.
  В угловом Z/2 x Z/4 НЕ универсально, критерии (y^2=X(X-r1)(X-r2): точка (r,0) делится
  на 2 <=> обе разности r-0, r-r' — квадраты):
    E_df, корни 0, r_b=S^2-16m^2n^2, r_c=-16m^2n^2:
        (0,0)   делится <=> -r_b = 16m^2n^2-S^2 — квадрат (т.к. -r_c=(4mn)^2 всегда квадрат)
        (r_b,0) делится <=> r_b = S^2-16m^2n^2 — квадрат (т.к. r_b-r_c=S^2 всегда квадрат)
        (r_c,0) не делится никогда (r_c<0)
      => Z/2 x Z/4  <=>  |S^2 - 16 m^2 n^2| = |(m^2+n^2-4mn)(m^2+n^2+4mn)| — квадрат.
    E_bh, корни 0, S^2, 4(m^2-n^2)^2:
      => Z/2 x Z/4  <=>  +-A*B — квадрат, A=3n^2-m^2, B=3m^2-n^2.
  Проверка этих критериев перебором по всем 127 парам:""")
crit_df = [(m0, n0) for (m0, n0) in PAIRS
           if ZZ(abs((m0**2+n0**2)**2 - 16*m0**2*n0**2)).is_square()]
crit_bh = [(m0, n0) for (m0, n0) in PAIRS
           if (lambda A_, B_: (A_*B_ > 0 and ZZ(A_*B_).is_square())
               or (-A_*B_ > 0 and ZZ(-A_*B_).is_square()))(3*n0**2-m0**2, 3*m0**2-n0**2)]
print("    |S^2-16m^2n^2| — квадрат у пар:", crit_df)
print("    пары с кручением E_df = Z/2xZ/4:", [k for k in PAIRS if tors_df[k] == (2, 4)])
print("    критерий совпал с расчётом Sage:",
      sorted(crit_df) == sorted([k for k in PAIRS if tors_df[k] == (2, 4)]))
print("    +-A*B — квадрат у пар:", crit_bh)
print("    пары с кручением E_bh = Z/2xZ/4:", [k for k in PAIRS if tors_bh[k] == (2, 4)])
print("    критерий совпал с расчётом Sage:",
      sorted(crit_bh) == sorted([k for k in PAIRS if tors_bh[k] == (2, 4)]))
print("""    ИТОГ: у обоих угловых факторов кручение равно Z/2 x Z/2 у всех 127 пар,
    а НЕ Z/2 x Z/4, как в рёберном случае. Это делает аргумент «ранг 0 + перечисление
    кручения» в угловом случае ещё слабее: перечислять пришлось бы всего 4 точки,
    а рациональных точек заведомо не меньше 12 (§4c).""")

sub("§4c. РЕШАЮЩИЙ ПОДСЧЁТ: сколько рациональных точек заведомо есть на угловых квартиках")
def pts_df(m0, n0):
    S = m0**2 + n0**2; out = []
    for t0 in [QQ(0), QQ(m0-n0)/QQ(m0+n0), -QQ(m0-n0)/QQ(m0+n0),
               QQ(m0+n0)/QQ(m0-n0), -QQ(m0+n0)/QQ(m0-n0)]:
        v = S**2*t0**4 + (2*S**2 - 64*m0**2*n0**2)*t0**2 + S**2
        if v > 0 and v.is_square():
            out.append(t0)
    return sorted(set(out))
def pts_bh(m0, n0):
    A_, B_ = QQ(3*n0**2 - m0**2), QQ(3*m0**2 - n0**2); out = []
    for t0 in [QQ(1), QQ(-1), QQ(m0-n0)/QQ(m0+n0), -QQ(m0-n0)/QQ(m0+n0),
               QQ(m0+n0)/QQ(m0-n0), -QQ(m0+n0)/QQ(m0-n0)]:
        v = (A_*t0**2 + B_)*(B_*t0**2 + A_)
        if v > 0 and v.is_square():
            out.append(t0)
    return sorted(set(out))
n_df = {k: len(pts_df(*k)) for k in PAIRS}
n_bh = {k: len(pts_bh(*k)) for k in PAIRS}
print("  B_df: различных t с квадратным c3*c5 (проверено пять кандидатов):", dist(n_df))
print("  B_bh: различных t с квадратным c1*c7 (проверено шесть кандидатов):", dist(n_bh))
print("  старший коэффициент B_df равен S^2 — квадрат всегда => ещё 2 точки на бесконечности")
print("  старший коэффициент B_bh равен A*B; квадрат у пар:",
      sum(1 for (m0, n0) in PAIRS if ZZ((3*n0**2-m0**2)*(3*m0**2-n0**2)) > 0
          and ZZ((3*n0**2-m0**2)*(3*m0**2-n0**2)).is_square()))
low_df = {k: 2*n_df[k] + 2 for k in PAIRS}
low_bh = {k: 2*n_bh[k] for k in PAIRS}
tord_df = {k: prod(tors_df[k]) for k in PAIRS}
tord_bh = {k: prod(tors_bh[k]) for k in PAIRS}
print("  нижняя оценка числа рациональных точек: B_df:", dist(low_df), " B_bh:", dist(low_bh))
print("  порядок кручения якобиана:            E_df:", dist(tord_df), " E_bh:", dist(tord_bh))
print("""  Квартика с рациональной точкой изоморфна своему якобиану, значит число её
  рациональных точек = |E(Q)|. Если точек больше, чем кручения, то rank >= 1.""")
print("     B_df: нижняя оценка > |кручение| у %d пар из 127" %
      sum(1 for k in PAIRS if low_df[k] > tord_df[k]))
print("     B_bh: нижняя оценка > |кручение| у %d пар из 127" %
      sum(1 for k in PAIRS if low_bh[k] > tord_bh[k]))
print("  [доказано перебором] => rank >= 1 у ОБОИХ угловых факторов при всех 127 парах,")
print("  и этот вывод не использует PARI: только счёт точек и кручение.")

sub("§4d. Проверка ранга напрямую (PARI ellrank)")
def rk(E):
    try:
        r = pari(E).ellrank()
        return (int(r[0]), int(r[1]))
    except Exception:
        return None
subset = PAIRS[:25] if FAST else PAIRS
t0_ = time.time()
rk_df = {}; rk_bh = {}; rk_edge = {}
for (m0, n0) in subset:
    rk_df[(m0,n0)] = rk(corner_E_df(m0, n0))
    rk_bh[(m0,n0)] = rk(corner_E_bh(m0, n0))
    rk_edge[(m0,n0)] = rk(edge_E_B(m0, n0))
print("  посчитано пар: %d за %.0f с" % (len(subset), time.time()-t0_))
print("  РЁБЕРНОЕ E_B : границы ранга", dist(rk_edge),
      "; ранг 0 у", sum(1 for v in rk_edge.values() if v == (0, 0)), "пар")
print("  УГЛОВОЕ E_df : границы ранга", dist(rk_df),
      "; ранг 0 у", sum(1 for v in rk_df.values() if v == (0, 0)), "пар")
print("  УГЛОВОЕ E_bh : границы ранга", dist(rk_bh),
      "; ранг 0 у", sum(1 for v in rk_bh.values() if v == (0, 0)), "пар")
if not FAST:
    exp = {(0, 0): 62, (1, 1): 60, (2, 2): 4, (0, 2): 1}
    print("  КОНТРОЛЬ: RESULTS_FOR_CLAUDE_2026-09-12_G1.md §3 даёт для E_B ровно")
    print("            62 раза [0,0], 60 раз [1,1], 4 раза [2,2], один раз [0,2]. Совпало:",
          dist(rk_edge) == exp)
    assert dist(rk_edge) == exp

sub("§4e. Явные неторсионные точки угловых факторов (генераторы от PARI) и их высоты")
for (m0, n0) in [(3, 2), (5, 2), (13, 8), (16, 5), (19, 5), (20, 19)]:
    for nmE, E in [("E_df", corner_E_df(m0, n0)), ("E_bh", corner_E_bh(m0, n0))]:
        try:
            r = pari(E).ellrank()
            gens = [E(QQ(g[0]), QQ(g[1])) for g in r[3]]
            nt = [P for P in gens if not P.has_finite_order()]
            hs = [RR(P.height()).n(24) for P in nt[:2]]
            print("  (%2d,%2d) %s: кручение %-6s ранг %s, неторсионных генераторов %d, высоты %s"
                  % (m0, n0, nmE, str(tuple(E.torsion_subgroup().invariants())),
                     (int(r[0]), int(r[1])), len(nt), hs))
        except Exception as ex:
            print("  (%2d,%2d) %s: не удалось (%s)" % (m0, n0, nmE, type(ex).__name__))

# ============================================================================
# §5. СПОСОБ 3 — НЕОБХОДИМЫЙ КЛАСС КУММЕРА
# ============================================================================
hdr("§5. СПОСОБ 3. Необходимый класс Куммера: угловой аналог условия (1,s,s)")

print("""Рёберный случай (REPORT_G1_SIX_KUMMER_2026-09-12.md, «Модели и доказательство»):
   C:  u0^2=F0, u4^2=F4, u8^2=F8            (главная диагональ a-e-i)
   E:  V^2=(X+b)(X+s m^4)(X+s n^4),  b = s m^2 n^2,  X = b t^2
   X+b = (mn)^2 u4^2,   X+s m^4 = s m^2 u0^2,   X+s n^4 = s n^2 u8^2
   => необходим класс delta = (1, s, s) в образе 2-спуска на E(Q)/2E(Q).

Угловой аналог. Свободные клетки G2 образуют крест; в нём ровно две магические линии:
средний столбец c1-c4-c7 и средняя строка c3-c4-c5. Берём средний столбец:
   C_c: u1^2=c1, u4^2=c4, u7^2=c7,   c1=alpha t^2+beta, c7=beta t^2+alpha, c4=s(1+t^2)
Все три ЛИНЕЙНЫ по T=t^2, поэтому Y^2=c1 c4 c7 — кубика по T, то есть эллиптическая кривая:
   E_c: Y^2 = (X + s*alpha*beta)(X + s*beta^2)(X + s*alpha^2),   X = s*alpha*beta*t^2
   X + s alpha beta = alpha*beta * c4 = alpha*beta * u4^2
   X + s beta^2     = s*beta  * c1  = s*beta  * u1^2
   X + s alpha^2    = s*alpha * c7  = s*alpha * u7^2
   => необходим класс delta_c = (alpha*beta, s*beta, s*alpha).
Это ТОЧНО тот же вид, что рёберный, с заменой (m^2, n^2) -> (beta, alpha); в рёберном случае
alpha=n^2, beta=m^2 — квадраты, и класс сворачивается в (1, s, s).""")

sub("§5a. Проверка тождеств карты C_c -> E_c символически над Q(m,n,t)")
Xsym = s_sym*alpha*beta*tt**2
id1 = Fld(Xsym + s_sym*alpha*beta - alpha*beta*g4)
id2 = Fld(Xsym + s_sym*beta**2 - s_sym*beta*g1)
id3 = Fld(Xsym + s_sym*alpha**2 - s_sym*alpha*g7)
print("  X + s*alpha*beta - alpha*beta*c4 =", id1)
print("  X + s*beta^2     - s*beta*c1     =", id2)
print("  X + s*alpha^2    - s*alpha*c7    =", id3)
assert id1 == 0 and id2 == 0 and id3 == 0
print("  все три тождества верны => класс delta_c = (alpha*beta, s*beta, s*alpha)")
sub("§5a'. То же самое для средней СТРОКИ c3-c4-c5 (второй угловой вариант)")
print("""  c3 = s(1+t^2)-4mnt, c5 = s(1+t^2)+4mnt, c4 = s(1+t^2) — по T=t^2 НЕ линейны
  (есть член t), поэтому строка даёт не кубику, а квартику: Y^2=c3 c4 c5 — кривая рода 1
  степени 6 по t. Аналог 2-спуска строится, но через якобиан B_df из §4, а не напрямую.

  Замечание о роде. Род кривой трёх условий в G-семействах равен 5 — это счёт проекта
  (BREMNER16_2026-09-12.md §2.2, [численно], в том числе для G2 при (m,n)=(5,2)); здесь он
  НЕ перепроверялся. Новое структурное наблюдение: в угловом случае все три формы среднего
  СТОЛБЦА чётны по t, поэтому C_c пропускается через кривую над T=t^2, и эллиптический
  фактор E_c получается прямой кубикой по T — ровно как F0,F4,F8 в рёберном случае.""")

def sqfree(x):
    x = QQ(x)
    if x == 0:
        return 0
    num = x.numerator()*x.denominator()
    return sign(num)*prod(p**(e % 2) for p, e in ZZ(num).abs().factor())

sub("§5b. КОНТРОЛЬ: общая формула на РЁБЕРНОМ случае обязана дать (1,s,s)")
for (m0, n0) in [(11, 4), (19, 5), (15, 1), (13, 8)]:
    s0 = QQ(m0**2 + n0**2)/2
    a_e, b_e = QQ(n0**2), QQ(m0**2)
    cls = tuple(sqfree(x) for x in (a_e*b_e, s0*b_e, s0*a_e))
    print("  (%2d,%2d): s=%-7s общая формула даёт %s ; в отчёте Кодекса было (1,s,s) = %s"
          % (m0, n0, s0, cls, (1, sqfree(s0), sqfree(s0))))
    assert cls == (1, sqfree(s0), sqfree(s0))

sub("§5c. Угловой класс delta_c: ВСЕГДА ли он в образе (все 127 пар)")
present = 0; mismatch = []
for (m0, n0) in PAIRS:
    s0 = QQ(m0**2+n0**2)/2
    a_ = QQ(3*n0**2 - m0**2)/2; b_ = QQ(3*m0**2 - n0**2)/2
    need = tuple(sqfree(x) for x in (a_*b_, s0*b_, s0*a_))
    t0 = QQ(m0-n0)/QQ(m0+n0)
    c1v = a_*t0**2 + b_; c4v = s0*(1+t0**2); c7v = b_*t0**2 + a_
    on_curve = (c1v > 0 and c1v.is_square() and c4v > 0 and c4v.is_square()
                and c7v > 0 and c7v.is_square())
    X0 = s0*a_*b_*t0**2
    got = tuple(sqfree(x) for x in (X0 + s0*a_*b_, X0 + s0*b_**2, X0 + s0*a_**2))
    if on_curve and got == need:
        present += 1
    else:
        mismatch.append((m0, n0, on_curve, need, got))
print("  пар, где латинская точка лежит на C_c И её класс РАВЕН требуемому delta_c:",
      present, "из", len(PAIRS))
print("  расхождений:", mismatch)
assert present == len(PAIRS)
print("""
  [доказано] Тест «нужного класса нет в образе 2-спуска» в угловом случае не исключает
  НИЧЕГО: класс delta_c = (alpha*beta, s*beta, s*alpha) реализуется рациональной точкой
  кривой C_c при любых (m,n). В рёберном случае он отсутствовал у (11,4), (19,5), (15,1)
  именно потому, что там рациональной точки нет вовсе.""")

sub("§5d. Явные угловые классы и модель E_c на нескольких парах")
for (m0, n0) in [(3, 2), (5, 2), (13, 8), (16, 5), (19, 5)]:
    s0 = QQ(m0**2+n0**2)/2
    a_ = QQ(3*n0**2 - m0**2)/2; b_ = QQ(3*m0**2 - n0**2)/2
    e1, e2, e3 = -s0*a_*b_, -s0*b_**2, -s0*a_**2
    u = lcm([QQ(x).denominator() for x in (e1, e2, e3)])
    E = Eroots(e1*u**2, e2*u**2, e3*u**2)
    need = tuple(sqfree(x) for x in (a_*b_, s0*b_, s0*a_))
    print("  (%2d,%2d): alpha=%-8s beta=%-8s s=%-7s" % (m0, n0, a_, b_, s0))
    print("           E_c = %s, кручение %s, ранг(PARI) %s"
          % (E.ainvs(), tuple(E.torsion_subgroup().invariants()), rk(E)))
    print("           требуемый delta_c = %s" % (need,))

# ============================================================================
# §6. УГЛОВОЕ СЕЧЕНИЕ
# ============================================================================
hdr("§6. Угловое СЕЧЕНИЕ: почему ни одно отношение угловой пары нельзя исключить «по аналогии»")

print("""Рёберная теорема фиксирует отношение РЁБЕРНОЙ пары (b:h) при центре n0:
   b=b0^2, h=h0^2, e=n0^2 с b0^2+h0^2=2n0^2; вторая рёберная пара — на конике u^2+v^2=2;
   четыре свободные клетки (УГЛЫ) = попарные средние клеток двух рёберных пар.
Угловое сечение фиксирует УГЛОВУЮ пару: a=b0^2, i=h0^2, e=n0^2; вторая угловая — на конике;
   свободны четыре РЁБЕРНЫЕ клетки.""")

def corner_section_cells(b0, h0, n0, t0):
    Bf = (1 + 2*t0 - t0**2)/(1 + t0**2)
    Hf = (1 - 2*t0 - t0**2)/(1 + t0**2)
    a = QQ(b0)**2; i = QQ(h0)**2; e = QQ(n0)**2
    c = e*Bf**2; g = e*Hf**2
    b = 3*e - a - c; h = 2*e - b
    d = 3*e - a - g; f = 2*e - d
    cs = [a, b, c, d, e, f, g, h, i]
    assert all(sum(cs[i] for i in L) == 3*e for L in LINES)
    return cs

SECS = [(17,7,13),(23,7,17),(71,49,61),(7,1,5),(31,17,25),(41,1,29),(47,23,37),
        (49,31,41),(73,17,53),(89,23,65),(79,47,65)]
print("\n  Три запрещённых для РЁБЕР отношения и остальные восемь сечений проекта —")
print("  поставленные в УГЛОВУЮ позицию, при t=0:")
for (b0, h0, n0) in SECS:
    assert b0**2 + h0**2 == 2*n0**2
    cs = corner_section_cells(b0, h0, n0, QQ(0))
    allsq = all(v > 0 and QQ(v).is_square() for v in cs)
    mark = " <-- ЗАПРЕЩЕНО ДЛЯ РЁБЕР" if (b0, h0, n0) in [(17,7,13),(23,7,17),(71,49,61)] else ""
    print("    (%2d:%2d:%2d): все девять клеток квадраты: %-5s ; различных значений: %d%s"
          % (b0, h0, n0, allsq, len(set(cs)), mark))
    assert allsq
cs = corner_section_cells(17, 7, 13, QQ(0))
print("\n  например (17:7:13):", [str(v) for v in cs], "= латинский квадрат из 289,49,169")
print("""
  [доказано] При t=0 угловое сечение даёт латинский квадрат
        b0^2 h0^2 n0^2 / h0^2 n0^2 b0^2 / n0^2 b0^2 h0^2,
  магический РОВНО потому, что b0^2+h0^2=2n0^2 — то есть по самому определению сечения.
  Значит для ЛЮБОГО отношения (b0:h0:n0) угловая пара реализуется в магическом квадрате
  из девяти квадратов (с повторениями). Формулировка «такого отношения у угловой пары не
  бывает» ЛОЖНА без оговорки «в невырожденном квадрате».""")

sub("§6a. Важная оговорка: у РЁБЕРНОГО сечения латинская точка тоже есть, только при другом t")
print("  условие Bf(t) = +-b0/n0 сводится к квадратному уравнению с дискриминантом 2n0^2-b0^2 = h0^2:")
for (b0, h0, n0) in SECS[:6]:
    roots = []
    for sgn in (1, -1):
        A_ = -(n0 + sgn*b0); B_ = 2*n0; C_ = n0 - sgn*b0
        if A_ == 0:
            if B_ != 0:
                roots.append(QQ(-C_)/QQ(B_))
            continue
        D = B_**2 - 4*A_*C_
        if D >= 0 and ZZ(D).is_square():
            sq = ZZ(D).sqrt()
            roots += [QQ(-B_+sq)/QQ(2*A_), QQ(-B_-sq)/QQ(2*A_)]
    print("    (%2d:%2d:%2d): 2n0^2-b0^2 = %6d = %d^2 ; t с латинским квадратом: %s"
          % (b0, h0, n0, 2*n0**2 - b0**2, h0, roots))
print("""  [наблюдение, требует сверки с текстом теоремы] Отсюда следует, что рёберная теорема
  об исключении отношений НЕ может доказывать «рациональных точек нет»: они есть у любого
  сечения. Она должна доказывать «все рациональные точки вырождены». Это не ослабляет её,
  но меняет то, какой угловой аналог нужно искать: не «точек нет», а «точки только вырождены».""")

sub("§6c. Локальный скан обоих СЕЧЕНИЙ (p<200): считаем, а не рассуждаем")
Rs = PolynomialRing(QQ, 'w'); w = Rs.gen()
def section_forms(kind, b0, h0, n0):
    """однородные бинарные формы степени 4 для четырёх свободных клеток сечения.
       kind='edge': фиксирована рёберная пара {1,7}; kind='corner': угловая пара {0,8}."""
    Bf = (1 + 2*w - w**2); Hf = (1 - 2*w - w**2); D = (1 + w**2)     # Bf^2+Hf^2 = 2 D^2
    assert Rs(Bf**2 + Hf**2 - 2*D**2) == 0
    e = QQ(n0)**2*D**2                      # все клетки умножены на квадрат D^2
    if kind == 'edge':
        known = {1: QQ(b0)**2*D**2, 7: QQ(h0)**2*D**2,
                 3: QQ(n0)**2*Bf**2, 5: QQ(n0)**2*Hf**2}
    else:
        known = {0: QQ(b0)**2*D**2, 8: QQ(h0)**2*D**2,
                 2: QQ(n0)**2*Bf**2, 6: QQ(n0)**2*Hf**2}
    free = [x for x in range(9) if x not in known and x != 4]
    rows = []; rhs = []
    for L in LINES:
        row = [1 if x in L else 0 for x in free]
        if any(row):
            rows.append(row)
            rhs.append(3*e - sum(known.get(x, Rs(0)) for x in L) - (e if 4 in L else 0))
    sol = Matrix(Rs.fraction_field(), rows).solve_right(vector(Rs.fraction_field(), rhs))
    cells = [None]*9; cells[4] = e
    for x, v in known.items():
        cells[x] = Rs(v)
    for j, x in enumerate(free):
        cells[x] = Rs(sol[j])
    assert all(Rs(sum(cells[i] for i in L) - 3*e) == 0 for L in LINES)
    out = []
    for i in free:
        cf = cells[i].coefficients(sparse=False) + [0]*5
        # однородизация t=u/v, умножение на v^4 (квадрат): sum cf[k] u^k v^(4-k)
        out.append([QQ(cf[k]) for k in range(5)])
    dn = lcm([QQ(x).denominator() for f in out for x in f])
    return [[ZZ(x*dn) for x in f] for f in out], free
def section_scan(kind, b0, h0, n0, p):
    fs, free = section_forms(kind, b0, h0, n0)
    QR = qrset(p); cnt = 0
    for (a, b) in [(x, 1) for x in range(p)] + [(1, 0)]:
        vals = [sum(int(f[k]) * pow(a, k, p) * pow(b, 4-k, p) for k in range(5)) % p for f in fs]
        if all(v in QR for v in vals):
            cnt += 1
    return cnt
print("  (для каждой тройки: минимум по p<200 числа выживших точек P^1(F_p))")
for (b0, h0, n0) in SECS:
    me = min(section_scan('edge', b0, h0, n0, p) for p in PRIMES200 if n0 % p)
    mc = min(section_scan('corner', b0, h0, n0, p) for p in PRIMES200 if n0 % p)
    print("    (%2d:%2d:%2d): рёберное сечение min выживших = %2d ; угловое = %2d"
          % (b0, h0, n0, me, mc))
print("""  [проверено перебором] Ни у рёберных, ни у угловых сечений локального препятствия нет:
  у обоих есть латинская точка. Значит теорема об исключении отношений (и рёберная, и любой
  её угловой аналог) обязана опираться на «все точки вырождены», а не на «точек нет».""")

sub("§6b. Где именно расходятся рёберное и угловое сечение")
print("""  Рёберное сечение, t=0: четыре свободные клетки (УГЛЫ) равны A,A,C,C с
      A=(h0^2+n0^2)/2, C=(b0^2+n0^2)/2 — квадраты лишь иногда.
  Угловое сечение, t=0: четыре свободные клетки (РЁБРА) равны b0^2,b0^2,h0^2,h0^2 —
      квадраты ВСЕГДА, по определению сечения.
  Проверка на всех 11 сечениях проекта:""")
for (b0, h0, n0) in SECS:
    A_ = QQ(h0**2 + n0**2)/2; C_ = QQ(b0**2 + n0**2)/2
    print("    (%2d:%2d:%2d): рёберное A=%-8s квадрат %-5s ; C=%-8s квадрат %-5s || угловое: b0^2,h0^2 квадраты True"
          % (b0, h0, n0, A_, A_.is_square(), C_, C_.is_square()))

# ============================================================================
# §7. ИТОГ
# ============================================================================
hdr("§7. ИТОГ и что делать дальше")
print("""
ОТВЕТ НА ЗАДАНИЕ — все три дешёвых способа в угловом случае неприменимы В ПРИНЦИПЕ,
и причина у всех трёх одна, посчитанная в §2:

  [доказано, тождество над Q(m,n)] Угловое семейство G2 при ЛЮБЫХ (m,n) имеет четыре
  рациональные точки t = +-(m-n)/(m+n), +-(m+n)/(m-n), в которых ВСЕ ДЕВЯТЬ клеток —
  квадраты (вырожденный латинский квадрат из трёх различных значений).

1. ЛОКАЛЬНЫЕ ПРЕПЯТСТВИЯ — невозможны. Кривая с рациональной точкой локально разрешима
   всюду. Посчитано: перебор P^1(F_p) для всех 127 пар и всех нечётных p < 200 — ни одного
   пустого случая (§3c), тогда как тот же код на рёберном семействе воспроизводит
   известные препятствия (13,8) mod 17 и (16,5) mod 41 и ловит десятки пар (§3a, §3b).

2. ФАКТОР РАНГА 0 — невозможен: rank >= 1 всегда. Обе угловые квартики имеют >= 12
   рациональных точек, а кручение их якобианов равно Z/2 x Z/2, порядок всего 4 (§4b, §4c);
   квартика с рациональной точкой изоморфна якобиану, значит ранг не нуль. Подтверждено
   PARI ellrank (ни одного ранга 0 из 127x2) и явными генераторами (§4d, §4e).
   ОТВЕТ НА ВОПРОС ПРО КРУЧЕНИЕ: оно НЕ такое, как в рёберном случае. Там Z/2 x Z/4
   универсально (три разности корней — квадраты (m^2-n^2)^2, (m^2+n^2)^2, (2mn)^2).
   В угловом Z/2 x Z/4 требует |S^2-16m^2n^2| = квадрат (для E_df) или +-A*B = квадрат
   (для E_bh); ни у одной из 127 пар это не выполнено, поэтому всюду Z/2 x Z/2.

3. НЕОБХОДИМЫЙ КЛАСС КУММЕРА — существует и выписан:
   E_c: Y^2=(X+s*alpha*beta)(X+s*beta^2)(X+s*alpha^2), X = s*alpha*beta*t^2,
   delta_c = (alpha*beta, s*beta, s*alpha),  alpha=(3n^2-m^2)/2, beta=(3m^2-n^2)/2.
   Это буквально рёберная формула с заменой (m^2,n^2) -> (beta,alpha); при alpha,beta
   квадратах сворачивается в (1,s,s) (§5b). Но исключить не может ничего: класс
   реализуется латинской точкой у всех 127 пар (§5c).

ЧТО ЭТО ЗНАЧИТ ДЛЯ ПРОЕКТА.
 - Угловой случай труден КАЧЕСТВЕННО иначе: там надо доказывать не «точек нет», а
   «все точки вырождены». Это Шаботи / MW-сито, а не локальные сравнения.
 - Дыра, которую отмечал Grok (теорема покрывает две пары из четырёх), дёшево не закрывается.
 - §6 показывает больше: она не закрывается и формулировкой «отношение запрещено» —
   любое отношение угловой пары реализуется в вырожденном квадрате. Правильная угловая
   формулировка обязана содержать слово «невырожденный».
 - §6a: это, по-видимому, верно и для РЁБЕРНОЙ теоремы. Проверить по её тексту.

ПОБОЧНАЯ НАХОДКА ПО РЁБЕРНОМУ СЛУЧАЮ (§3b'), не относится к заданию, но дешёвая.
 Девятиклеточная локальная система G1 пуста над F_p при некотором p < 200 не только у
 (13,8) и (16,5), но и у (11,4), (15,1), (19,5), (19,16). Из семи открытых пар
 остаётся одна — (15,8), и ровно потому, что 15^2+8^2=17^2 даёт рациональную точку.
 Это НЕ заменяет доказанное C(Q)=пусто на трёхклеточной кривой (то сильнее), но для
 вопроса «есть ли в слое магический квадрат из девяти квадратов» достаточно и дёшево.
 Проверено двумя независимыми выводами клеток и сверено с local_13_8_p17_прямая.log §7.

КОНТРОЛИ, НА КОТОРЫХ ДЕРЖИТСЯ ЭТОТ РАСЧЁТ (все сошлись).
 - Девять клеток G2 выведены заново из магичности и совпали с bremner16/engine2.sage (§1a).
 - Тот же локальный код на G1 даёт (13,8)->[17,47,73,79,103,191] и (16,5)->41, как в
   local_13_8_p17_прямая.log §7 (§3a, §3b').
 - Тот же код якобиана на рёберной квартике z^2=F0F8 даёт в точности E_B из
   REPORT_G1_SIX_KUMMER (§4a), а распределение рангов E_B — ровно 62/60/4/1, как в
   RESULTS_FOR_CLAUDE §3 (§4d).
 - Общая формула класса Куммера на рёберных параметрах даёт (1,s,s) и совпадает с
   числами отчёта Кодекса для (11,4), (19,5), (15,1) (§5b).
 - Знак якобиана квартики проверен счётом точек #C(F_p)=#E(F_p) (§4a'); с неверным знаком
   получается квадратичный твист на -1, у которого кручение ложно выглядит как Z/2 x Z/4.

ЧТО В УГЛОВОМ СЛУЧАЕ ВСЁ-ТАКИ МОЖЕТ СРАБОТАТЬ (не дёшево, но реально).
 - Ранг E_df равен 1 у 53 пар из 127 (§4d). При ранге 1 и известных вырожденных точках
   естественно считать эллиптический Шаботи / MW-сито на угловой кривой рода 5: задача
   не «есть ли точка», а «есть ли точка вне известного вырожденного множества».
 - Вырожденное множество углового слоя выписано явно: t in {0, +-1, oo} (повторы клеток)
   плюс четыре латинские точки +-(m-n)/(m+n), +-(m+n)/(m-n). Это конечный явный список,
   его можно закладывать в сито как известные точки.

ЧТО НЕ ДОКАЗАНО И НЕ ЗАЯВЛЯЕТСЯ.
 - НЕ доказано, что в угловом семействе нет невырожденного квадрата. Показано лишь, что
   три конкретных дешёвых метода его не найдут и не исключат.
 - НЕ проверены p >= 200 и (m,n) вне 1<=n<m<=20 — но для §3 это и не нужно: отсутствие
   локального препятствия следует из наличия рациональной точки, а она есть тождественно.
 - Ранг углового фактора посчитан PARI (программный аргумент). Независимое подтверждение
   rank>=1 — §4c, счёт точек — от PARI не зависит.

НЕ ДЕЛАТЬ.
 - Не запускать сканы локальных препятствий по угловому семейству на больших p: результат
   известен заранее и он пустой (§2 => §3).
 - Не искать у углового семейства фактор ранга 0: его нет ни у одной из двух квартик,
   и причина структурная, а не вычислительная.
 - Не переносить теорему об исключении отношений на угловые пары «по аналогии»: §6.
 - Не путать латинскую точку с решением: в ней три различных значения, а не девять.
""")
print("время работы: %.0f с" % (time.time() - T00))
