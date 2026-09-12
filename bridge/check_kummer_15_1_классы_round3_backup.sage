# -*- coding: utf-8 -*-
# check_kummer_15_1_классы.sage
#
# НЕЗАВИСИМАЯ ПОПЫТКА СЛОМАТЬ заявление Codex:  C_{15,1}(Q) = пусто.
# Источник заявления: RESULT_G1_15_1_KUMMER_EXCLUDED_2026-09-12.md
#
# Установка: не подтверждать, а ломать.  Вся арифметика классов — своя:
#   * свой квадратсвободный представитель со знаком (факторизация целого),
#   * свой групповой закон на y^2 = x^3 + a2 x^2 + a4 x + a6,
#   * своя delta с ЯВНОЙ регуляризацией на 2-кручении И независимая перепроверка
#     этой регуляризации через точки, которые 2-кручением НЕ являются,
#   * СВОЙ неторсионный генератор, найденный перебором x = p/q^2 (isqrt),
#     без чужой точки и без E.gens(),
#   * свой критерий квадрата в Q_p.
#
# Sage используется как независимый источник ГРАНИЦ РАНГА (eclib / PARI ellrank)
# и как свидетель для самотестов; арифметика delta от него не зависит.
#
# Тяжёлые переборы вынесены в exec(raw-строка) — там чистый Python без препарсера.
#
# Автор: Claude (аудит-атака, раунд 3), 2026-09-12.

import itertools, random, sys
from math import isqrt

def hdr(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)

FAIL = []
def check(name, cond):
    ok = bool(cond)
    print(("  [OK]   " if ok else "  [FAIL] ") + name)
    if not ok:
        FAIL.append(name)
    return ok

# =========================================================== 0. параметры
hdr("0. ПАРАМЕТРЫ СЕМЕЙСТВА G1 ДЛЯ (m,n) = (15,1) — СЧИТАЮ САМ")
m = 15
n = 1
check("gcd(m,n)=1", gcd(m, n) == 1)
s = QQ(m**2 + n**2) / 2
b = s * m**2 * n**2
e1 = -b
e2 = -s * m**4
e3 = -s * n**4
print("  s  =", s)
print("  b  =", b)
print("  e1 = %s,  e2 = %s,  e3 = %s" % (e1, e2, e3))
check("s == 113 (как у Codex)", s == 113)
check("b == 25425 (как у Codex)", b == 25425)
check("E: V^2=(X+25425)(X+5720625)(X+113) (как у Codex)",
      (-e1, -e2, -e3) == (25425, 5720625, 113))
check("корни попарно различны", len(set([e1, e2, e3])) == 3)

# ============================================ 1. тождества C -> E (символьно)
hdr("1. ТОЖДЕСТВА C -> E: СИМВОЛЬНАЯ ПРОВЕРКА В Q[t]")
Rt = PolynomialRing(QQ, 't')
tv = Rt.gen()
F0 = m**2 + n**2 * tv**2
F4 = s * (1 + tv**2)
F8 = n**2 + m**2 * tv**2
Xt = b * tv**2
print("  F0 =", F0)
print("  F4 =", F4)
print("  F8 =", F8)
print("  X  =", Xt)
check("X - e1 == (mn)^2 * F4", (Xt - e1) - (m * n) ** 2 * F4 == 0)
check("X - e2 == s*m^2 * F0", (Xt - e2) - s * m**2 * F0 == 0)
check("X - e3 == s*n^2 * F8", (Xt - e3) - s * n**2 * F8 == 0)
check("(X-e1)(X-e2)(X-e3) == b^2*F0*F4*F8   (V = b u0u4u8 лежит на E)",
      (Xt - e1) * (Xt - e2) * (Xt - e3) - b**2 * F0 * F4 * F8 == 0)
print("  рациональные корни F0,F4,F8:", F0.roots(QQ), F4.roots(QQ), F8.roots(QQ))
check("ни одна F_i не зануляется при рациональном t",
      F0.roots(QQ) == [] and F4.roots(QQ) == [] and F8.roots(QQ) == [])
check("X = b t^2 >= 0 > e_i  =>  X-e_i > 0: образ не 2-кручение и не O",
      b > 0 and e1 < 0 and e2 < 0 and e3 < 0)

# ======================================== 2. свой квадратсвободный представитель
hdr("2. СВОЙ КВАДРАТСВОБОДНЫЙ ПРЕДСТАВИТЕЛЬ КЛАССА В Q*/Q*^2")

def sqfree(q):
    """Квадратсвободный представитель класса q in Q*/Q*^2, со знаком."""
    q = QQ(q)
    if q == 0:
        raise ValueError("нулевой аргумент sqfree")
    z = ZZ(q.numerator()) * ZZ(q.denominator())
    sgn = 1 if z > 0 else -1
    z = z.abs()
    r = ZZ(1)
    for (p, ex) in ZZ(z).factor():
        if ex % 2 == 1:
            r *= p
    return sgn * r

_st = [(QQ(1), 1), (QQ(4), 1), (QQ(-9), -1), (QQ(50400), 14), (QQ(-50400), -14),
       (QQ(3) / QQ(2), 6), (QQ(-8) / QQ(18), -1), (QQ(113) * 4, 113),
       (QQ(5720625) - QQ(75825), 2)]
check("самотест sqfree на 9 значениях", all(sqfree(a) == ZZ(v) for (a, v) in _st))
for (a, v) in _st:
    print("     sqfree(%s) = %s   (ожидалось %s)" % (a, sqfree(a), v))
random.seed(20260912)
ok2 = True
for _ in range(400):
    qq = QQ(random.randint(-10 ** 7, 10 ** 7)) / QQ(random.randint(1, 10 ** 5))
    if qq == 0:
        continue
    r = QQ(sqfree(qq)) / qq
    if not (r > 0 and r.is_square()):
        ok2 = False
        print("     ПЛОХО:", qq, sqfree(qq))
check("на 400 случайных q: sqfree(q)/q — положительный квадрат", ok2)

required = (sqfree((m * n) ** 2), sqfree(s * m**2), sqfree(s * n**2))
print("  ТРЕБУЕМЫЙ класс delta(P) для любой конечной точки C:", required)
check("требуемый класс = (1, s, s) = (1,113,113)", required == (ZZ(1), ZZ(113), ZZ(113)))

# ======================================== 3. кривая, кручение, СВОЙ поиск точки
hdr("3. КРИВАЯ E, КРУЧЕНИЕ, СВОЙ ПОИСК ТОЧЕК ПЕРЕБОРОМ x = p/q^2")
a2 = ZZ(-(e1 + e2 + e3))
a4 = ZZ(e1 * e2 + e1 * e3 + e2 * e3)
a6 = ZZ(-e1 * e2 * e3)
print("  y^2 = x^3 + (%s) x^2 + (%s) x + (%s)" % (a2, a4, a6))

def f(x):
    x = QQ(x)
    return (x - e1) * (x - e2) * (x - e3)

Rx = PolynomialRing(QQ, 'x')
xx = Rx.gen()
check("f(x) == x^3 + a2 x^2 + a4 x + a6",
      (xx - e1) * (xx - e2) * (xx - e3) == xx**3 + a2 * xx**2 + a4 * xx + a6)

E = EllipticCurve(QQ, [0, a2, 0, a4, a6])
print("  E =", E)
print("  disc =", E.discriminant())
print("  conductor =", E.conductor())
T = E.torsion_subgroup()
print("  E(Q)_tors =", T.invariants(), " порядок =", T.order())
check("кручение ровно (Z/2)^2", sorted([ZZ(i) for i in T.invariants()]) == [2, 2])
print("  (даже при большем кручении |T/2T| = |T[2]| = 4, так что |E(Q)/2E(Q)| = 2^(r+2).)")

# ---- СВОЙ перебор. x = p/q^2 (знаменатель обязан быть точным квадратом) ----
SCAN = r'''
from math import isqrt, gcd
A, B, Cc = 25425, 5720625, 113          # f(x) = (x+A)(x+B)(x+Cc)
res = []
# q = 1: целые x. f(x) >= 0 на [-B,-A] U [-Cc, oo)
for x in range(-B, -A + 1):
    v = (x + A) * (x + B) * (x + Cc)
    if v < 0:
        continue
    r = isqrt(v)
    if r * r == v:
        res.append((x, 1, r))
for x in range(-Cc, INT_HI + 1):
    v = (x + A) * (x + B) * (x + Cc)
    r = isqrt(v)
    if r * r == v:
        res.append((x, 1, r))
# q >= 2: x = p/q^2, числитель (p+A q^2)(p+B q^2)(p+Cc q^2) должен быть квадратом
for q in range(2, QMAX + 1):
    q2 = q * q
    lo, hi = -B * q2, PHI * q2
    for p in range(lo, hi + 1):
        if gcd(p, q) != 1:
            continue
        v = (p + A * q2) * (p + B * q2) * (p + Cc * q2)
        if v < 0:
            continue
        r = isqrt(v)
        if r * r == v:
            res.append((p, q, r))
FOUND = res
'''
ns = {"INT_HI": 3000000, "QMAX": 5, "PHI": 30000}
print("  перебор: q=1 по x in [-5720625,-25425] U [-113, 3*10^6];  q=2..5, x=p/q^2, |x|<=3*10^4")
exec(SCAN, ns)
raw = ns["FOUND"]
pts = []
for (p, q, r) in raw:
    x = QQ(p) / QQ(q * q)
    y = QQ(r) / QQ(q ** 3)
    assert y * y == f(x), ("плохая точка", p, q, r)
    pts.append((x, y))
print("  найдено точек (x, |y|):", len(pts))
for (x, y) in sorted(pts, key=lambda P: (P[0].denominator(), abs(P[0]))):
    print("     x = %-16s  y = %s" % (x, y))
check("свой перебор нашёл неторсионного кандидата (y != 0)", any(y != 0 for (x, y) in pts))

cands = [(x, y) for (x, y) in pts if y != 0]
cands.sort(key=lambda P: (P[0].denominator(), abs(P[0])))
MY_G = cands[0]
print("  МОЙ кандидат в генераторы: G =", MY_G)
check("G лежит на E: y^2 == f(x)", MY_G[1] ** 2 == f(MY_G[0]))

# ================================================= 4. СВОЙ групповой закон
hdr("4. СВОЙ ГРУППОВОЙ ЗАКОН НА y^2 = x^3 + a2 x^2 + a4 x + a6")
O = None

def neg(P):
    return O if P is O else (P[0], -P[1])

def add(P, Q):
    if P is O:
        return Q
    if Q is O:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0:
            return O
        lam = (3 * x1**2 + 2 * a2 * x1 + a4) / (2 * y1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    nu = y1 - lam * x1
    x3 = lam**2 - a2 - x1 - x2
    y3 = -(lam * x3 + nu)
    return (x3, y3)

def mul(k, P):
    k = ZZ(k)
    if k < 0:
        return mul(-k, neg(P))
    R = O
    Q = P
    while k > 0:
        if k % 2 == 1:
            R = add(R, Q)
        Q = add(Q, Q)
        k = k // 2
    return R

def on_curve(P):
    return P is O or P[1] ** 2 == f(P[0])

def to_sage(P):
    return E(0) if P is O else E(P[0], P[1])

tor = [O, (QQ(e1), QQ(0)), (QQ(e2), QQ(0)), (QQ(e3), QQ(0))]
random.seed(7)
ok = True
for _ in range(50):
    P = add(add(mul(random.randint(-6, 6), MY_G), mul(random.randint(0, 1), tor[1])),
            mul(random.randint(0, 1), tor[2]))
    Q = add(add(mul(random.randint(-6, 6), MY_G), mul(random.randint(0, 1), tor[1])),
            mul(random.randint(0, 1), tor[2]))
    S = add(P, Q)
    if not (on_curve(P) and on_curve(Q) and on_curve(S)):
        ok = False
    if to_sage(P) + to_sage(Q) != to_sage(S):
        ok = False
        print("     расхождение с Sage:", P, Q, S)
check("мой групповой закон = Sage на 50 случайных парах, всё на кривой", ok)
check("2T_i = O для i=1,2,3", all(mul(2, tor[i]) is O for i in (1, 2, 3)))
check("T1 + T2 = T3", add(tor[1], tor[2]) == tor[3])
mults = [mul(k, MY_G) for k in range(1, 15)]
check("kG != O при k=1..14 => G бесконечного порядка", all(P is not O for P in mults))
print("  x(kG), k=1..5:", [P[0] for P in mults[:5]])

# ========================================= 5. СВОЯ delta и её регуляризация
hdr("5. DELTA: СВОЁ ОПРЕДЕЛЕНИЕ + НЕЗАВИСИМАЯ ПРОВЕРКА РЕГУЛЯРИЗАЦИИ")
EE = [e1, e2, e3]

def delta(P):
    if P is O:
        return (ZZ(1), ZZ(1), ZZ(1))
    x, y = P
    out = []
    for i in range(3):
        d = x - EE[i]
        if d == 0:
            j, k = [u for u in range(3) if u != i]
            d = (EE[i] - EE[j]) * (EE[i] - EE[k])
        out.append(sqfree(d))
    return tuple(out)

def prod_cls(A, B):
    return tuple(sqfree(QQ(A[i]) * QQ(B[i])) for i in range(3))

print("  delta(O)  =", delta(O))
print("  delta(T1) =", delta(tor[1]), "   T1 = (%s, 0)" % e1)
print("  delta(T2) =", delta(tor[2]), "   T2 = (%s, 0)" % e2)
print("  delta(T3) =", delta(tor[3]), "   T3 = (%s, 0)" % e3)
print("  delta(G)  =", delta(MY_G))

print("\n  Перепроверка регуляризации через точки, НЕ являющиеся 2-кручением:")
ok = True
for i in (1, 2, 3):
    GT = add(MY_G, tor[i])
    if GT is O or GT[1] == 0:
        print("     G+T%d вырождена — пропуск" % i)
        continue
    derived = prod_cls(delta(MY_G), delta(GT))
    direct = delta(tor[i])
    print("     T%d: регуляризация %s   vs   delta(G)*delta(G+T%d) %s   %s"
          % (i, direct, i, derived, "совпало" if derived == direct else "РАСХОЖДЕНИЕ"))
    if derived != direct:
        ok = False
check("регуляризация delta на 2-кручении подтверждена независимо", ok)

def rnd_point():
    return add(add(mul(random.randint(-8, 8), MY_G),
                   mul(random.randint(0, 1), tor[1])),
               mul(random.randint(0, 1), tor[2]))

random.seed(1234)
bad = 0
for _ in range(200):
    P = rnd_point()
    Q = rnd_point()
    if delta(add(P, Q)) != prod_cls(delta(P), delta(Q)):
        bad += 1
        print("     НАРУШЕНИЕ ГОМОМОРФНОСТИ:", P, Q)
check("delta(P+Q) = delta(P)delta(Q) на 200 случайных парах (нарушений %d)" % bad, bad == 0)

bad = 0
for _ in range(100):
    d = delta(rnd_point())
    if sqfree(QQ(d[0]) * QQ(d[1]) * QQ(d[2])) != 1:
        bad += 1
check("произведение трёх координат delta — всегда квадрат (100 точек)", bad == 0)
print("  контроль: sqfree(1*113*113) =", sqfree(QQ(113) * QQ(113)),
      "-> требуемый класс НЕ отсекается этим тривиальным условием")

# ========================================= 6. ГРАНИЦЫ РАНГА (независимые)
hdr("6. ГРАНИЦЫ РАНГА E(Q) — НЕЗАВИСИМЫЕ ВЫЗОВЫ")
rank_upper = None
Emin = E.minimal_model()
print("  минимальная модель:", Emin.a_invariants())

try:
    rb = E.rank_bounds()
    print("  Sage E.rank_bounds()  =", rb)
    rank_upper = ZZ(rb[1])
except Exception as ex:
    print("  E.rank_bounds() упал:", ex)

try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    mw = mwrank_EllipticCurve(list(Emin.a_invariants()))
    mw.two_descent(verbose=0)
    print("  eclib rank            =", mw.rank())
    print("  eclib rank_bound      =", mw.rank_bound())
    print("  eclib certain         =", mw.certain())
    print("  eclib selmer_rank     =", mw.selmer_rank())
    if mw.certain():
        rb2 = ZZ(mw.rank_bound())
        rank_upper = rb2 if rank_upper is None else min(rank_upper, rb2)
except Exception as ex:
    print("  eclib mwrank упал:", ex)

try:
    pr = Emin.pari_curve().ellrank()
    print("  PARI ellrank          =", pr, " (формат [lower, upper, sha_exp, points])")
    rb3 = ZZ(pr[1])
    rank_upper = rb3 if rank_upper is None else min(rank_upper, rb3)
except Exception as ex:
    print("  PARI ellrank упал:", ex)

try:
    print("  E.selmer_rank()       =", E.selmer_rank())
except Exception as ex:
    print("  E.selmer_rank() упал:", ex)
try:
    print("  E.analytic_rank()     =", E.analytic_rank(), " (справочно, не в доказательстве)")
except Exception as ex:
    print("  analytic_rank упал:", ex)

print("\n  ИТОГОВАЯ верхняя граница ранга, которую я принимаю:", rank_upper)
check("верхняя граница rank E(Q) <= 1 подтверждена независимо",
      rank_upper is not None and rank_upper <= 1)

# ========================================= 7. ВОСЕМЬ КЛАССОВ — МОЙ СЧЁТ
hdr("7. ТАБЛИЦА КЛАССОВ delta(a*G + c*T1 + d*T2) — МОЙ СЧЁТ, МОЙ ГЕНЕРАТОР")
rows = []
for (a, c, d) in itertools.product([0, 1], repeat=3):
    P = add(add(mul(a, MY_G), mul(c, tor[1])), mul(d, tor[2]))
    rows.append(((a, c, d), P, delta(P)))

print("  мой генератор G =", MY_G)
print()
print("  | (a,c,d) |        x(P)        |    d1   |   d2   |   d3   |")
print("  |---------|--------------------|---------|--------|--------|")
for (acd, P, dd) in rows:
    xs = "O" if P is O else str(P[0])
    print("  | %s | %-18s | %7s | %6s | %6s |" % (str(acd), xs, dd[0], dd[1], dd[2]))

classes = [dd for (_, _, dd) in rows]
check("все 8 классов ПОПАРНО РАЗЛИЧНЫ (иначе образ неполон)", len(set(classes)) == 8)
if len(set(classes)) != 8:
    from collections import Counter
    print("     повторы:", [k for k, v in Counter(classes).items() if v > 1])
check("ТРЕБУЕМОГО класса (1,113,113) СРЕДИ НИХ НЕТ", required not in set(classes))
if required in set(classes):
    print("     !!! НАЙДЕН — заявление Codex ОПРОВЕРГНУТО !!!")
check("8 классов образуют подгруппу (замкнутость на 64 произведениях)",
      all(prod_cls(A, B) in set(classes) for A in classes for B in classes))

# ============================ 8. АТАКА НА ПОЛНОТУ
hdr("8. АТАКА: ищу точку E(Q) с классом ВНЕ этих восьми")
extra = set()
for (x, y) in pts:
    extra.add(delta((x, y)))
    if y != 0:
        extra.add(delta((x, -y)))
print("  классы ВСЕХ найденных перебором точек (%d шт.):" % len(pts))
for c in sorted(extra):
    print("     ", c)
check("все найденные перебором точки дают классы из моих восьми", extra <= set(classes))

big = set()
for k in range(-60, 61):
    for c in (0, 1):
        for d in (0, 1):
            big.add(delta(add(add(mul(k, MY_G), mul(c, tor[1])), mul(d, tor[2]))))
print("  классы для k in [-60,60] со всеми сдвигами кручения:", len(big), "штук")
check("они лежат в моих восьми и их ровно 8",
      big <= set(classes) and len(big) == 8)

# ================= 9. ПРЯМАЯ АТАКА НА C: поиск контрпримера
hdr("9. ПРЯМАЯ АТАКА НА C: поиск рационального t с тремя квадратами")
CSCAN = r'''
from math import isqrt, gcd
m2, n2, sv = 225, 1, 113
hits = []
two = 0
for q in range(1, NB + 1):
    qq = q * q
    for p in range(0, NB + 1):
        if gcd(p, q) != 1:
            continue
        pp = p * p
        A = m2 * qq + n2 * pp
        Bv = n2 * qq + m2 * pp
        ra = isqrt(A)
        if ra * ra != A:
            continue
        rb = isqrt(Bv)
        if rb * rb != Bv:
            continue
        two += 1
        D = sv * (pp + qq)
        rd = isqrt(D)
        if rd * rd == D:
            hits.append((p, q))
HITS = hits
TWO = two
'''
ns2 = {"NB": 1200}
exec(CSCAN, ns2)
print("  перебор t=p/q, 0<=p<=1200, 1<=q<=1200, gcd(p,q)=1")
print("  пар, где F0 и F8 обе квадраты (контроль живости перебора):", ns2["TWO"])
print("  полных точек C(Q):", len(ns2["HITS"]))
if ns2["HITS"]:
    print("  !!! КОНТРПРИМЕР:", ns2["HITS"][:10])
check("контрпримеров не найдено (отсутствие находки — НЕ доказательство)",
      len(ns2["HITS"]) == 0)

# ================================ 10. БЕСКОНЕЧНОСТЬ
hdr("10. ТОЧКИ НАД t = БЕСКОНЕЧНОСТЬ")
Rz = PolynomialRing(QQ, 'z')
zv = Rz.gen()
G0 = n**2 + m**2 * zv**2
G4 = s * (1 + zv**2)
G8 = m**2 + n**2 * zv**2
print("  z=1/t, U_i=u_i/t;  при z=0:  U0^2=%s, U4^2=%s, U8^2=%s"
      % (G0(0), G4(0), G8(0)))
check("U0^2 = n^2 — квадрат", QQ(G0(0)).is_square())
check("U8^2 = m^2 — квадрат", QQ(G8(0)).is_square())
check("U4^2 = s = 113 НЕ квадрат => над t=oo нет рациональных точек",
      not QQ(G4(0)).is_square())
print("  поле вычетов всех мест над t=oo содержит Q(sqrt 113) => степень >= 2.")

# ====================================== 11. ТОРСОР D И ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ
hdr("11. ТОРСОР D ДЛЯ КЛАССА (1,113,113) И ЛОКАЛЬНЫЕ ТОЧКИ C")
c12 = e2 - e1
c13 = e3 - e1
print("  D:  z1^2 - 113 z2^2 = %s" % c12)
print("      z1^2 - 113 z3^2 = %s" % c13)
check("совпадает с D у Codex (-5695200, 25312)", (c12, c13) == (QQ(-5695200), QQ(25312)))
check("(15 u4)^2 = X-e1", (m ** 2) * F4 - (Xt - e1) == 0)
check("113*(15 u0)^2 = X-e2", 113 * (m ** 2) * F0 - (Xt - e2) == 0)
check("113*u8^2 = X-e3", 113 * F8 - (Xt - e3) == 0)
print("  => C(Q_v) непусто влечёт D(Q_v) непусто для каждого места v.")

def is_sq_Qp(a, p):
    """Мой критерий квадрата в Q_p* для рационального a != 0."""
    a = QQ(a)
    if a == 0:
        raise ValueError
    v = a.valuation(p)
    u = a / QQ(p) ** v            # единица
    if v % 2 != 0:
        return False
    num = ZZ(u.numerator())
    den = ZZ(u.denominator())
    if p == 2:
        # u = num/den, оба нечётны; квадрат <=> num*den = 1 mod 8
        return (num * den) % 8 == 1
    inv = den.inverse_mod(p)
    r = (num % p) * inv % p
    if r == 0:
        raise ValueError("единица делится на p")
    return kronecker(r, p) == 1

# самотест критерия против Sage
random.seed(99)
ok = True
for p in [2, 3, 5, 7, 11, 13, 113]:
    K = Qp(p, 40)
    for _ in range(60):
        a = QQ(random.randint(-500, 500)) / QQ(random.randint(1, 60))
        if a == 0:
            continue
        if is_sq_Qp(a, p) != K(a).is_square():
            ok = False
            print("     расхождение критерия:", a, p)
check("мой критерий квадрата в Q_p совпал с Sage на 420 случаях", ok)

print("\n  локальные точки на C: ищу целое t с тремя квадратами в Q_p")
PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61,
          67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139]
missing = []
for p in PRIMES:
    hit = None
    for tt in range(0, 400):
        v0 = QQ(m**2 + n**2 * tt**2)
        v4 = QQ(s) * QQ(1 + tt**2)
        v8 = QQ(n**2 + m**2 * tt**2)
        if is_sq_Qp(v0, p) and is_sq_Qp(v4, p) and is_sq_Qp(v8, p):
            hit = tt
            break
    print("     p=%-4s C(Q_p): %s" % (p, ("t=%s" % hit) if hit is not None else "не найдено при 0<=t<400"))
    if hit is None:
        missing.append(p)
check("для всех проверенных p найдена явная точка C(Q_p)", not missing)
print("     C(R): t=0 даёт 225, 113, 1 — все > 0 => непусто")
print("  (замечание: это подтверждает локальную разрешимость только для перечисленных p;")
print("   общий аргумент Codex для p>=101 через Хассе-Вейля я НЕ пересчитывал полностью.)")

# ========================================= 12. СВЕРКА С ТАБЛИЦЕЙ CODEX
hdr("12. СВЕРКА МОЕЙ ТАБЛИЦЫ С ТАБЛИЦЕЙ CODEX")
codex_rows = {
    (0, 0, 0): (1, 1, 1),
    (0, 0, 1): (-1582, 226, -7),
    (0, 1, 0): (-1, 1582, -1582),
    (0, 1, 1): (1582, 7, 226),
    (1, 0, 0): (-14, 2, -7),
    (1, 0, 1): (113, 113, 1),
    (1, 1, 0): (14, 791, 226),
    (1, 1, 1): (-113, 14, -1582),
}
codex_set = set(tuple(ZZ(u) for u in v) for v in codex_rows.values())
mine_set = set(classes)
print("  моё множество классов  :", sorted(mine_set))
print("  Codex множество классов:", sorted(codex_set))
check("множества классов совпадают", mine_set == codex_set)
check("КЛЮЧЕВОЕ: (1,113,113) отсутствует и у меня, и у Codex",
      required not in mine_set and required not in codex_set)
check("все 8 строк Codex квадратсвободны (как он заявляет)",
      all(sqfree(QQ(u)) == ZZ(u) for v in codex_rows.values() for u in v))
check("у Codex тоже 8 различных строк", len(codex_set) == 8)
check("множество Codex замкнуто относительно умножения",
      all(prod_cls(A, B) in codex_set for A in codex_set for B in codex_set))
same_rows = all(tuple(ZZ(u) for u in codex_rows[acd]) == dd for (acd, P, dd) in rows)
print("  построчное совпадение (a,c,d) -> класс:", same_rows)

GC = (QQ(-75825), QQ(146764800))
check("точка Codex G=(-75825,146764800) действительно на E (мой счёт)",
      GC[1] ** 2 == f(GC[0]))
print("  delta(G_codex) =", delta(GC), "  (Codex пишет (-14, 2, -7))")
check("delta(G_codex) = (-14,2,-7) по моему счёту",
      delta(GC) == (ZZ(-14), ZZ(2), ZZ(-7)))
rel = None
for k in range(-8, 9):
    for c in (0, 1):
        for d in (0, 1):
            P = add(add(mul(k, MY_G), mul(c, tor[1])), mul(d, tor[2]))
            if P is not O and P == GC:
                rel = (k, c, d)
print("  связь точек: G_codex = %s  (в базисе моего G и кручения)" % (rel,))

# ЛОВУШКА С ПОРЯДКОМ КООРДИНАТ
hdr("13. ЛОВУШКА: ЧУВСТВИТЕЛЬНОСТЬ К ПОРЯДКУ КОРНЕЙ (e1,e2,e3)")
from itertools import permutations
for perm in permutations([0, 1, 2]):
    tgt = tuple(required[i] for i in perm)
    print("  перестановка %s -> класс %s : %s"
          % (perm, tgt, "ЕСТЬ В ОБРАЗЕ (!!)" if tgt in mine_set else "нет"))
print("  Порядок (e1,e2,e3) = (-b, -s m^4, -s n^4) критичен: (113,113,1) В ОБРАЗЕ,")
print("  а (1,113,113) — нет.  Порядок зафиксирован тождествами раздела 1 (символьно).")

# ============================================================ ИТОГ
hdr("ИТОГ")
if FAIL:
    print("  ПРОВАЛЕННЫЕ ПРОВЕРКИ (%d):" % len(FAIL))
    for x in FAIL:
        print("    -", x)
else:
    print("  Все проверки пройдены.")
print("""
  Логика закрытия (моя формулировка):
    delta: E(Q) -> (Q*/Q*^2)^3 — гомоморфизм, 2E(Q) содержится в ядре,
    значит образ есть образ конечной группы E(Q)/2E(Q) порядка 2^(r+2).
    Верхняя граница rank <= 1  =>  |образ| <= 8.
    Я предъявил 8 РАЗЛИЧНЫХ классов от явных рациональных точек
    =>  образ = ровно эти восемь.
    Требуемый класс (1,113,113) среди них отсутствует
    =>  конечных точек C(Q) нет.  Над t=oo нужен sqrt(113) в Q — нет.
    =>  C(Q) пусто.
  Единственное звено, опирающееся на ПО: ВЕРХНЯЯ ГРАНИЦА РАНГА.
""")
