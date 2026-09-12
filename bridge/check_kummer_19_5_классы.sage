# -*- coding: utf-8 -*-
# ПОПЫТКА СЛОМАТЬ заявление Codex: C_{19,5}(Q) = пусто (необходимый класс Куммера (1,193,193)).
# Всё считается с нуля: s, b, корни, тождества, СОБСТВЕННАЯ неторсионная точка,
# СОБСТВЕННАЯ арифметика на кривой, СОБСТВЕННОЕ перечисление классов delta,
# и отдельная атака на единственное узкое место — верхнюю границу ранга.
# Claude (независимая проверка), 2026-09-12.   Запуск: sage check_kummer_19_5_классы.sage

from sage.all import *
import sys

def hdr(t):
    print("\n" + "="*78); print(t); print("="*78); sys.stdout.flush()

# ------------------------------------------------------------------ утилиты
def sqclass(x):
    """Квадратсвободный представитель СО ЗНАКОМ класса x в Q*/Q*^2. Своя реализация."""
    x = QQ(x)
    if x == 0:
        raise ValueError("нулевой аргумент квадратного класса")
    z = ZZ(x.numerator()) * ZZ(x.denominator())      # p/q ~ p*q
    out = 1
    for (p, e) in factor(abs(z)):
        if e % 2 == 1:
            out *= p
    return out if z > 0 else -out

def tmul(A, B):
    return tuple(sqclass(QQ(a)*QQ(c)) for a, c in zip(A, B))

# ================================================================== [1] МОДЕЛЬ
hdr("[1] МОДЕЛЬ — считаю сам, чужие числа только сверяю")
m, n = 19, 5
assert gcd(m, n) == 1
s = QQ(m**2 + n**2)/2
b = s*m**2*n**2
print("  s = (m^2+n^2)/2 =", s, "   Codex 193 ->", s == 193)
print("  b = s*m^2*n^2   =", b, "   Codex 1741825 ->", b == 1741825)
e1, e2, e3 = -b, -s*m**4, -s*n**4
print("  e1,e2,e3 =", e1, e2, e3)
print("  Codex (-1741825,-25151953,-120625) ->",
      (e1, e2, e3) == (-1741825, -25151953, -120625))
assert len({e1, e2, e3}) == 3
print("  s квадратсвободно:", s == sqclass(s), "  s простое:", ZZ(s).is_prime())

# ================================================================== [2] ТОЖДЕСТВА
hdr("[2] ТОЧНЫЕ ТОЖДЕСТВА (полиномиально, X = b t^2)")
R = PolynomialRing(QQ, 4, names=('t', 'u0', 'u4', 'u8'))
(t, u0, u4, u8) = R.gens()
F0 = m**2 + n**2*t**2
F4 = s*(1 + t**2)
F8 = n**2 + m**2*t**2
X = b*t**2
print("  X-e1 - (mn)^2*F4 =", (X-e1) - (m*n)**2*F4)
print("  X-e2 - s*m^2*F0  =", (X-e2) - s*m**2*F0)
print("  X-e3 - s*n^2*F8  =", (X-e3) - s*n**2*F8)
assert (X-e1) - (m*n)**2*F4 == 0
assert (X-e2) - s*m**2*F0 == 0
assert (X-e3) - s*n**2*F8 == 0
print("  V=b*u0u4u8:  произведение трёх разностей - V^2 =",
      ((m*n)**2*u4**2)*(s*m**2*u0**2)*(s*n**2*u8**2) - (b*u0*u4*u8)**2)
req = (sqclass((m*n)**2), sqclass(s*m**2), sqclass(s*n**2))
print("  НЕОБХОДИМЫЙ класс любой конечной точки C:", req,
      "  Codex (1,193,193) ->", req == (1, 193, 193))
print("  произведение координат req квадрат:", sqclass(req[0]*req[1]*req[2]) == 1)
print("  F0,F4,F8 > 0 при любом вещественном t  ->  u_i != 0, X != e_i  ->  delta определён")

# ================================================================== [3] БЕСКОНЕЧНОСТЬ
hdr("[3] ТОЧКИ НАД t = бесконечность")
print("  z=1/t, U_i=u_i/t:  U0^2=%s, U4^2=%s, U8^2=%s при z=0" % (n**2, s, m**2))
print("  U4^2 = s = %s квадрат в Q? %s" % (s, QQ(s).is_square()))
assert not QQ(s).is_square()
print("  -> точек над бесконечностью НЕТ")

# ================================================================== [4] КРИВАЯ
hdr("[4] КРИВАЯ E")
a2 = -(e1+e2+e3); a4 = e1*e2+e1*e3+e2*e3; a6 = -e1*e2*e3
E = EllipticCurve(QQ, [0, a2, 0, a4, a6])
print("  E:", E)
print("  минимальная модель:", E.minimal_model().ainvs())
print("  disc =", factor(E.discriminant()))
print("  conductor =", E.conductor(), "=", factor(E.conductor()))
print("  кручение:", E.torsion_subgroup().invariants())
assert E.torsion_subgroup().invariants() == (2, 2)
print("  |E(Q)[2]| = 4  ->  |E(Q)/2E(Q)| = 2^(r+2)")

# ================================================================== [5] РАНГ
hdr("[5] РАНГ — ЕДИНСТВЕННОЕ УЗКОЕ МЕСТО, БЬЮ ПО НЕМУ")
print("  -- маршрут A: 2-спуск --")
try:
    print("  Sage E.rank_bounds() (Simon)      ->", E.rank_bounds())
except Exception as ex:
    print("  Sage rank_bounds ОШИБКА:", ex)
try:
    print("  Sage E.rank() (eclib/mwrank)      -> ", E.rank())
except Exception as ex:
    print("  Sage E.rank() НЕ СМОГ:", str(ex).splitlines()[0])
print("  Sage E.selmer_rank()              ->", E.selmer_rank(),
      " (= rank + dim Sha[2] + 2)")
mypts = []
for eff in [0, 1, 2, 4]:
    try:
        out = pari(E).ellrank(eff)
        print("  PARI ellrank(effort=%s)            -> %s" % (eff, out))
        for pt in out[3]:
            mypts.append((QQ(pt[0]), QQ(pt[1])))
    except Exception as ex:
        print("  PARI ellrank effort %s ОШИБКА: %s" % (eff, ex))

print("\n  -- маршрут B: аналитический (независим от 2-спуска) --")
print("  root number (знак фунц. уравнения) =", E.root_number(), " -> аналит. ранг НЕЧЁТЕН")
try:
    ub = E.analytic_rank_upper_bound(max_Delta=2.0, adaptive=True, root_number=-1)
    print("  analytic_rank_upper_bound         ->", ub, " (УСЛОВНО ПО GRH, см. docstring Sage)")
except Exception as ex:
    print("  analytic_rank_upper_bound ОШИБКА:", ex)
print("  analytic_rank (численно)          ->", E.analytic_rank())
print("  Gross-Zagier + Kolyvagin: аналит. ранг <=1  =>  алг. ранг = аналит. ранг (безусловная теорема)")

print("\n  -- маршрут C: спуск по 2-изогении (другая реализация, Sage/R.Miller) --")
try:
    from sage.schemes.elliptic_curves.descent_two_isogeny import two_descent_by_two_isogeny
    n1, n2, n1p, n2p = two_descent_by_two_isogeny(E)
    print("  #Sel^phi' = %s, #Sel^phi = %s  ->  rank <= log2(%s)+log2(%s)-2 = %s"
          % (n2, n2p, n2, n2p, ZZ(n2).exact_log(2)+ZZ(n2p).exact_log(2)-2))
    print("  (эта реализация даёт только <=3 — она НЕ подтверждает границу 1)")
except Exception as ex:
    print("  ОШИБКА:", ex)

print("\n  -- маршрут D: матрица спаривания Касселса-Тейта, которую печатает PARI (debug 3) --")
CT = matrix(GF(2), [[0,1,1,1,1],[1,0,0,1,0],[1,0,0,1,0],[1,1,1,0,1],[1,0,0,1,0]])
print(CT)
print("  ранг над F2 =", CT.rank(), " (альтернированная форма -> ранг чётен: ",
      CT.rank() % 2 == 0, ")")
print("  радикал спаривания имеет размерность", 5 - CT.rank(), "-> |радикал| =", 2**(5-CT.rank()))
print("  ТЕОРЕМА: образ delta(E(Q)) лежит в радикале CT. Значит |образ| <= 8 и rank <= 1.")

# ================================================================== [6] МОЯ ТОЧКА
hdr("[6] СОБСТВЕННАЯ НЕТОРСИОННАЯ ТОЧКА")
mypts = list(dict.fromkeys(mypts))
print("  мои прогоны PARI дали точки:")
for P in mypts:
    print("    ", P, "  h =", E(P[0], P[1]).height())
CODEX = (QQ(-28561963657)/QQ(1369), QQ(-2089086742828800)/QQ(50653))
MY = None
for P in mypts:
    if P[0] != CODEX[0]:
        MY = P; break
if MY is None:
    MY = mypts[0]
print("  БЕРУ СВОЮ ТОЧКУ  G =", MY)
print("  она отличается от точки Codex:", MY[0] != CODEX[0])
if len(mypts) >= 2:
    HM = E.height_pairing_matrix([E(p[0], p[1]) for p in mypts])
    print("  матрица высот моих точек, det =", HM.determinant(),
          " -> они ЗАВИСИМЫ, лишнего ранга нет" if abs(HM.determinant()) < 1e-6 else " -> НЕЗАВИСИМЫ!")

# ================================================================== [7] СВОЯ АРИФМЕТИКА
hdr("[7] СОБСТВЕННАЯ АРИФМЕТИКА НА КРИВОЙ (не Sage-овская)")
A2, A4, A6 = QQ(a2), QQ(a4), QQ(a6)

def on_curve(P):
    if P is None: return True
    x, y = P; return y*y == x**3 + A2*x*x + A4*x + A6
def neg(P): return None if P is None else (P[0], -P[1])
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if y1 + y2 == 0: return None
        lam = (3*x1*x1 + 2*A2*x1 + A4)/(2*y1)
    else:
        lam = (y2 - y1)/(x2 - x1)
    x3 = lam*lam - A2 - x1 - x2
    return (x3, lam*(x1 - x3) - y1)
def mul(k, P):
    k = ZZ(k)
    if k < 0: return mul(-k, neg(P))
    Rr, Q = None, P
    while k > 0:
        if k & 1: Rr = add(Rr, Q)
        Q = add(Q, Q); k >>= 1
    return Rr

T1 = (e1, QQ(0)); T2 = (e2, QQ(0)); T3 = (e3, QQ(0))
for nm, P in [("G", MY), ("T1", T1), ("T2", T2), ("T3", T3)]:
    print("  %-3s подстановкой на кривой: %s" % (nm, on_curve(P)))
    assert on_curve(P)
assert add(T1, T2) == T3
print("  T1+T2 == T3:", True, "   2*T1 == O:", mul(2, T1) is None)
ok = True
for k in [2, 3, 5, 7, 11]:
    sg = k*E(MY[0], MY[1])
    if (QQ(sg[0]), QQ(sg[1])) != mul(k, MY): ok = False
print("  моё сложение == Sage на k*G, k=2,3,5,7,11:", ok)
assert ok

# ================================================================== [8] delta
hdr("[8] СОБСТВЕННОЕ ПЕРЕЧИСЛЕНИЕ КЛАССОВ delta")
EE = [e1, e2, e3]
def delta(P):
    if P is None: return (1, 1, 1)
    x, y = P; out = []
    for i in range(3):
        if x == EE[i]:
            j, k = [q for q in range(3) if q != i]
            out.append(sqclass((EE[i]-EE[j])*(EE[i]-EE[k])))
        else:
            out.append(sqclass(x - EE[i]))
    return tuple(out)

rows = []
for a in [0, 1]:
    for c in [0, 1]:
        for d in [0, 1]:
            P = add(mul(a, MY), add(mul(c, T1), mul(d, T2)))
            assert on_curve(P)
            rows.append(((a, c, d), P, delta(P)))
print("  | (a,c,d) |     d1 |     d2 |     d3 | d1*d2*d3 квадрат |")
print("  |---------|--------|--------|--------|------------------|")
for (k, P, D) in rows:
    print("  | %s | %6s | %6s | %6s | %s |"
          % (k, D[0], D[1], D[2], sqclass(D[0]*D[1]*D[2]) == 1))
cls = [D for (_, _, D) in rows]
print("\n  различных классов: %d из %d" % (len(set(cls)), len(cls)))
print("  требуемый класс", req, "присутствует:", req in set(cls))

print("\n  ФАКТОРИЗАЦИИ для ручной проверки строки (1,0,0)  (P = G):")
for i in range(3):
    v = MY[0] - EE[i]
    print("    x-e%d = %s = %s   -> класс %s" % (i+1, v, factor(v), sqclass(v)))
print("  ФАКТОРИЗАЦИИ кручения:")
print("    T1: (e1-e2)(e1-e3)=%s, e1-e2=%s, e1-e3=%s"
      % (factor((e1-e2)*(e1-e3)), factor(e1-e2), factor(e1-e3)))
print("    T2: e2-e1=%s, (e2-e1)(e2-e3)=%s, e2-e3=%s"
      % (factor(e2-e1), factor((e2-e1)*(e2-e3)), factor(e2-e3)))

# ================================================================== [9] ГОМОМОРФНОСТЬ
hdr("[9] ГОМОМОРФНОСТЬ delta НА ВСЕХ 64 ПАРАХ + ГРУППОВАЯ СТРУКТУРА")
P8 = {k: P for (k, P, _) in rows}; D8 = {k: D for (k, _, D) in rows}
bad = 0
for k1 in P8:
    for k2 in P8:
        if delta(add(P8[k1], P8[k2])) != tmul(D8[k1], D8[k2]):
            bad += 1; print("   НАРУШЕНИЕ", k1, k2)
print("  нарушений:", bad, "из 64")
print("  delta(2P) == (1,1,1) для всех восьми P:",
      all(delta(mul(2, P)) == (1, 1, 1) for (_, P, _) in rows))
G8 = set(cls)
print("  множество замкнуто относительно умножения:",
      all(tmul(x, y) in G8 for x in G8 for y in G8))
print("  элементарная 2-группа:", all(tmul(x, x) == (1, 1, 1) for x in G8))
print("  req в образе:", req in G8)

# ================================================================== [10] СВЕРКА С CODEX
hdr("[10] СВЕРКА С ТАБЛИЦЕЙ CODEX")
print("  точка Codex на кривой (подстановкой):", on_curve(CODEX))
rel = None
for k in [-1, 1]:
    for T in [None, T1, T2, T3]:
        if add(mul(k, MY), T) == CODEX:
            rel = (k, T)
print("  связь моей точки с точкой Codex:  Gc = %s" % (rel,))
cls_c = []
for a in [0, 1]:
    for c in [0, 1]:
        for d in [0, 1]:
            cls_c.append(delta(add(mul(a, CODEX), add(mul(c, T1), mul(d, T2)))))
codex_table = [(1,1,1), (-1,4053,-4053), (-4053,386,-42), (4053,42,386),
               (-3,386,-1158), (3,42,14), (1351,1,1351), (-1351,4053,-3)]
print("  мой пересчёт по точке Codex совпал с его таблицей как множество:",
      set(cls_c) == set(codex_table))
print("  множество МОИХ классов == множеству Codex:", set(cls) == set(codex_table))
print("  (построчный порядок отличается — это лишь другой выбор генератора)")

# ================================================================== [11] SELMER
hdr("[11] ЛЕЖИТ ЛИ req В 2-СЕЛМЕРЕ? (если да — ранговая граница НЕСУЩАЯ)")
print("  торсор D:  z1^2 - 193 z2^2 = %s w^2 ;  z1^2 - 193 z3^2 = %s w^2"
      % (e2-e1, e3-e1))
print("  e2-e1 =", factor(e2-e1), "   e3-e1 =", factor(e3-e1))
def loc_C(p, bound=4000):
    K = Qp(p, prec=60)
    for tt in range(0, bound):
        f0 = K(m**2 + n**2*tt**2); f4 = K(s*(1+tt**2)); f8 = K(n**2 + m**2*tt**2)
        if f0 == 0 or f4 == 0 or f8 == 0: continue
        if f0.is_square() and f4.is_square() and f8.is_square():
            return tt
    return None
print("  явные целые t с квадратными F0,F4,F8 в Q_p (даёт точку C(Q_p) и потому D(Q_p)):")
for p in [2, 3, 5, 7, 19, 193, 11, 13, 17, 23, 29, 31, 37, 41, 43]:
    print("    p=%-4s t=%s" % (p, loc_C(p)))
print("  над R: t=0 даёт F0=361>0, F4=193>0, F8=25>0 -> точки есть")
print("  ВЫВОД: req ВСЮДУ ЛОКАЛЬНО РАЗРЕШИМ -> req лежит в Sel_2 -> без границы ранга НИЧЕГО НЕ СЛЕДУЕТ.")

# ================================================================== [12] ПРЯМОЙ ПОИСК
hdr("[12] ПРЯМОЙ ПОИСК ТОЧЕК НА ТОРСОРЕ req (попытка сломать вывод в лоб)")
print("  193 | (e2-e1) и 193 | (e3-e1)  ->  193 | z1, z1 = 193*a:")
print("     193 a^2 + %s w^2 = z2^2   и   193 a^2 - %s w^2 = z3^2"
      % (-(e2-e1)//193, (e3-e1)//193))
hits = []
for w in range(1, 601):
    c1 = (-(e2-e1)//193)*w*w; c2 = ((e3-e1)//193)*w*w
    for aa in range(0, 601):
        if gcd(aa, w) != 1: continue
        tq = 193*aa*aa
        if tq - c2 < 0: continue
        if ZZ(tq - c2).is_square() and ZZ(tq + c1).is_square():
            hits.append((aa, w)); print("   НАЙДЕНА ТОЧКА:", aa, w)
print("  a,w <= 600: найдено", len(hits), "(отсутствие НЕ есть доказательство отсутствия)")

# ================================================================== ИТОГ
hdr("ИТОГ")
print("  s =", s, "  b =", b, "  требуемый класс:", req)
print("  моих различных классов: %d/8 ; требуемый среди них: %s" % (len(set(cls)), req in set(cls)))
print("  нарушений гомоморфности:", bad)
print("  множества классов (моё и Codex) совпадают:", set(cls) == set(codex_table))
print("  граница ранга <=1: PARI ellrank + матрица Касселса ранга 2 (радикал dim 3)")
print("  вторая, независимая линия: root number -1 + GRH-оценка аналит.ранга 1 + Gross-Zagier-Kolyvagin")
print("  НЕ подтвердили границу 1: eclib/mwrank (<=3), Simon (<=3), спуск по 2-изогении (<=3)")
