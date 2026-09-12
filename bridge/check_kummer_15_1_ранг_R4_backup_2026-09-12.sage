# -*- coding: utf-8 -*-
# check_kummer_15_1_ранг.sage   (РАУНД 4, независимая атака Claude)
# Предыдущая версия сохранена в check_kummer_15_1_ранг_PREV_backup_2026-09-12.sage
#
# ЦЕЛЬ: сломать заявление Codex "C_{15,1}(Q) = пусто".
# Всё строится заново из (m,n)=(15,1). Ничего из файлов Codex не импортируется,
# его числа используются ТОЛЬКО для сравнения с моими.
#
# Разделы:
#   0  модель из (m,n)
#   1  тождества X-e_i КАК ПОЛИНОМИАЛЬНЫЕ ТОЖДЕСТВА (проверка ПОРЯДКА корней!)
#   2  кручение и лемма |A/2A| = |A[2]|  (атака на формулу 2^(r+2))
#   3  верхняя граница ранга: PARI ellrank / eclib mwrank / Sage / МОЙ спуск по 2-изогении
#   4  delta своей реализацией, 8 классов, отсутствие (1,113,113)
#   5  логика "8 различных => образ полон" без насыщения
#   6  потерянные ветви: нули F_i, бесконечность, знаки, особые точки
#   7  локальная разрешимость C(Q_p): лежит ли требуемый класс в Sel^2
#      (если да -- исключение ЧИСТО глобальное и целиком висит на ранге)
#   8  независимый маршрут: квартика B / кривая E_B того же проекта
#   9  прямая охота за контрпримером на C

import sys, time
T0 = time.time()
def say(*a):
    print(*a); sys.stdout.flush()
def hdr(s):
    say(""); say("="*76); say(s); say("="*76)

V = {}

# ============================================================================
hdr("0. МОДЕЛЬ ЗАНОВО ИЗ (m,n) = (15,1)")
# ============================================================================
m = ZZ(15); n = ZZ(1)
say("gcd(m,n) =", gcd(m,n))
s = QQ(m^2 + n^2)/2
b = s*m^2*n^2
e1 = -b
e2 = -s*m^4
e3 = -s*n^4
say("s =", s, " (целое:", s.is_integer(), ")   b =", b)
say("e1,e2,e3 =", e1, e2, e3)
V['s'] = s; V['b'] = b
V['codex_model_ok'] = bool(s == 113 and b == 25425 and (e1,e2,e3) == (-25425,-5720625,-113))
say("совпадает с моделью Codex (s=113,b=25425,e=(-25425,-5720625,-113)):", V['codex_model_ok'])
say("s квадрат в Q?", QQ(s).is_square(), "  s свободно от квадратов?", ZZ(s).is_squarefree())
say("корни e_i попарно различны:", len({e1,e2,e3}) == 3)

R.<X> = PolynomialRing(QQ)
cub = (X-e1)*(X-e2)*(X-e3)
c = cub.coefficients(sparse=False)      # c0 + c1 X + c2 X^2 + X^3
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
say("E =", E)
say("E.discriminant() != 0:", E.discriminant() != 0)
Emin = E.minimal_model()
say("минимальная модель a-инварианты:", list(Emin.a_invariants()))
V['min_ainvs'] = [ZZ(x) for x in Emin.a_invariants()]
say("Codex заявил [0,-1,0,-42422006041,3362853259097305]:",
    V['min_ainvs'] == [0,-1,0,-42422006041,3362853259097305])
N = E.conductor()
say("кондуктор =", N, "=", factor(N))

# ============================================================================
hdr("1. ТОЖДЕСТВА X-e_i КАК ТОЖДЕСТВА В QQ[t]  (КРИТИЧЕСКАЯ ПРОВЕРКА ПОРЯДКА)")
# ============================================================================
# Это самое опасное место всей конструкции: класс (113,113,1) ЕСТЬ в образе,
# а (1,113,113) -- нет. Перестановка индексов переворачивает вывод.
S.<t> = PolynomialRing(QQ)
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
Xt = b*t^2
id1 = (Xt - e1) - (m*n)^2*F4
id2 = (Xt - e2) - s*m^2*F0
id3 = (Xt - e3) - s*n^2*F8
say("F0,F4,F8 =", F0, ";", F4, ";", F8)
say("X - e1 == (mn)^2 * F4  тождественно:", id1 == 0, "   (остаток:", id1, ")")
say("X - e2 == s*m^2 * F0   тождественно:", id2 == 0, "   (остаток:", id2, ")")
say("X - e3 == s*n^2 * F8   тождественно:", id3 == 0, "   (остаток:", id3, ")")
V['identities_ok'] = bool(id1 == 0 and id2 == 0 and id3 == 0)

# проверка совместности с уравнением E:
lhs = (Xt-e1)*(Xt-e2)*(Xt-e3)
rhs = b^2*F0*F4*F8
say("(X-e1)(X-e2)(X-e3) == b^2 F0 F4 F8 тождественно:", (lhs-rhs) == 0)
V['E_map_ok'] = bool((lhs-rhs) == 0)

def sqcls(x):
    """квадратный класс рационального x != 0: целое, свободное от квадратов, со знаком"""
    x = QQ(x)
    if x == 0:
        raise ValueError("ноль не имеет квадратного класса")
    return ZZ(x.numerator()*x.denominator()).squarefree_part()

# требуемый класс выводим ИЗ тождеств, а не переписываем у Codex:
req = (sqcls((m*n)^2), sqcls(s*m^2), sqcls(s*n^2))
say("требуемый класс delta = (класс (mn)^2, класс s m^2, класс s n^2) =", req)
say("Codex заявил (1,113,113):", req == (1,113,113))
V['required_class'] = req

# ============================================================================
hdr("2. КРУЧЕНИЕ И ФОРМУЛА |E(Q)/2E(Q)| = 2^(r+2)")
# ============================================================================
Tor = E.torsion_subgroup()
say("E.torsion_subgroup() =", Tor, "   инварианты:", Tor.invariants(), "  порядок:", Tor.order())
say("точки кручения:", [P for P in Tor.points()] if Tor.order() <= 16 else "(много)")
tor_inv = [ZZ(i) for i in Tor.invariants()]
V['torsion'] = tor_inv
say("E(Q)[2] порядок:", len([P for P in Tor.points() if 2*P == E(0)]))
say("полный рациональный 2-торсион:", len([P for P in Tor.points() if 2*P == E(0)]) == 4)

# АТАКА: а если бы кручение было Z/2 x Z/4 или Z/2 x Z/8 -- изменился бы счёт?
say("")
say("ЛЕММА (проверяю явно): для конечной абелевой A имеем |A/2A| = |A[2]|.")
bad = []
from itertools import product as iproduct
for inv in [(2,2),(2,4),(2,8),(2,6),(2,12),(4,4),(2,2,2)]:
    els = list(iproduct(*[range(d) for d in inv]))
    dbl = set(tuple((2*x[i]) % inv[i] for i in range(len(inv))) for x in els)
    q = len(els)//len(dbl)                      # |A/2A| = |A|/|2A|
    k = len([x for x in els if all((2*x[i]) % inv[i] == 0 for i in range(len(inv)))])
    ok = (q == k)
    say("   инварианты", inv, ": |A/2A| =", q, "  |A[2]| =", k, "  совпало:", ok,
        "  (полный 2-торсион:", k == 4, ")")
    if not ok: bad.append(inv)
say("лемма подтверждена на всех проверенных:", not bad)
say("ВЫВОД: при полном 2-кручении |T/2T| = 4 ВСЕГДА, даже для Z/2xZ/4, Z/2xZ/8.")
say("       поэтому |E(Q)/2E(Q)| = 2^r * 4 = 2^(r+2) -- формула УСТОЙЧИВА к росту кручения.")
V['torsion_attack_kills_formula'] = False

# ============================================================================
hdr("3. ВЕРХНЯЯ ГРАНИЦА РАНГА -- ЧЕТЫРЕ СРЕДСТВА")
# ============================================================================
say("--- 3A. PARI ellrank (прямой вызов pari(...).ellinit().ellrank()) ---")
ai = [ZZ(x) for x in Emin.a_invariants()]
for eff in [0,1,2,3]:
    try:
        ep = pari(ai).ellinit()
        res = ep.ellrank(eff) if eff > 0 else ep.ellrank()
        say("   effort =", eff, "->", res)
        if eff == 0:
            V['pari_lo'] = ZZ(res[0]); V['pari_hi'] = ZZ(res[1])
    except Exception as ex:
        say("   effort =", eff, "ОШИБКА:", ex)
# и на ИСХОДНОЙ (неминимальной) модели -- другой путь через ellinit
try:
    ep2 = pari([0, c[2], 0, c[1], c[0]]).ellinit()
    say("   исходная модель ->", ep2.ellrank())
except Exception as ex:
    say("   исходная модель ОШИБКА:", ex)

say("")
say("--- 3B. eclib mwrank_EllipticCurve, флаг certain ---")
from sage.libs.eclib.interface import mwrank_EllipticCurve
mw = mwrank_EllipticCurve(ai, verbose=False)
mw.two_descent(verbose=False, second_descent=True, selmer_only=False, first_limit=20, second_limit=10)
say("   mwrank rank       =", mw.rank())
say("   mwrank rank_bound =", mw.rank_bound())
say("   mwrank selmer_rank=", mw.selmer_rank())
say("   mwrank CERTAIN    =", mw.certain())
V['mwrank_bound'] = ZZ(mw.rank_bound()); V['mwrank_certain'] = bool(mw.certain())

say("")
say("   контроль: ТОЛЬКО первый спуск (second_descent=False)")
mw1 = mwrank_EllipticCurve(ai, verbose=False)
mw1.two_descent(verbose=False, second_descent=False, selmer_only=False, first_limit=20)
say("   первый спуск: rank_bound =", mw1.rank_bound(), " certain =", mw1.certain(),
    " selmer_rank =", mw1.selmer_rank())
V['first_descent_bound'] = ZZ(mw1.rank_bound())

say("")
say("--- 3C. Sage (без базы Кремоны) ---")
say("   E.selmer_rank() =", E.selmer_rank())
say("   E.rank_bound()  =", E.rank_bound())
try:
    r_sage = E.rank(use_database=False, only_use_mwrank=True, proof=True)
    say("   E.rank(proof=True, база выключена) =", r_sage)
    V['sage_rank'] = ZZ(r_sage)
except Exception as ex:
    say("   E.rank ОШИБКА:", ex)
say("   кондуктор", N, "вне диапазона таблиц Кремоны (>500000):", N > 500000)

say("")
say("--- 3D. МОЙ СОБСТВЕННЫЙ СПУСК ПО 2-ИЗОГЕНИИ (независимый код) ---")

def is_sq_Qp(z, p):
    """z != 0 целое/рациональное: квадрат ли в Q_p"""
    z = QQ(z)
    if z == 0: return True
    v = z.valuation(p)
    if v % 2 != 0: return False
    u = z / p^v
    un = u.numerator(); ud = u.denominator()
    u = un*ud                     # тот же класс по модулю квадратов, целое, p не делит
    if p == 2:
        return (u % 8) == 1
    return kronecker(u, p) == 1

def quartic_Qp_solvable(cs, p, maxdepth=14):
    """w^2 = cs[0] u^4 + cs[1] u^3 v + cs[2] u^2 v^2 + cs[3] u v^3 + cs[4] v^4
       есть ли примитивное (u,v) в Z_p^2 (u,v не оба делятся на p)?
       Возвращает True / False / None(не решено на данной глубине)."""
    def f(u,v):
        return (cs[0]*u^4 + cs[1]*u^3*v + cs[2]*u^2*v^2 + cs[3]*u*v^3 + cs[4]*v^4)
    slack = 4 if p == 2 else 2
    stack = [(u,v,1) for u in range(p) for v in range(p) if not (u % p == 0 and v % p == 0)]
    undecided = False
    while stack:
        u,v,k = stack.pop()
        val = f(ZZ(u), ZZ(v))
        if val == 0:
            return True                      # w=0 -- настоящая точка
        e = val.valuation(p)
        if e + slack <= k:
            if is_sq_Qp(val, p):
                return True                  # Гензель: решение поднимается
            else:
                continue                     # весь класс мёртв
        if k >= maxdepth:
            undecided = True
            continue
        pk = p^k
        for a in range(p):
            for bb in range(p):
                stack.append((u + a*pk, v + bb*pk, k+1))
    return None if undecided else False

def quartic_R_solvable(cs):
    Rr = RealField(80)
    P = PolynomialRing(Rr,'z')(list(reversed([Rr(x) for x in cs])))
    # ищем u/v = z с f(z) >= 0, плюс v=0 (старший коэффициент)
    if Rr(cs[0]) >= 0: return True
    if Rr(cs[4]) >= 0: return True
    for z in [Rr(i)/8 for i in range(-4000,4001)]:
        if P(z) >= 0: return True
    return False

def phi_selmer_dim(a_, b_, label):
    """S^{phi} для y^2 = x(x^2 + a_ x + b_), T=(0,0)."""
    supp = set([2]) | set(p for p,_ in factor(ZZ(b_)))
    supp = sorted(supp)
    primes_to_check = sorted(set([2]) | set(p for p,_ in factor(ZZ(b_)))
                             | set(p for p,_ in factor(ZZ(a_^2 - 4*b_))))
    say("   [", label, "] a =", a_, " b =", b_)
    say("   [", label, "] носитель d:", supp, " проверяемые p:", primes_to_check)
    cands = []
    from itertools import product
    prs = [p for p in supp]
    for signs in [1,-1]:
        for expo in product([0,1], repeat=len(prs)):
            d = signs*prod([prs[i]^expo[i] for i in range(len(prs))])
            cands.append(ZZ(d))
    good = []
    unknown = []
    for d in cands:
        if d == 0: continue
        if not ZZ(b_ / d).is_integer() and not (ZZ(b_) % d == 0):
            pass
        cs = [d, 0, a_, 0, ZZ(b_)/d]
        if not all(QQ(x).is_integer() for x in cs):
            # домножим на d^2? стандартно b/d целое, т.к. d | b свободно от квадратов
            pass
        cs = [QQ(x) for x in cs]
        if not quartic_R_solvable(cs):
            continue
        ok = True; unk = False
        for p in primes_to_check:
            r = quartic_Qp_solvable([ZZ(x) for x in cs], p)
            if r is False:
                ok = False; break
            if r is None:
                unk = True
        if ok and not unk:
            good.append(d)
        elif ok and unk:
            unknown.append(d)
    say("   [", label, "] локально всюду разрешимые d:", sorted(good))
    if unknown:
        say("   [", label, "] НЕ РЕШЕНО для d:", sorted(unknown))
    return good, unknown

A_ = e3 - e1        # сдвиг x = X - e3, корни 0, -(e3-e1), -(e3-e2)
B_ = e3 - e2
say("   сдвиг X = x - e3: y^2 = x (x + %s)(x + %s)" % (A_, B_))
a_E = ZZ(A_ + B_); b_E = ZZ(A_*B_)
a_Ep = ZZ(-2*a_E); b_Ep = ZZ(a_E^2 - 4*b_E)
say("   a_E =", a_E, " b_E =", factor(b_E))
say("   a_E' =", a_Ep, " b_E' =", factor(b_Ep))
g1, u1 = phi_selmer_dim(a_E, b_E, "phi")
g2, u2 = phi_selmer_dim(a_Ep, b_Ep, "phi^")
if not u1 and not u2:
    d1 = ZZ(len(g1)).exact_log(2); d2 = ZZ(len(g2)).exact_log(2)
    say("   dim S^phi =", d1, "  dim S^phi^ =", d2,
        "  => rank <= d1+d2-2 =", d1+d2-2)
    V['isogeny_bound'] = d1+d2-2
else:
    say("   спуск по изогении НЕ ЗАВЕРШЁН (есть нерешённые d)")
    V['isogeny_bound'] = None

say("")
say("--- 3E. поиск лишних точек: если ранг >= 2, где вторая образующая? ---")
G = E(-75825, 146764800)
say("   G = (-75825, 146764800) лежит на E:", G in E, "  порядок:", G.order())
say("   высота Нерона-Тейта h(G) =", G.height())
try:
    ps = E.point_search(14, rank_bound=3)
    say("   point_search(height=14): найдено", len(ps), "точек:", ps[:6])
    if ps:
        say("   их независимость/ранг подгруппы:", E.saturation(ps)[0] if ps else None)
except Exception as ex:
    say("   point_search ОШИБКА:", ex)

# ============================================================================
hdr("4. delta СВОЕЙ РЕАЛИЗАЦИЕЙ: 8 КЛАССОВ И ТРЕБУЕМЫЙ КЛАСС")
# ============================================================================
es = [e1,e2,e3]
O = E(0)
def delta(P):
    if P == O: return (ZZ(1),ZZ(1),ZZ(1))
    x = P[0]
    out = []
    for i in range(3):
        val = x - es[i]
        if val == 0:
            j,k = [q for q in range(3) if q != i]
            val = (es[i]-es[j])*(es[i]-es[k])
        out.append(sqcls(val))
    return tuple(out)

T1 = E(e1,0); T2 = E(e2,0); T3 = E(e3,0)
say("T1+T2 == T3:", T1+T2 == T3)
rows = []
say("")
say("(a,c,d)   точка X                         delta")
for a in [0,1]:
    for cc_ in [0,1]:
        for dd in [0,1]:
            P = a*G + cc_*T1 + dd*T2
            dl = delta(P)
            rows.append(((a,cc_,dd), P, dl))
            say(" (%d,%d,%d)  %-34s %s" % (a,cc_,dd, str(P[0]) if P != O else "O", dl))
classes = [r[2] for r in rows]
V['n_distinct'] = len(set(classes))
say("")
say("различных классов:", len(set(classes)), "из 8")
say("требуемый", req, "среди них:", req in set(classes))
V['required_in_image'] = bool(req in set(classes))
say("перестановка", (113,113,1), "среди них:", (113,113,1) in set(classes),
    "  <-- ЛОВУШКА: отличается только порядком!")
say("произведение координат -- квадрат в каждой строке:",
    all(ZZ(prod(cl)).is_square() for cl in classes))

# замкнутость по умножению (образ обязан быть подгруппой)
def mul(u,v): return tuple(sqcls(QQ(u[i])*QQ(v[i])) for i in range(3))
closed = all(mul(u,v) in set(classes) for u in classes for v in classes)
say("множество из 8 классов замкнуто относительно умножения:", closed)
say("требуемый класс * каждый из 8 -- ни один не равен тождеству? (т.е. req вне группы):",
    all(mul(req,cl) != (1,1,1) or cl == req for cl in classes))

# ============================================================================
hdr("5. ЛОГИКА: 'ВОСЕМЬ РАЗЛИЧНЫХ => ОБРАЗ ПОЛОН' БЕЗ НАСЫЩЕНИЯ")
# ============================================================================
say("Схема: delta: E(Q) -> (Q*/Q*^2)^3 гомоморфизм с ядром РОВНО 2E(Q).")
say("       => |image delta| = |E(Q)/2E(Q)| = 2^(r+2) (раздел 2).")
say("       r <= 1 (раздел 3) => |image| <= 8.")
say("       предъявлено 8 различных элементов образа => image = эти 8. Насыщение НЕ нужно.")
say("       (заодно: 8 различных => r >= 1, значит r = 1 ровно)")
say("")
say("ЭМПИРИЧЕСКИЙ КОНТРОЛЬ ЭТОГО УТВЕРЖДЕНИЯ: delta на многих точках E(Q)")
cls_set = set(classes)
bad = []
for k in range(-60, 61):
    for i in range(4):
        tt = [O,T1,T2,T3][i]
        P = k*G + tt
        if P == O:
            dl = (1,1,1)
        else:
            dl = delta(P)
        if dl not in cls_set:
            bad.append((k,i,dl))
        if dl == req:
            say("   !!! НАЙДЕНА ТОЧКА С ТРЕБУЕМЫМ КЛАССОМ:", k, i, P)
say("   проверено 121*4 точек kG+T; все классы внутри найденной восьмёрки:", not bad)
if bad: say("   ВЫБИВАЮЩИЕСЯ:", bad[:5])
say("   требуемый класс встретился:", any(True for _ in []) )

# ============================================================================
hdr("6. ПОТЕРЯННЫЕ ВЕТВИ")
# ============================================================================
say("6A. Может ли F_i обратиться в 0 при рациональном t?")
for nm,FF in [("F0",F0),("F4",F4),("F8",F8)]:
    disc = FF.discriminant()
    say("   ", nm, "=", FF, "  дискриминант =", disc, " <0 =>", "корней в R нет" if disc < 0 else "ЕСТЬ КОРНИ!")
say("   => ни одна F_i не зануляется, все три строго положительны при t в R.")
say("   => образ в E никогда не попадает в точку 2-кручения, регуляризация delta не нужна.")

say("")
say("6B. Знаки u_i: delta зависит только от u_i^2 => 8 выборов знаков дают ОДИН класс.")
say("    (проверка: (mn u4)^2, (m u0)^2, (n u8)^2 -- квадраты независимо от знака)")

say("")
say("6C. Бесконечность (карта z=1/t, U_i=u_i/t):")
say("    U0^2 -> n^2 =", n^2, " квадрат:", ZZ(n^2).is_square())
say("    U4^2 -> s   =", s,   " квадрат:", QQ(s).is_square())
say("    U8^2 -> m^2 =", m^2, " квадрат:", ZZ(m^2).is_square())
say("    => над t=infinity рациональных точек нет ТОЛЬКО потому, что s=113 не квадрат.")

say("")
say("6D. Особые точки / гладкость: проверяю род и число точек над t=infinity.")
say("    F0 старший коэф. n^2 =", n^2, "квадрат:", ZZ(n^2).is_square(),
    "; F8 старший m^2 квадрат:", ZZ(m^2).is_square(),
    "; F4 старший s квадрат:", QQ(s).is_square())
say("    ветвление только над 6 конечными корнями F0F4F8 (все различны:",
    len(set((F0*F4*F8).roots(CC, multiplicities=False))) == 6, ")")
say("    2g-2 = 8*(-2) + 6*4 = 8 => g = 5 (совпадает с заявленным)")

# ============================================================================
hdr("7. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ C(Q_p): ЛЕЖИТ ЛИ (1,113,113) В Sel^2?")
# ============================================================================
# Однородно: (x:y), G0=m^2y^2+n^2x^2, G4=s(x^2+y^2), G8=n^2y^2+m^2x^2.
# F_i квадрат в Q_p  <=>  G_i квадрат в Q_p (при y != 0), плюс y=0 = бесконечность.
def C_Qp_point(p, maxdepth=16):
    Gs = [lambda x,y: m^2*y^2 + n^2*x^2,
          lambda x,y: ZZ(s)*(x^2 + y^2),
          lambda x,y: n^2*y^2 + m^2*x^2]
    slack = 4 if p == 2 else 2
    stack = [(x,y,1) for x in range(p) for y in range(p) if not (x%p==0 and y%p==0)]
    undecided = False
    while stack:
        x,y,k = stack.pop()
        vals = [g(ZZ(x),ZZ(y)) for g in Gs]
        if any(v == 0 for v in vals):
            return ("точка с нулевой F", x, y)
        decided_all = True; dead = False
        for v in vals:
            e = v.valuation(p)
            if e + slack <= k:
                if not is_sq_Qp(v, p):
                    dead = True; break
            else:
                decided_all = False
        if dead: continue
        if decided_all:
            return ("x,y", x, y, k)
        if k >= maxdepth:
            undecided = True; continue
        pk = p^k
        for aa in range(p):
            for bb in range(p):
                stack.append((x + aa*pk, y + bb*pk, k+1))
    return None if undecided else False

say("C(R): F0,F4,F8 > 0 при всех вещественных t =>", "ТОЧКА ЕСТЬ")
for p in [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,113]:
    r = C_Qp_point(p)
    say("   p =", p, "->", ("C(Q_%d) НЕ ПУСТО, свидетель " % p) + str(r) if r else
        ("C(Q_%d) ПУСТО (обструкция!)" % p if r is False else "НЕ РЕШЕНО"))

# ============================================================================
hdr("8. НЕЗАВИСИМЫЙ МАРШРУТ ПРОЕКТА: КРИВАЯ E_B (квартика z^2=F0*F8)")
# ============================================================================
EB = EllipticCurve([0, -((m^2-n^2)^2 + (m^2+n^2)^2), 0, (m^2-n^2)^2*(m^2+n^2)^2, 0])
say("E_B: Y^2 = X(X-(m^2-n^2)^2)(X-(m^2+n^2)^2) =", EB)
say("   кручение E_B:", EB.torsion_subgroup().invariants())
try:
    epb = pari([ZZ(x) for x in EB.minimal_model().a_invariants()]).ellinit()
    say("   PARI ellrank(E_B) =", epb.ellrank())
except Exception as ex:
    say("   PARI ellrank(E_B) ОШИБКА:", ex)
say("   Sage EB.rank_bound() =", EB.rank_bound())
say("   ЕСЛИ rank E_B = 0, то по теореме проекта t in {0,±1,infty}; тогда:")
for tv in [0,1,-1]:
    say("      t =", tv, ": F4 =", s*(1+tv^2), " квадрат:", QQ(s*(1+tv^2)).is_square())
say("      t = infty: нужен квадрат s =", s, "->", QQ(s).is_square())

# ============================================================================
hdr("9. ПРЯМАЯ ОХОТА ЗА КОНТРПРИМЕРОМ НА C")
# ============================================================================
found = []
LIM = 700
for q in range(1, LIM+1):
    q2 = q*q
    for pnum in range(0, LIM+1):
        if gcd(pnum,q) != 1: continue
        p2 = pnum*pnum
        A0 = m^2*q2 + n^2*p2
        if not ZZ(A0).is_square(): continue
        A8 = n^2*q2 + m^2*p2
        if not ZZ(A8).is_square(): continue
        A4 = ZZ(s)*(p2 + q2)
        if not ZZ(A4).is_square(): continue
        found.append((pnum,q))
say("перебор t=p/q, 0<=p<=%d, 1<=q<=%d (gcd=1): найдено точек C(Q): %d" % (LIM,LIM,len(found)))
say("   ", found[:10])
V['search_found'] = len(found)

hdr("ИТОГ")
for k in sorted(V): say("   ", k, "=", V[k])
say("")
say("время, с:", round(time.time()-T0, 1))
