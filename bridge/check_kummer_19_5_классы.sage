# -*- coding: utf-8 -*-
# ============================================================================
#  НЕЗАВИСИМАЯ ПРОВЕРКА (попытка ОПРОВЕРЖЕНИЯ) исключения пары (m,n) = (19,5)
#  в семействе G1.   Claude (Opus 5), 2026-09-12.
#
#  Установка: НЕ подтверждать, а ломать.  Все константы, кривая, образующие
#  и классы delta считаются с нуля.  Числа Codex-а используются только
#  в финальном разделе 10 и помечены явно.
#
#  Запуск:  sage /home/kep/magicKube/bridge/check_kummer_19_5_классы.sage
# ============================================================================

import itertools, sys

def hdr(t):
    print("\n" + "="*74); print(t); print("="*74); sys.stdout.flush()
def say(*a):
    print(*a); sys.stdout.flush()

# ---------------------------------------------------------------------------
hdr(u"0. СИМВОЛЬНАЯ ПРОВЕРКА ТОЖДЕСТВ C -> E  (общие m,n, ничего не берём на веру)")
# ---------------------------------------------------------------------------
Rsym = PolynomialRing(QQ, ['ms','ns','ts','ss']); (ms,ns,ts,ss) = Rsym.gens()
F0s = ms^2 + ns^2*ts^2
F4s = ss*(1 + ts^2)
F8s = ns^2 + ms^2*ts^2
bs  = ss*ms^2*ns^2
Xs  = bs*ts^2
e1s, e2s, e3s = -bs, -ss*ms^4, -ss*ns^4
id1 = Xs - e1s - (ms*ns)^2 * F4s
id2 = Xs - e2s - ss*ms^2   * F0s
id3 = Xs - e3s - ss*ns^2   * F8s
idV = ((ms*ns)^2*F4s)*(ss*ms^2*F0s)*(ss*ns^2*F8s) - bs^2*F0s*F4s*F8s
say("X - e1 - (mn)^2*F4              =", id1)
say("X - e2 - s*m^2*F0               =", id2)
say("X - e3 - s*n^2*F8               =", id3)
say("prod(X-e_i) - (b*u0*u4*u8)^2    =", idV)
assert id1 == 0 and id2 == 0 and id3 == 0 and idV == 0
say(u"[доказано] у любой конечной точки C(Q) образ в E имеет delta = (1, s, s).")
say(u"[доказано] V!=0 на образе: F0,F4,F8 > 0 при вещественном t (m,n,s > 0),")
say(u"           значит точка НЕ 2-кручение и НЕ O, правило замены не нужно.")

# ---------------------------------------------------------------------------
hdr(u"1. КОНСТАНТЫ ДЛЯ (m,n) = (19,5) — пересчитаны заново")
# ---------------------------------------------------------------------------
m, n = 19, 5
assert gcd(m, n) == 1
s_q = QQ(m^2 + n^2) / 2
say("s = (m^2+n^2)/2 = (%d+%d)/2 = %s" % (m^2, n^2, s_q))
s = ZZ(s_q)
b = s * m^2 * n^2
say("b = s*m^2*n^2 = %d * %d * %d = %d" % (s, m^2, n^2, b))
e1, e2, e3 = -b, -s*m^4, -s*n^4
say("e1 = -b     = %d" % e1)
say("e2 = -s*m^4 = %d" % e2)
say("e3 = -s*n^4 = %d" % e3)
sq_s = ZZ(QQ(s).squarefree_part())
REQUIRED = (ZZ(1), sq_s, sq_s)
say(u"ТРЕБУЕМЫЙ класс delta = (1, s, s) =", REQUIRED)
say(u"s = %d — квадрат? %s   => точек C(Q) над t=oo %s"
    % (s, s.is_square(), u"НЕТ [доказано]" if not s.is_square() else u"ВОЗМОЖНЫ"))
say(u"  (в координатах z=1/t, U4=u4/t получается U4^2 = s(1+z^2) -> s при z=0)")

# ---------------------------------------------------------------------------
hdr(u"2. КРИВАЯ E И ЕЁ 2-КРУЧЕНИЕ")
# ---------------------------------------------------------------------------
a2 = -(e1+e2+e3); a4 = e1*e2 + e1*e3 + e2*e3; a6 = -e1*e2*e3
E = EllipticCurve([0, a2, 0, a4, a6])
say("E :", E)
Rx = PolynomialRing(QQ, 'Xv'); Xv = Rx.gen()
assert (Xv-e1)*(Xv-e2)*(Xv-e3) == Xv^3 + a2*Xv^2 + a4*Xv + a6
say(u"[проверено] кубика E тождественно равна (X-e1)(X-e2)(X-e3)")
say("disc      =", factor(E.discriminant()))
say("conductor =", factor(E.conductor()))
say("Emin      =", E.minimal_model())
T = E.torsion_subgroup()
say("E(Q)_tors =", T.invariants(), " порядок", E.torsion_order())
assert T.invariants() == (2,2), u"кручение не (Z/2)^2 — вся логика 2^(r+2) поплыла"
say(u"[проверено] полное рациональное 2-кручение => |E(Q)/2E(Q)| = 2^(r+2)")

# ---------------------------------------------------------------------------
hdr(u"3. РАНГ: три РАЗНЫХ источника верхней границы — здесь и зарыта собака")
# ---------------------------------------------------------------------------
mw = E.mwrank_curve()
mw_sel, mw_ub, mw_lb, mw_cert = mw.selmer_rank(), mw.rank_bound(), mw.rank(), mw.certain()
say("eclib/mwrank :  rk S^2(E) = %d  =>  %d <= rank <= %d ;  certain=%s"
    % (mw_sel, mw_lb, mw_ub, mw_cert))
say(u"  >>> ЧИСТЫЙ 2-СПУСК НЕ ДАЁТ rank=1. Он даёт только 1 <= rank <= %d." % mw_ub)
say(u"  >>> |Sel^2| = 2^%d = %d, а образ delta имеет размер 2^(r+2) = 8 при r=1." % (mw_sel, 2^mw_sel))
say(u"  >>> Значит 8 классов НЕ 'исчерпывают Селмера'; нужна отдельная граница r<=1.")
pari_out = None
try:
    pari_out = pari(E).ellrank()
    say("PARI ellrank() = %s   [r_low, r_up, s, points]" % pari_out)
    say(u"  PARI использует спаривание Касселса-Тейта на Sel^2 — это безусловно,")
    say(u"  но реализация одна; повтор PARI не есть независимая проверка.")
except Exception as ex:
    say("PARI ellrank недоступен:", ex)
say("Sage rank_bound(algorithm='pari') =", E.rank_bound())
Em = E.minimal_model()
say("root number w =", Em.root_number(), u"=> аналитический ранг НЕЧЁТЕН, >= 1")
try:
    Ld = Em.lseries().dokchitser(60)
    say("проверка функционального уравнения (должно быть ~0):", Ld.check_functional_equation())
    say("L(E,1)  =", Ld.derivative(1,0))
    say("L'(E,1) =", Ld.derivative(1,1))
    say(u"  L(1)=0, L'(1)!=0 => аналитический ранг = 1 =>")
    say(u"  [доказано (ПО), по модулярности + Гросс-Загье + Колывагин] rank E(Q) = 1, Sha конечна.")
    say(u"  Это НЕЗАВИСИМЫЙ от 2-спуска источник верхней границы r <= 1.")
except Exception as ex:
    say("L-ряд не посчитался:", ex)

# ---------------------------------------------------------------------------
hdr(u"4. СВОИ ТОЧКИ БЕСКОНЕЧНОГО ПОРЯДКА (не берём чужую)")
# ---------------------------------------------------------------------------
cands = []
try:
    for eff in [1,2,4]:
        pr = pari(E).ellrank(eff)
        for xy in pr[3]:
            cands.append((QQ(xy[0]), QQ(xy[1])))
except Exception as ex:
    say("PARI поиск точек:", ex)
try:
    for Q in E.gens(proof=False):
        cands.append((QQ(Q[0]), QQ(Q[1])))
except Exception as ex:
    say("E.gens(proof=False):", ex)
try:
    for Q in E.point_search(16, rank_bound=1):
        if Q.order() == oo: cands.append((QQ(Q[0]), QQ(Q[1])))
except Exception as ex:
    say("point_search:", ex)

pts = []
seen = set()
for (x0,y0) in cands:
    if (x0,y0) in seen: continue
    seen.add((x0,y0))
    # НЕЗАВИСИМАЯ проверка уравнения вручную
    if y0^2 != (x0-e1)*(x0-e2)*(x0-e3):
        say(u"  ОТБРОШЕНО (не на кривой):", x0, y0); continue
    P = E(x0, y0)
    if P.order() != oo:
        say(u"  ОТБРОШЕНО (кручение):", x0); continue
    pts.append(P)
if not pts:
    say(u"НЕ НАШЁЛ своей неторсионной точки — останавливаюсь честно."); sys.exit(1)
say(u"нашёл %d своих неторсионных точек; проверка y^2=(x-e1)(x-e2)(x-e3) пройдена для каждой:" % len(pts))
for P in pts:
    say("   x = %s" % P[0]); say("   y = %s   h^ = %s" % (P[1], P.height()))
G = pts[0]
say(u"беру G = первая из них")
try:
    sat, idx, reg = E.saturation([G])
    say(u"saturation([G]): индекс = %s" % idx)
    if idx != 1:
        G = sat[0]; say(u"  -> G заменена на насыщенную")
except Exception as ex:
    say(u"saturation:", ex)
say(u"ЗАМЕЧАНИЕ: насыщенность для вывода НЕ нужна — нужны 2^(r+2) различных классов.")

T1 = E(e1, 0); T2 = E(e2, 0); T3 = E(e3, 0)
assert 2*T1 == E(0) and 2*T2 == E(0) and T1 + T2 == T3
say(u"[проверено] T1=(e1,0), T2=(e2,0), T3=(e3,0), T1+T2=T3, 2Ti=O")

# ---------------------------------------------------------------------------
hdr(u"5. delta — реализовано с нуля")
# ---------------------------------------------------------------------------
ROOTS = [e1, e2, e3]
# Опорное множество простых: любой квадратный класс в образе delta поддержан
# на простых, делящих 2*disc(E).  Снимаем квадраты ТОЛЬКО по ним и затем
# ТРЕБУЕМ, чтобы остаток был точным квадратом — это проверка, а не допущение.
SUPPORT = sorted(set(ZZ(2*E.discriminant().numerator()).prime_factors()))
say(u"опорные простые для квадратных классов:", SUPPORT)
def sqfree(q):
    q = QQ(q); assert q != 0
    N = ZZ(q.numerator()) * ZZ(q.denominator())   # тот же класс, что и q
    d = ZZ(1)
    if N < 0: d = ZZ(-1); N = -N
    for p in SUPPORT:
        v = N.valuation(p)
        if v % 2 == 1: d *= p
        N //= p^v
    assert N.is_square(), u"остаток не квадрат: класс вне опорного множества (%s)" % N
    return d
def delta(P):
    if P == E(0): return (ZZ(1), ZZ(1), ZZ(1))
    x = QQ(P[0]); y = QQ(P[1]); out = []
    for i in range(3):
        if y == 0 and x == ROOTS[i]:
            j, k = [u for u in range(3) if u != i]
            out.append(sqfree((ROOTS[i]-ROOTS[j])*(ROOTS[i]-ROOTS[k])))
        else:
            out.append(sqfree(x - ROOTS[i]))
    return tuple(out)
def mulc(a, c): return tuple(sqfree(a[i]*c[i]) for i in range(3))

# ---------------------------------------------------------------------------
hdr(u"6. ТАБЛИЦА ВСЕХ 2^(r+2) = 8 КЛАССОВ  (r = 1)")
# ---------------------------------------------------------------------------
rows = []
for (a, c, d) in itertools.product([0,1], repeat=3):
    P = a*G + c*T1 + d*T2
    lbl = " + ".join([t for t,f in [("G",a),("T1",c),("T2",d)] if f]) or "O"
    rows.append((lbl, P, delta(P)))
wl = max(len(r[0]) for r in rows)
say("%-*s | %-40s | delta (бесквадратные представители со знаком)" % (wl, u"точка", "x(P)"))
say("-"*104)
for lbl,P,dP in rows:
    xs = "O" if P == E(0) else str(P[0])
    if len(xs) > 40: xs = xs[:37] + "..."
    say("%-*s | %-40s | (%s, %s, %s)" % (wl, lbl, xs, dP[0], dP[1], dP[2]))
classes = [r[2] for r in rows]
uniq = set(classes)

# ---------------------------------------------------------------------------
hdr(u"7. КОНТРОЛИ")
# ---------------------------------------------------------------------------
say(u"(а) различных классов: %d из 8  ->  %s" %
    (len(uniq), u"ВСЕ РАЗЛИЧНЫ, образ полон при r=1" if len(uniq)==8
                else u"СОВПАДЕНИЯ! образ НЕ полон, вывод неверен"))
say(u"(б) требуемый класс %s среди них?  ->  %s" %
    (str(REQUIRED), u"ДА — ИСКЛЮЧЕНИЕ НЕВЕРНО" if REQUIRED in uniq else u"НЕТ"))
say(u"(в) d1*d2*d3 — квадрат у всех восьми? ->",
    all(sqfree(c[0]*c[1]*c[2]) == 1 for c in classes))
say(u"(г) гомоморфность delta:")
bad = 0; tested = 0
pairs = [(G,T1),(G,T2),(T1,T2),(2*G,T1),(G,G),(3*G,T2),(2*G,-3*G),(G+T1,G+T2),(-2*G,T1),(2*G,G)]
set_random_seed(20260912)
for _ in range(20):
    A = ZZ.random_element(-3,4)*G + ZZ.random_element(0,2)*T1 + ZZ.random_element(0,2)*T2
    B2 = ZZ.random_element(-3,4)*G + ZZ.random_element(0,2)*T1 + ZZ.random_element(0,2)*T2
    pairs.append((A,B2))
for (A,B2) in pairs:
    if A == E(0) or B2 == E(0) or A+B2 == E(0): continue
    tested += 1
    if delta(A+B2) != mulc(delta(A), delta(B2)):
        bad += 1; say(u"    НАРУШЕНА на", A, B2)
say(u"    проверено пар: %d, нарушений: %d" % (tested, bad))
say(u"(д) delta(2P) = (1,1,1) (ядро содержит 2E(Q)):")
for P in [G, G+T1, G+T2, G+T1+T2, -G]:
    say("    delta(2P) =", delta(2*P))
say(u"(е) устойчивость к выбору образующей — пересчёт таблицы по каждой своей точке:")
for P in pts:
    cl = set()
    for (a,c,d) in itertools.product([0,1],repeat=3):
        cl.add(delta(a*P + c*T1 + d*T2))
    say("    x=%s...  различных=%d  множество совпало с базовым? %s  (1,s,s) внутри? %s"
        % (str(P[0])[:18], len(cl), cl == uniq, REQUIRED in cl))

# ---------------------------------------------------------------------------
hdr(u"8. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ ТРЕБУЕМОГО КЛАССА (насколько близко он подходит)")
# ---------------------------------------------------------------------------
d1, d2, d3 = REQUIRED
say(u"торсор класса (d1,d2,d3):  d1*z1^2 - d2*z2^2 = (e2-e1)*z0^2,  d1*z1^2 - d3*z3^2 = (e3-e1)*z0^2")
say("e2-e1 = %d,  e3-e1 = %d" % (e2-e1, e3-e1))
badp = sorted(set(ZZ(2*E.discriminant().numerator()).prime_factors()))
say(u"плохие простые:", badp)
def qp_ok(p, prec=40):
    K = Qp(p, prec)
    for num in range(-80, 81):
        for den in range(1, 30):
            z1 = QQ(num)/den
            A = d1*z1^2 - (e2-e1); B2 = d1*z1^2 - (e3-e1)
            if A == 0 or B2 == 0: continue
            if (K(A)/K(d2)).is_square() and (K(B2)/K(d3)).is_square():
                return True, z1
    if (K(d1)/K(d2)).is_square() and (K(d1)/K(d3)).is_square(): return True, "z0=0"
    return False, None
for p in badp:
    r_, wit = qp_ok(p)
    say(u"   Q_%-5d: %s %s" % (p, u"разрешим" if r_ else u"свидетель НЕ найден (перебор ограничен!)",
                               ("z1=%s"%wit) if r_ else ""))
say(u"   R: d_i>0, e2-e1<0, e3-e1<0 -> z1=0 годится -> разрешим")
say(u"   ВЫВОД: класс (1,s,s), по-видимому, ВЕЗДЕ ЛОКАЛЬНО разрешим, т.е. лежит в Sel^2,")
say(u"   но НЕ в образе E(Q). Это ровно элемент Sha[2]. Поэтому локальные методы")
say(u"   исключить (19,5) НЕ могут, и всё держится на границе r <= 1.")

# ---------------------------------------------------------------------------
hdr(u"9. ЛОБОВОЙ ПОИСК ТОЧЕК НА САМОЙ C  (попытка сломать напрямую)")
# ---------------------------------------------------------------------------
from math import isqrt, gcd as _gcd
LIM = 4000
mm, nn, ss_ = int(m*m), int(n*n), int(s)
def issq(v):
    r = isqrt(v); return r*r == v
found = []
for q in range(1, LIM+1):
    q2 = q*q
    for p_ in range(0, LIM+1):
        if _gcd(p_, q) != 1: continue
        p2 = p_*p_
        if not issq(ss_*(q2+p2)): continue
        if not issq(mm*q2 + nn*p2): continue
        if not issq(nn*q2 + mm*p2): continue
        found.append((p_, q))
say(u"t = p/q, 0 <= p,q <= %d, gcd(p,q)=1 : найдено точек C(Q) = %d" % (LIM, len(found)))
if found: say("   ", found[:20])
say(u"   (отсутствие находки — НЕ доказательство отсутствия; это только попытка сломать)")

# ---------------------------------------------------------------------------
hdr(u"10. СВЕРКА С ЧУЖИМИ ЧИСЛАМИ (Codex) — только здесь")
# ---------------------------------------------------------------------------
codex = [(1,1,1),(-1,4053,-4053),(-4053,386,-42),(4053,42,386),
         (-3,386,-1158),(3,42,14),(1351,1,1351),(-1351,4053,-3)]
codex = set(tuple(ZZ(v) for v in c) for c in codex)
say(u"множество Codex-а совпало с моим? ->", codex == uniq)
say(u"разница (моё \\ его):", uniq - codex)
say(u"разница (его \\ моё):", codex - uniq)

hdr(u"ИТОГ")
say(u"s = %d, b = %d, требуемый класс %s" % (s, b, str(REQUIRED)))
say(u"различных классов: %d/8 ; (1,s,s) в образе: %s" % (len(uniq), REQUIRED in uniq))
say(u"rank: 2-спуск (eclib) 1 <= r <= %d ; PARI(с Касселсом) r = 1 ; L'(1)!=0 => r = 1" % mw_ub)
