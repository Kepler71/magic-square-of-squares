# -*- coding: utf-8 -*-
# check_kummer_15_1_ранг.sage
# АТАКА на заявление Codex: C_{15,1}(Q) = пусто.
# Углы атаки:
#   (A) верхняя граница rank E(Q) двумя независимыми средствами (PARI ellrank, eclib mwrank)
#   (B) структура кручения и формула |E(Q)/2E(Q)| = 2^(r+2)
#   (C) законность вывода "2^(r+2) различных классов => образ полон" без насыщения
#   (D) потерянные ветви: порядок корней e_i, знаки u_i, нули u_i, бесконечность, особые точки
#   (E) прямой поиск точек на C (попытка предъявить контрпример)
# Всё строится ЗАНОВО из (m,n)=(15,1). Файлы Codex не импортируются.

import sys, time, json

T0 = time.time()
def say(*a):
    print(*a); sys.stdout.flush()
def hdr(s):
    say(""); say("="*78); say(s); say("="*78)

VERDICT = {}

# ===========================================================================
hdr("0. ПОСТРОЕНИЕ СЕМЕЙСТВА ЗАНОВО ИЗ (m,n)=(15,1)")
# ===========================================================================
m = QQ(15); n = QQ(1)
say("gcd(m,n) =", gcd(ZZ(m), ZZ(n)))
s = (m^2 + n^2)/2
b = s*m^2*n^2
e1 = -b
e2 = -s*m^4
e3 = -s*n^4
say("s =", s, "  знаменатель:", s.denominator())
say("b = s*m^2*n^2 =", b)
say("e1,e2,e3 =", e1, e2, e3)
say("Codex заявил s=113, b=25425, e=(-25425,-5720625,-113):",
    (s == 113 and b == 25425 and (e1,e2,e3) == (-25425,-5720625,-113)))
VERDICT['model_matches_codex'] = bool(s == 113 and b == 25425 and (e1,e2,e3)==(-25425,-5720625,-113))

say("s квадрат в Q?", QQ(s).is_square())
say("s свободно от квадратов?", ZZ(s).is_squarefree(), " факторизация s:", factor(ZZ(s)))

R.<X> = PolynomialRing(QQ)
cub = (X-e1)*(X-e2)*(X-e3)
cc = cub.coefficients(sparse=False)
E = EllipticCurve([0, cc[2], 0, cc[1], cc[0]])
say("E:", E)
say("disc(E) =", factor(E.discriminant()))
say("conductor =", E.conductor(), " =", factor(E.conductor()))
Emin = E.minimal_model()
say("минимальная модель:", Emin.a_invariants())
say("Codex min_ainvs [0,-1,0,-42422006041,3362853259097305]:",
    list(Emin.a_invariants()) == [0,-1,0,-42422006041,3362853259097305])

# ===========================================================================
hdr("1. ТОЖДЕСТВА X-e_i СИМВОЛИЧЕСКИ (проверка ПОРЯДКА корней!)")
# ===========================================================================
# Это самое опасное место: если индексация e_i перепутана, требуемый класс
# окажется другим и вывод развалится. В образе Codex есть (113,113,1) --
# перестановка требуемого (1,113,113). Поэтому проверяем ОЧЕНЬ аккуратно.
S.<t> = PolynomialRing(QQ)
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
say("F0 =", F0, "   F4 =", F4, "   F8 =", F8)
Xt = b*t^2
id1 = Xt - e1 - (m*n)^2 * F4
id2 = Xt - e2 - s*m^2 * F0
id3 = Xt - e3 - s*n^2 * F8
say("X-e1 - (mn)^2*F4 == 0 :", id1 == 0, "   [X-e1 =", Xt-e1, "= 225*F4 ]")
say("X-e2 - s*m^2*F0 == 0 :", id2 == 0, "   [X-e2 =", Xt-e2, "= 25425*F0]")
say("X-e3 - s*n^2*F8 == 0 :", id3 == 0, "   [X-e3 =", Xt-e3, "= 113*F8 ]")
VERDICT['identities_symbolic'] = bool(id1==0 and id2==0 and id3==0)

# Проверка V^2 = произведение
u0v, u4v, u8v = var_names = None, None, None
Pr.<tt,U0,U4,U8> = PolynomialRing(QQ)
Xp = b*tt^2
lhs = (b*U0*U4*U8)^2
rhs = (Xp-e1)*(Xp-e2)*(Xp-e3)
# подставляем U_i^2 -> F_i(tt)
rel = lhs - rhs
red = rel.reduce([U0^2-(m^2+n^2*tt^2), U4^2-s*(1+tt^2), U8^2-(n^2+m^2*tt^2)])
say("V=b*u0*u4*u8 удовлетворяет V^2=(X-e1)(X-e2)(X-e3) по модулю уравнений C:", red == 0)
VERDICT['V_identity'] = bool(red == 0)

# требуемый квадратный класс
def sqfree(x):
    x = QQ(x)
    if x == 0:
        return ZZ(0)
    num = x.numerator(); den = x.denominator()
    z = num*den                      # тот же квадратный класс
    sgn = 1 if z > 0 else -1
    z = abs(z)
    r = ZZ(1)
    for p,ex in factor(z):
        if ex % 2 == 1:
            r *= p
    return ZZ(sgn*r)

req = (sqfree((m*n)^2), sqfree(s*m^2), sqfree(s*n^2))
say("ТРЕБУЕМЫЙ класс delta = (класс (mn)^2, класс s*m^2, класс s*n^2) =", req)
say("Codex заявил (1,113,113):", req == (1,113,113))
VERDICT['required_class'] = [int(x) for x in req]

# ===========================================================================
hdr("2. КРУЧЕНИЕ: явно, и влияние на формулу |E(Q)/2E(Q)|")
# ===========================================================================
Tor = E.torsion_subgroup()
say("E(Q)_tors =", Tor, "  порядок =", Tor.order())
say("структура:", Tor.invariants())
T2pts = [P for P in E.torsion_points() if 2*P == E(0) and P != E(0)]
say("точки порядка 2:", T2pts)
say("|E(Q)[2]| =", len(T2pts)+1)
VERDICT['torsion_order'] = int(Tor.order())
VERDICT['torsion_invariants'] = [int(x) for x in Tor.invariants()]
VERDICT['E2_order'] = int(len(T2pts)+1)

say("")
say("--- контроль формулы |E(Q)/2E(Q)| = 2^r * |T[2]| для разных T ---")
# E(Q) = Z^r x T. |E/2E| = 2^r * |T/2T|, и |T/2T| = |T[2]| для конечной абелевой группы.
for inv in [(2,2),(2,4),(2,8),(2,2,3),(4,4),(2,)]:
    G0 = AbelianGroup(len(inv), inv)
    # считаем |T/2T| и |T[2]| перебором
    els = list(G0)
    twoT = set()
    for g in els:
        twoT.add(g^2)
    quot = len(els)//len(twoT)
    t2 = len([g for g in els if g^2 == G0.one()])
    say("  T =", inv, " |T| =", len(els), " |T/2T| =", quot, " |T[2]| =", t2,
        " совпало:", quot == t2)
say("ВЫВОД: |T/2T| = |T[2]| всегда; при полном рациональном 2-кручении |T[2]| = 4")
say("       ЛЮБОЕ увеличение кручения (Z/2xZ/4, Z/2xZ/8, добавление нечётной части)")
say("       НЕ меняет 2^(r+2): |E/2E| = 2^r*|T[2]| = 2^r*4 = 2^(r+2).")

# ===========================================================================
hdr("3A. ВЕРХНЯЯ ГРАНИЦА РАНГА: PARI ellrank (независимо)")
# ===========================================================================
ainv = [ZZ(x) for x in Emin.a_invariants()]
say("ainvs (минимальная модель) =", ainv)
try:
    ep = pari(ainv).ellinit()
    rk = ep.ellrank()
    say("PARI ellrank ->", rk)
    lo = ZZ(rk[0]); hi = ZZ(rk[1])
    say("  нижняя граница ранга:", lo, "  верхняя граница ранга:", hi)
    say("  дополнительный параметр (ранг спаривания Кассельса / s):", rk[2])
    say("  найденные точки:", rk[3])
    VERDICT['pari_rank_lo'] = int(lo); VERDICT['pari_rank_hi'] = int(hi)
except Exception as ex:
    say("PARI ellrank НЕ ОТРАБОТАЛ:", ex)
    VERDICT['pari_rank_lo'] = None; VERDICT['pari_rank_hi'] = None

# то же на НЕминимальной исходной модели -- контроль
try:
    ep2 = pari([ZZ(x) for x in E.a_invariants()]).ellinit()
    rk2 = ep2.ellrank()
    say("PARI ellrank на исходной (неминимальной) модели ->", rk2[0], rk2[1])
except Exception as ex:
    say("  (исходная модель: ellrank не отработал:", ex, ")")

# ===========================================================================
hdr("3B. ВЕРХНЯЯ ГРАНИЦА РАНГА: eclib mwrank (независимо), certain?")
# ===========================================================================
try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    mw = mwrank_EllipticCurve([int(x) for x in ainv])
    mw.two_descent(verbose=False, selmer_only=False, second_descent=True)
    say("mwrank rank        =", mw.rank())
    say("mwrank rank_bound  =", mw.rank_bound())
    say("mwrank selmer_rank =", mw.selmer_rank())
    say("mwrank CERTAIN     =", mw.certain())
    VERDICT['mwrank_rank'] = int(mw.rank())
    VERDICT['mwrank_bound'] = int(mw.rank_bound())
    VERDICT['mwrank_selmer'] = int(mw.selmer_rank())
    VERDICT['mwrank_certain'] = bool(mw.certain())
except Exception as ex:
    say("eclib mwrank НЕ ОТРАБОТАЛ:", ex)
    VERDICT['mwrank_certain'] = None

# ===========================================================================
hdr("3C. ТРЕТЬЕ СРЕДСТВО: 2-Selmer через Sage, analytic rank, Sha")
# ===========================================================================
try:
    sr = E.selmer_rank()
    say("E.selmer_rank() =", sr, " (это dim_F2 Sel^2, включая вклад кручения)")
    VERDICT['sage_selmer_rank'] = int(sr)
except Exception as ex:
    say("selmer_rank не отработал:", ex)
try:
    rb = E.rank_bound()
    say("E.rank_bound() =", rb)
    VERDICT['sage_rank_bound'] = int(rb)
except Exception as ex:
    say("rank_bound не отработал:", ex)
try:
    r_sage = E.rank(proof=True)
    say("E.rank(proof=True) =", r_sage)
    VERDICT['sage_rank'] = int(r_sage)
except Exception as ex:
    say("E.rank(proof=True) не отработал:", ex)
try:
    ar = E.analytic_rank()
    say("аналитический ранг (только как СОГЛАСОВАНИЕ, не доказательство) =", ar)
    VERDICT['analytic_rank'] = int(ar)
except Exception as ex:
    say("analytic_rank не отработал:", ex)
try:
    sha = E.sha().an_numerical()
    say("численное #Sha (BSD, условно) =", sha)
except Exception as ex:
    say("sha не отработал:", ex)

# ===========================================================================
hdr("4. ГЕНЕРАТОРЫ И ТОЧКА G ИЗ ЗАЯВЛЕНИЯ CODEX")
# ===========================================================================
G = E(-75825, 146764800)
say("G =", G, " лежит на E: да (иначе Sage бы бросил исключение)")
say("порядок G:", G.order())
say("высота Нерона-Тейта G:", G.height())
try:
    gens = E.gens(proof=True)
    say("E.gens(proof=True) =", gens)
    VERDICT['gens'] = [str(g) for g in gens]
    if len(gens) == 1:
        g0 = gens[0]
        # выразим G через g0 по модулю кручения: сравним высоты
        say("h(G)/h(gen) =", (G.height()/g0.height()).n())
except Exception as ex:
    say("E.gens не отработал:", ex)

# ===========================================================================
hdr("5. ОТОБРАЖЕНИЕ 2-СПУСКА delta -- СВОЯ РЕАЛИЗАЦИЯ")
# ===========================================================================
roots = [e1, e2, e3]
def delta(P):
    """полный 2-спуск: E(Q) -> (Q*/Q*^2)^3, ядро 2E(Q)."""
    if P == E(0):
        return (ZZ(1), ZZ(1), ZZ(1))
    x = P[0]; y = P[1]
    out = []
    for i in range(3):
        if x == roots[i]:
            j,k = [a for a in range(3) if a != i]
            out.append(sqfree((roots[i]-roots[j])*(roots[i]-roots[k])))
        else:
            out.append(sqfree(x - roots[i]))
    return tuple(out)

T1 = E(e1, 0); T2 = E(e2, 0); T3 = E(e3, 0)
say("T1,T2,T3 =", T1, T2, T3, "  T1+T2 == T3:", T1+T2 == T3)

basis8 = []
labels = []
for a in [0,1]:
    for c in [0,1]:
        for d in [0,1]:
            P = a*G + c*T1 + d*T2
            basis8.append(P); labels.append((a,c,d))

say("")
say("%-10s %-34s %s" % ("(a,c,d)", "точка", "delta"))
classes = []
for lab, P in zip(labels, basis8):
    dl = delta(P)
    classes.append(dl)
    say("%-10s %-34s %s" % (str(lab), str(P)[:34], str(dl)))

say("")
say("все 8 классов различны:", len(set(classes)) == 8)
VERDICT['eight_distinct'] = bool(len(set(classes)) == 8)
say("ТРЕБУЕМЫЙ класс", req, "среди них:", tuple(req) in set(classes))
VERDICT['required_in_image'] = bool(tuple(req) in set(classes))

# сверка с таблицей Codex
codex_tbl = {
 (0,0,0):(1,1,1), (0,0,1):(-1582,226,-7), (0,1,0):(-1,1582,-1582),
 (0,1,1):(1582,7,226), (1,0,0):(-14,2,-7), (1,0,1):(113,113,1),
 (1,1,0):(14,791,226), (1,1,1):(-113,14,-1582)}
ok = all(tuple(map(ZZ, codex_tbl[l])) == c for l,c in zip(labels, classes))
say("моя таблица совпала с таблицей Codex построчно:", ok)
VERDICT['table_matches_codex'] = bool(ok)
if not ok:
    for l,c in zip(labels, classes):
        if tuple(map(ZZ, codex_tbl[l])) != c:
            say("  РАСХОЖДЕНИЕ", l, "мой:", c, " Codex:", codex_tbl[l])

# каждый класс должен иметь произведение = квадрат
say("")
for l,c in zip(labels, classes):
    pr = sqfree(c[0]*c[1]*c[2])
    if pr != 1:
        say("  ОШИБКА: произведение координат не квадрат для", l, c, "->", pr)
say("произведение трёх координат = квадрат во всех строках:",
    all(sqfree(c[0]*c[1]*c[2]) == 1 for c in classes))

# ===========================================================================
hdr("6. ОБРАЗ -- ГРУППА? ГОМОМОРФНОСТЬ? ИНЪЕКТИВНОСТЬ delta?")
# ===========================================================================
def mul(c1,c2):
    return tuple(sqfree(a*bq) for a,bq in zip(c1,c2))
Sset = set(classes)
closed = all(mul(c1,c2) in Sset for c1 in classes for c2 in classes)
say("множество из 8 классов замкнуто относительно умножения (подгруппа):", closed)
VERDICT['image_is_group'] = bool(closed)

# гомоморфность delta на большом наборе точек
say("")
say("проверка delta(P+Q)=delta(P)delta(Q) на случайных точках E(Q):")
pool = []
for a in range(-6, 7):
    for c in [0,1]:
        for d in [0,1]:
            pool.append(a*G + c*T1 + d*T2)
bad = 0; tested = 0
for P in pool:
    for Q in pool:
        tested += 1
        if delta(P+Q) != mul(delta(P), delta(Q)):
            bad += 1
            if bad <= 5:
                say("  НАРУШЕНИЕ:", P, Q)
say("  протестировано пар:", tested, "  нарушений:", bad)
VERDICT['hom_violations'] = int(bad)

# все точки из пула попадают в 8 классов (согласованность с полнотой образа)
outside = [P for P in pool if delta(P) not in Sset]
say("точек пула с классом ВНЕ найденных восьми:", len(outside))
VERDICT['pool_outside'] = int(len(outside))

# независимый поиск точек на E большего размера: не вылезет ли 9-й класс?
say("")
say("поиск дополнительных точек E(Q) (height bound) -- не даст ли 9-й класс:")
try:
    extra = E.point_search(12, verbose=False)
    say("  найдено точек:", len(extra))
    newcl = set()
    for P in extra:
        for c in [0,1]:
            for d in [0,1]:
                for a in range(-3,4):
                    Q = a*P + c*T1 + d*T2
                    newcl.add(delta(Q))
    say("  всего различных классов от них:", len(newcl))
    say("  все внутри найденных восьми:", newcl.issubset(Sset))
    say("  требуемый класс появился:", tuple(req) in newcl)
    VERDICT['point_search_new_classes'] = bool(not newcl.issubset(Sset))
    VERDICT['point_search_req'] = bool(tuple(req) in newcl)
except Exception as ex:
    say("  point_search не отработал:", ex)

# ЯДРО delta: проверим, что delta(2P)=(1,1,1) на пуле (необходимое условие ker>=2E)
say("")
k2 = all(delta(2*P) == (1,1,1) for P in pool)
say("delta(2P)=(1,1,1) для всех P пула (ker delta содержит 2E(Q)):", k2)
VERDICT['ker_contains_2E'] = bool(k2)

# ===========================================================================
hdr("7. ЛОГИКА ЗАКРЫТИЯ -- ПРОВЕРКА САМОГО РАССУЖДЕНИЯ")
# ===========================================================================
r_hi = VERDICT.get('pari_rank_hi', None)
say("Аргумент Codex:")
say("  (i)   delta инъективна на E(Q)/2E(Q)     [стандартная теорема, полное 2-кручение]")
say("  (ii)  |E(Q)/2E(Q)| = 2^rank * |E(Q)[2]| = 2^(rank+2)")
say("  (iii) rank <= 1  =>  |E(Q)/2E(Q)| <= 8  =>  |Im delta| <= 8")
say("  (iv)  предъявлены 8 РАЗЛИЧНЫХ значений delta на рациональных точках => |Im delta| >= 8")
say("  (v)   значит Im delta = ровно эти 8; требуемого среди них нет => C(Q) конечных точек нет")
say("")
say("КРИТИЧЕСКИЙ ВОПРОС: нужна ли насыщенность G?")
say("  НЕТ. Насыщенность нужна была бы, если бы мы выводили полноту из того, что")
say("  G порождает E(Q) mod tors. Здесь используется только СЧЁТ: 8 <= |Im| и |Im| <= 8.")
say("  Если бы G = k*G0 с k>1, то 8 различных классов всё равно лежат в Im delta,")
say("  и оценка |Im delta| <= 2^(r+2) = 8 всё равно принудительно даёт равенство.")
say("  => аргумент КОРРЕКТЕН при условии верхней границы ранга. Насыщение не требуется.")
say("")
say("ГДЕ АРГУМЕНТ МОГ БЫ СЛОМАТЬСЯ:")
say("  - если rank >= 2 (тогда |E/2E| >= 16 и 8 классов не исчерпывают образ);")
say("  - если E(Q)[2] был бы меньше (Z/2)^2 -- но тогда 2^(r+2) неверно;")
say("  - если бы какая-то из восьми 'точек' не была рациональной точкой E;")
say("  - если бы порядок корней e_i был перепутан (требуемый класс был бы другим).")
say("Все четыре пункта проверены выше по отдельности.")

# ===========================================================================
hdr("8. ПОТЕРЯННЫЕ ВЕТВИ: нули u_i, знаки, особые точки, бесконечность")
# ===========================================================================
say("(a) нули F_i при рациональном (даже вещественном) t:")
for nm,F in [("F0",F0),("F4",F4),("F8",F8)]:
    rts = F.roots(QQ)
    say("   ", nm, "=", F, " рациональных корней:", rts,
        " дискриминант<0 / положительно определён:", all(cf > 0 for cf in F.coefficients()))
say("    => ни один u_i не обращается в нуль при рациональном t; особых точек аффинной")
say("       модели над Q нет (якобиан падает в ранге только при u_i=0 и F_i'=0).")
VERDICT['no_zero_ui'] = True

say("")
say("(b) знаки корней u_i: замена u_i -> -u_i меняет только знак V=b*u0*u4*u8,")
say("    т.е. P -> -P. delta(-P)=delta(P), так как x-координата та же.")
chk = all(delta(-P) == delta(P) for P in pool)
say("    проверено на пуле: delta(-P)=delta(P) всюду:", chk)
VERDICT['sign_invariance'] = bool(chk)

say("")
say("(c) бесконечность t=oo: z=1/t, U_i=u_i/t.")
say("    U0^2 = m^2 z^2 + n^2 -> n^2 =", n^2, " (квадрат, ок)")
say("    U4^2 = s(z^2+1)     -> s   =", s, " квадрат в Q?", QQ(s).is_square())
say("    U8^2 = n^2 z^2 + m^2 -> m^2 =", m^2, " (квадрат, ок)")
say("    Расширение Q(z)[U4]/(U4^2-s(z^2+1)) в точке z=0: s(z^2+1) -- единица со значением")
say("    s=113, неквадрат => место ИНЕРТНО, единственное место над z=0 имеет поле вычетов")
say("    Q(sqrt(113)), степень 2 => рациональных точек над t=oo НЕТ.")
Kq = QQ
say("    113 квадрат в Q:", QQ(113).is_square())
VERDICT['infinity_excluded'] = bool(not QQ(s).is_square())

say("")
say("(d) может ли точка C(Q) отобразиться в 2-кручение E?")
say("    X = b*t^2 >= 0 при вещественном t; x-координаты 2-кручения:", e1, e2, e3, "-- все < 0.")
say("    Кроме того X-e_i = (квадрат)*u_i^2 != 0. Так что образ не 2-кручение.")
say("    Но даже если бы был: delta(T_i) входит в те же 8 классов, требуемого среди них нет.")
say("    точки кручения с X>=0 в E(Q):", [P for P in E.torsion_points() if P != E(0) and P[0] >= 0])

# ===========================================================================
hdr("9. ПРЯМОЙ ПОИСК ТОЧЕК НА C -- ПОПЫТКА ПРЕДЪЯВИТЬ КОНТРПРИМЕР")
# ===========================================================================
# t = p/q, gcd(p,q)=1, q>0, p>=0 (F_i чётны по t).
# однородно: q^2*F0 = 225 q^2 + p^2, q^2*F4 = 113(q^2+p^2), q^2*F8 = q^2 + 225 p^2
# все три должны быть квадратами.
import math
N = 4000
found = []
cnt = 0
for q in range(1, N+1):
    q2 = q*q
    for p in range(0, N+1):
        pp = p*p
        ssum = pp + q2
        if ssum % 113:
            continue
        cnt += 1
        v = 113*ssum
        r = math.isqrt(v)
        if r*r != v:
            continue
        A = 225*q2 + pp
        r = math.isqrt(A)
        if r*r != A:
            continue
        B = q2 + 225*pp
        r = math.isqrt(B)
        if r*r != B:
            continue
        if math.gcd(p,q) == 1:
            found.append((p,q))
say("перебор t=p/q, 0<=p<=%d, 1<=q<=%d; прошло фильтр 113|p^2+q^2: %d пар" % (N,N,cnt))
say("НАЙДЕННЫЕ точки C(Q) (кроме тривиальных):", found)
VERDICT['direct_search_found'] = [list(x) for x in found]
VERDICT['direct_search_N'] = N

# ===========================================================================
hdr("10. КОНТРОЛЬ: тот же тест на паре, где точка ЗАВЕДОМО есть")
# ===========================================================================
# (15,8) -- Codex говорит, что там требуемый класс ПРИСУТСТВУЕТ и известна
# вырожденная точка C. Если наш delta/логика ломаются, контроль это покажет.
def build(mm, nn):
    mm = QQ(mm); nn = QQ(nn)
    ss = (mm^2+nn^2)/2; bb = ss*mm^2*nn^2
    ee = [-bb, -ss*mm^4, -ss*nn^4]
    ccub = (X-ee[0])*(X-ee[1])*(X-ee[2])
    co = ccub.coefficients(sparse=False)
    EE = EllipticCurve([0, co[2], 0, co[1], co[0]])
    return ss, bb, ee, EE

for (mm,nn) in [(15,8),(15,1)]:
    ss,bb,ee,EE = build(mm,nn)
    rq = (sqfree((QQ(mm)*QQ(nn))^2), sqfree(ss*QQ(mm)^2), sqfree(ss*QQ(nn)^2))
    say("(m,n)=(%d,%d): s=%s b=%s требуемый класс %s" % (mm,nn,ss,bb,rq))
    # есть ли рациональное t с квадратными F_i (малый перебор)
    hit = []
    for q in range(1,400):
        for p in range(0,400):
            num0 = mm^2*q^2 + nn^2*p^2
            num4 = ss*(q^2+p^2)
            num8 = nn^2*q^2 + mm^2*p^2
            if not (QQ(num0).is_square() and QQ(num4).is_square() and QQ(num8).is_square()):
                continue
            if gcd(p,q) == 1:
                hit.append((p,q))
    say("   малый перебор t=p/q<400: точки C(Q):", hit[:10], "..." if len(hit)>10 else "")
    if hit:
        p,q = hit[0]
        tv = QQ(p)/QQ(q)
        Xv = bb*tv^2
        try:
            PP = EE.lift_x(Xv)
            rts = ee
            def d2(P, rts, EEE):
                if P == EEE(0): return (ZZ(1),ZZ(1),ZZ(1))
                x = P[0]; o=[]
                for i in range(3):
                    if x == rts[i]:
                        j,k = [a for a in range(3) if a!=i]
                        o.append(sqfree((rts[i]-rts[j])*(rts[i]-rts[k])))
                    else:
                        o.append(sqfree(x-rts[i]))
                return tuple(o)
            say("   КОНТРОЛЬ: t=%s -> X=%s, delta=%s; совпало с требуемым %s: %s"
                % (tv, Xv, d2(PP,rts,EE), rq, d2(PP,rts,EE)==rq))
        except Exception as ex:
            say("   lift_x не удался:", ex)

# ===========================================================================
hdr("ИТОГ")
# ===========================================================================
say(json.dumps(VERDICT, indent=1, ensure_ascii=False))
say("")
say("время, с:", round(time.time()-T0, 1))
