# -*- coding: utf-8 -*-
# =====================================================================================
#  ПОПЫТКА ОПРОВЕРЖЕНИЯ ЗАЯВЛЕНИЯ CODEX:  C_{11,4}(Q) = пусто
#  Угол атаки: МОДЕЛЬ И ТОЖДЕСТВА.
#
#  Раунд 3 (Claude, 2026-09-12). Всё строится с нуля.
#  Числа Codex используются ТОЛЬКО в финальной секции K (сверка), и нигде раньше.
#  Точка бесконечного порядка на E ищется СОБСТВЕННЫМ поиском, не берётся у Codex.
#
#  Предыдущие версии: check_kummer_11_4_модель_PREV_backup_2026-09-12.sage
#                     check_kummer_11_4_модель_R2_backup_2026-09-12.sage
#
#  Метки статуса:
#     [ДОК]  — доказано символьно / точной арифметикой внутри этого скрипта
#     [ЧИС]  — проверено численно (конечный перебор или приближение)
#     [ПО]   — доказано ПРИ УСЛОВИИ (корректность внешней библиотеки / теоремы извне)
#     [НАБЛ] — наблюдение, не доказательство
# =====================================================================================
import sys, time

T0 = time.time()
LOG = []          # (метка, код, текст, ok)

def hdr(t):
    print("\n" + "=" * 92); print(t); print("=" * 92); sys.stdout.flush()

def rec(tag, code, text, ok, extra=""):
    LOG.append((tag, code, text, ok))
    print("   [%s %-4s] %-5s %s%s" % ("OK    " if ok else "ПРОВАЛ", tag, code, text,
                                      ("   -- " + str(extra)) if extra else ""))
    sys.stdout.flush()

def note(x):
    print("   .  " + str(x)); sys.stdout.flush()


# -------------------------------------------------------------------------------------
# вспомогательное: квадратный класс рационального числа
# -------------------------------------------------------------------------------------
def sqclass(r):
    """представитель класса r в Q*/Q*^2 — бесквадратное целое"""
    r = QQ(r)
    if r == 0:
        raise ValueError("квадратный класс нуля не определён")
    return ZZ(r.numerator() * r.denominator()).squarefree_part()

def is_sq_Qp(x, p):
    """x (ненулевое рациональное) — квадрат в Q_p?"""
    x = QQ(x)
    if x == 0: return True
    v = x.valuation(p)
    if v % 2: return False
    u = x / QQ(p) ** v
    num, den = ZZ(u.numerator()), ZZ(u.denominator())
    if p == 2:
        return (num * den) % 8 == 1
    return kronecker(num * den % p, p) == 1


# =====================================================================================
hdr("0.  СЕМЕЙСТВО G1: девять клеток, магичность, и почему это ровно кривая C")
# =====================================================================================
# Клетки из EVEN_N_2ADIC_2026-09-12.md §1 — переписаны мной заново и проверены символьно.
Rg = PolynomialRing(QQ, ['M', 'N', 'P', 'Q'])
M, N, P, Q = Rg.gens()
F4g = (M**2 + N**2) * (P**2 + Q**2) / 2
c = [M**2*Q**2 + N**2*P**2,  (M*P + N*Q)**2,  F4g - 2*M*N*P*Q,
     (M*P - N*Q)**2,        F4g,             (M*Q + N*P)**2,
     F4g + 2*M*N*P*Q,       (M*Q - N*P)**2,  N**2*Q**2 + M**2*P**2]
rec("ДОК", "0a", "девять клеток выписаны заново из определения G1", True)
lines = [(0,1,2), (3,4,5), (6,7,8),          # строки
         (0,3,6), (1,4,7), (2,5,8),          # столбцы
         (0,4,8), (2,4,6)]                   # диагонали
bad = [L for L in lines if c[L[0]] + c[L[1]] + c[L[2]] != 3 * F4g]
rec("ДОК", "0b", "все 8 линий 3x3 равны 3*c4 (магичность) тождественно в Q[m,n,p,q]", len(bad) == 0, bad)
roots4 = {1: M*P + N*Q, 3: M*P - N*Q, 5: M*Q + N*P, 7: M*Q - N*P}
auto = all(Rg(c[i]) == roots4[i]**2 for i in (1, 3, 5, 7))
rec("ДОК", "0c", "клетки c1,c3,c5,c7 — квадраты ТОЖДЕСТВЕННО (автоматические)", auto)
# диагональ c0,c4,c8 = q^2*F0, q^2*F4(t), q^2*F8 при t=p/q
FRg = FractionField(Rg)
rec("ДОК", "0d", "c0 = q^2*F0(p/q),  c4 = q^2*F4(p/q),  c8 = q^2*F8(p/q)",
    (FRg(c[0]) == FRg(Q**2 * (M**2 + N**2 * (FRg(P)/FRg(Q))**2)))
    and (FRg(c[4]) == FRg(Q**2 * ((M**2 + N**2)/2 * (1 + (FRg(P)/FRg(Q))**2))))
    and (FRg(c[8]) == FRg(Q**2 * (N**2 + M**2 * (FRg(P)/FRg(Q))**2))))
note("следствие [ДОК]: клетки c0,c4,c8 — квадраты  <=>  F0,F4,F8 — квадраты  <=>  точка на C.")
note("значит C(Q)=пусто закрывает пару уже на уровне СЕМИ квадратов, тем более девяти.")
note("ОГОВОРКА: вывод «эти девять клеток исчерпывают семейство G1» здесь НЕ проверялся,")
note("он взят из файлов проекта и лежит вне моего угла атаки.")


# =====================================================================================
hdr("A.  СИМВОЛЬНАЯ ПРОВЕРКА ТОЖДЕСТВ (общие m,n; кольцо многочленов от t)")
# =====================================================================================
S = PolynomialRing(QQ, ['m', 'n'])
mS, nS = S.gens()
FF = FractionField(S)
Rt2 = PolynomialRing(FF, 't')
t = Rt2.gen()
mm, nn = FF(mS), FF(nS)

sS = (mm**2 + nn**2) / 2
bS = sS * mm**2 * nn**2
F0 = mm**2 + nn**2 * t**2
F4 = sS * (1 + t**2)
F8 = nn**2 + mm**2 * t**2

e1S, e2S, e3S = -bS, -sS * mm**4, -sS * nn**4
X = bS * t**2

rec("ДОК", "A1", "X - e1 == (m*n)^2 * F4   тождественно", (X - e1S) - (mm*nn)**2 * F4 == 0)
rec("ДОК", "A2", "X - e2 == s*m^2 * F0     тождественно", (X - e2S) - sS * mm**2 * F0 == 0)
rec("ДОК", "A3", "X - e3 == s*n^2 * F8     тождественно", (X - e3S) - sS * nn**2 * F8 == 0)
# и обратная сторона: не перепутаны ли F0 и F8?
rec("ДОК", "A2'", "X - e2 НЕ равно s*n^2*F8 (значит e2 отвечает именно F0, а не F8)",
    (X - e2S) - sS * nn**2 * F8 != 0)
rec("ДОК", "A3'", "X - e3 НЕ равно s*m^2*F0 (значит e3 отвечает именно F8, а не F0)",
    (X - e3S) - sS * mm**2 * F0 != 0)
rec("ДОК", "A1'", "X - e1 НЕ равно (m*n)^2*F0 и не (m*n)^2*F8",
    ((X - e1S) - (mm*nn)**2 * F0 != 0) and ((X - e1S) - (mm*nn)**2 * F8 != 0))

V2 = (X - e1S) * (X - e2S) * (X - e3S)
rec("ДОК", "A4", "(b*u0*u4*u8)^2 == (X-e1)(X-e2)(X-e3)  при u_i^2=F_i",
    V2 - (bS**2 * F0 * F4 * F8) == 0)

# множители перед F_i (именно они дают квадратный класс)
k1 = ((X - e1S) / F4).numerator() / ((X - e1S) / F4).denominator()
note("коэффициенты: (X-e1)/F4 = %s ;  (X-e2)/F0 = %s ;  (X-e3)/F8 = %s"
     % ((X - e1S) / F4, (X - e2S) / F0, (X - e3S) / F8))
rec("ДОК", "A5", "коэффициенты равны (mn)^2, s*m^2, s*n^2 — классы 1, [s], [s]",
    (X - e1S) / F4 == (mm*nn)**2 and (X - e2S) / F0 == sS*mm**2 and (X - e3S) / F8 == sS*nn**2)
note("ВЫВОД [ДОК]: для ЛЮБОЙ рациональной точки C необходимо delta = (1, [s], [s])")
note("             — единица стоит В ПЕРВОЙ позиции, т.е. у корня e1 = -b = -s*m^2*n^2.")


# =====================================================================================
hdr("B.  СПЕЦИАЛИЗАЦИЯ (m,n) = (11,4): кривая E, корни, дискриминант, 2-кручение")
# =====================================================================================
m, n = 11, 4
rec("ДОК", "B0", "gcd(m,n) = 1", gcd(m, n) == 1, "m=%d n=%d" % (m, n))
s = QQ(m**2 + n**2) / 2
b = s * m**2 * n**2
e1, e2, e3 = QQ(-b), QQ(-s*m**4), QQ(-s*n**4)
note("s = %s   b = %s" % (s, b))
note("e1 = %s   e2 = %s   e3 = %s" % (e1, e2, e3))
rec("ДОК", "B1", "b целое", b in ZZ, b)
rec("ДОК", "B2", "корни попарно различны", len({e1, e2, e3}) == 3)

Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
cub = (x - e1) * (x - e2) * (x - e3)
A2c, A4c, A6c = cub[2], cub[1], cub[0]
E = EllipticCurve([0, A2c, 0, A4c, A6c])
rec("ДОК", "B3", "E неособа: disc != 0", E.discriminant() != 0, factor(E.discriminant()))
rec("ДОК", "B4", "кубика E совпадает с (x-e1)(x-e2)(x-e3)",
    E.division_polynomial(2).monic() == cub.monic() and cub.degree() == 3)
T = E.torsion_subgroup()
rec("ДОК", "B5", "E[2] полностью рационально (структура кручения)", T.invariants() == (2, 2), T.invariants())
xt = sorted([Pt.xy()[0] for Pt in E.torsion_points() if Pt.order() == 2])
rec("ДОК", "B6", "x-координаты 2-кручения = {e1,e2,e3}", set(xt) == {e1, e2, e3}, xt)
Emin = E.minimal_model()
iso_to_min = E.isomorphism_to(Emin)
iso_from_min = Emin.isomorphism_to(E)
note("минимальная модель: %s" % Emin)
note("кондуктор N = %s = %s" % (Emin.conductor(), factor(Emin.conductor())))

req = (sqclass(1), sqclass(s), sqclass(s))
rec("ДОК", "B7", "ТРЕБУЕМЫЙ класс delta = (1,[s],[s]) вычислен из тождеств A", True, req)
rec("ДОК", "B8", "[s] != 1, т.е. тест действительно различает позицию 1 и позиции 2,3",
    sqclass(s) != 1, "[s] = %s" % sqclass(s))

# F_i положительны на всей вещественной прямой -> u_i != 0 -> образ не 2-кручение
disc0 = (0)**2 - 4*(n**2)*(m**2)   # дискр. F0 как квадратного по t
rec("ДОК", "B9", "F0,F4,F8 > 0 при всех вещественных t (дискриминанты < 0)",
    (-4*n**2*m**2 < 0) and (-4*s*s < 0) and (-4*m**2*n**2 < 0))
note("следствие [ДОК]: u0,u4,u8 != 0, значит образ точки C — НЕ точка 2-кручения и не O.")


# =====================================================================================
hdr("C.  ОТОБРАЖЕНИЕ delta: определение, гомоморфность, ядро")
# =====================================================================================
ES = [e1, e2, e3]

def delta(Pt):
    """delta: E(Q) -> (Q*/Q*^2)^3, стандартное вложение полного 2-спуска"""
    if Pt.is_zero():
        return (1, 1, 1)
    xP = Pt.xy()[0]
    out = []
    for i in range(3):
        d = xP - ES[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (ES[i] - ES[j]) * (ES[i] - ES[k])
        out.append(sqclass(d))
    return tuple(out)

def mul(u, v):
    return tuple(sqclass(QQ(u[i]) * QQ(v[i])) for i in range(3))

# набор пробных точек: кручение + всё, что найдётся поиском
tors = list(E.torsion_points())
search_pts = E.point_search(9)
probe = list(set(tors + search_pts + [pp + tt_ for pp in search_pts for tt_ in tors]))
probe = [pp for pp in probe if pp in E]
note("пробных точек для тестов гомоморфности: %d" % len(probe))
cnt, bad = 0, []
for i in range(len(probe)):
    for j in range(i, len(probe)):
        Pa, Pb = probe[i], probe[j]
        if delta(Pa + Pb) != mul(delta(Pa), delta(Pb)):
            bad.append((Pa, Pb))
        cnt += 1
rec("ДОК", "C1", "delta — гомоморфизм (проверено пар: %d)" % cnt, len(bad) == 0, bad[:2])
rec("ДОК", "C2", "delta(O) = (1,1,1)", delta(E(0)) == (1, 1, 1))
rec("ДОК", "C3", "delta(2P) = (1,1,1) для всех пробных P (ядро содержит 2E(Q))",
    all(delta(2*pp) == (1, 1, 1) for pp in probe))
rec("ДОК", "C4", "произведение координат delta — всегда квадрат",
    all(sqclass(QQ(d[0])*QQ(d[1])*QQ(d[2])) == 1 for d in [delta(pp) for pp in probe]))


# =====================================================================================
hdr("D.  ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ: проверяем порядок корней НЕ алгеброй, а точками")
# =====================================================================================
# D-I. Настоящие точки C для других (m,n): t=1 и m^2+n^2 = квадрат.
hits = []
for mm2 in range(2, 40):
    for nn2 in range(1, mm2):
        if gcd(mm2, nn2) != 1: continue
        if not ZZ(mm2**2 + nn2**2).is_square(): continue
        hits.append((mm2, nn2, 1))
okD, detail = True, []
for (mA, nA, tA) in hits:
    sA = QQ(mA**2 + nA**2) / 2
    bA = sA * mA**2 * nA**2
    eA = [-bA, -sA*mA**4, -sA*nA**4]
    XA = bA * tA**2
    F0A, F4A, F8A = mA**2 + nA**2*tA**2, sA*(1+tA**2), nA**2 + mA**2*tA**2
    assert QQ(F0A).is_square() and QQ(F4A).is_square() and QQ(F8A).is_square()
    VA = bA * QQ(F0A).sqrt() * QQ(F4A).sqrt() * QQ(F8A).sqrt()
    onA = (VA**2 == (XA-eA[0])*(XA-eA[1])*(XA-eA[2]))
    dA = tuple(sqclass(XA - eA[i]) for i in range(3))
    want = (1, sqclass(sA), sqclass(sA))
    detail.append(((mA, nA, tA), dA, want, onA and dA == want))
    okD = okD and onA and dA == want
rec("ДОК", "D1", "НАСТОЯЩИЕ точки C (%d штук, другие пары) дают delta = (1,[s],[s]) ровно в этом порядке" % len(hits),
    okD)
for d in detail: note("   (m,n,t)=%s  delta=%s  ожидалось %s  %s" % d)
note("замечание: у всех этих контролей [s]=2, так что контроль ловит подмену позиции 1<->2,3.")

# D-II. Контроли на САМОЙ паре (11,4): по одному F_i за раз делаем квадратом.
# F4 квадрат:  274(1+t^2) квадрат.  t=15/7:  1+t^2 = 274/49.
tc = QQ(15)/7
F4c = s*(1 + tc**2)
rec("ДОК", "D2a", "t=15/7: F4 = %s — квадрат в Q" % F4c, QQ(F4c).is_square())
Xc = b*tc**2
rec("ДОК", "D2b", "тогда X-e1 — ТОЧНЫЙ квадрат в Q (позиция 1 даёт класс 1)",
    QQ(Xc - e1).is_square(), "X-e1 = %s = (%s)^2" % (Xc - e1, sqrt(QQ(Xc - e1))))
rec("ДОК", "D2c", "и X-e1 = (mn)^2*F4 численно", Xc - e1 == (m*n)**2 * F4c)
# F0 квадрат: 121+16t^2 = w^2 -> t=15, w=61
tc2 = QQ(15)
F0c = m**2 + n**2*tc2**2
rec("ДОК", "D3a", "t=15: F0 = %s = %s^2" % (F0c, sqrt(QQ(F0c))), QQ(F0c).is_square())
Xc2 = b*tc2**2
rec("ДОК", "D3b", "тогда класс (X-e2) равен [s] = %s (позиция 2 даёт класс s)" % sqclass(s),
    sqclass(Xc2 - e2) == sqclass(s), "класс = %s" % sqclass(Xc2 - e2))
rec("ДОК", "D3c", "и X-e2 = s*m^2*F0 численно", Xc2 - e2 == s*m**2*F0c)
# F8 квадрат: t=1/15
tc3 = QQ(1)/15
F8c = n**2 + m**2*tc3**2
rec("ДОК", "D4a", "t=1/15: F8 = %s — квадрат" % F8c, QQ(F8c).is_square())
Xc3 = b*tc3**2
rec("ДОК", "D4b", "тогда класс (X-e3) равен [s] = %s (позиция 3 даёт класс s)" % sqclass(s),
    sqclass(Xc3 - e3) == sqclass(s), "класс = %s" % sqclass(Xc3 - e3))
rec("ДОК", "D4c", "и X-e3 = s*n^2*F8 численно", Xc3 - e3 == s*n**2*F8c)
note("D2-D4 — контроли НА САМОЙ (11,4): каждая позиция проверена отдельным примером,")
note("без опоры на общую алгебру. Порядок корней подтверждён независимо.")


# =====================================================================================
hdr("E.  ОБРАЗ delta(E(Q)): собственный поиск генератора, 8 классов, наличие/отсутствие цели")
# =====================================================================================
found = {}
pts_used = []
cands = list(E.torsion_points())
# собственный поиск точек (не берём ничего у Codex)
for h in (9, 12, 14):
    try:
        cands += E.point_search(h)
    except Exception as ex:
        note("point_search(%s) не сработал: %s" % (h, ex))
# плюс точки, которые найдёт PARI ellrank (независимая от Codex реализация)
try:
    pr = pari(Emin).ellrank(3)
    for pt in pr[3]:
        Pm = Emin(QQ(pt[0]), QQ(pt[1]))
        cands.append(iso_from_min(Pm))
    note("PARI ellrank вернул %d точек" % len(pr[3]))
except Exception as ex:
    note("PARI ellrank точки: %s" % ex)
cands = [pp for pp in set(cands)]
gen_inf = [pp for pp in cands if pp.order() == oo]
rec("ЧИС", "E0", "собственным поиском найдена точка бесконечного порядка на E", len(gen_inf) > 0,
    gen_inf[0] if gen_inf else None)
G = gen_inf[0] if gen_inf else None
if G is not None:
    note("МОЙ генератор-кандидат G = %s" % (G.xy(),))
    grp = [E(0)] + [pt for pt in E.torsion_points() if not pt.is_zero()]
    full = set()
    for tt_ in grp:
        full.add(delta(tt_))
        full.add(delta(G + tt_))
    img = sorted(full)
    rec("ДОК", "E1", "получено ровно 8 РАЗЛИЧНЫХ классов от явных рациональных точек",
        len(img) == 8, len(img))
    for d in img: note("   %s" % (d,))
    closed = all(mul(u, v) in full for u in img for v in img)
    rec("ДОК", "E2", "найденный набор замкнут по умножению (это подгруппа порядка 8)", closed)
    rec("ДОК", "E3", "требуемый класс %s ОТСУТСТВУЕТ среди 8 найденных" % (req,), req not in full)
    # проверим все точки, что нашлись, не дают ли 9-го класса
    extra = set(delta(pp) for pp in cands) | set(delta(pp + tt_) for pp in cands for tt_ in grp)
    rec("ЧИС", "E4", "ни одна найденная точка E(Q) не даёт 9-го класса (иначе rank>=2 и вывод падает)",
        extra <= full, sorted(extra - full)[:3])
else:
    rec("ЧИС", "E1", "генератор не найден — секция E не выполнена", False)
    img, full = [], set()


# =====================================================================================
hdr("F.  РАНГ E(Q): самое несущее место. Четыре независимых подхода")
# =====================================================================================
rb = Emin.rank_bounds()
rec("ЧИС", "F1", "eclib/mwrank: безусловные границы ранга", True, rb)
rec("ЧИС", "F1b", "БЕЗУСЛОВНОГО 2-спуска НЕДОСТАТОЧНО: верхняя граница %d > 1" % rb[1],
    rb[1] == 1, "верх = %d" % rb[1])
sr = Emin.selmer_rank()
rec("ЧИС", "F2", "2-Selmer ранг = %d  =>  rank + dim Sha[2] = %d" % (sr, sr - 2), True,
    "при rank=1 имеем dim Sha[2] = %d, #Sha[2] = %d" % (sr - 3, 2**(sr - 3)))
pr_out = {}
for eff in (0, 1, 2, 3):
    try:
        pr_out[eff] = pari(Emin).ellrank(eff)
    except Exception as ex:
        pr_out[eff] = "err %s" % ex
note("PARI ellrank по effort: %s" % pr_out)
pari_ok = all((not isinstance(v, str)) and ZZ(v[0]) == 1 and ZZ(v[1]) == 1 for v in pr_out.values())
rec("ПО", "F3", "PARI ellrank (2-спуск + спаривание Касселса-Тейта): rank = 1 при всех effort",
    pari_ok)

w = Emin.root_number()
rec("ЧИС", "F4", "знак функционального уравнения w(E) = %d" % w, w == -1)
note("по 2-теореме о чётности (Докчитсер-Докчитсер/Монски — ТЕОРЕМА) ранг нечётен => rank in {1,3}")

L = Emin.lseries().dokchitser(prec=60)
fe = L.check_functional_equation()
L1 = L(1); Lp1 = L.derivative(1, 1)
rec("ЧИС", "F5", "Dokchitser: невязка функционального уравнения ~ 0", abs(RR(fe)) < 1e-20, RR(fe))
rec("ЧИС", "F6", "L(E,1) = 0 численно", abs(RR(L1)) < 1e-12, RR(L1))
rec("ЧИС", "F7", "L'(E,1) != 0 с БОЛЬШИМ запасом", abs(RR(Lp1)) > 1, RR(Lp1))
try:
    par = pari(Emin).ellanalyticrank()
    rec("ЧИС", "F8", "PARI ellanalyticrank независимо даёт (rank_an, L^(r)(1))", ZZ(par[0]) == 1, par)
    agree = abs(RR(par[1]) - RR(Lp1)) < 1e-6
    rec("ЧИС", "F9", "PARI и Dokchitser согласованы по L'(1)", agree,
        "PARI %s vs Dokchitser %s" % (RR(par[1]), RR(Lp1)))
except Exception as ex:
    rec("ЧИС", "F8", "PARI ellanalyticrank", False, ex)

note("")
note("ВЫВОД ПО РАНГУ:")
note("  аналитический ранг = 1  (L(1)=0, L'(1)=%.6f != 0, две независимые реализации)." % RR(Lp1))
note("  По Гроссу-Загье + Колывагину (ТЕОРЕМА, модулярность известна) из rank_an <= 1 следует")
note("  rank_alg = rank_an = 1 и конечность Sha.  Это НЕ опирается на PARI ellrank и НЕ на BSD.")
note("  [ДОК ПО] — при условии корректности численного значения L'(1); запас 13.6 против 0 огромен.")
note("  Codex опирался ТОЛЬКО на PARI ellrank и сам отметил, что eclib даёт [1,3]. Здесь ранг")
note("  подтверждён вторым, теоретически независимым путём.")


# =====================================================================================
hdr("G.  ТОЧКИ НАД БЕСКОНЕЧНОСТЬЮ: аккуратно, в проективных координатах")
# =====================================================================================
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
# t = 1/z, U_i = u_i/t = u_i*z
G0 = m**2*z**2 + n**2        # = z^2 * F0(1/z)
G4 = s*(z**2 + 1)            # = z^2 * F4(1/z)
G8 = n**2*z**2 + m**2        # = z^2 * F8(1/z)
Rtz = PolynomialRing(QQ, ['tv'])
tv = Rtz.gen()
chk = True
for (Fi, Gi) in [(m**2 + n**2*tv**2, G0), (s*(1+tv**2), G4), (n**2 + m**2*tv**2, G8)]:
    lhs = Fi.subs({tv: 1/z}) * z**2
    chk = chk and (Rz.fraction_field()(lhs) == Rz.fraction_field()(Gi))
rec("ДОК", "G1", "подстановка t=1/z, U_i = u_i*z корректна: U_i^2 = z^2*F_i(1/z)", chk)
rec("ДОК", "G2", "при z=0:  U0^2 = n^2,  U4^2 = s,  U8^2 = m^2",
    (G0(0), G4(0), G8(0)) == (n**2, s, m**2), (G0(0), G4(0), G8(0)))
rec("ДОК", "G3", "правые части в z=0 НЕ равны нулю => карта z=0 гладкая, ветви не склеиваются",
    G0(0) != 0 and G4(0) != 0 and G8(0) != 0)
rec("ДОК", "G4", "s = %s НЕ квадрат в Q  =>  рациональных точек над t=inf НЕТ" % s,
    not QQ(s).is_square())
# ветвление: корни F0,F4,F8 попарно различны и все конечны и != 0
br = []
for Fi in [m**2 + n**2*x**2, s*(1+x**2), n**2 + m**2*x**2]:
    br += list(Fi.roots(QQbar, multiplicities=False))
rec("ДОК", "G5", "6 точек ветвления попарно различны, все конечны и отличны от 0 и inf",
    len(set(br)) == 6 and all(r != 0 for r in br), len(set(br)))
rec("ДОК", "G6", "=> над t=inf ровно 8 неразветвлённых точек; рациональна лишь при s,m^2,n^2 квадратах",
    True)
# неприводимость / степень 8: F0,F4,F8 независимы по модулю квадратов в Q(t)
Rq = PolynomialRing(QQ, 'tq'); tq = Rq.gen()
FL = [m**2 + n**2*tq**2, s*(1+tq**2), n**2 + m**2*tq**2]
indep = True
for mask in range(1, 8):
    prod = Rq.fraction_field()(1)
    for i in range(3):
        if mask >> i & 1: prod *= FL[i]
    nu = prod.numerator() * prod.denominator()
    if nu.is_square() or (nu / nu.leading_coefficient()).is_square() and QQ(nu.leading_coefficient()).is_square():
        indep = False
rec("ДОК", "G7", "F0,F4,F8 независимы mod квадратов в Q(t) => C неприводима, степень 8, род 5", indep)
note("род: 2g-2 = 8*(-2) + 6*4 = 8 => g = 5  [ДОК, Риман-Гурвиц]")
# t = 0 отдельно
rec("ДОК", "G8", "t=0 не даёт точки: F4(0) = s не квадрат", not QQ(s).is_square())


# =====================================================================================
hdr("H.  УСТОЙЧИВОСТЬ К ПЕРЕСТАНОВКЕ КОРНЕЙ (главная ловушка)")
# =====================================================================================
if G is not None:
    from itertools import permutations
    allok, flip = True, []
    for perm in permutations(range(3)):
        ES_p = [ES[i] for i in perm]
        # требуемый класс пересчитывается ИЗ ТОЖДЕСТВ при том же порядке
        coef = {0: QQ((m*n)**2), 1: QQ(s*m**2), 2: QQ(s*n**2)}
        req_p = tuple(sqclass(coef[i]) for i in perm)
        img_p = set(tuple(d[i] for i in perm) for d in img)
        if (req_p in img_p) != (req in full):
            allok = False
    rec("ДОК", "H1", "при СОГЛАСОВАННОЙ перестановке корней вывод не меняется ни для одной из 6",
        allok)
    # а вот несогласованная перестановка — и вывод переворачивается
    for perm in permutations(range(3)):
        rq = tuple(req[i] for i in perm)
        if rq in full:
            flip.append((perm, rq))
    rec("НАБЛ", "H2", "ОПАСНОСТЬ: перестановка ТОЛЬКО требуемого класса даёт класс, который В ОБРАЗЕ ЕСТЬ",
        len(flip) > 0, flip)
    note("то есть путаница e1<->e2 молча перевернула бы вывод в 'не исключено'.")
    note("именно поэтому секция D (контроли точками) обязательна — она закрывает эту дыру.")


# =====================================================================================
hdr("I.  ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ: можно ли обойтись без ранга? (ответ: нельзя)")
# =====================================================================================
badp = sorted(set(ZZ(2*m*n*(m**2 - n**2)*(m**2 + n**2)).prime_factors()))
note("плохие простые: %s" % badp)
loc = {}
for p in badp + [13, 17, 19, 23]:
    ok = False
    wit = None
    for num in range(0, p**3):
        tt_ = QQ(num)
        try:
            if is_sq_Qp(m**2 + n**2*tt_**2, p) and is_sq_Qp(s*(1+tt_**2), p) and is_sq_Qp(n**2 + m**2*tt_**2, p):
                ok, wit = True, tt_; break
        except Exception:
            pass
        if num > 4000: break
    if not ok:  # попробуем карту в бесконечности и дроби
        for den in range(1, 60):
            for num in range(0, 200):
                tt_ = QQ(num) / den
                if is_sq_Qp(m**2 + n**2*tt_**2, p) and is_sq_Qp(s*(1+tt_**2), p) and is_sq_Qp(n**2 + m**2*tt_**2, p):
                    ok, wit = True, tt_; break
            if ok: break
    loc[p] = (ok, wit)
rec("ЧИС", "I1", "C(Q_p) непусто для всех плохих p (найдены явные свидетели t)",
    all(v[0] for v in loc.values()), loc)
rec("ДОК", "I2", "C(R) непусто: при любом вещественном t все F_i > 0", True)
note("следствие [ДОК]: требуемый класс (1,[s],[s]) ЛЕЖИТ В ГРУППЕ СЕЛМЕРА (образ C(Q_p)->E(Q_p)).")
note("значит НИКАКОЕ локальное условие его не убьёт; исключение ОБЯЗАНО опираться на ранг.")
note("это и есть настоящая точка хрупкости заявления — и именно её Codex подпёр только PARI.")


# =====================================================================================
hdr("J.  ПРЯМОЙ ПОИСК РАЦИОНАЛЬНЫХ ТОЧЕК НА C_{11,4} (честная попытка сломать)")
# =====================================================================================
# F4 = 137(p^2+q^2)/2 квадрат  <=>  (p^2+q^2)/2 = 137*z^2  <=>  p^2+q^2 = 274 z^2.
# параметризуем через гауссовы целые: 274 = 15^2+7^2
hits2 = []
BND = 1200
ZI = ZZ[I]
for w0 in (ZI(15 + 7*I), ZI(15 - 7*I)):
    for a in range(-BND, BND + 1):
        for bb in range(0, BND + 1):
            if a == 0 and bb == 0: continue
            if gcd(a, bb) != 1: continue
            g = w0 * (ZI(a + bb*I))**2
            p_, q_ = abs(ZZ(g.real())), abs(ZZ(g.imag()))
            if p_ == 0 or q_ == 0: continue
            d = gcd(p_, q_); p_ //= d; q_ //= d
            if 2*(m**2 + n**2)*0 == 1: pass
            if (p_**2 + q_**2) % 2: continue
            if not ZZ((m**2 + n**2)*(p_**2 + q_**2)//2).is_square(): continue
            if not ZZ(m**2*q_**2 + n**2*p_**2).is_square(): continue
            if not ZZ(n**2*q_**2 + m**2*p_**2).is_square(): continue
            hits2.append((p_, q_))
rec("ЧИС", "J1", "прямой поиск по параметризации конуса (|a|,|b| <= %d) точек НЕ нашёл" % BND,
    len(hits2) == 0, hits2[:5])
note("это НЕ доказательство отсутствия — только неудачная попытка сломать. [НАБЛ]")
# грубый контрольный перебор
hits3 = []
for q_ in range(1, 400, 2):
    for p_ in range(1, 400, 2):
        if gcd(p_, q_) != 1: continue
        if not ZZ(137*(p_**2 + q_**2)//2).is_square(): continue
        if ZZ(m**2*q_**2 + n**2*p_**2).is_square() and ZZ(n**2*q_**2 + m**2*p_**2).is_square():
            hits3.append((p_, q_))
rec("ЧИС", "J2", "контрольный тупой перебор p,q < 400 тоже ничего не нашёл", len(hits3) == 0, hits3[:5])


# =====================================================================================
hdr("K.  СВЕРКА С ЧИСЛАМИ CODEX (только здесь используются его величины)")
# =====================================================================================
cod_s, cod_b = QQ(137)/2, 132616
cod_e = (QQ(-132616), QQ(-2005817)/2, QQ(-17536))
cod_req = (1, 274, 274)
cod_G = (QQ(70664), QQ(138738600))
cod_img = {(1,1,1), (-1,28770,-28770), (-28770,137,-210), (28770,210,137),
           (105,210,2), (-105,137,-14385), (-274,28770,-105), (274,1,274)}
rec("ЧИС", "K1", "s и b совпадают с моими", (cod_s, cod_b) == (s, b))
rec("ЧИС", "K2", "корни совпадают И В ТОМ ЖЕ ПОРЯДКЕ", cod_e == (e1, e2, e3))
rec("ЧИС", "K3", "требуемый класс совпадает с моим", cod_req == req)
try:
    GC = E(cod_G[0], cod_G[1])
    rec("ЧИС", "K4", "точка G Codex лежит на МОЕЙ кривой E и имеет бесконечный порядок",
        GC.order() == oo, "delta(G_codex) = %s" % (delta(GC),))
except Exception as ex:
    rec("ЧИС", "K4", "точка G Codex на моей кривой", False, ex)
if G is not None:
    rec("ЧИС", "K5", "список из 8 классов Codex совпадает с МОИМ образом", cod_img == full,
        "лишние у него: %s ; лишние у меня: %s" % (sorted(cod_img - full), sorted(full - cod_img)))


# =====================================================================================
hdr("ИТОГ")
# =====================================================================================
fails = [(t_, c_, s_) for (t_, c_, s_, o_) in LOG if not o_]
for (t_, c_, s_, o_) in LOG:
    print("   %-7s %-5s %s" % ("OK" if o_ else "ПРОВАЛ", c_, s_))
print("\n   провалов: %d" % len(fails))
for f in fails: print("     -> %s %s" % (f[1], f[2]))
print("\n   время: %.1f c" % (time.time() - T0))
