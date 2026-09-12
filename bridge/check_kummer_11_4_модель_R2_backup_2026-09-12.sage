# -*- coding: utf-8 -*-
# =====================================================================================
#  ПРОВЕРКА/ПОПЫТКА ОПРОВЕРЖЕНИЯ ЗАЯВЛЕНИЯ CODEX:  C_{11,4}(Q) = пусто
#  Угол атаки: МОДЕЛЬ И ТОЖДЕСТВА.  Всё строится с нуля, числа Codex используются
#  ТОЛЬКО в финальной секции K (сверка).
#
#  Автор прогона: Claude (раунд «модель»), 2026-09-12.
#  Предыдущая версия этого файла сохранена как check_kummer_11_4_модель_PREV_backup_2026-09-12.sage
#
#  Метки статуса печатаются у каждого пункта:
#     [ДОК]  — доказано символьно/точной арифметикой в этом скрипте
#     [ЧИС]  — проверено численно (конечный перебор, приближение)
#     [ПО]   — доказано ПРИ УСЛОВИИ (GRH / корректность внешней библиотеки)
#     [НАБЛ] — наблюдение, не доказательство
# =====================================================================================
import sys, os, time
from math import gcd as _gcd, isqrt as _isqrt

T_START = time.time()
RESULTS = []          # (метка, код, текст, ok)


def hdr(t):
    print("\n" + "=" * 86)
    print(t)
    print("=" * 86)
    sys.stdout.flush()


def rec(tag, code, text, ok, extra=""):
    RESULTS.append((tag, code, text, ok))
    print("   [%s %s] %-4s %s%s" % ("OK  " if ok else "ПРОВАЛ", tag, code, text,
                                    ("  -- " + str(extra)) if extra else ""))
    sys.stdout.flush()


def note(s):
    print("   . " + s)
    sys.stdout.flush()


BADP = [2, 3, 5, 7, 11, 137]     # плохие простые для (11,4); проверяются в секции B


def sqcls(q, primes=None):
    """класс квадратов рационального числа q != 0 как бесквадратное целое (со знаком).
    Если задан список primes, класс вычисляется быстро: снимаем валюации по этим
    простым, остаток ОБЯЗАН быть точным квадратом (это ещё и контроль)."""
    q = QQ(q)
    if q == 0:
        raise ValueError("sqcls(0)")
    v = ZZ(q.numerator() * q.denominator())
    if primes is None:
        return v.squarefree_part()
    sgn = -1 if v < 0 else 1
    v = abs(v)
    c = 1
    for p in primes:
        k = v.valuation(p)
        if k % 2:
            c *= p
        v = v // p**k
    if not v.is_square():
        # выходит за группу, порождённую primes — считаем честно (медленно)
        return ZZ(q.numerator() * q.denominator()).squarefree_part()
    return sgn * c


# =====================================================================================
hdr("A.  СИМВОЛЬНЫЕ ТОЖДЕСТВА В Q(m,n)[t].  НИКАКИХ ЧИСЕЛ.")
# =====================================================================================
# Всё определяется заново по условию задачи; ничего не копируется.
Fmn.<M, N> = FractionField(PolynomialRing(QQ, 'M,N'))
Rt.<t> = PolynomialRing(Fmn)

s_sym = (M**2 + N**2) / 2
F0 = M**2 + N**2 * t**2
F4 = s_sym * (1 + t**2)
F8 = N**2 + M**2 * t**2
b_sym = s_sym * M**2 * N**2
e_sym = [-b_sym, -s_sym * M**4, -s_sym * N**4]        # порядок e1,e2,e3 ИЗ УСЛОВИЯ
X_sym = b_sym * t**2

note("s  = %s" % s_sym)
note("b  = s*m^2*n^2 = %s" % b_sym)
note("e1 = -b       = %s" % e_sym[0])
note("e2 = -s*m^4   = %s" % e_sym[1])
note("e3 = -s*n^4   = %s" % e_sym[2])
note("X  = b*t^2")

# --- A1..A3: три тождества, КАЖДОЕ с явным множителем
id_data = [
    ("A1", 0, "F4", F4, M**2 * N**2,                 "(mn)^2 * F4  = (mn)^2 * u4^2"),
    ("A2", 1, "F0", F0, s_sym * M**2,                "s*m^2 * F0   = s*m^2 * u0^2"),
    ("A3", 2, "F8", F8, s_sym * N**2,                "s*n^2 * F8   = s*n^2 * u8^2"),
]
for code, i, fname, Fpoly, coef, descr in id_data:
    lhs = X_sym - e_sym[i]
    rhs = coef * Fpoly
    rec("ДОК", code, "X - e%d  ==  %s" % (i + 1, descr), bool(lhs - rhs == 0),
        "разность = %s" % (lhs - rhs))

# --- A4: множитель восстанавливается однозначно делением (контроль, что он не «подогнан»)
for code, i, fname, Fpoly, coef, descr in id_data:
    q = (X_sym - e_sym[i]) / Fpoly
    rec("ДОК", "A4" + code[-1], "(X-e%d)/%s — константа по t и равна %s" % (i + 1, fname, coef),
        bool(Rt(q.numerator()).degree() == 0 if q in Rt else True) and bool(q == coef), "q = %s" % q)

# --- A5: НЕГАТИВНЫЙ контроль: «перепутанные» пары НЕ дают тождества
bad = 0
for i in range(3):
    for (fname, Fpoly) in [("F4", F4), ("F0", F0), ("F8", F8)]:
        q = (X_sym - e_sym[i]) / Fpoly
        is_const = (q.denominator().degree() == 0 and q.numerator().degree() == 0) if hasattr(q, "denominator") else False
        right = (i == 0 and fname == "F4") or (i == 1 and fname == "F0") or (i == 2 and fname == "F8")
        if is_const != right:
            bad += 1
            note("   НЕОЖИДАННО: (X-e%d)/%s const=%s" % (i + 1, fname, is_const))
rec("ДОК", "A5", "соответствие e1<->F4, e2<->F0, e3<->F8 ОДНОЗНАЧНО (все прочие пары не дают константы)",
    bad == 0, "нарушений %d" % bad)

# --- A6: кубика и V
V_sym = b_sym * t**0  # заглушка; проверяем тождество для V^2
lhs = (b_sym)**2 * F0 * F4 * F8
rhs = (X_sym - e_sym[0]) * (X_sym - e_sym[1]) * (X_sym - e_sym[2])
rec("ДОК", "A6", "(b*u0*u4*u8)^2 == (X-e1)(X-e2)(X-e3)  тождественно", bool(lhs - rhs == 0),
    "разность = %s" % (lhs - rhs))

# --- A7: следствие про классы квадратов, символьно
note("Следствие (символьно): для ЛЮБОЙ точки C(Q) с конечным t")
note("   X-e1 = (mn)^2 u4^2  -> класс 1        (множитель (mn)^2 — точный квадрат)")
note("   X-e2 = s*m^2 * u0^2 -> класс [s]      (множитель m^2 — точный квадрат)")
note("   X-e3 = s*n^2 * u8^2 -> класс [s]      (множитель n^2 — точный квадрат)")
ok = (M**2 * N**2).is_square() if hasattr(M**2 * N**2, "is_square") else True
rec("ДОК", "A7", "множители (mn)^2, m^2, n^2 — квадраты в Q(m,n) => класс = (1,[s],[s])", True,
    "именно в этом порядке при e=(−b, −s m^4, −s n^4)")

# --- A8: дискриминанты F_i по t (для последующего «u_i != 0»)
for nm, Fp in [("F0", F0), ("F4", F4), ("F8", F8)]:
    P = Rt(Fp)
    note("%s: deg=%d, ст.коэф=%s, св.член=%s, disc_t=%s" % (nm, P.degree(), P.leading_coefficient(),
                                                            P[0], P.discriminant()))
rec("ДОК", "A8", "disc_t(F0)=-4m^2n^2, disc_t(F4)=-(m^2+n^2)^2, disc_t(F8)=-4m^2n^2 < 0 при m,n>0",
    bool(Rt(F0).discriminant() == -4 * M**2 * N**2) and bool(Rt(F8).discriminant() == -4 * M**2 * N**2)
    and bool(Rt(F4).discriminant() == -(M**2 + N**2)**2))

# --- A9: шесть точек ветвления различны  <=>  m != n
res_pairs = [("F0,F4", Rt(F0).resultant(Rt(F4))), ("F0,F8", Rt(F0).resultant(Rt(F8))),
             ("F4,F8", Rt(F4).resultant(Rt(F8)))]
for nm, r in res_pairs:
    note("Res(%s) = %s" % (nm, factor(r.numerator()) if r != 0 else 0))
rec("ДОК", "A9", "корни F0,F4,F8 попарно различны <=> m^2 != n^2 (символьно: резольвенты != 0)",
    all(r != 0 for _, r in res_pairs))


# =====================================================================================
hdr("B.  СПЕЦИАЛИЗАЦИЯ (m,n) = (11,4).  СОБСТВЕННЫЕ ЧИСЛА.")
# =====================================================================================
m, n = 11, 4
rec("ДОК", "B0", "gcd(m,n)=1", gcd(m, n) == 1)
s = QQ(m**2 + n**2) / 2
b = s * m**2 * n**2
e = [-b, -s * m**4, -s * n**4]
RX.<Xv> = PolynomialRing(QQ)
cub = (Xv - e[0]) * (Xv - e[1]) * (Xv - e[2])
E = EllipticCurve([0, cub[2], 0, cub[1], cub[0]])
note("s = %s,   b = %s" % (s, b))
note("e1 = %s ,  e2 = %s ,  e3 = %s" % tuple(e))
note("E: %s" % E)
rec("ДОК", "B1", "b целое", b in ZZ, "b=%s" % b)
rec("ДОК", "B2", "e1,e2,e3 попарно различны", len(set(e)) == 3)
rec("ДОК", "B3", "E неособа (disc != 0)", E.discriminant() != 0, "disc = %s" % factor(E.discriminant()))
rec("ДОК", "B4", "правая часть E совпадает с (X-e1)(X-e2)(X-e3)",
    E.defining_polynomial()(Xv, 0, 1) == -cub or (Xv**3 + E.a2() * Xv**2 + E.a4() * Xv + E.a6()) == cub)
tors = E.torsion_subgroup()
rec("ДОК", "B5", "E(Q)_tors = (Z/2)^2, т.е. 2-кручение всё рационально",
    tuple(tors.invariants()) == (2, 2), "инварианты %s" % (tors.invariants(),))
two_tors_x = sorted([P[0] for P in E.torsion_points() if not P.is_zero() and 2 * P == E(0)])
rec("ДОК", "B5b", "x-координаты 2-кручения = {e1,e2,e3}", set(two_tors_x) == set(e), two_tors_x)
note("|E(Q)/2E(Q)| = |E(Q)[2]| * 2^r = 4*2^r = 2^(r+2)  (верно для любой конечной T: |T/2T|=|T[2]|)")
sq_s = sqcls(s)
REQ = (1, sq_s, sq_s)
note("ТРЕБУЕМЫЙ КЛАСС (мой расчёт): %s     [s]=sqfree(%s)=%s" % (REQ.__str__(), s, sq_s))
rec("ДОК", "B6", "s НЕ является квадратом в Q", not QQ(s).is_square(), "s=%s" % s)
# положительность F_i при вещественном t
rec("ДОК", "B7", "F0,F4,F8 > 0 для всех вещественных t  =>  u0,u4,u8 != 0 на C(R)",
    all(RX(f).discriminant() < 0 and RX(f).leading_coefficient() > 0
        for f in [m**2 + n**2 * Xv**2, s * (1 + Xv**2), n**2 + m**2 * Xv**2]))
rec("ДОК", "B8", "X = b t^2 >= 0 > e_i  =>  образ точки C НИКОГДА не 2-кручение и не O",
    b > 0 and max(e) < 0, "max e_i = %s" % max(e))


# =====================================================================================
hdr("C.  СОБСТВЕННАЯ РЕАЛИЗАЦИЯ delta И КОНТРОЛЬ ЕЁ СВОЙСТВ")
# =====================================================================================
def delta(P, roots, primes=None):
    if P.is_zero():
        return (1, 1, 1)
    x = P[0]
    out = []
    for i in range(3):
        v = x - roots[i]
        if v == 0:
            j, k = [z for z in range(3) if z != i]
            v = (roots[i] - roots[j]) * (roots[i] - roots[k])
        out.append(sqcls(v, primes))
    return tuple(out)


def mulcls(a, c, primes=None):
    return tuple(sqcls(QQ(a[i]) * QQ(c[i]), primes) for i in range(3))


# генератор (получаю САМ, через PARI ellrank + насыщение)
Emin = E.minimal_model()
iso = Emin.isomorphism_to(E)
t0 = time.time()
pe = pari(Emin).ellrank(3)
gens_min = [Emin([QQ(P[0]), QQ(P[1])]) for P in pe[3]]
note("PARI ellrank(effort=3) на Emin: [r_lo, r_hi, s, points] = %s   (%.1f c)" % (pe, time.time() - t0))
sat = Emin.saturation(gens_min) if gens_min else ([], 1, 0)
note("saturation: индекс = %s, точки = %s" % (sat[1], sat[0]))
G = iso(sat[0][0]) if sat[0] else None
note("образующая свободной части на МОЕЙ модели E:  G = %s" % G)
rec("ДОК", "C0", "G лежит на E и имеет бесконечный порядок",
    G is not None and G in E and G.order() == Infinity)

pool = list(E.torsion_points())
if G is not None:
    for k in range(-4, 5):
        pool.append(k * G)
        for Tt in E.torsion_points():
            pool.append(k * G + Tt)
pool = list(set(pool))

badhom = None
cnt = 0
for P in pool[:24]:
    for Q in pool[:24]:
        if P.is_zero() or Q.is_zero() or (P + Q).is_zero():
            continue
        if P[0] in e or Q[0] in e or (P + Q)[0] in e:
            continue   # формула для 2-кручения проверяется отдельно (C1b)
        cnt += 1
        if mulcls(delta(P, e), delta(Q, e)) != delta(P + Q, e):
            badhom = (P, Q)
            break
    if badhom:
        break
rec("ДОК", "C1", "delta — гомоморфизм на проверенных парах (нет 2-кручения)", badhom is None,
    "пар проверено %d, контрпример %s" % (cnt, badhom))

badhom2 = None
for P in pool:
    for Tt in E.torsion_points():
        if mulcls(delta(P, e), delta(Tt, e)) != delta(P + Tt, e):
            badhom2 = (P, Tt)
rec("ДОК", "C1b", "delta — гомоморфизм и с участием 2-кручения (формула (ei-ej)(ei-ek))",
    badhom2 is None, "контрпример %s" % (badhom2,))
rec("ДОК", "C2", "произведение трёх координат delta всегда квадрат",
    all(QQ(prod(delta(P, e))).is_square() for P in pool))
rec("ДОК", "C3", "delta(O) = (1,1,1)", delta(E(0), e) == (1, 1, 1))
rec("ДОК", "C4", "delta(2P) = (1,1,1) для всех точек пула",
    all(delta(2 * P, e) == (1, 1, 1) for P in pool if not (2 * P).is_zero()))


# =====================================================================================
hdr("D.  ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ: пары (m,n), у которых точка на C ИЗВЕСТНА")
# =====================================================================================
note("Если порядок корней перепутан, требуемый класс окажется перестановкой и")
note("вывод перевернётся.  Нужны случаи с ЯВНОЙ точкой C(Q) и [s] != 1.")


def build(mm, nn):
    ss = QQ(mm**2 + nn**2) / 2
    bb = ss * mm**2 * nn**2
    rr = [-bb, -ss * mm**4, -ss * nn**4]
    cc = (Xv - rr[0]) * (Xv - rr[1]) * (Xv - rr[2])
    return ss, bb, rr, EllipticCurve([0, cc[2], 0, cc[1], cc[0]])


def cpoint_to_E(mm, nn, p, q):
    """t = p/q, возвращает (X,V) и u-координаты, если точка на C существует"""
    A0 = mm**2 * q**2 + nn**2 * p**2
    A8 = nn**2 * q**2 + mm**2 * p**2
    A4n = (mm**2 + nn**2) * (p**2 + q**2)          # = 2*s*(p^2+q^2)
    if not (ZZ(A0).is_square() and ZZ(A8).is_square()):
        return None
    if not QQ(QQ(A4n) / 2).is_square():
        return None
    u0 = QQ(ZZ(A0).sqrt()) / q
    u8 = QQ(ZZ(A8).sqrt()) / q
    u4 = QQ(QQ(A4n) / 2).sqrt() / q
    return (u0, u4, u8)


ctrl = []
found_ctrl = 0
# (a) пифагоровы (m,n): t=1 даёт F0=F4=F8=m^2+n^2
for (mm, nn) in [(15, 8), (20, 21), (7, 24), (119, 120), (44, 117), (3, 4), (5, 12), (9, 40)]:
    if gcd(mm, nn) != 1:
        continue
    u = cpoint_to_E(mm, nn, 1, 1)
    if u is None:
        continue
    ss, bb, rr, EE = build(mm, nn)
    tval = QQ(1)
    Xp = bb * tval**2
    Vp = bb * u[0] * u[1] * u[2]
    P = EE(Xp, Vp)
    d = delta(P, rr)
    want = (1, sqcls(ss), sqcls(ss))
    ctrl.append((mm, nn, tval, d, want, d == want))
    found_ctrl += 1
    note("(m,n)=(%3d,%3d) t=%s  u=%s  delta=%s  требуемый=%s  %s"
         % (mm, nn, tval, tuple(u).__str__(), d.__str__(), want.__str__(),
            "СОВПАЛО" if d == want else "!!! РАСХОЖДЕНИЕ !!!"))

# (b) ПОИСК контроля с t != +-1 — гораздо более жёсткий тест ориентации
note("")
note("поиск контрольной точки с t != +-1 (перебор m,n <= 60, |p|,q <= 120) ...")
t0 = time.time()
extra_ctrl = []
for mm in range(1, 61):
    for nn in range(1, mm):
        if gcd(mm, nn) != 1:
            continue
        for q in range(1, 121):
            for p in range(1, 121):
                if p == q or _gcd(p, q) != 1:
                    continue
                A0 = mm**2 * q**2 + nn**2 * p**2
                r0 = _isqrt(A0)
                if r0 * r0 != A0:
                    continue
                A8 = nn**2 * q**2 + mm**2 * p**2
                r8 = _isqrt(A8)
                if r8 * r8 != A8:
                    continue
                A4 = (mm**2 + nn**2) * (p**2 + q**2)
                if A4 % 2:
                    continue
                r4 = _isqrt(A4 // 2)
                if r4 * r4 != A4 // 2:
                    continue
                extra_ctrl.append((mm, nn, p, q))
note("найдено кандидатов с t != ±1: %d  (%.1f c)   %s" % (len(extra_ctrl), time.time() - t0,
                                                          extra_ctrl[:6]))
for (mm, nn, p, q) in extra_ctrl[:6]:
    u = cpoint_to_E(mm, nn, p, q)
    ss, bb, rr, EE = build(mm, nn)
    tval = QQ(p) / q
    P = EE(bb * tval**2, bb * u[0] * u[1] * u[2])
    d = delta(P, rr)
    want = (1, sqcls(ss), sqcls(ss))
    ctrl.append((mm, nn, tval, d, want, d == want))
    note("(m,n)=(%d,%d) t=%s  delta=%s  требуемый=%s  %s" % (mm, nn, tval, d.__str__(),
         want.__str__(), "СОВПАЛО" if d == want else "!!! РАСХОЖДЕНИЕ !!!"))

rec("ДОК", "D1", "все контрольные точки дают delta = (1,[s],[s]) ИМЕННО в этом порядке",
    len(ctrl) > 0 and all(c[5] for c in ctrl), "контролей %d" % len(ctrl))
rec("ДОК", "D2", "у контролей [s] != 1 => тест РАЗЛИЧАЕТ позицию 1 и позиции 2,3",
    all(c[4][1] != 1 for c in ctrl))


# =====================================================================================
hdr("E.  ОБРАЗ delta ДЛЯ (11,4) И ОТСУТСТВИЕ ТРЕБУЕМОГО КЛАССА")
# =====================================================================================
T1, T2 = [P for P in E.torsion_points() if not P.is_zero()][:2]
basis = [G, T1, T2]
img = set()
for c1 in range(2):
    for c2 in range(2):
        for c3 in range(2):
            P = c1 * G + c2 * T1 + c3 * T2
            img.add(delta(P, e))
img = sorted(img)
note("образ delta, порождённый G и двумя точками 2-кручения (%d классов):" % len(img))
for c in img:
    note("     %s" % (c.__str__()))
rec("ДОК", "E1", "ровно 8 РАЗЛИЧНЫХ классов от явных точек E(Q)", len(img) == 8)
grp_ok = all(mulcls(a, c) in img for a in img for c in img)
rec("ДОК", "E2", "найденные классы образуют подгруппу (замкнуты по умножению)", grp_ok)
rec("ДОК", "E3", "требуемый класс %s ОТСУТСТВУЕТ среди них" % (REQ.__str__()), REQ not in img)
note("ЛОГИКА: если rank E(Q) = 1, то |E(Q)/2E(Q)| = 2^3 = 8, delta инъективен на E(Q)/2E(Q),")
note("        8 найденных классов = ВЕСЬ образ, и требуемого класса в нём нет =>")
note("        ни одна точка C(Q) с конечным t не существует.")
note("ВНИМАНИЕ: весь вывод держится на равенстве rank = 1 (см. секцию F).")


# =====================================================================================
hdr("F.  РАНГ E(Q): ТРИ НЕЗАВИСИМЫЕ ДОРОЖКИ, ЧЕСТНЫЕ МЕТКИ")
# =====================================================================================
note("N = %s = %s" % (Emin.conductor(), factor(Emin.conductor())))
note("root number w = %s" % Emin.root_number())
sr = Emin.selmer_rank()
note("dim_F2 Sel_2(E/Q) = %s   (полный 2-спуск, eclib)" % sr)
note("для кривой с полным 2-кручением: dim Sel_2 = rank + 2 + dim Sha[2]")
note("=> rank + dim Sha[2] = %s" % (sr - 2))
rec("ДОК", "F1", "ЧИСТЫЙ 2-спуск даёт лишь rank <= %s, а НЕ rank <= 1" % (sr - 2), True,
    "dim Sel_2 = %s" % sr)
rec("ДОК", "F2", "rank >= 1 ДОКАЗАН предъявлением точки бесконечного порядка", G is not None)
note("")
note("ДОРОЖКА 1 (PARI ellrank): [%s, %s], dim Sha[2] = %s" % (pe[0], pe[1], pe[2]))
note("  верхняя граница 1 получена НЕ 2-спуском, а спариванием Кассельса-Тейта на Sel_2")
note("  (оно альтернирующее, ранг чётен; здесь ранг 2 => dim ker = 3 => rank <= 1).")
rec("ПО", "F3", "rank <= 1 — ПРИ УСЛОВИИ корректности реализации CT-спаривания в PARI",
    ZZ(pe[1]) == 1, "PARI r_hi = %s" % pe[1])
note("")
note("ДОРОЖКА 2 (нулевые суммы Бобера, analytic_rank_upper_bound):")
ub = None
for D in [1.0, 1.3, 1.6]:
    try:
        t0 = time.time()
        u_ = Emin.analytic_rank_upper_bound(max_Delta=D, adaptive=False, root_number=-1)
        note("   Delta=%.1f -> верхняя граница аналитического ранга = %s  (%.1f c)" % (D, u_, time.time() - t0))
        ub = u_ if ub is None else min(ub, u_)
    except Exception as ex:
        note("   Delta=%.1f сбой: %s" % (D, ex))
note("  ВАЖНО: docstring Sage: 'conditional on the Generalized Riemann Hypothesis'.")
note("  Значит это ПО GRH, а не безусловно.  Дальше: w=-1 => порядок нуля нечётен =>")
note("  ord_{s=1} L(E,s) = 1 => (Гросс-Загир + Колывагин) rank = 1 и Sha конечна.")
rec("ПО", "F4", "rank = 1 — ПО GRH (нулевые суммы) + модулярность + GZ/Колывагин",
    ub is not None and ZZ(ub) == 1, "верхняя граница = %s" % ub)
note("")
note("ДОРОЖКА 3 (аналитический ранг численно, БЕЗ строгой оценки ошибки):")
try:
    t0 = time.time()
    ar = Emin.analytic_rank()
    note("   Sage analytic_rank = %s  (%.1f c)" % (ar, time.time() - t0))
    rec("ЧИС", "F5", "численный аналитический ранг = 1 (без сертификата ошибки)", ZZ(ar) == 1, ar)
except Exception as ex:
    rec("ЧИС", "F5", "численный аналитический ранг", False, ex)
try:
    shan = Emin.sha().an_numerical()
    note("   аналитический порядок Sha (BSD, численно) ~ %s ; 2-спуск дал dim Sha[2] = %s => |Sha[2]| = %s"
         % (shan, pe[2], 2**ZZ(pe[2])))
    rec("ЧИС", "F6", "BSD-согласованность: #Sha_an ~ 4 = |Sha[2]|", abs(RR(shan) - 4) < 0.1, shan)
except Exception as ex:
    rec("ЧИС", "F6", "BSD-согласованность", False, ex)


# =====================================================================================
hdr("G.  ГЕОМЕТРИЯ C: ГЛАДКОСТЬ ОБЕИХ КАРТ, БЕСКОНЕЧНОСТЬ, РОД")
# =====================================================================================
# карта 1: (t,u0,u4,u8)
S1 = PolynomialRing(QQ, ['tt', 'a0', 'a4', 'a8'])
tt, a0, a4, a8 = S1.gens()
eqs1 = [a0**2 - (m**2 + n**2 * tt**2), a4**2 - s * (1 + tt**2), a8**2 - (n**2 + m**2 * tt**2)]
Jac1 = matrix(S1, [[f.derivative(v) for v in S1.gens()] for f in eqs1])
sing1 = S1.ideal(eqs1 + Jac1.minors(3))
d1 = sing1.dimension()
rec("ДОК", "G1", "аффинная карта t-конечно ГЛАДКА над Q~ (особое множество пусто)", d1 == -1,
    "dim особого множества = %s" % d1)

# карта 2: (z=1/t, U_i=u_i/t)
S2 = PolynomialRing(QQ, ['zz', 'A0', 'A4', 'A8'])
zz, A0, A4, A8 = S2.gens()
eqs2 = [A0**2 - (m**2 * zz**2 + n**2), A4**2 - s * (zz**2 + 1), A8**2 - (n**2 * zz**2 + m**2)]
Jac2 = matrix(S2, [[f.derivative(v) for v in S2.gens()] for f in eqs2])
sing2 = S2.ideal(eqs2 + Jac2.minors(3))
d2 = sing2.dimension()
rec("ДОК", "G2", "аффинная карта z=1/t ГЛАДКА над Q~", d2 == -1, "dim особого множества = %s" % d2)
note("склейка карт на t != 0:  z = 1/t,  U_i = u_i/t  — изоморфизм (обратное t=1/z, u_i=U_i/z).")
note("Проверка склейки подстановкой: u0^2=m^2+n^2t^2  =>  (u0/t)^2 = m^2/t^2+n^2 = m^2 z^2 + n^2 . OK")
note("Поэтому ГЛАДКАЯ ПРОЕКТИВНАЯ МОДЕЛЬ = объединение этих двух гладких аффинных карт,")
note("и точки 'над t=inf' — это в точности точки карты 2 с z=0.")
note("При z=0:  A0^2 = n^2 = %s,  A4^2 = s = %s,  A8^2 = m^2 = %s" % (n**2, s, m**2))
rec("ДОК", "G3", "над t=inf обязательно A4^2 = s; s не квадрат => рациональных точек НЕТ",
    not QQ(s).is_square())
note("число точек над z=0 над Q~: 2*2*2 = 8 (все три значения != 0 => нет ветвления),")
note("все они определены над Q(sqrt(s)) = Q(sqrt(%s)) и разбиты на 4 сопряжённые пары." % sq_s)
# подтверждение: sqrt(F0), sqrt(F8) лежат в Q((1/t)), sqrt(F4) — нет
Lser.<zs> = LaurentSeriesRing(QQ, default_prec=12)
r0 = (n**2 + m**2 * zs**2).sqrt()
r8 = (m**2 + n**2 * zs**2).sqrt()
rec("ДОК", "G4", "sqrt(F0)/t и sqrt(F8)/t лежат в Q[[1/t]] (ряды с рациональными коэф.)",
    (r0**2 - (n**2 + m**2 * zs**2)).is_zero() and (r8**2 - (m**2 + n**2 * zs**2)).is_zero(),
    "r0 = %s ..." % str(r0)[:40])
note("а sqrt(F4)/t = sqrt(s)*sqrt(1+z^2) требует sqrt(s) => расширение степени 2. Отсюда G3.")

# род: Риман-Гурвиц + перекрёстная проверка разложением (Кани-Розен)
note("")
note("РОД.  Накрытие C -> P^1_t группой (Z/2)^3, степень 8.")
Fs = {"F0": RX(m**2 + n**2 * Xv**2), "F4": RX(s * (1 + Xv**2)), "F8": RX(n**2 + m**2 * Xv**2)}
allroots = []
for nm, f in Fs.items():
    allroots += list(f.roots(CC, multiplicities=False))
distinct = len(set([tuple([RR(z.real()).n(30), RR(z.imag()).n(30)]) for z in allroots])) == 6
rec("ДОК", "G5", "6 точек ветвления попарно различны", distinct, [CC(z).n(20) for z in allroots])
note("над каждой точкой ветвления: 4 точки с e=2  => вклад 6*4*(2-1) = 24")
note("над t=inf ветвления нет (см. G3/G4)")
note("2g-2 = 8*(-2) + 24 = 8  =>  g = 5")
# Кани-Розен: сумма родов семи квадратичных подполей
subs = {"F0": Fs["F0"], "F4": Fs["F4"], "F8": Fs["F8"],
        "F0F4": Fs["F0"] * Fs["F4"], "F0F8": Fs["F0"] * Fs["F8"],
        "F4F8": Fs["F4"] * Fs["F8"], "F0F4F8": Fs["F0"] * Fs["F4"] * Fs["F8"]}
tot = 0
for nm, f in subs.items():
    d = f.degree()
    sqfree = f.is_squarefree()
    g_ = (d - 2) // 2 if d % 2 == 0 else (d - 1) // 2
    tot += g_
    note("   подполе sqrt(%-7s): deg=%d, squarefree=%s, род=%d" % (nm, d, sqfree, g_))
rec("ДОК", "G6", "сумма родов 7 квадратичных подполей = 5 (Кани-Розен) — совпало с Риман-Гурвицем",
    tot == 5, "сумма = %d" % tot)
indep = all(not (subs[k]).is_square() for k in subs)
rec("ДОК", "G7", "ни одно произведение непустого подмножества {F0,F4,F8} не квадрат => C ГЕОМЕТРИЧЕСКИ НЕПРИВОДИМА, степень 8",
    indep)


# =====================================================================================
hdr("H.  ПЕРЕСТАНОВКИ КОРНЕЙ: ломается ли вывод?")
# =====================================================================================
from itertools import permutations
allsafe = True
dangerous = []
for sg in permutations(range(3)):
    roots_p = [e[i] for i in sg]
    # требуемый класс пересчитывается В ТОМ ЖЕ порядке
    base_req = {0: 1, 1: sq_s, 2: sq_s}
    req_p = tuple(base_req[i] for i in sg)
    img_p = set()
    for c1 in range(2):
        for c2 in range(2):
            for c3 in range(2):
                P = c1 * G + c2 * T1 + c3 * T2
                img_p.add(delta(P, roots_p))
    inimg = req_p in img_p
    note("sg=%s корни=%s требуемый=%s  -> в образе: %s"
         % (sg.__str__(), [str(x) for x in roots_p], req_p.__str__(), inimg))
    if inimg:
        allsafe = False
rec("ДОК", "H1", "ни при какой СОГЛАСОВАННОЙ перестановке требуемый класс не попадает в образ", allsafe)
# несогласованный (ошибочный) вариант
for sg in permutations(range(3)):
    req_bad = tuple({0: 1, 1: sq_s, 2: sq_s}[i] for i in sg)
    if req_bad != REQ and req_bad in set(img):
        dangerous.append((sg, req_bad))
note("")
note("НЕСОГЛАСОВАННЫЙ случай (ровно та ошибка, которую мы ищем):")
for sg, rb in dangerous:
    note("   *** если требуемый класс ОШИБОЧНО записать как %s (перестановка %s при НЕпереставленных корнях)"
         % (rb.__str__(), sg.__str__()))
    note("       — он В ОБРАЗЕ ЕСТЬ, и всё 'доказательство' молча рушится.")
rec("НАБЛ", "H2", "тест НЕТРИВИАЛЕН: перестановка требуемого класса попадает в образ",
    len(dangerous) > 0, "опасных перестановок %d" % len(dangerous))
note("Но секции A1-A5 и D1 фиксируют соответствие e1<->F4, e2<->F0, e3<->F8 однозначно,")
note("а контрольные точки D подтверждают его на РЕАЛЬНЫХ точках C(Q). Ошибки порядка нет.")


# =====================================================================================
hdr("I.  ПРЯМАЯ ПОПЫТКА СЛОМАТЬ: ПОИСК ТОЧЕК C_{11,4}(Q) В ОГРОМНОМ ДИАПАЗОНЕ")
# =====================================================================================
note("t=p/q, gcd(p,q)=1.  Условия:")
note("   (i)   %d q^2 + %d p^2 = квадрат        (u0)" % (m**2, n**2))
note("   (ii)  %d q^2 + %d p^2 = квадрат        (u8)" % (n**2, m**2))
note("   (iii) s(p^2+q^2) квадрат  <=>  %s (p^2+q^2) квадрат  <=>  p^2+q^2 = %s k^2" % (2 * sq_s // 2, sq_s))
note("(iii) — коника; параметризуем ПОЛНОСТЬЮ через Z[i] (Гильберт 90):")
note("   %s = 15^2 + 7^2 ;  p+qi = (15+7i)(x+yi)^2 / (x^2+y^2) * (x^2+y^2)" % sq_s)
note("   p = 15(x^2-y^2) - 14xy,   q = 7(x^2-y^2) + 30xy,   k = x^2+y^2")
# контроль полноты параметризации грубым перебором
Bchk = 900
bf = set()
for p_ in range(0, Bchk + 1):
    for q_ in range(1, Bchk + 1):
        if _gcd(p_, q_) != 1:
            continue
        S_ = p_ * p_ + q_ * q_
        if S_ % sq_s:
            continue
        r_ = _isqrt(S_ // sq_s)
        if r_ * r_ == S_ // sq_s:
            bf.add((p_, q_))
par = set()
Achk = 160
for x_ in range(-Achk, Achk + 1):
    for y_ in range(0, Achk + 1):
        if (x_ == 0 and y_ == 0) or _gcd(x_, y_) != 1:
            continue
        p_ = 15 * (x_ * x_ - y_ * y_) - 14 * x_ * y_
        q_ = 7 * (x_ * x_ - y_ * y_) + 30 * x_ * y_
        g_ = _gcd(abs(p_), abs(q_))
        if g_ == 0:
            continue
        p_ //= g_
        q_ //= g_
        for (pp, qq) in [(p_, q_), (p_, -q_), (-p_, q_), (-p_, -q_)]:
            if qq > 0 and 0 <= pp <= Bchk and qq <= Bchk:
                par.add((pp, qq))
rec("ДОК", "I1", "параметризация коники ПОЛНА (совпала с грубым перебором в коробке %d)" % Bchk,
    bf == par, "перебор %d решений, параметризация %d" % (len(bf), len(par)))


def mkmask(Md):
    ss_ = [False] * Md
    for i in range(Md):
        ss_[(i * i) % Md] = True
    return ss_


MODS = [64, 63, 65, 11, 13]
MASKS = [mkmask(k) for k in MODS]


def issq_fast(v):
    if v < 0:
        return False
    for k, mk in zip(MODS, MASKS):
        if not mk[v % k]:
            return False
    r = _isqrt(v)
    return r * r == v


A_SEARCH = int(os.environ.get("A_SEARCH", "6000"))
t0 = time.time()
hits = []
tested = 0
maxpq = 0
for x_ in range(-A_SEARCH, A_SEARCH + 1):
    for y_ in range(0, A_SEARCH + 1):
        if (x_ == 0 and y_ == 0) or _gcd(x_, y_) != 1:
            continue
        d_ = x_ * x_ - y_ * y_
        ee_ = x_ * y_
        p_ = abs(15 * d_ - 14 * ee_)
        q_ = abs(7 * d_ + 30 * ee_)
        if p_ == 0 or q_ == 0:
            continue
        g_ = _gcd(p_, q_)
        p_ //= g_
        q_ //= g_
        tested += 1
        if p_ > maxpq:
            maxpq = p_
        if q_ > maxpq:
            maxpq = q_
        if issq_fast(121 * q_ * q_ + 16 * p_ * p_) and issq_fast(16 * q_ * q_ + 121 * p_ * p_):
            hits.append((p_, q_))
note("проверено параметров (x,y): %d, кандидатов t: %d, макс. |p|,q = %d, время %.1f c"
     % ((2 * A_SEARCH + 1) * (A_SEARCH + 1), tested, maxpq, time.time() - t0))
rec("ЧИС", "I2", "точек C_{11,4}(Q) НЕ найдено при |p|,q до %d (это НЕ доказательство отсутствия)" % maxpq,
    len(hits) == 0, "найдено: %s" % hits[:5])


# =====================================================================================
hdr("J.  ТОРСОР КРИТИЧЕСКОГО КЛАССА (1,%s,%s): ПОПЫТКА НАЙТИ НА НЁМ ТОЧКУ" % (sq_s, sq_s))
# =====================================================================================
note("2-накрытие для класса (d1,d2,d3):  x-e_i = d_i z_i^2")
note("   d1 z1^2 - d2 z2^2 = e2-e1 ,   d1 z1^2 - d3 z3^2 = e3-e1")
note("масштабируем корни на 4 (классы квадратов не меняются), чтобы всё стало целым.")
e4 = [4 * x for x in e]
Aq = e4[1] - e4[0]
Bq = e4[2] - e4[0]
note("4e = %s ;  A = e2-e1 = %s ;  B = e3-e1 = %s" % ([str(x) for x in e4], Aq, Bq))
Pxy.<xq, yq> = PolynomialRing(QQ)


def quartic_for(d):
    d1, d2, d3 = d
    C1 = Conic(QQ, [d1, -d2, -Aq])
    if not C1.has_rational_point():
        return None, None
    dp = C1.parametrization()[0].defining_polynomials()
    z1 = Pxy(dp[0]); z2 = Pxy(dp[1]); w = Pxy(dp[2])
    assert d1 * z1**2 - d2 * z2**2 - Aq * w**2 == 0
    Q = d3 * (d1 * z1**2 - Bq * w**2)          # = (d3 z3)^2
    return Q, (z1, z2, w)


def reduce_quartic(Q):
    co = [ZZ(Q.coefficient({xq: 4 - j, yq: j})) for j in range(5)]
    g_ = gcd(co)
    sqp = ZZ(g_).squarefree_part()
    co = [ZZ(c * sqp / g_) for c in co]
    F = sum(co[j] * xq**(4 - j) * yq**j for j in range(5))
    try:
        from sage.rings.polynomial.binary_form_reduce import smallest_poly
        Fr, Mtr = smallest_poly(F, norm_type='height', prec=200)
        co = [ZZ(Fr.coefficient({xq: 4 - j, yq: j})) for j in range(5)]
    except Exception as ex:
        note("   (редукция бинарной формы недоступна: %s)" % ex)
    return co


def quartic_invariants(co):
    a, bq, cq, dq, eq = co
    I = 12 * a * eq - 3 * bq * dq + cq**2
    J = 72 * a * cq * eq + 9 * bq * cq * dq - 27 * a * dq**2 - 27 * eq * bq**2 - 2 * cq**3
    return I, J


def search_quartic(co, H):
    hs = []
    for a_ in range(-H, H + 1):
        for c_ in range(1, H + 1):
            if _gcd(a_, c_) != 1:
                continue
            v = co[0] * a_**4 + co[1] * a_**3 * c_ + co[2] * a_**2 * c_**2 + co[3] * a_ * c_**3 + co[4] * c_**4
            if issq_fast(int(v)):
                hs.append((a_, c_))
    return hs


H_Q = int(os.environ.get("H_QUARTIC", "600"))
for dcls, label in [((105, 210, 2), "КОНТРОЛЬ (класс delta(G), точка ОБЯЗАНА быть)"),
                    (REQ, "КРИТИЧЕСКИЙ класс")]:
    dcls = tuple(dcls)
    note("")
    note("--- класс %s : %s" % (dcls.__str__(), label))
    Q, par3 = quartic_for(dcls)
    if Q is None:
        note("    коника d1 z1^2 - d2 z2^2 = A w^2 НЕ имеет рациональных точек => класс НЕ в образе")
        continue
    co = reduce_quartic(Q)
    note("    приведённая квартика: %s" % co)
    I, J = quartic_invariants(co)
    Jac = EllipticCurve([0, 0, 0, -27 * I, -27 * J])
    ok_j = (Jac.j_invariant() == E.j_invariant())
    quad_tw = Jac.is_isomorphic(Emin) or any(Jac.quadratic_twist(dd).is_isomorphic(Emin)
                                             for dd in [-1, 2, -2, sq_s, -sq_s, 3, -3, 6, -6])
    note("    инварианты I,J -> якобиан Y^2=X^3-27IX-27J ; j совпал с j(E): %s ; изоморфен E (возможно, после квадр. кручения): %s"
         % (ok_j, quad_tw))
    t0 = time.time()
    hs = search_quartic(co, H_Q)
    note("    поиск точек на квартике, |a|,c <= %d : %d штук, %.1f c ; %s"
         % (H_Q, len(hs), time.time() - t0, hs[:6]))
    if dcls == (105, 210, 2):
        rec("ДОК", "J1", "КОНТРОЛЬ: на торсоре известного класса точки НАХОДЯТСЯ (машина торсоров работает)",
            len(hs) > 0, hs[:3])
        rec("ДОК", "J2", "j-инвариант якобиана квартики = j(E) (торсор построен верно)", ok_j)
    else:
        rec("ЧИС", "J3", "на торсоре КРИТИЧЕСКОГО класса точек НЕ найдено при |a|,c <= %d (НЕ доказательство)" % H_Q,
            len(hs) == 0, hs[:3])
        rec("ДОК", "J4", "j-инвариант якобиана критической квартики = j(E)", ok_j)
        crit_co = co

# локальная разрешимость критической квартики — независимая проверка того,
# что класс лежит в группе Селмера (значит локального препятствия НЕТ)
note("")
note("локальная разрешимость критической квартики Y^2 = F(a,c):")


def has_Qp_point(co, p, kmax=4):
    Qp = Qp_ = None
    from sage.rings.padics.factory import Qp as _Qp
    K = _Qp(p, 40)
    for k in range(1, kmax + 1):
        M = p**k
        for a_ in range(M):
            for c_ in range(M):
                if a_ % p == 0 and c_ % p == 0:
                    continue
                v = co[0] * a_**4 + co[1] * a_**3 * c_ + co[2] * a_**2 * c_**2 + co[3] * a_ * c_**3 + co[4] * c_**4
                if v == 0:
                    return True, (a_, c_)
                if K(v).is_square():
                    return True, (a_, c_)
    return False, None


badp = sorted(set(ZZ(Emin.conductor()).prime_factors() + [2]))
loc = {}
for p in badp:
    okp, wit = has_Qp_point(crit_co, p, kmax=3 if p > 5 else 4)
    loc[p] = okp
    note("   p=%-4s Q_p-точка: %s  %s" % (p, okp, ("(a,c)=%s" % (wit,)) if wit else ""))
# вещественное место
Rpoly = RX(sum(crit_co[j] * Xv**(4 - j) for j in range(5)))
real_ok = any(Rpoly(RR(v)) > 0 for v in [-10**6, -100, -10, -1, 0, 1, 10, 100, 10**6]) or crit_co[0] > 0
loc['R'] = real_ok
note("   место R: F принимает положительные значения: %s" % real_ok)
rec("ЧИС", "J5", "критический класс ЛОКАЛЬНО РАЗРЕШИМ ВЕЗДЕ => он в Sel_2, локального препятствия НЕТ",
    all(loc.values()), loc)
note("СЛЕДСТВИЕ: исключение (11,4) НЕ может опираться на локальные соображения.")
note("Оно опирается ровно на rank E(Q) = 1, т.е. на то, что критический класс лежит в Sha[2].")


# =====================================================================================
hdr("K.  СВЕРКА С ЧИСЛАМИ CODEX (только здесь используются его величины)")
# =====================================================================================
CODEX = dict(s=QQ(137) / 2, b=132616, req=(1, 274, 274), rank=1, nclasses=8)
rec("ДОК", "K1", "s совпадает", s == CODEX['s'], "мой s=%s" % s)
rec("ДОК", "K2", "b совпадает", b == CODEX['b'], "мой b=%s" % b)
rec("ДОК", "K3", "требуемый класс совпадает", REQ == CODEX['req'], "мой %s" % (REQ.__str__()))
rec("ДОК", "K4", "число классов совпадает (8)", len(img) == CODEX['nclasses'])
rec("ПО", "K5", "rank = 1 совпадает (но см. F3/F4: это НЕ следствие 2-спуска)", ZZ(pe[1]) == CODEX['rank'])


# =====================================================================================
hdr("ИТОГ")
# =====================================================================================
fails = [(t_, c_, x_) for (t_, c_, x_, o_) in RESULTS if not o_]
for (t_, c_, x_, o_) in RESULTS:
    print("   %-6s %-6s %-4s %s" % ("OK" if o_ else "ПРОВАЛ", t_, c_, x_))
print("\n   провалов: %d" % len(fails))
for f in fails:
    print("     - %s %s %s" % f)
print("\n   общее время: %.1f c" % (time.time() - T_START))
