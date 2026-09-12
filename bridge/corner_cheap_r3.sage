# -*- coding: utf-8 -*-
# СТАТУС: расчёт
# ДЛЯ: Claude / Codex / Grok / пользователь; продолжение bridge/corner_cheap.sage (раунд 1)
#       и bridge/corner_cheap_r2.sage (раунд 2)
# ИТОГ: дешёвые победы в УГЛОВОМ случае. Раунды 1–2 проверили 2 и 10 эллиптических факторов
#       из 21; здесь перебраны ВСЕ 21 фактора у обоих семейств на всех 127 парах, плюс
#       строгий (валюационный) локальный тест вместо наивного mod p.
# ОТМЕНЯЕТ: ничего; уточняет формулировку «способ 2 неприменим В ПРИНЦИПЕ» из раунда 1
# ПРОВЕРЕНО: клетки обоих семейств выведены здесь заново; контроли — (13,8) p=17, (16,5) p=41,
#            E_B = 62/60/4/1, класс (1,s,s), знак якобиана по Ленгу
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужен полный ответ «закрывается ли угловой случай дёшево»
#
# Запуск: sage /home/kep/magicKube/bridge/corner_cheap_r3.sage
#         CORNER_FAST=1 — только 12 пар вместо 127

import functools, itertools, os, time
from collections import Counter
print = functools.partial(print, flush=True)
T0 = time.time()
FAST = os.environ.get('CORNER_FAST', '') == '1'

def hdr(s):
    print("\n" + "=" * 100); print(s); print("=" * 100)
def sub(s):
    print("\n" + "-" * 100); print(s); print("-" * 100)

PAIRS = [(m, n) for m in range(2, 21) for n in range(1, m) if gcd(m, n) == 1]
assert len(PAIRS) == 127
SMALL = [(3,2),(5,2),(11,4),(13,8),(15,1),(15,8),(16,5),(19,5),(19,16),(7,1),(12,5),(20,19)]
SCAN = SMALL if FAST else PAIRS
LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
OPEN7 = [(13,8),(15,8),(16,5)]                 # открытые пары семиклеточной задачи (HANDOFF §2)
OPEN9 = [(15,8)]                               # единственная открытая для девятиклеточной

Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
RT = PolynomialRing(QQ, 'T'); Tv = RT.gen()

# ============================================================================================
hdr("§1. Клетки обоих семейств выведены заново (точные рациональные, без общего множителя)")
# ============================================================================================

def cells_G1(m, n):
    """РЁБЕРНОЕ семейство: полные пары — (c1,c7) и (c3,c5); свободен «икс» {0,2,4,6,8}."""
    s = QQ(m**2 + n**2) / 2
    F0 = m**2 + n**2*t**2
    F4 = s*(1 + t**2)
    F8 = n**2 + m**2*t**2
    return [F0, (m*t + n)**2, F4 - 2*m*n*t, (m*t - n)**2, F4,
            (m + n*t)**2, F4 + 2*m*n*t, (m - n*t)**2, F8]

def cells_G2(m, n):
    """УГЛОВОЕ семейство: полные пары — (c0,c8) и (c2,c6); свободен «крест» {1,3,4,5,7}."""
    s = QQ(m**2 + n**2) / 2
    al = QQ(3*n**2 - m**2) / 2
    be = QQ(3*m**2 - n**2) / 2
    c4 = s*(1 + t**2)
    return [(m*t + n)**2, al*t**2 + be, (m*t - n)**2, c4 - 4*m*n*t, c4,
            c4 + 4*m*n*t, (m + n*t)**2, be*t**2 + al, (m - n*t)**2]

FREE = {'G1': (0, 2, 4, 6, 8), 'G2': (1, 3, 4, 5, 7)}
CELLS = {'G1': cells_G1, 'G2': cells_G2}

for tag in ('G1', 'G2'):
    for (m, n) in [(3,2), (13,8), (16,5), (19,5)]:
        cc = CELLS[tag](m, n)
        assert all(sum(cc[i] for i in L) - 3*cc[4] == 0 for L in LINES), (tag, m, n)
        for i in range(9):
            if i not in FREE[tag]:
                assert cc[i].is_square(), (tag, m, n, i)   # автоматические квадраты
print("  [символьно] у обоих семейств все 8 линий = 3*центр тождественно по t;")
print("  четыре клетки вне «свободного» множества — квадраты многочленов тождественно.")
print("  G1 свободен  {0,2,4,6,8} (икс: два УГЛА + центр);  неквадратны могут быть углы.")
print("  G2 свободен  {1,3,4,5,7} (крест: два РЕБРА + центр); неквадратны могут быть рёбра.")

sub("§1a. чётность по t: от неё зависит, сколько эллиптических факторов вообще существует")
for tag in ('G1', 'G2'):
    cc = CELLS[tag](5, 2)
    ev = [i for i in FREE[tag] if cc[i](-t) == cc[i]]
    od = [i for i in FREE[tag] if cc[i](-t) != cc[i]]
    print(f"  {tag}: чётные по t свободные клетки {ev}; нечётные {od}")
print("  в обоих семействах картина ОДНА И ТА ЖЕ: три чётные клетки и пара нечётных,")
print("  переставляемых инволюцией t -> -t. Значит и набор факторов одинаков (§2).")

# ============================================================================================
hdr("§2. ПОЛНЫЙ СПИСОК эллиптических факторов девятиклеточной кривой (их 21, а не 2 и не 10)")
# ============================================================================================
print("""Девятиклеточная кривая X: все пять свободных клеток — квадраты. Это (Z/2)^5-накрытие
прямой, род 49. Её промежуточные факторы — кривые z^2 = произведение подмножества S клеток:
   |S| = 1  коника (род 0)               — информации нет
   |S| = 2  квартика (род 1)             — 10 штук, ЭЛЛИПТИЧЕСКИЕ
   |S| = 3  секстика (род 2)             — 10 штук
   |S| = 4  степень 8 (род 3)            —  5 штук
   |S| = 5  степень 10 (род 4)           —  1 штука
контроль рода: 10*1 + 10*2 + 5*3 + 1*4 = 49.
Но у ЧЁТНЫХ по t подмножеств произведение — многочлен от T = t^2, и кривая ещё раз
распадается инволюцией t -> -t:  z^2 = G(T)  и  w^2 = T*G(T)  (w = t*z).
   |S|=3 чётное: G кубика -> род 1 и T*G квартика -> род 1        (2 новых эллиптических)
   |S|=4 чётное: G квартика -> род 1 и T*G квинтика -> род 2      (1 новая эллиптическая)
   |S|=5 чётное: G квинтика -> род 2 и T*G секстика -> род 2      (новых нет)
   |S|=2 чётное: G квадратика -> род 0, T*G кубика -> род 1, но она 2-ИЗОГЕННА якобиану
                 самой квартики (накрытие неразветвлено), ранг тот же — не новая.
Итого РАЗНЫХ дешёвых мест, где может встретиться ранг 0:
   10 (пары) + 4*2 (чётные тройки) + 3 (чётные четвёрки) = 21.
Раунд 1 проверил 2 из 21 у G2 (и 1 у G1 — ту, на которой стоит весь проект);
раунд 2 проверил 10 из 21. Здесь — все 21.""")

def sq_free_scale(f):
    """умножить многочлен на КВАДРАТ так, чтобы коэффициенты стали целыми (кривая не меняется)"""
    d = lcm([QQ(c).denominator() for c in f.coefficients()])
    return (d**2) * f

def ell_from_cubic(g):
    """E: z^2 = g(T), deg g = 3 -> Вейерштрасс"""
    c = g.coefficients(sparse=False)
    g0, g1, g2, g3 = c[0], c[1], c[2], c[3]
    return EllipticCurve([0, g2*g3, 0, g1*g3**2, g0*g3**3])

def ell_from_quartic(q):
    """якобиан z^2 = q(x), deg q = 4:  y^2 = x^3 - 27 I x - 27 J"""
    c = q.coefficients(sparse=False); c = c + [0]*(5 - len(c))
    a, b, cc, d, e = c[4], c[3], c[2], c[1], c[0]
    I = 12*a*e - 3*b*d + cc**2
    J = 72*a*cc*e + 9*b*cc*d - 27*a*d**2 - 27*e*b**2 - 2*cc**3
    return EllipticCurve([0, 0, 0, -27*I, -27*J])

def even_in_t(f):
    return f(-t) == f

def to_T(f):
    """чётный многочлен от t -> многочлен от T = t^2"""
    c = f.coefficients(sparse=False)
    return RT([c[i] for i in range(0, len(c), 2)])

def factors_of(tag, m, n):
    """все эллиптические факторы: список (имя, кривая)"""
    cc = CELLS[tag](m, n); fr = FREE[tag]
    out = []
    for S in itertools.combinations(fr, 2):
        q = sq_free_scale(cc[S[0]] * cc[S[1]])
        if q.degree() == 4 and q.discriminant() != 0:
            out.append(("c%d c%d" % S, ell_from_quartic(q)))
    for k in (3, 4, 5):
        for S in itertools.combinations(fr, k):
            f = prod(cc[i] for i in S)
            if not even_in_t(f):
                continue
            G = sq_free_scale(to_T(f))
            if G.degree() == 3 and G.discriminant() != 0:
                out.append(("G" + "".join(str(i) for i in S), ell_from_cubic(G)))
            if G.degree() == 4 and G.discriminant() != 0:
                out.append(("G" + "".join(str(i) for i in S), ell_from_quartic(G)))
            TG = sq_free_scale(Tv * G)
            if TG.degree() == 4 and TG.discriminant() != 0:
                out.append(("TG" + "".join(str(i) for i in S), ell_from_quartic(TG)))
    return out

sub("§2a. сколько факторов реально получается (контроль ожидаемого числа 21)")
for tag in ('G1', 'G2'):
    names = [nm for nm, _ in factors_of(tag, 5, 2)]
    print(f"  {tag}: {len(names)} факторов: {names}")

sub("§2b. КОНТРОЛЬ якобиана по Ленгу: #C(F_p) обязано совпасть с #Jac(F_p), а не с твистом")
Rx = PolynomialRing(QQ, 'x'); xx = Rx.gen()
for tag in ('G1', 'G2'):
    for (m, n) in [(5, 2), (13, 8)]:
        cc = CELLS[tag](m, n); fr = FREE[tag]
        q = sq_free_scale(cc[fr[0]] * cc[fr[1]])
        E = ell_from_quartic(q); Etw = E.quadratic_twist(-1)
        ok, row = True, []
        for p in [11, 19, 23, 29, 31, 37, 43, 53]:
            if ZZ(E.discriminant()) % p == 0 or ZZ(q.discriminant()) % p == 0:
                continue
            Fp = GF(p); qp = q.change_ring(Fp)
            nC = 0
            for a in Fp:
                v = qp(a)
                nC += 1 if v == 0 else (2 if v.is_square() else 0)
            lead = Fp(q.coefficients(sparse=False)[4])
            nC += 2 if (lead != 0 and lead.is_square()) else (1 if lead == 0 else 0)
            nE = E.change_ring(Fp).cardinality()
            row.append((p, nC, nE, Etw.change_ring(Fp).cardinality()))
            ok = ok and (nC == nE)
        print(f"  {tag} ({m:2d},{n}) (p,#C,#Jac,#твист)={row[:4]} ... #C==#Jac всюду: {ok}")
        assert ok

sub("§2c. КОНТРОЛЬ на эталоне проекта: фактор c0*c8 семейства G1 = E_B, ранги 62/60/4/1")
cnt_EB = Counter(); mism = 0
for (m, n) in PAIRS:
    cc = cells_G1(m, n)
    E = ell_from_quartic(sq_free_scale(cc[0]*cc[8]))
    EB = EllipticCurve([0, -((m**2-n**2)**2 + (m**2+n**2)**2), 0,
                        (m**2-n**2)**2*(m**2+n**2)**2, 0])
    if not E.is_isomorphic(EB):
        mism += 1
    r = pari(EB).ellrank(); cnt_EB[(ZZ(r[0]), ZZ(r[1]))] += 1
print("  расхождений Jac(c0*c8) с E_B из REPORT_G1_SIX_KUMMER:", mism)
print("  распределение границ ранга E_B:", dict(cnt_EB))
print("  в RESULTS_FOR_CLAUDE_2026-09-12_G1.md §3 было {(0,0):62,(1,1):60,(2,2):4,(0,2):1}")
assert mism == 0
print("  кручение E_B по всем 127 парам:",
      dict(Counter(tuple(EllipticCurve([0, -((m**2-n**2)**2 + (m**2+n**2)**2), 0,
                                        (m**2-n**2)**2*(m**2+n**2)**2, 0]).torsion_subgroup().invariants())
                   for (m, n) in PAIRS)))

# ============================================================================================
hdr("§3. ГЛАВНЫЙ ПЕРЕБОР: есть ли у УГЛОВОГО семейства хоть один фактор ранга 0")
# ============================================================================================
def rank0_scan(tag, pairs):
    zero, bounds, tors, timeout_like = [], Counter(), Counter(), 0
    for (m, n) in pairs:
        for nm, E in factors_of(tag, m, n):
            tors[(nm, tuple(E.torsion_subgroup().invariants()))] += 1
            try:
                r = pari(E).ellrank()
                lo, hi = ZZ(r[0]), ZZ(r[1])
            except Exception:
                timeout_like += 1
                continue
            bounds[(lo, hi)] += 1
            if hi == 0:
                zero.append((m, n, nm))
    return zero, bounds, tors, timeout_like

for tag in ('G2', 'G1'):
    sub(f"§3{'a' if tag=='G2' else 'b'}. {tag}: все 21 фактор x {len(SCAN)} пар")
    t1 = time.time()
    zero, bounds, tors, fails = rank0_scan(tag, SCAN)
    print(f"  посчитано кривых: {sum(bounds.values())} за {time.time()-t1:.0f} c; сбоев ellrank: {fails}")
    print(f"  распределение границ ранга: {dict(sorted(bounds.items()))}")
    print(f"  ФАКТОРЫ РАНГА 0: {len(zero)}")
    if zero:
        by = Counter(nm for (_, _, nm) in zero)
        print(f"    по имени фактора: {dict(by)}")
        prs = sorted(set((m, n) for (m, n, _) in zero))
        print(f"    пары, у которых есть фактор ранга 0: {len(prs)} из {len(SCAN)}")
        print(f"    {prs}")
        for pr in (OPEN7 if tag == 'G1' else []):
            got = [nm for (m, n, nm) in zero if (m, n) == pr]
            print(f"    ОТКРЫТАЯ пара {pr}: факторы ранга 0 = {got}")
    tt = Counter(k[1] for k in tors.elements())
    print(f"  кручение по всем факторам: {dict(tt)}")
    globals()['ZERO_' + tag] = zero

# ============================================================================================
hdr("§4. ЧТО ДАЁТ РАНГ 0: перечисление точек и проверка, вырождены ли они")
# ============================================================================================
print("""Если у фактора E ранг 0, то E(Q) = кручение, кривая-фактор имеет конечное число
рациональных точек, а девятиклеточная кривая X отображается в неё с конечными слоями.
Значит X(Q) конечно И ПЕРЕЧИСЛИМО: надо поднять каждую точку кручения обратно в t.
Это и есть «дешёвая победа»: она НЕ требует ни ранга, ни Шаботи.""")

def t_values_from_rank0(tag, m, n, name, E):
    """вернуть множество рациональных t, совместимых с точкой кручения фактора"""
    cc = CELLS[tag](m, n)
    if name.startswith('TG') or name.startswith('G'):
        return None                 # подъём через T=t^2 — делаем только для пар-факторов
    i, j = [int(u) for u in name.replace('c', '').split()]
    q = sq_free_scale(cc[i]*cc[j])
    # все рациональные t с квадратным q(t) лежат на квартике; при rank(Jac)=0 их конечное число.
    return ('quartic', q)

sub("§4a. для найденных факторов ранга 0 — конкретный перечень и вырожденность")
found_any = False
for tag in ('G2', 'G1'):
    zero = globals()['ZERO_' + tag]
    if not zero:
        print(f"  {tag}: факторов ранга 0 нет — перечислять нечего.")
        continue
    found_any = True
    shown = 0
    for (m, n, nm) in zero:
        if shown >= 12:
            print("  ... (показаны первые 12)"); break
        cc = CELLS[tag](m, n)
        ts = [QQ(0), QQ(1), QQ(-1)]
        if m != n:
            ts += [QQ(m-n)/(m+n), -QQ(m-n)/(m+n), QQ(m+n)/(m-n), -QQ(m+n)/(m-n)]
        good = [t0 for t0 in ts if all(QQ(cc[i](t0)).is_square() for i in range(9))]
        ndist = [len(set(QQ(cc[i](t0)) for i in range(9))) for t0 in good]
        print(f"  {tag} ({m:2d},{n:2d}) фактор {nm}: известные рациональные t с девятью "
              f"квадратами {good}, различных значений {ndist}")
        shown += 1
if not found_any:
    print("  ни одного фактора ранга 0 — ни у углового, ни у рёберного семейства в этом наборе.")

# ============================================================================================
hdr("§5. СПОСОБ 1 СТРОГО: локальная разрешимость над Z_p, а не наивное сравнение mod p")
# ============================================================================================
print("""Наивный тест «значение клетки — ненулевой квадрат в F_p» НЕ является локальным
препятствием: рациональная точка может редуцироваться в точку, где клетка делится на p
(в EVEN_N_2ADIC_2026-09-12.md:396 это уже отмечалось как незакрытая оговорка).
Строгий тест: перебрать (a:b) в P^1(Z/p^N) и требовать, чтобы КАЖДАЯ клетка была квадратом
в Q_p; если валюация >= N — точка «неопределённая» и требует большего N.
Пусто при некотором N  =>  точек в Q_p нет  =>  препятствие настоящее.""")

def cells_hom(tag, m, n):
    """девять клеток как ЦЕЛЫЕ однородные формы степени 2 от (a,b), t = a/b;
       общий множитель — квадрат, поэтому квадратность каждой клетки не меняется."""
    cc = CELLS[tag](m, n)
    out = []
    for f in cc:
        c = f.coefficients(sparse=False) + [0, 0, 0]
        d = lcm([QQ(v).denominator() for v in c[:3]])
        # f(t)*b^2*d^2  = d^2*(c2 a^2 + c1 a b + c0 b^2)
        out.append(tuple(ZZ(d**2 * c[k]) for k in (2, 1, 0)))   # (A2 a^2 + A1 a b + A0 b^2)
    return out

def local_solvable(tag, m, n, p, N=5):
    """строгий перебор P^1(Z/p^N); возвращает (есть_определённо_точка, число_неопределённых)"""
    forms = cells_hom(tag, m, n)
    free = FREE[tag]
    pN = p**N
    QR = set((x*x) % p for x in range(1, p))
    reps = [(a, 1) for a in range(pN)] + [(1, p*b) for b in range(p**(N-1))]
    definite, undet = 0, 0
    for (a, b) in reps:
        ok, und = True, False
        for i in free:
            A2, A1, A0 = forms[i]
            v = (A2*a*a + A1*a*b + A0*b*b) % pN
            if v == 0:
                und = True; continue
            e = 0
            while v % p == 0:
                v //= p; e += 1
            if e % 2 == 1:
                ok = False; break
            if (v % p) not in QR:
                ok = False; break
        if ok and not und:
            definite += 1
        elif ok and und:
            undet += 1
    return definite, undet

sub("§5a. КОНТРОЛЬ: известные препятствия G1 обязаны устоять и при строгом тесте")
for (m, n, p) in [(13, 8, 17), (16, 5, 41), (13, 8, 47), (16, 5, 47), (11, 4, 23), (19, 16, 17)]:
    d3, u3 = local_solvable('G1', m, n, p, N=3)
    print(f"  G1 ({m:2d},{n:2d}) p={p:3d}, N=3: определённых точек {d3}, неопределённых {u3}"
          f"   {'ПУСТО СТРОГО — препятствие настоящее' if d3 == 0 and u3 == 0 else ''}")
d, u = local_solvable('G1', 15, 8, 17, N=3)
print(f"  контроль-наоборот: G1 (15,8) p=17: определённых {d}, неопределённых {u} (обязано быть > 0)")
assert d > 0

sub("§5b. G2: строгий тест на тех же p (ожидание — точка есть всегда, из-за латинской точки)")
for (m, n) in [(13, 8), (16, 5), (11, 4), (5, 2), (15, 8)]:
    row = []
    for p in [3, 5, 7, 17, 41, 47]:
        d, u = local_solvable('G2', m, n, p, N=3)
        row.append((p, d))
    print(f"  G2 ({m:2d},{n:2d}) (p, определённых точек) = {row}")

sub("§5c. полный наивный скан G2 (p < 200) — для сверки с раундом 1")
PR = [int(q) for q in prime_range(3, 200)]
def naive_survivors(tag, m, n, p):
    forms = cells_hom(tag, m, n); free = FREE[tag]
    QR = set((x*x) % p for x in range(1, p)); cnt = 0
    for (a, b) in [(x, 1) for x in range(p)] + [(1, 0)]:
        if all(((A2*a*a + A1*a*b + A0*b*b) % p) in QR for (A2, A1, A0) in [forms[i] for i in free]):
            cnt += 1
    return cnt
kill2 = []
for (m, n) in SCAN:
    for p in PR:
        if naive_survivors('G2', m, n, p) == 0:
            kill2.append((m, n, p)); break
print("  G2: пар, исключённых наивным сканом при p<200:", len(kill2), kill2[:10])
kill1 = []
for (m, n) in SCAN:
    for p in PR:
        if naive_survivors('G1', m, n, p) == 0:
            kill1.append((m, n, p)); break
print("  G1 (контроль): пар, исключённых наивным сканом при p<200:", len(kill1), "из", len(SCAN))
print("  G1 НЕ исключены:", sorted(set(SCAN) - set((m, n) for (m, n, _) in kill1)))

# ============================================================================================
hdr("§6. СПОСОБ 3: угловой класс Куммера и почему он ничего не исключает")
# ============================================================================================
def kummer_triple(tag, m, n):
    """средний столбец (G2) / главная диагональ (G1): три ЧЁТНЫЕ клетки -> кубика по T"""
    cc = CELLS[tag](m, n)
    ev = [i for i in FREE[tag] if cc[i](-t) == cc[i]]
    assert len(ev) == 3, ev
    return ev

def sqfree(v):
    v = QQ(v)
    if v == 0: return 0
    num = ZZ(v.numerator()*v.denominator())
    return sign(num)*prod(pp**(e % 2) for pp, e in num.abs().factor())

sub("§6a. общая формула класса и КОНТРОЛЬ: в рёберном случае обязано выйти (1,s,s)")
for (m, n) in [(11, 4), (19, 5), (15, 1), (13, 8), (16, 5)]:
    s = QQ(m**2 + n**2)/2
    al1, be1 = QQ(n**2), QQ(m**2)                       # G1: alpha,beta — КВАДРАТЫ
    al2, be2 = QQ(3*n**2 - m**2)/2, QQ(3*m**2 - n**2)/2  # G2
    d1 = (sqfree(al1*be1), sqfree(s*be1), sqfree(s*al1))
    d2 = (sqfree(al2*be2), sqfree(s*be2), sqfree(s*al2))
    print(f"  ({m:2d},{n:2d}) s={s}:  рёберный класс {d1}  |  угловой класс {d2}")
print("  рёберные значения совпадают с (1,s,s) из REPORT_G1_SIX_KUMMER (11,4)->(1,274,274) и т.д.")

sub("§6b. лежит ли угловой класс в образе: проверка через ЛАТИНСКУЮ точку, все 127 пар")
inimg, bad = 0, []
for (m, n) in PAIRS:
    s = QQ(m**2 + n**2)/2
    al, be = QQ(3*n**2 - m**2)/2, QQ(3*m**2 - n**2)/2
    r1, r2, r3 = s*al*be, s*be**2, s*al**2
    E = EllipticCurve([0, r1 + r2 + r3, 0, r1*r2 + r1*r3 + r2*r3, r1*r2*r3])
    tst = QQ(m - n)/(m + n)
    X = s*al*be*tst**2
    yy = (X + r1)*(X + r2)*(X + r3)
    if not QQ(yy).is_square():
        bad.append((m, n)); continue
    P = E(X, QQ(yy).sqrt())
    cls = (sqfree(X + r1), sqfree(X + r2), sqfree(X + r3))
    need = (sqfree(al*be), sqfree(s*be), sqfree(s*al))
    if cls == need:
        inimg += 1
    else:
        bad.append((m, n, cls, need))
print(f"  пар, где латинская точка лежит на E_c и даёт РОВНО требуемый класс: {inimg} из 127")
print(f"  расхождений: {bad[:5]}")
assert inimg == 127

# ============================================================================================
hdr("§7. СТРУКТУРА: род угловых кривых и разложение якобиана")
# ============================================================================================
sub("§7a. род кривых (проверено построением функционального поля над F_p, а не формулой)")
def genus_of(tag, m, n, S, p=10007):
    """род кривой {u_i^2 = c_i, i in S} над F_p — последовательные квадратичные расширения"""
    Fp = GF(p)
    K = FunctionField(Fp, 'x'); xg = K.gen()
    cur = K
    for k, i in enumerate(S):
        f = CELLS[tag](m, n)[i]
        co = [Fp(QQ(v)) for v in f.coefficients(sparse=False)]
        poly = sum(co[j]*xg**j for j in range(len(co)))
        Ry = PolynomialRing(cur, 'y%d' % k); yg = Ry.gen()
        cur = cur.extension(yg**2 - cur(poly), 'w%d' % k)
    return cur.genus()

for tag in ('G1', 'G2'):
    fr = FREE[tag]
    g2_ = genus_of(tag, 5, 2, fr[:2])
    g3_ = genus_of(tag, 5, 2, fr[:3])
    print(f"  {tag} (5,2): род при 2 условиях = {g2_} (ожидание 1); при 3 условиях = {g3_} (ожидание 5)")

sub("§7b. три ЧЁТНЫЕ свободные клетки дают кривую рода 5, распадающуюся полностью")
print("""  У обоих семейств тройка ЧЁТНЫХ свободных клеток — это
     G1: {c0, c4, c8} = главная диагональ a-e-i  (кривая C проекта, класс (1,s,s));
     G2: {c1, c4, c7} = средний СТОЛБЕЦ b-e-h.
  Кривая трёх условий имеет род 5 и накрывается (Z/2)^3; её якобиан распадается на
  якобианы семи промежуточных двойных накрытий:
     три квартики z^2 = c_i c_j        -> три эллиптические кривые;
     одна секстика z^2 = c_i c_j c_k   -> род 2, а так как произведение ЧЁТНО по t,
        она сама распадается: z^2 = G(T) и w^2 = T*G(T), обе рода 1.
  Итого Jac ~ E_ij x E_ik x E_jk x E_G x E_TG — ПЯТЬ эллиптических кривых, 1+1+1+1+1 = 5.
  Это и объясняет, почему «дешёвые» методы вообще работали: искать ранг 0 можно в пяти местах.""")
for tag in ('G1', 'G2'):
    for (m, n) in [(5, 2), (13, 8)]:
        cc = CELLS[tag](m, n); ev = kummer_triple(tag, m, n)
        Es = []
        for S in itertools.combinations(ev, 2):
            Es.append(ell_from_quartic(sq_free_scale(cc[S[0]]*cc[S[1]])))
        G = sq_free_scale(to_T(prod(cc[i] for i in ev)))
        Es.append(ell_from_cubic(G))
        Es.append(ell_from_quartic(sq_free_scale(Tv*G)))
        # контроль разложения: #X(F_p) против суммы следов Фробениуса
        okp = []
        for p in [101, 103, 107, 109, 113]:
            if any(ZZ(E.discriminant()) % p == 0 for E in Es):
                continue
            Fp = GF(p); QR = set((x*x) % p for x in range(1, p))
            nX = 0
            for a in range(p):
                vals = [ZZ(Fp(QQ(cc[i](a)))) for i in ev]
                if all(v != 0 for v in vals) and all(v in QR for v in vals):
                    nX += 8
            aps = sum(p + 1 - E.change_ring(Fp).cardinality() for E in Es)
            okp.append((p, nX, aps))
        print(f"  {tag} ({m:2d},{n:2d}) тройка {ev}: (p, #точек аффинной части X, сумма a_p пяти кривых) = {okp}")

sub("§7c. ранги пяти факторов этой тройки на открытых парах")
for tag in ('G1', 'G2'):
    for (m, n) in ([(13,8),(15,8),(16,5),(11,4),(19,5),(15,1),(19,16)] if tag == 'G1' else [(13,8),(15,8),(16,5),(5,2),(3,2)]):
        cc = CELLS[tag](m, n); ev = kummer_triple(tag, m, n)
        row = []
        for S in itertools.combinations(ev, 2):
            E = ell_from_quartic(sq_free_scale(cc[S[0]]*cc[S[1]]))
            r = pari(E).ellrank(); row.append(("c%dc%d" % S, (ZZ(r[0]), ZZ(r[1])),
                                               tuple(E.torsion_subgroup().invariants())))
        G = sq_free_scale(to_T(prod(cc[i] for i in ev)))
        for nm, E in [("G", ell_from_cubic(G)), ("TG", ell_from_quartic(sq_free_scale(Tv*G)))]:
            r = pari(E).ellrank(); row.append((nm, (ZZ(r[0]), ZZ(r[1])),
                                               tuple(E.torsion_subgroup().invariants())))
        print(f"  {tag} ({m:2d},{n:2d}): " + "; ".join(f"{a} ранг {b} кручение {c}" for a, b, c in row))

hdr("§8. КОНЕЦ")
print(f"время работы: {time.time()-T0:.0f} c")
