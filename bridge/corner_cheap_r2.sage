# -*- coding: utf-8 -*-
# corner_cheap_r2.sage — НЕЗАВИСИМАЯ ПЕРЕПРОВЕРКА и ДОСЧЁТ к bridge/corner_cheap.sage
#
# СТАТУС: расчёт
# ДЛЯ: Claude/Codex, тема «дешёвые победы в угловом случае»
# ИТОГ: три дешёвых способа в углах проверены заново с нуля; закрыты дыры прогона 22:50–23:01
#       (там посчитаны 2 произведения свободных клеток из 10, E_c только на 5 парах из 127,
#        а угловое СЕЧЕНИЕ на ранг 0 не проверялось вовсе)
# ОТМЕНЯЕТ: ничего
# ПРОВЕРЕНО: девять клеток G2 и четыре клетки углового сечения выведены заново из магичности
# ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужен ответ «закрываются ли углы дёшево» с полным перебором вариантов
#
# Запуск: sage bridge/corner_cheap_r2.sage
import functools, itertools, sys, time
from collections import Counter
print = functools.partial(print, flush=True)
T0 = time.time()

def hdr(s):
    print("\n" + "=" * 100); print(s); print("=" * 100)
def sub(s):
    print("\n" + "-" * 100); print(s); print("-" * 100)

PAIRS = [(m, n) for m in range(2, 21) for n in range(1, m) if gcd(m, n) == 1]
assert len(PAIRS) == 127, len(PAIRS)
LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]

def sqf(v):
    v = QQ(v); num = ZZ(v.numerator() * v.denominator())
    return sign(num) * prod(p^(e % 2) for p, e in num.abs().factor())

# ============================================================================================
hdr("§0. Девять клеток УГЛОВОГО семейства G2 выведены заново из аксиом магичности")
# ============================================================================================
R.<m, n, t> = QQ[]
# Разметка c0..c8 = a b c / d e f / g h i.  ПОЛНЫЕ (квадратные) пары — УГЛОВЫЕ (c0,c8),(c2,c6).
P, Q, Rr, S = m*t + n, m - n*t, m*t - n, m + n*t          # параметризация Сегре
assert P^2 + Q^2 == Rr^2 + S^2
c0, c8, c2, c6 = P^2, Q^2, Rr^2, S^2
c4 = (P^2 + Q^2) / 2
c1 = 3*c4 - c0 - c2          # только магичность, руками ничего не подставляется
c3 = 3*c4 - c0 - c6
c5 = 2*c4 - c3
c7 = 2*c4 - c1
cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]
ok = all(sum(cells[i] for i in L) == 3*c4 for L in LINES)
print("  все 8 линий равны 3*c4 тождественно над Q(m,n,t):", ok); assert ok
sS, al, be = (m^2 + n^2)/2, (3*n^2 - m^2)/2, (3*m^2 - n^2)/2
for name, d in [("c4 = s(1+t^2)", c4 - sS*(1 + t^2)),
                ("c1 = alpha t^2 + beta", c1 - (al*t^2 + be)),
                ("c7 = beta t^2 + alpha", c7 - (be*t^2 + al)),
                ("c3 = c4 - 4mn t", c3 - (sS*(1 + t^2) - 4*m*n*t)),
                ("c5 = c4 + 4mn t", c5 - (sS*(1 + t^2) + 4*m*n*t))]:
    print(f"  {name:24s} разность = {d}"); assert d == 0
print("  alpha + beta = m^2+n^2 = 2s:", (al + be - 2*sS) == 0)
print("  автоматические квадраты — ЧЕТЫРЕ УГЛА c0,c2,c6,c8; свободен КРЕСТ c1,c3,c4,c5,c7.")
print("  [доказано символьно] вывод не опирается на файлы проекта.")

# --- клетки в масштабе x2 (без знаменателей); квадратность ПРОИЗВЕДЕНИЙ от этого не меняется
def coeffs_G2(mm, nn):
    """девять клеток G2 как тройки (a2,a1,a0) для a2 t^2 + a1 t + a0, всё умножено на 2"""
    Sq, A1, B1 = mm^2 + nn^2, 3*nn^2 - mm^2, 3*mm^2 - nn^2
    return [(2*mm^2, 4*mm*nn, 2*nn^2), (A1, 0, B1), (2*mm^2, -4*mm*nn, 2*nn^2),
            (Sq, -8*mm*nn, Sq),        (Sq, 0, Sq),  (Sq, 8*mm*nn, Sq),
            (2*nn^2, 4*mm*nn, 2*mm^2), (B1, 0, A1),  (2*nn^2, -4*mm*nn, 2*mm^2)]
def coeffs_G1(mm, nn):
    """девять клеток РЁБЕРНОГО G1 в том же масштабе x2"""
    Sq = mm^2 + nn^2
    return [(2*nn^2, 0, 2*mm^2), (2*mm^2, 4*mm*nn, 2*nn^2), (Sq, -4*mm*nn, Sq),
            (2*mm^2, -4*mm*nn, 2*nn^2), (Sq, 0, Sq), (2*nn^2, 4*mm*nn, 2*mm^2),
            (Sq, 4*mm*nn, Sq), (2*nn^2, -4*mm*nn, 2*mm^2), (2*mm^2, 0, 2*nn^2)]
def val(co, t0):
    return [a2*t0^2 + a1*t0 + a0 for (a2, a1, a0) in co]
for fn in (coeffs_G1, coeffs_G2):
    for (mm, nn) in [(3,2),(5,2),(13,8),(16,5),(19,5)]:
        for t0 in [QQ(3)/7, QQ(-2)/5, QQ(11)/4]:
            cc = val(fn(mm, nn), t0)
            assert all(sum(cc[i] for i in L) == 3*cc[4] for L in LINES), (fn, mm, nn, t0)
print("  контроль магичности численно (оба семейства, 5 пар x 3 значения t, 8 линий): True")

# ============================================================================================
hdr("§1. ЛАТИНСКАЯ ТОЧКА: у КАЖДОГО углового слоя есть точка с девятью квадратами")
# ============================================================================================
sub("§1a. символьно над Q(m,n)")
ts = (m - n)/(m + n)
sqroots = {1: m^2 + 2*m*n - n^2, 3: m^2 - 2*m*n - n^2, 4: m^2 + n^2,
           5: m^2 + 2*m*n - n^2, 7: m^2 - 2*m*n - n^2}
cpol = {1: al*ts^2 + be, 3: sS*(1 + ts^2) - 4*m*n*ts, 4: sS*(1 + ts^2),
        5: sS*(1 + ts^2) + 4*m*n*ts, 7: be*ts^2 + al}
for k in sorted(cpol):
    d = (cpol[k]*(m + n)^2 - sqroots[k]^2).numerator()
    print(f"  c{k}*(m+n)^2 - ({sqroots[k]})^2 = {d}"); assert d == 0
print("  углы: c0*(m+n)^2 = ((m t*+n)(m+n))^2 — квадрат тождественно.")
print("  [доказано символьно] при t* = (m-n)/(m+n) все девять клеток — квадраты, но")
print("  различных значений ТРИ: латинский вырожденный квадрат, а не решение задачи.")

sub("§1b. те же четыре t численно на всех 127 парах (перебор)")
bad, dist = [], set()
for (mm, nn) in PAIRS:
    co = coeffs_G2(mm, nn)
    for t0 in [QQ(mm-nn)/(mm+nn), -QQ(mm-nn)/(mm+nn), QQ(mm+nn)/(mm-nn), -QQ(mm+nn)/(mm-nn)]:
        cc = val(co, t0)
        if not all(QQ(v).is_square() for v in cc): bad.append((mm, nn, t0))
        dist.add(len(set(cc)))
print("  пар/точек без девяти квадратов:", len(bad), bad[:5], "; различных значений:", sorted(dist))
assert bad == [] and dist == {3}

sub("§1c. КОНТРАСТ: у рёберного G1 такой точки почти никогда нет")
g1_ok = sorted((mm, nn) for (mm, nn) in PAIRS
               if any(all(QQ(v).is_square() for v in val(coeffs_G1(mm, nn), t0)) for t0 in [0, 1, -1]))
crit = sorted((mm, nn) for (mm, nn) in PAIRS
              if ZZ(mm^2 + nn^2).is_square() or QQ((mm^2 + nn^2)/2).is_square())
print("  G1 даёт девять квадратов в t in {0,1,-1} у пар:", g1_ok)
print("  критерий m^2+n^2=кв. или (m^2+n^2)/2=кв. даёт:", crit, "; совпало:", g1_ok == crit)
print(f"  ИТОГ: G2 — 127/127, G1 — {len(g1_ok)}/127. Здесь и сидит вся разница дешёвых методов.")

# ============================================================================================
hdr("§2. СПОСОБ 1 — ЛОКАЛЬНЫЕ ПРЕПЯТСТВИЯ: перебор P^1(F_p), нечётные p < 200")
# ============================================================================================
PRIMES = list(primes(3, 200))
SQS = {p: set(ZZ((GF(p)(i))^2) for i in range(1, p)) for p in PRIMES}
def survivors(co, p, sq=None):
    sq = sq if sq is not None else SQS[p]
    out = []
    for t0 in range(p):
        good = True
        for (a2, a1, a0) in co:
            v = (a2*t0*t0 + a1*t0 + a0) % p
            if v == 0 or v not in sq: good = False; break
        if good: out.append(t0)
    if all((a2 % p) != 0 and (a2 % p) in sq for (a2, a1, a0) in co): out.append('oo')
    return out

sub("§2a. КОНТРОЛЬ кода на РЁБЕРНОМ случае: (13,8) mod 17 и (16,5) mod 41 обязаны быть пусты")
for (mm, nn, p) in [(13, 8, 17), (16, 5, 41), (15, 8, 17), (15, 8, 41)]:
    s = survivors(coeffs_G1(mm, nn), p)
    print(f"  G1 ({mm:2d},{nn}) p={p:3d}: выживших = {len(s):2d}  "
          f"{'ПУСТО — известное препятствие воспроизведено' if not s else ''}")
assert survivors(coeffs_G1(13, 8), 17) == [] and survivors(coeffs_G1(16, 5), 41) == []

sub("§2b. полный скан обоих семейств: 127 пар x 45 простых")
for tag, fn in [("G1 (рёберное)", coeffs_G1), ("G2 (УГЛОВОЕ)", coeffs_G2)]:
    killed, firstp, mins = [], Counter(), []
    for (mm, nn) in PAIRS:
        co = fn(mm, nn); f = None; mi = 10^9
        for p in PRIMES:
            s = survivors(co, p); mi = min(mi, len(s))
            if not s and f is None: f = p
        mins.append(mi)
        if f: killed.append((mm, nn)); firstp[f] += 1
    print(f"  {tag}: локально исключено пар из 127: {len(killed)}; "
          f"минимум выживших по всем (m,n,p): {min(mins)}")
    if killed: print("     первое пустое p:", sorted(firstp.items()))
    if tag.startswith("G1"):
        left = sorted(set(PAIRS) - set(killed))
        print("     локально НЕ исключены:", left, "; совпадает со списком §1c:", left == g1_ok)
    else:
        assert not killed, killed
        print("     [перебор] НИ ОДНА угловая пара не исключается локально.")

sub("§2c. и дальше сканировать бессмысленно: редукция латинской точки выживает всегда")
cnt = 0
for (mm, nn) in PAIRS:
    co = coeffs_G2(mm, nn); t0 = QQ(mm - nn)/(mm + nn)
    for p in PRIMES:
        if (mm + nn) % p == 0: continue
        vals = [ZZ(GF(p)(v)) for v in val(co, t0)]
        if any(v == 0 for v in vals): continue
        assert all(v in SQS[p] for v in vals), (mm, nn, p)
        cnt += 1
print(f"  проверено пар-простых с ненулевой редукцией: {cnt}; девять квадратов в F_p всюду: True")
print("  [доказано] у G2 локального препятствия не бывает НИ ПРИ КАКОМ p: кривая с рациональной")
print("  точкой всюду локально разрешима. СПОСОБ 1 в углах мёртв структурно, а не «до 200».")

sub("§2d. НОВОЕ: сколько в F_p ЛИШНИХ точек сверх вырожденных (шансы решета)")
for (mm, nn) in [(3, 2), (5, 2), (13, 8), (16, 5), (19, 5)]:
    co = coeffs_G2(mm, nn); tst = QQ(mm - nn)/(mm + nn)
    Dg = [QQ(0), QQ(1), QQ(-1), tst, -tst, 1/tst, -1/tst]
    row = []
    for p in [101, 199, 401, 701, 997]:
        sq = set(ZZ((GF(p)(i))^2) for i in range(1, p))
        s = set(survivors(co, p, sq))
        d = set(ZZ(GF(p)(x)) for x in Dg) | {'oo'}
        row.append((p, len(s), len(s & d), len(s - d)))
    print(f"  ({mm:2d},{nn}): (p, |S_p|, вырожденных, ЛИШНИХ) = {row}")
print("  [проверено численно] лишние точки есть при каждом p => решето по одному p ничего")
print("  не даёт; нужен настоящий MW-sieve с группой, а это уже не «дёшево».")

# ============================================================================================
hdr("§3. СПОСОБ 2 — ФАКТОР РАНГА 0: ВСЕ 10 произведений свободных клеток (раньше считали 2)")
# ============================================================================================
Rx.<x> = QQ[]
def jac_quartic(q):
    """якобиан z^2=q(x), deg q = 4:  y^2 = x^3 - 27 I x - 27 J"""
    co = q.coefficients(sparse=False); co = co + [0]*(5 - len(co))
    a, b, c, d, e = co[4], co[3], co[2], co[1], co[0]
    I = 12*a*e - 3*b*d + c^2
    J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return EllipticCurve([0, 0, 0, -27*I, -27*J])

sub("§3a. КОНТРОЛЬ формулы якобиана на рёберной квартике z^2=F0*F8 (обязан выйти E_B)")
for (mm, nn) in [(3, 2), (13, 8), (16, 5), (19, 5), (20, 19)]:
    q = (mm^2 + nn^2*x^2)*(nn^2 + mm^2*x^2)
    EB = EllipticCurve([0, -((mm^2-nn^2)^2 + (mm^2+nn^2)^2), 0, (mm^2-nn^2)^2*(mm^2+nn^2)^2, 0])
    same = jac_quartic(q).is_isomorphic(EB)
    print(f"  ({mm:2d},{nn}) Jac(F0F8) ~= E_B: {same}; кручение E_B = {EB.torsion_subgroup().invariants()}")
    assert same

sub("§3b. КОНТРОЛЬ знака (Ленг): #C(F_p) обязано равняться #Jac(F_p), а не #твиста")
def quartic_of(mm, nn, pr):
    co = coeffs_G2(mm, nn)
    f = lambda k: co[k][0]*x^2 + co[k][1]*x + co[k][2]
    return f(pr[0])*f(pr[1])
for (mm, nn) in [(2, 1), (5, 2), (13, 8)]:
    q = quartic_of(mm, nn, (3, 5)); E = jac_quartic(q); Et = E.quadratic_twist(-1)
    row, okk = [], True
    for p in [7, 11, 19, 23, 29, 31, 37, 43]:
        if ZZ(E.discriminant()) % p == 0 or ZZ(q.discriminant()) % p == 0: continue
        qp = q.change_ring(GF(p)); lead = GF(p)(q.coefficients(sparse=False)[4])
        nC = sum(0 if qp(t0) == 0 else (2 if qp(t0).is_square() else 0) for t0 in GF(p)) \
             + sum(1 for t0 in GF(p) if qp(t0) == 0) \
             + (2 if lead != 0 and lead.is_square() else (1 if lead == 0 else 0))
        nE = E.change_ring(GF(p)).cardinality(); nT = Et.change_ring(GF(p)).cardinality()
        row.append((p, nC, nE, nT)); okk = okk and (nC == nE)
    print(f"  ({mm:2d},{nn}) (p,#C,#Jac,#твист) = {row} -> #C == #Jac всюду: {okk}")
    assert okk

sub("§3c. ВСЕ 10 произведений: кручение, нижняя оценка числа точек, есть ли ранг 0")
FREE = [1, 3, 4, 5, 7]; PRODS = list(itertools.combinations(FREE, 2))
print("  произведения c_i*c_j, i<j из креста:", PRODS)
tor_cnt = {pr: Counter() for pr in PRODS}; low = {pr: Counter() for pr in PRODS}
cert = Counter(); rank0 = []
for (mm, nn) in PAIRS:
    tst = QQ(mm - nn)/(mm + nn)
    Dg = [QQ(0), QQ(1), QQ(-1), tst, -tst, 1/tst, -1/tst]
    for pr in PRODS:
        q = quartic_of(mm, nn, pr)
        assert q.degree() == 4 and q.discriminant() != 0, (mm, nn, pr)
        E = jac_quartic(q); tor = tuple(E.torsion_subgroup().invariants())
        tor_cnt[pr][tor] += 1
        nt = len(set(t0 for t0 in Dg if q(t0) != 0 and QQ(q(t0)).is_square()))
        lead = q.coefficients(sparse=False)[4]
        npts = 2*nt + (2 if QQ(lead).is_square() else 0)
        low[pr][npts] += 1
        if npts > E.torsion_order(): cert['rank>=1 доказан счётом точек'] += 1
        else:
            cert['счётом НЕ доказан -> ellrank'] += 1
            r = pari(E).ellrank()
            if ZZ(r[1]) == 0: rank0.append(((mm, nn), pr))
for pr in PRODS:
    print(f"  c{pr[0]}*c{pr[1]}: кручение {dict(tor_cnt[pr])}; нижняя оценка #точек {dict(low[pr])}")
print("  сертификаты:", dict(cert))
print("  ПАРЫ С ФАКТОРОМ РАНГА 0 среди 127 x 10 = 1270 квартик:", rank0)

sub("§3d. независимый контроль: PARI ellrank на всех 1270 квартиках")
t1 = time.time(); rb = {pr: Counter() for pr in PRODS}; zero = []
for (mm, nn) in PAIRS:
    for pr in PRODS:
        r = pari(jac_quartic(quartic_of(mm, nn, pr))).ellrank()
        rb[pr][(ZZ(r[0]), ZZ(r[1]))] += 1
        if ZZ(r[1]) == 0: zero.append(((mm, nn), pr))
for pr in PRODS: print(f"  c{pr[0]}*c{pr[1]}: границы ранга {dict(rb[pr])}")
print(f"  ранг 0 встретился: {zero}   ({time.time()-t1:.0f} c)")
cb = Counter()
for (mm, nn) in PAIRS:
    EB = EllipticCurve([0, -((mm^2-nn^2)^2 + (mm^2+nn^2)^2), 0, (mm^2-nn^2)^2*(mm^2+nn^2)^2, 0])
    r = pari(EB).ellrank(); cb[(ZZ(r[0]), ZZ(r[1]))] += 1
print("  КОНТРОЛЬ рёберного эталона E_B:", dict(cb), "— в RESULTS_FOR_CLAUDE §3 было 62/60/4/1")

# ============================================================================================
hdr("§4. СПОСОБ 3 — КЛАСС КУММЕРА и четвёртый фактор E_c (средний столбец, T=t^2), 127 пар")
# ============================================================================================
sub("§4a. вывод E_c заново, символьно")
A1s, B1s, Ss = 3*n^2 - m^2, 3*m^2 - n^2, m^2 + n^2
X = A1s*B1s*Ss*t^2
for name, lhs, rhs in [("X + S*A1*B1 = A1*B1*C4", X + Ss*A1s*B1s, A1s*B1s*(Ss*(t^2 + 1))),
                       ("X + S*B1^2  = S*B1*C1 ", X + Ss*B1s^2,    Ss*B1s*(A1s*t^2 + B1s)),
                       ("X + S*A1^2  = S*A1*C7 ", X + Ss*A1s^2,    Ss*A1s*(B1s*t^2 + A1s))]:
    print(f"  {name}: разность = {lhs - rhs}"); assert lhs - rhs == 0
print("  => E_c: W^2 = (X+S A1 B1)(X+S B1^2)(X+S A1^2), X = S A1 B1 t^2;")
print("     необходимый УГЛОВОЙ класс Куммера delta_c = (A1*B1, S*B1, S*A1).")

sub("§4b. КОНТРОЛЬ: та же общая формула на РЁБЕРНЫХ параметрах обязана дать (1,s,s)")
for (mm, nn) in [(11, 4), (19, 5), (15, 1), (13, 8), (16, 5)]:
    a_e, b_e, s_e = nn^2, mm^2, QQ(mm^2 + nn^2)/2
    print(f"  ({mm:2d},{nn}): (alpha*beta, s*beta, s*alpha) = "
          f"({sqf(a_e*b_e)}, {sqf(s_e*b_e)}, {sqf(s_e*a_e)})   [у Кодекса (1, {sqf(2*s_e)}, {sqf(2*s_e)})]")

sub("§4c. ранг и кручение E_c на ВСЕХ 127 парах; торсионна ли латинская точка")
cE, cT, zero_c, lat = Counter(), Counter(), [], Counter()
for (mm, nn) in PAIRS:
    A1, B1, Sq = 3*nn^2 - mm^2, 3*mm^2 - nn^2, mm^2 + nn^2
    r1, r2, r3 = Sq*A1*B1, Sq*B1^2, Sq*A1^2
    E = EllipticCurve([0, r1 + r2 + r3, 0, r1*r2 + r1*r3 + r2*r3, r1*r2*r3])
    cT[tuple(E.torsion_subgroup().invariants())] += 1
    r = pari(E).ellrank(); cE[(ZZ(r[0]), ZZ(r[1]))] += 1
    if ZZ(r[1]) == 0: zero_c.append((mm, nn))
    Xl = A1*B1*Sq*(QQ(mm - nn)/(mm + nn))^2
    yy = (Xl + r1)*(Xl + r2)*(Xl + r3)
    assert QQ(yy).is_square(), (mm, nn)
    Pl = E(Xl, QQ(yy).sqrt())
    lat['кручение' if Pl.order() != +Infinity else 'бесконечный порядок'] += 1
print("  кручение E_c:", dict(cT))
print("  границы ранга E_c:", dict(cE))
print("  пары с rank E_c = 0:", zero_c)
print("  латинская точка на E_c:", dict(lat))
print("  [перебор] нужный класс delta_c реализуется латинской точкой => СПОСОБ 3 в углах")
print("  не исключает ничего: класс всегда в образе.")

# ============================================================================================
hdr("§5. УГЛОВОЕ СЕЧЕНИЕ: тот же конвейер, что рёберный, но A=h0^2, C=b0^2, beta=4n0^2")
# ============================================================================================
SECTIONS = [(17,7,13),(23,7,17),(71,49,61),(7,1,5),(31,17,25),(41,1,29),
            (47,23,37),(49,31,41),(73,17,53),(89,23,65),(79,47,65)]
sub("§5a. четыре свободные клетки углового сечения выведены заново (символьно)")
Rs.<b0, h0, n0, u> = QQ[]
Rc, Sc, den = n0*(1 - 2*u - u^2), n0*(1 + 2*u - u^2), (1 + u^2)
print("  R^2+S^2-2n0^2(1+u^2)^2 =", Rc^2 + Sc^2 - 2*n0^2*den^2)
a_, i_, e_ = b0^2*den^2, h0^2*den^2, n0^2*den^2
c_, g_ = Rc^2, Sc^2
b_ = 3*e_ - a_ - c_; d_ = 3*e_ - a_ - g_; f_ = 2*e_ - d_; h_ = 2*e_ - b_
sq9 = [a_, b_, c_, d_, e_, f_, g_, h_, i_]
Iq = Rs.ideal(b0^2 + h0^2 - 2*n0^2)
okl = all(Iq.reduce(sum(sq9[j] for j in L) - 3*e_) == 0 for L in LINES)
print("  все 8 линий = 3*центр по модулю b0^2+h0^2=2n0^2:", okl); assert okl
W = {"C*P+B*Q": b0^2*den^2 + 4*n0^2*u*(u^2 - 1), "A*P+B*Q": h0^2*den^2 + 4*n0^2*u*(u^2 - 1),
     "C*P-B*Q": b0^2*den^2 - 4*n0^2*u*(u^2 - 1), "A*P-B*Q": h0^2*den^2 - 4*n0^2*u*(u^2 - 1)}
got = {}
for nm, expr in zip(["b", "d", "f", "h"], [b_, d_, f_, h_]):
    for k2, w in W.items():
        if Iq.reduce(expr - w) == 0: got[nm] = k2
print("  свободные (РЁБЕРНЫЕ) клетки <-> формы:", got); assert len(got) == 4
print("  => УГЛОВОЕ сечение = семейство family/pipeline.sage с A = h0^2, C = b0^2, beta = 4n0^2.")
print("     рёберное: A=(h0^2+n0^2)/2, C=(b0^2+n0^2)/2, beta=2n0^2 (там beta=A+C, здесь beta=2(A+C)).")
print("     Ядро конвейера (Q, Q', C_J2) зависит только от (A:C) — beta входит лишь в возврат к t.")

sub("§5b. инварианты: b' = (C-A)/(C+A) и поле k = Q(sqrt(C^2-A^2))")
print(f"  {'сечение':>16} | {'kk':>6} | {'b_рёб':>12} | {'b_угл':>12} | {'k рёберное':>12} | {'k УГЛОВОЕ':>12}")
for (bb, hh, nnn) in SECTIONS:
    assert bb^2 + hh^2 == 2*nnn^2
    kk = ZZ((bb^2 - hh^2)/2)
    Ae, Ce = QQ(hh^2 + nnn^2)/2, QQ(bb^2 + nnn^2)/2
    Ac, Cc = QQ(hh^2), QQ(bb^2)
    print(f"  ({bb:3d},{hh:2d},{nnn:2d})   | {kk:6d} | {str((Ce-Ae)/(Ce+Ae)):>12} | "
          f"{str((Cc-Ac)/(Cc+Ac)):>12} | Q(v{sqf(Ce^2-Ae^2):5d}) | Q(v{sqf(Cc^2-Ac^2):5d})")
print("  [доказано] b'_угл = 2*b'_рёб;  k_угл = Q(sqrt(kk)),  k_рёб = Q(sqrt(2kk)), kk=(b0^2-h0^2)/2.")

sub("§5c. ранг 0 у квартик УГЛОВОГО сечения? (в corner_cheap.log это не считалось)")
zeros_sec = []
print(f"  {'сечение':>16} | {'Jac(A-квартика)':>28} | {'Jac(C-квартика)':>28} | изв.точек")
for (bb, hh, nnn) in SECTIONS:
    out = []
    for A_ in (hh^2, bb^2):
        q = A_*(x^2 + 1)^2 + 4*nnn^2*x*(x^2 - 1)
        E = jac_quartic(q); tor = tuple(E.torsion_subgroup().invariants())
        r = pari(E).ellrank()
        nt = len(set(t0 for t0 in [QQ(0), QQ(1), QQ(-1)] if q(t0) != 0 and QQ(q(t0)).is_square()))
        npts = 2*nt + (2 if QQ(q.coefficients(sparse=False)[4]).is_square() else 0)
        out.append((tor, (ZZ(r[0]), ZZ(r[1])), npts, E.torsion_order()))
        if ZZ(r[1]) == 0: zeros_sec.append((bb, hh, nnn, A_))
    print(f"  ({bb:3d},{hh:2d},{nnn:2d})   | tors {str(out[0][0]):>10} ранг {str(out[0][1]):>6} | "
          f"tors {str(out[1][0]):>10} ранг {str(out[1][1]):>6} | {out[0][2]}/{out[1][2]} "
          f"(|tors| {out[0][3]}/{out[1][3]})")
print("  сечения с фактором ранга 0:", zeros_sec)

sub("§5d. КОНТРОЛЬ: те же квартики в РЁБЕРНОЙ позиции — там точки t=0 нет")
for (bb, hh, nnn) in SECTIONS:
    Ae, Ce = QQ(hh^2 + nnn^2)/2, QQ(bb^2 + nnn^2)/2
    qs = []
    for A_ in (Ae, Ce):
        q = A_*(x^2 + 1)^2 + 4*nnn^2*x*(x^2 - 1)      # beta углового не важна для наличия t=0
        qs.append(A_.is_square())
    print(f"  ({bb:3d},{hh:2d},{nnn:2d}) рёберное: A={Ae} кв.{qs[0]}, C={Ce} кв.{qs[1]} "
          f"|| угловое: A={hh^2} кв.True, C={bb^2} кв.True")

sub("§5e. НОВОЕ: угловое сечение как задача о конгруэнтных числах")
print("  b=S^2-kk, h=R^2+kk, d=R^2-kk, f=S^2+kk при R^2+S^2=2n0^2, kk=(b0^2-h0^2)/2:")
print("  нужны ДВА x=R,S с квадратными x^2 +- kk и R^2+S^2=2n0^2. Вырожденное решение R=S=n0.")
for (bb, hh, nnn) in SECTIONS:
    kk = ZZ((bb^2 - hh^2)/2); ks = ZZ(sqf(kk))
    Ek = EllipticCurve([0, 0, 0, -ks^2, 0]); r = pari(Ek).ellrank()
    assert nnn^2 - kk == hh^2 and nnn^2 + kk == bb^2
    print(f"  ({bb:3d},{hh:2d},{nnn:2d}): kk={kk:6d}, бескв.часть {ks:5d}, конгруэнтная кривая "
          f"y^2=x^3-{ks}^2x ранг {(ZZ(r[0]), ZZ(r[1]))}; n0^2-kk={hh}^2, n0^2+kk={bb}^2")
print("  [доказано] kk конгруэнтно ВСЕГДА: h0^2,n0^2,b0^2 — прогрессия квадратов с разностью kk.")
print("  => дешёвого препятствия нет и здесь: ранг конгруэнтной кривой >= 1 по построению.")

sub("§5f. локальный скан сечений (p<200): угловое против рёберного")
for (bb, hh, nnn) in SECTIONS:
    mi = {0: 10^9, 1: 10^9}
    for p in PRIMES:
        F = GF(p)
        if (bb*hh*nnn*(bb^2 - hh^2)) % p == 0: continue
        for tag in (0, 1):
            if tag == 0: A_, C_, B_ = F(hh^2), F(bb^2), F(4*nnn^2)
            else:        A_, C_, B_ = F(hh^2 + nnn^2)/2, F(bb^2 + nnn^2)/2, F(2*nnn^2)
            c2 = 0
            for t0 in F:
                Pp, Qq = (t0^2 + 1)^2, t0*(t0^2 - 1)
                vv = [A_*Pp + B_*Qq, A_*Pp - B_*Qq, C_*Pp + B_*Qq, C_*Pp - B_*Qq]
                if all(v != 0 and v.is_square() for v in vv): c2 += 1
            mi[tag] = min(mi[tag], c2)
    print(f"  ({bb:3d},{hh:2d},{nnn:2d}): min выживших t при p<200 — угловое {mi[0]}, рёберное {mi[1]}")

hdr("§6. КОНЕЦ")
print(f"время работы: {time.time()-T0:.0f} c")
