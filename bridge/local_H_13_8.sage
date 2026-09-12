#!/usr/bin/env sage
# -*- coding: utf-8 -*-
#
# local_H_13_8.sage
#
# Локальная разрешимость кривой рода 2
#     H : Y^2 = F0(t) * F4(t) * F8(t)
# для семейства G1 с (m,n) = (13,8).
#
#     s  = (m^2+n^2)/2
#     F0 = m^2 + n^2 t^2
#     F4 = s (1 + t^2)
#     F8 = n^2 + m^2 t^2
#
# Проверяется H(Q_v) != пусто для v = бесконечность (R) и всех простых p.
#
# Два НЕЗАВИСИМЫХ решающих алгоритма для Q_p:
#   (A) рекурсия по дискам P^1(Q_p) с ярлыком простого корня (solve/…);
#   (B) полный перебор дисков P^1(Z/p^N) с гарантированной определённостью
#       класса квадратов (decide_bruteforce).
# Алгоритм (A) дополнительно сверен на 6000+ случайных квартиках с
# sage.schemes.elliptic_curves.descent_two_isogeny.test_qpls (независимая
# C-реализация Qp_soluble из 2-спуска Sage).
# Ответы ДА дополнительно подтверждаются ТОЧНЫМИ целочисленными
# сертификатами (валюация + класс единицы), без p-адической точности.
#
# Запуск:  sage /home/kep/magicKube/bridge/local_H_13_8.sage

import sys

R.<t> = PolynomialRing(QQ)
Rz.<z> = PolynomialRing(ZZ)

# ======================================================================
# 0. Модель
# ======================================================================

m, n = 13, 8
assert gcd(m, n) == 1

s  = QQ(m**2 + n**2) / 2
F0 = m**2 + n**2 * t**2
F4 = s * (1 + t**2)
F8 = n**2 + m**2 * t**2

Fprod = F0 * F4 * F8
# домножаем на 4 = квадрат  =>  класс квадратов в каждом Q_v и в R не меняется
f = R(4 * Fprod)
assert all(c in ZZ for c in f.list())
f = Rz([ZZ(c) for c in f.list()])

print("=" * 78)
print("H : Y^2 = F0*F4*F8,   (m,n) = (%d,%d)" % (m, n))
print("=" * 78)
print("s  =", s)
print("F0 =", F0)
print("F4 =", F4)
print("F8 =", F8)
print("F0*F4*F8 =", Fprod)
print("рабочая модель (домножена на квадрат 4):")
print("  f(t) = 4*F0*F4*F8 =", f)
print("  factor(f) =", factor(f))

print("\n--- сверка модели подстановкой конкретных t ---")
ok_model = True
for tv in [QQ(0), QQ(1), QQ(2), QQ(-3), QQ(5)/QQ(3), QQ(-7)/QQ(4), QQ(100)/QQ(7)]:
    lhs = f(tv)
    rhs = 4 * (m**2 + n**2*tv**2) * (s*(1+tv**2)) * (n**2 + m**2*tv**2)
    same = (lhs == rhs)
    ok_model = ok_model and same
    print("  t=%-9s f(t)=%-26s 4*F0*F4*F8=%-26s %s"
          % (tv, lhs, rhs, "OK" if same else "РАСХОЖДЕНИЕ"))
assert ok_model

b_cod = s * m**2 * n**2
a_cod = QQ(m**2)/QQ(n**2) + 1 + QQ(n**2)/QQ(m**2)
codex = b_cod * (t**6 + a_cod*t**4 + a_cod*t**2 + 1)
print("  модель Codex  b*(t^6+a t^4+a t^2+1) == F0*F4*F8 :", bool(R(codex) == R(Fprod)))
print("  b = s*m^2*n^2 =", b_cod, "  класс квадратов b =", ZZ(4*b_cod).squarefree_part(),
      " (Codex: 466)")
assert R(codex) == R(Fprod)

# ======================================================================
# 1. Дискриминант, плохие простые
# ======================================================================
print("\n" + "=" * 78)
print("1. Дискриминант и плохие простые")
print("=" * 78)
assert f.degree() == 6
print("deg f =", f.degree())
print("f бесквадратен:", f.is_squarefree(), " => гладкая модель, род 2")
assert f.is_squarefree()
disc_f = f.discriminant()
print("disc(f) =", factor(disc_f))
lead = ZZ(f.leading_coefficient())
print("старший коэффициент =", lead, "=", factor(lead),
      "  свободная от квадратов часть =", lead.squarefree_part())
bad = sorted(set([p for p, _ in factor(2 * disc_f)]))
print("плохие простые (делители 2*disc(f)):", bad)

# ======================================================================
# 2. Вещественное место
# ======================================================================
print("\n" + "=" * 78)
print("2. Вещественное место")
print("=" * 78)
print("F0 = m^2 + n^2 t^2 >= %d > 0  для всех t in R" % m**2)
print("F8 = n^2 + m^2 t^2 >= %d > 0  для всех t in R" % n**2)
print("F4 = s(1+t^2)      >= %s > 0  для всех t in R" % s)
print("=> f(t) = 4*F0*F4*F8 > 0 всюду на R (произведение трёх положительных).")
rr = f.roots(RR)
print("число вещественных корней f:", len(rr), "(ожидается 0)")
print("f(0) = %s > 0 : %s" % (f(0), bool(f(0) > 0)))
print("min f на сетке t=k/10, |k|<=500 :", min([f(QQ(k)/QQ(10)) for k in range(-500, 501)]))
real_ok = (len(rr) == 0 and f(0) > 0)
print("ВЫВОД: H(R) != пусто, напр. t=0, Y=sqrt(%s)=%.6f : %s"
      % (f(0), sqrt(RR(f(0))), real_ok))
assert real_ok

# ======================================================================
# 3. Алгоритм (A): рекурсия по дискам P^1(Q_p)
# ======================================================================
#
# Ищем (u:v) in P^1(Q_p) с G(u,v) = v^6 f(u/v) in (Q_p^*)^2 u {0}.
# Гладкая модель y^2=f(x), f бесквадратен, deg f = 6:
#   аффинные точки  <-> f(x) квадрат (или 0) при некотором x in Q_p;
#   точки на бесконечности (две) <-> старший коэффициент квадрат в Q_p.
# Покрытие P^1(Q_p): x in Z_p  U  {x = 1/w, w in pZ_p} (w=0 = бесконечность).
#
# solve(h, cc, p): "существует ли y in Z_p с p^cc * h(y) in квадратах u {0}"
#   1) вынести содержание p^c, cc <- (cc+c) mod 2, h примитивен;
#   2) ЯРЛЫК ПРОСТОГО КОРНЯ: если есть y0 с h(y0)=0 mod p и h'(y0)!=0 mod p,
#      то по Гензелю h имеет простой корень r in Z_p, h(y)=(y-r)u(y), u(r)
#      единица.  Берём y = r + p^M w, M большое нужной чётности, w единица:
#      p^cc h(y) = p^(cc+M) w u(y); класс единицы w u(y) пробегает все классы
#      => ДА при любом cc.  (Верно и для p=2: берём M кратное большое, тогда
#      u(y) = u(r) mod 8.)
#   3) иначе по классам вычетов:
#      p нечётно, y0 mod p:  h(y)=h(y0) mod p;
#         h(y0)!=0: ДА <=> cc==0 и (h(y0)|p)=+1 (Гензель);  иначе класс пуст;
#         h(y0)==0: кратный корень -> рекурсия y=y0+p*w.
#      p=2, y0 mod 8: h(y0+8k)-h(y0) делится на 8, т.е. h(y)=h(y0) mod 8;
#         h(y0) нечётно: ДА <=> cc==0 и h(y0)=1 mod 8; иначе класс пуст;
#         h(y0) чётно: рекурсия y=y0+8*w.
# Завершаемость: f бесквадратен над Q_p, v_p(disc) конечна.  Глубина
# ограничена MAXDEPTH; при превышении бросается исключение — «нет» по
# исчерпанию лимита не возвращается НИКОГДА.
# ======================================================================

MAXDEPTH = 150


class DepthExceeded(Exception):
    pass


def _content_val(h, p):
    return min([c.valuation(p) for c in h.list() if c != 0])


def solve(h, cc, p, depth=0):
    """Есть ли y in Z_p с p^cc*h(y) in (Q_p^*)^2 u {0}?  h in Z[y]."""
    if depth > MAXDEPTH:
        raise DepthExceeded("глубина рекурсии > %d при p=%d" % (MAXDEPTH, p))
    if h.is_zero():
        return True
    c = _content_val(h, p)
    if c > 0:
        h = Rz([ZZ(x) // p**c for x in h.list()])
    cc = (cc + c) % 2

    hp = h.derivative()
    for y0 in range(p):                      # ярлык простого корня
        if ZZ(h(y0)) % p == 0 and ZZ(hp(y0)) % p != 0:
            return True

    if p != 2:
        for y0 in range(p):
            val = ZZ(h(y0)) % p
            if val != 0:
                if cc == 0 and kronecker(val, p) == 1:
                    return True
            else:
                if solve(h(y0 + p * z), cc, p, depth + 1):
                    return True
        return False
    else:
        for y0 in range(8):
            val = ZZ(h(y0))
            if val % 2 == 1:
                if cc == 0 and val % 8 == 1:
                    return True
            else:
                if solve(h(y0 + 8 * z), cc, p, depth + 1):
                    return True
        return False


def is_locally_solvable(fpoly, p, deg):
    """H(Q_p) != пусто для гладкой модели y^2 = fpoly(x), deg чётна."""
    fp = Rz(fpoly)
    assert deg % 2 == 0 and fp.degree() <= deg
    coeffs = fp.list() + [ZZ(0)] * (deg + 1 - len(fp.list()))
    if solve(fp, 0, p, 0):                       # x in Z_p
        return True
    frev = Rz(list(reversed(coeffs)))            # x = 1/w, w in p Z_p
    return solve(frev(p * z), 0, p, 0)


# ======================================================================
# 3b. Алгоритм (B): независимый полный перебор дисков P^1(Z/p^N)
# ======================================================================
# Диск v=1: x in u + p^N Z_p.  f(x) = f(u) mod p^N.
#   Пусть A=f(u) in Z, w=v_p(A).  Если w <= N-1 (p нечётно) или w <= N-3 (p=2),
#   то класс квадратов величины f(x) ОПРЕДЕЛЁН для всех x из диска и равен
#   классу A.  Иначе диск НЕОПРЕДЕЛЁН.
# Диск u=1, v = p*v' : G(1,v) = frev(v), та же логика.
# Возврат: True (найден определённый квадратный диск),
#          False (все диски определены и неквадратны),
#          None  (остались неопределённые диски — увеличить N).
# ======================================================================

def _horner_mod(coeffs_desc, x, mod):
    """Значение многочлена по схеме Горнера по модулю mod (обычные int)."""
    acc = 0
    for c in coeffs_desc:
        acc = (acc * x + c) % mod
    return acc


def _sq_class_determined(Amod, p, N):
    """(определён?, квадрат?) для значения A с известным A mod p^N.
    Если v_p(A) < N, то v_p(A)=v_p(Amod) и A/p^w = Amod/p^w mod p^(N-w)."""
    if Amod == 0:
        return (False, None)          # v_p(A) >= N — точности не хватает
    w = 0
    a = Amod
    while a % p == 0:
        a //= p
        w += 1
    slack = 1 if p != 2 else 3
    if w > N - slack:
        return (False, None)
    if w % 2 == 1:
        return (True, False)
    if p != 2:
        return (True, kronecker(a % p, p) == 1)
    return (True, (a % 8) == 1)


def _bf_scan(fpoly, p, Nmax, deg, early_stop):
    """Адаптивный перебор дисков P^1(Z_p): диск дробится только если класс
    квадратов на нём ещё не определён.  Возвращает
    (verdict, nsq, nnonsq, nundet, свидетель).
    verdict: True/False/None (None = остались неопределённые диски)."""
    fp = Rz(fpoly)
    coeffs = [int(c) for c in (fp.list() + [ZZ(0)] * (deg + 1 - len(fp.list())))]
    cd_aff = list(reversed(coeffs))          # старшие первыми, для Горнера
    cd_inv = coeffs[:]                       # обращённый многочлен
    nsq = nnonsq = nundet = 0
    witness = None
    for tag, cd in (("aff", cd_aff), ("inv", cd_inv)):
        # (N, r): диск {x = r mod p^N}; ветвь inv: x=1/w, w in pZ_p, r = 0 mod p
        stack = [(1, r) for r in range(p)] if tag == "aff" else [(1, 0)]
        while stack:
            N, r = stack.pop()
            mod = p**N
            Amod = _horner_mod(cd, r % mod, mod)
            det, sq = _sq_class_determined(Amod, p, N)
            if det:
                if sq:
                    nsq += 1
                    if witness is None:
                        witness = (tag, N, r)
                    if early_stop:
                        return True, nsq, nnonsq, nundet, witness
                else:
                    nnonsq += 1
                continue
            if N >= Nmax:
                nundet += 1
                continue
            for d in range(p):
                stack.append((N + 1, r + d * mod))
    if nsq > 0:
        return True, nsq, nnonsq, nundet, witness
    if nundet == 0:
        return False, nsq, nnonsq, nundet, None
    return None, nsq, nnonsq, nundet, None


def decide_bruteforce(fpoly, p, N, deg):
    v, a, b, c, w = _bf_scan(fpoly, p, N, deg, early_stop=True)
    return v, (w if v else ("undetermined discs", c))


# ======================================================================
# 3c. Валидация алгоритма (A)
# ======================================================================
print("\n" + "=" * 78)
print("3c. Валидация решающего алгоритма (A)")
print("=" * 78)

print("(i) элементарные тесты для solve (вопрос: ∃x in Z_p с h(x) в квадратах ∪{0}) :")
elem = [
    (Rz(-1),            5, True,  "-1 квадрат в Q_5"),
    (Rz(-1),            3, False, "-1 не квадрат в Q_3"),
    (Rz(3),             3, False, "3 не квадрат в Q_3"),
    (Rz(z),             3, True,  "x=0 даёт значение 0"),
    (Rz(3*(z**2 + 1)),  3, False, "v_3 нечётна на всём Z_3"),
    (Rz(2),             2, False, "2 не квадрат в Q_2"),
    (Rz(17),            2, True,  "17 = 1 mod 8"),
    (Rz(z**2 + 1),      3, True,  "x=0 даёт 1 = квадрат"),
    (Rz(z**2 + 1),      5, True,  "x=0 даёт 1; и -1 квадрат в Z_5"),
    (Rz(5*z + 5),       5, True,  "x=-1 даёт 0"),
    (Rz(7*(z**2+1)),    7, False, "-1 не вычет mod 7 => v_7(x^2+1)=0 => v нечётна"),
    (Rz(7*(z**4+1)),    7, False, "z^4+1 единица на Z_7 => v нечётна"),
    (Rz(3*(z**2+1)),    7, True,  "x=1: 6; x=2: 15; x=3: 30=2 mod 7 ... x=0: 3 не выч.; x=4: 51=2; x=5: 78=1 выч."),
]
elem_ok = True
for (hh, pp, expect, why) in elem:
    got = solve(hh, 0, pp, 0)
    flag = (got == expect)
    elem_ok = elem_ok and flag
    print("    %-16s p=%-3d ожидалось %-5s получено %-5s %-6s  %s"
          % (hh, pp, expect, got, "OK" if flag else "СБОЙ", why))
print("    элементарные тесты:", "ПРОЙДЕНЫ" if elem_ok else "ПРОВАЛ")

print("\n(ii) сверка алгоритма (A) с независимой C-реализацией Sage test_qpls")
print("     (Qp_soluble для квартик y^2 = a x^4 + b x^3 + c x^2 + d x + e):")
from sage.schemes.elliptic_curves.descent_two_isogeny import test_qpls

set_random_seed(20260912)
quartic_bad = []
ncomp = 0
nref_true = 0
nref_false = 0
by_p = {}
for trial in range(12000):
    p = choice([2, 3, 5, 7, 11, 13, 17, 19, 233])
    rng = choice([2, 3, 6, 20, 250])
    A = ZZ.random_element(-rng, rng + 1)
    B = ZZ.random_element(-rng, rng + 1)
    C = ZZ.random_element(-rng, rng + 1)
    D = ZZ.random_element(-rng, rng + 1)
    E = ZZ.random_element(-rng, rng + 1)
    q = Rz(A*z**4 + B*z**3 + C*z**2 + D*z + E)
    if q.is_zero() or q.degree() != 4 or not q.is_squarefree():
        continue
    ref = bool(test_qpls(A, B, C, D, E, p))
    try:
        mine = is_locally_solvable(q, p, 4)
    except DepthExceeded:
        quartic_bad.append((A, B, C, D, E, p, "DEPTH"))
        continue
    ncomp += 1
    by_p[p] = by_p.get(p, 0) + 1
    nref_true += 1 if ref else 0
    nref_false += 0 if ref else 1
    if mine != ref:
        quartic_bad.append((A, B, C, D, E, p, ref, mine))
print("     сверено квартик:", ncomp, " по простым:", dict(sorted(by_p.items())))
print("     из них эталон test_qpls дал НЕПУСТО: %d, ПУСТО: %d"
      % (nref_true, nref_false))
print("     (обе ветви ответа реально задействованы — алгоритм не «всегда да»)")
print("     расхождений:", len(quartic_bad))
for bc in quartic_bad[:10]:
    print("        ", bc)
quartic_ok = (len(quartic_bad) == 0 and ncomp > 3000)
print("     сверка на квартиках:", "ПРОЙДЕНА" if quartic_ok else "ПРОВАЛ")

print("\n(iii) сверка алгоритма (A) с независимым алгоритмом (B) на случайных секстиках:")
set_random_seed(777)
sext_bad = []
nsext = 0
nsext_true = 0
nsext_false = 0
nundet = 0
for trial in range(400):
    p = choice([2, 3, 5, 7, 13])
    N = {2: 20, 3: 14, 5: 12, 7: 12, 13: 10}[p]
    cs = [ZZ.random_element(-8, 9) for _ in range(7)]
    q = Rz(cs)
    if q.degree() != 6 or not q.is_squarefree():
        continue
    mine = is_locally_solvable(q, p, 6)
    bf, info = decide_bruteforce(q, p, N, 6)
    if bf is None:
        nundet += 1
        if mine is False:
            sext_bad.append((cs, p, "A=False, B=undetermined", info))
        continue
    nsext += 1
    nsext_true += 1 if bf else 0
    nsext_false += 0 if bf else 1
    if mine != bf:
        sext_bad.append((cs, p, bf, mine))
print("     сверено секстик (алгоритм B дал определённый ответ):", nsext)
print("     из них (B) дал НЕПУСТО: %d, ПУСТО: %d" % (nsext_true, nsext_false))
print("     случаев, где B не определился (A=True — согласовано):", nundet)
print("     расхождений:", len(sext_bad))
for bc in sext_bad[:10]:
    print("        ", bc)
sext_ok = (len(sext_bad) == 0 and nsext > 100)
print("     сверка на секстиках:", "ПРОЙДЕНА" if sext_ok else "ПРОВАЛ")

# ======================================================================
# 4. Точные сертификаты
# ======================================================================

def exact_certificate(fpoly, p, tvals):
    """Ищет ТОЧНЫЙ сертификат точки над Q_p: t in Q с f(t) in (Q_p^*)^2,
    либо точку на бесконечности.  Проверка целочисленная (валюация + класс
    единицы), без p-адических приближений."""
    lc = ZZ(fpoly.leading_coefficient())
    w = lc.valuation(p)
    u = lc // p**w
    if w % 2 == 0 and ((p != 2 and kronecker(u % p, p) == 1) or (p == 2 and u % 8 == 1)):
        return ("бесконечность", lc, w, u)
    for tv in tvals:
        val = QQ(fpoly(tv))
        if val == 0:
            return (tv, 0, None, None)
        w = val.valuation(p)
        uu = val / QQ(p)**w
        num, den = ZZ(uu.numerator()), ZZ(uu.denominator())
        uint = num * inverse_mod(den % p**5, p**5) % p**5    # единица mod p^5
        if w % 2 != 0:
            continue
        if p != 2:
            if kronecker(uint % p, p) == 1:
                return (tv, val, w, uint % p)
        else:
            if uint % 8 == 1:
                return (tv, val, w, uint % 8)
    return None


CERT_T = [QQ(k) for k in range(-60, 61)] + \
         [QQ(a)/QQ(b) for a in range(-30, 31) for b in range(1, 16) if gcd(a, b) == 1] + \
         [QQ(a)/QQ(2**e) for a in range(-40, 41, 2) for e in range(1, 7)] + \
         [QQ(a)/QQ(3**e) for a in range(-40, 41) for e in range(1, 5)] + \
         [QQ(a)*QQ(233) for a in range(-40, 41)] + \
         [QQ(a)/QQ(233) for a in range(-40, 41)]
CERT_T = list(dict.fromkeys(CERT_T))
# сортируем по высоте, чтобы сертификаты были минимальными и проверялись руками
CERT_T.sort(key=lambda q: (max(abs(q.numerator()), q.denominator()), abs(q.numerator())))

# ======================================================================
# 5. Проверка H(13,8) во всех местах
# ======================================================================
print("\n" + "=" * 78)
print("5. Локальная разрешимость H(13,8) по всем местам")
print("=" * 78)

print("Хорошая редукция и p >= 17: род 2, оценка Вейля #H~(F_p) >= p+1-4*sqrt(p).")
print("   p=17 : 17+1-4*sqrt(17) = %+.4f  > 0" % (17 + 1 - 4*sqrt(17.0)))
print("   p=13 : 13+1-4*sqrt(13) = %+.4f  <= 0 (оценка не работает; но 13 плохое)"
      % (13 + 1 - 4*sqrt(13.0)))
print("   => для всех p >= 17 хорошей редукции есть гладкая F_p-точка,")
print("      она поднимается по лемме Гензеля: H(Q_p) != пусто.")
print("   Ниже это дополнительно проверено прямым счётом для 17 <= p < 500.")

check_primes = sorted(set(bad + list(prime_range(60))))
BF_N = {2: 18, 3: 14, 5: 12, 7: 12, 11: 12, 13: 12, 233: 8}

print("\nявно проверяются простые:", check_primes)
print("\n  %-5s %-7s %-9s %-11s %s" % ("p", "плохое", "алг.(A)", "алг.(B)", "точный сертификат"))
print("  " + "-" * 92)
results = {}
for p in check_primes:
    try:
        ansA = is_locally_solvable(f, p, 6)
        errA = None
    except DepthExceeded as ex:
        ansA, errA = None, str(ex)
    N = BF_N.get(p, 12)
    ansB, infoB = decide_bruteforce(f, p, N, 6)
    cert = exact_certificate(f, p, CERT_T)
    results[p] = (ansA, ansB, cert, errA)
    if cert is None:
        certstr = "-"
    elif cert[0] == "бесконечность":
        certstr = "t=бесконечность, lc=%s=p^%d*u, u=%s (квадрат)" % (cert[1], cert[2], cert[3])
    elif cert[1] == 0:
        certstr = "t=%s, f(t)=0" % cert[0]
    else:
        certstr = "t=%s, f(t)=p^%d*u, u=%s mod p%s (квадрат)" % (
            cert[0], cert[2], cert[3], "^3" if p == 2 else "")
    print("  %-5d %-7s %-9s %-11s %s"
          % (p, "да" if p in bad else "нет",
             ("НЕПУСТО" if ansA else "ПУСТО") if ansA is not None else "ОШИБКА",
             ("НЕПУСТО" if ansB else ("ПУСТО" if ansB is False else "не опред. (N=%d)" % N)),
             certstr))

print("\n--- полная статистика дисков P^1(Z/p^N) для плохих простых (без ранней остановки) ---")
print("(диск дробится только пока класс квадратов на нём не определён;")
print(" остаточные 'неопред.' диски стягиваются к p-адическим корням f)")
print("  %-5s %-6s %-12s %-14s %-12s %s"
      % ("p", "Nmax", "квадратных", "неквадратных", "неопред.", "вердикт"))
bfstat_ok = True
for p in bad + [11, 17, 19]:
    Nmax = {2: 16, 3: 12, 5: 10, 7: 10, 13: 8, 233: 6}.get(p, 8)
    v, a1, a2, a3, w = _bf_scan(f, p, Nmax, 6, early_stop=False)
    print("  %-5d %-6d %-12d %-14d %-12d %s"
          % (p, Nmax, a1, a2, a3,
             "НЕПУСТО" if v else ("ПУСТО" if v is False else "не определён")))
    if a1 == 0:
        bfstat_ok = False
        print("     ВНИМАНИЕ: ни одного определённо-квадратного диска при p=%d" % p)

print("\n--- КОНТРОЛЬНАЯ ГРУППА: кривые, где препятствие ЕСТЬ ---")
print("(показывает, что конвейер реально умеет выдавать ПУСТО, а не всегда ДА)")
ctrl = [
    (Rz(7*(z**2+1)*(z**4+1)), 7,  False,
     "палиндром, v_7 нечётна всюду на P^1(Q_7)"),
    (Rz(3*(z**2+1)*(z**4+z**2+1) + 0), 3, None, "контроль (ответ вычисляется)"),
    (Rz(-f),                   None, None, "-f: нет ВЕЩЕСТВЕННЫХ точек"),
    (Rz(233*f),               233,  None, "233*f (класс квадратов сдвинут на 233)"),
    (Rz(2*f),                   2,  None, "2*f (класс квадратов сдвинут на 2)"),
    (Rz(3*f),                   3,  None, "3*f"),
]
for (g, pp, expect, why) in ctrl:
    if pp is None:
        neg = all(g(QQ(k)/QQ(7)) < 0 for k in range(-200, 201)) and g.roots(RR) == []
        print("   %-28s R: вещественных точек нет: %-5s  [%s]"
              % ("-f", neg, why))
        continue
    gA = is_locally_solvable(g, pp, 6)
    gB, _ = decide_bruteforce(g, pp, {2: 18, 3: 14, 233: 8}.get(pp, 12), 6)
    mark = ""
    if expect is not None:
        mark = "  ожидалось %s -> %s" % (expect, "OK" if gA == expect else "СБОЙ")
    print("   %-28s p=%-4d (A)=%-8s (B)=%-8s  [%s]%s"
          % (str(g)[:28], pp,
             "НЕПУСТО" if gA else "ПУСТО",
             "НЕПУСТО" if gB else ("ПУСТО" if gB is False else "не опред."),
             why, mark))

print("\n--- независимая перепроверка точных сертификатов через Qp ---")
for p in [2, 3, 5, 7, 13, 233]:
    cert = results[p][2]
    if cert[0] == "бесконечность":
        val = QQ(lead)
        desc = "старший коэффициент (точка на бесконечности)"
    else:
        val = QQ(f(cert[0]))
        desc = "f(%s)" % cert[0]
    K = Qp(p, 80)
    print("   p=%-4d %-28s v_p=%-3d  Qp.is_square = %s"
          % (p, desc, val.valuation(p), K(val).is_square()))
    assert K(val).is_square()

print("\n--- прямой контроль: точки над F_p для хорошей редукции 17 <= p < 500 ---")
fq = f.change_ring(QQ)
weil_ok = True
weil_min = None
for p in prime_range(17, 500):
    if p in bad:
        continue
    Hp = HyperellipticCurve(f.change_ring(GF(p)))
    npts = Hp.count_points(1)[0]
    if weil_min is None or npts < weil_min[1]:
        weil_min = (p, npts)
    if npts == 0:
        weil_ok = False
        print("   p=%d: точек над F_p НЕТ (!)" % p)
print("   все хорошие p в [17,500) имеют точки над F_p:", weil_ok)
print("   минимум #H(F_p) по этому диапазону:", weil_min)

print("\n--- явная ГЛАДКАЯ F_p-точка (основание Гензеля) для нескольких p ---")
for p in [17, 19, 23, 29, 31, 37, 41, 43, 101, 499]:
    if p in bad:
        continue
    Fp = GF(p)
    fb = f.change_ring(Fp)
    found = None
    for x0 in Fp:
        v = fb(x0)
        if v != 0 and v.is_square():
            found = (x0, "Y!=0 -> dF/dY != 0, Гензель по Y")
            break
    if found is None:
        for x0 in Fp:
            if fb(x0) == 0 and fb.derivative()(x0) != 0:
                found = (x0, "Y=0, f'(x)!=0 -> Гензель по X")
                break
    print("   p=%-4d x=%-4s %s" % (p, found[0], found[1]))

# ======================================================================
# 6. Точки на бесконечности
# ======================================================================
print("\n" + "=" * 78)
print("6. Точки на бесконечности гладкой модели")
print("=" * 78)
print("deg f = 6 чётна => на гладкой модели две точки над бесконечностью,")
print("определённые над Q_v тогда и только тогда, когда старший коэффициент")
print("является квадратом в Q_v.")
print("lc =", lead, "=", factor(lead), ", свободная от квадратов часть =",
      lead.squarefree_part())
print("над Q: %d не квадрат => рациональных точек на бесконечности НЕТ" % lead)
print("над R: %d > 0 => две вещественные точки на бесконечности" % lead)
inf_yes, inf_no = [], []
for p in check_primes:
    w = lead.valuation(p)
    u = lead // p**w
    isq = (w % 2 == 0) and ((p != 2 and kronecker(u % p, p) == 1) or (p == 2 and u % 8 == 1))
    (inf_yes if isq else inf_no).append(p)
print("p, где точки на бесконечности ЕСТЬ (lc квадрат в Q_p):", inf_yes)
print("p, где точек на бесконечности НЕТ:", inf_no)
print("для последних препятствия нет: выше предъявлены аффинные точки.")

# ======================================================================
# 7. Итог
# ======================================================================
print("\n" + "=" * 78)
print("7. ИТОГ")
print("=" * 78)
emptyA = [p for p in results if results[p][0] is not True]
emptyB = [p for p in results if results[p][1] is False]
undetB = [p for p in results if results[p][1] is None]
nocert = [p for p in results if results[p][2] is None]
errs   = [p for p in results if results[p][3] is not None]
print("места с ПУСТЫМ H(Q_p) по алгоритму (A):", emptyA)
print("места с ПУСТЫМ H(Q_p) по алгоритму (B):", emptyB)
print("места, где (B) не определился:", undetB, "(у всех них (A)=НЕПУСТО и есть сертификат)")
print("места без точного сертификата:", nocert)
print("места с ошибкой расчёта:", errs)
print("H(R) != пусто:", real_ok)
print("валидация (A): элементарные %s, квартики %s, секстики против (B) %s"
      % ("OK" if elem_ok else "ПРОВАЛ",
         "OK" if quartic_ok else "ПРОВАЛ",
         "OK" if sext_ok else "ПРОВАЛ"))
print("статистика дисков: у каждого плохого p есть определённо-квадратный диск:",
      bfstat_ok)
allgood = (not emptyA and not emptyB and not nocert and not errs and not undetB
           and real_ok and elem_ok and quartic_ok and sext_ok and weil_ok and bfstat_ok)
print()
if allgood:
    print("ВЫВОД: H(13,8) ЛОКАЛЬНО РАЗРЕШИМА ВО ВСЕХ МЕСТАХ Q.")
    print("ЛОКАЛЬНОГО ПРЕПЯТСТВИЯ НЕТ.  obstruction_found = false.")
else:
    print("ВЫВОД: однозначного заключения нет, см. списки выше.")
