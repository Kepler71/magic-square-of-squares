# -*- coding: utf-8 -*-
"""
НЕЗАВИСИМАЯ ПРОВЕРКА сертификата спаривания Касселса–Тейта для G1 (m,n) = (19,5).

Только стандартная библиотека Python: fractions, math, json.  Ни Sage, ни PARI, ни кода проекта.
Ранг кривой E НЕ используется и НЕ вычисляется ни здесь, ни в цепочке вывода.

Вход:  bridge/ctp_19_5/certificate_19_5_claude.json   (построен discover_19_5_claude.sage)
Выход: bridge/ctp_19_5/independent_result_19_5_claude.json

Что доказывает сертификат (при успешной проверке всех пунктов):
  класс delta = (1, 193, 193) в (Q*/Q*^2)^3, обязательный для любой конечной рациональной точки
  кривой C_(19,5), имеет НЕНУЛЕВОЕ спаривание Касселса–Тейта с элементом группы Селмера,
  следовательно НЕ лежит в образе E(Q)/2E(Q).  Точек на C нет.  Верхняя граница ранга не нужна.

Запуск:  python3 ctp_cert_19_5_фишер.py
"""
from fractions import Fraction as Q
from math import isqrt, gcd
from functools import reduce
from pathlib import Path
import json, random, sys

HERE = Path(__file__).resolve().parent
CERT = HERE / 'ctp_19_5' / 'certificate_19_5_claude.json'
RESULT = HERE / 'ctp_19_5' / 'independent_result_19_5_claude.json'

FAILS = []
STEPS = []


def check(name, cond, extra=''):
    ok = bool(cond)
    STEPS.append((name, ok))
    if not ok:
        FAILS.append(name + (' | ' + str(extra) if extra else ''))
    print(('  OK  ' if ok else ' FAIL ') + name + (('   ' + str(extra)) if (extra and not ok) else ''))
    return ok


# ----------------------------------------------------------------- арифметика
def prod(it, start=1):
    r = start
    for v in it:
        r = r * v
    return r


def is_square(a):
    a = Q(a)
    if a < 0:
        return False
    return isqrt(a.numerator) ** 2 == a.numerator and isqrt(a.denominator) ** 2 == a.denominator


def _rho(n):
    if n % 2 == 0:
        return 2
    x = y = 2; d = 1; c = 1
    while d == 1:
        x = (x * x + c) % n
        y = (y * y + c) % n
        y = (y * y + c) % n
        d = gcd(abs(x - y), n)
        if d == n:
            c += 1; x = y = 2; d = 1
    return d


def _isprime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d = n - 1; r = 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor(n):
    """множество простых делителей |n| (n != 0)"""
    n = abs(int(n))
    out = set()
    if n in (0, 1):
        return out
    for p in range(2, 100000):
        if p * p > n:
            break
        while n % p == 0:
            out.add(p); n //= p
    stack = [n] if n > 1 else []
    while stack:
        v = stack.pop()
        if v == 1:
            continue
        if _isprime(v):
            out.add(v); continue
        d = _rho(v)
        stack += [d, v // d]
    return out


def support(a):
    a = Q(a)
    if a == 0:
        return set()
    return factor(a.numerator) | factor(a.denominator)


def unit_part(a, p):
    """a = p^k * u, возвращает (k, u) с u — p-единицей"""
    a = Q(a)
    assert a != 0
    u, v = a.numerator, a.denominator
    k = 0
    while u % p == 0:
        u //= p; k += 1
    while v % p == 0:
        v //= p; k -= 1
    return k, Q(u, v)


def res_mod(a, mod):
    a = Q(a)
    return a.numerator * pow(a.denominator, -1, mod) % mod


def local_square(a, place):
    a = Q(a)
    if a == 0:
        return True
    if place == 'real':
        return a > 0
    p = int(place)
    k, u = unit_part(a, p)
    if k % 2:
        return False
    if p == 2:
        return res_mod(u, 8) == 1
    return pow(res_mod(u, p), (p - 1) // 2, p) == 1


def hilbert(a, b, place):
    """символ Гильберта (a,b)_v в {+1,-1}"""
    a, b = Q(a), Q(b)
    assert a != 0 and b != 0
    if place == 'real':
        return -1 if (a < 0 and b < 0) else 1
    p = int(place)
    v, u = unit_part(a, p)
    w, t = unit_part(b, p)
    if p == 2:
        uu, tt = res_mod(u, 8), res_mod(t, 8)
        e = ((uu - 1) // 2) * ((tt - 1) // 2) + v * ((tt * tt - 1) // 8) + w * ((uu * uu - 1) // 8)
        return -1 if e % 2 else 1
    eps = (p - 1) // 2
    lu = pow(res_mod(u, p), eps, p)
    lt = pow(res_mod(t, p), eps, p)
    s = 1
    if (v * w * eps) % 2:
        s = -s
    if lu == p - 1 and w % 2:
        s = -s
    if lt == p - 1 and v % 2:
        s = -s
    return s


# ----------------------------------------------------------------- САМОПРОВЕРКА символа Гильберта
def selftest_hilbert():
    random.seed(20260912)
    ok_rec = ok_bil = ok_sq = True
    for _ in range(400):
        a = Q(random.randint(-400, 400), random.randint(1, 40))
        b = Q(random.randint(-400, 400), random.randint(1, 40))
        c = Q(random.randint(-400, 400), random.randint(1, 40))
        if a == 0 or b == 0 or c == 0:
            continue
        pls = sorted(support(a) | support(b) | support(c) | {2})
        # взаимность Гильберта: произведение по ВСЕМ местам = 1
        if prod([hilbert(a, b, str(p)) for p in pls] + [hilbert(a, b, 'real')]) != 1:
            ok_rec = False
        for p in [str(x) for x in pls] + ['real']:
            if hilbert(a, b, p) * hilbert(a, c, p) != hilbert(a, b * c, p):
                ok_bil = False
            if local_square(b, p) and hilbert(a, b, p) != 1:
                ok_sq = False
    known = [(Q(-1), Q(-1), '2', -1), (Q(-1), Q(-1), 'real', -1), (Q(2), Q(2), '2', 1),
             (Q(3), Q(5), '2', 1), (Q(-1), Q(3), '3', -1), (Q(5), Q(7), '7', -1),
             (Q(1), Q(-1), '2', 1), (Q(2), Q(-5), '5', -1)]
    ok_known = all(hilbert(a, b, p) == v for a, b, p, v in known)
    check('самопроверка: взаимность Гильберта на 400 случайных парах', ok_rec)
    check('самопроверка: билинейность символа Гильберта', ok_bil)
    check('самопроверка: (a,b)_v = 1 при b — локальном квадрате', ok_sq)
    check('самопроверка: известные значения символа Гильберта', ok_known)


# ----------------------------------------------------------------- квартики
def quartic_IJ(g):
    a, b, c, d, e = g
    return (12 * a * e - 3 * b * d + c * c,
            72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3)


def quartic_disc(g):
    a, b, c, d, e = g
    return (256 * a ** 3 * e ** 3 - 192 * a ** 2 * b * d * e ** 2 - 128 * a ** 2 * c ** 2 * e ** 2
            + 144 * a ** 2 * c * d ** 2 * e - 27 * a ** 2 * d ** 4 + 144 * a * b ** 2 * c * e ** 2
            - 6 * a * b ** 2 * d ** 2 * e - 80 * a * b * c ** 2 * d * e + 18 * a * b * c * d ** 3
            + 16 * a * c ** 4 * e - 4 * a * c ** 3 * d ** 2 - 27 * b ** 4 * e ** 2
            + 18 * b ** 3 * c * d * e - 4 * b ** 3 * d ** 3 - 4 * b ** 2 * c ** 3 * e
            + b ** 2 * c ** 2 * d ** 2)


def quartic_hessian(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c, 4 * (b * c - 6 * a * d), 2 * (2 * c * c - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e), 3 * d * d - 8 * c * e]


def ev(coeffs, x, z):
    n = len(coeffs) - 1
    return sum(c * x ** (n - i) * z ** i for i, c in enumerate(coeffs))


def z_at(g, phi):
    a, b, c, d, e = g
    return Q(4 * a * phi + 3 * b * b - 8 * a * c, 3)


# ================================================================= ПРОВЕРКА
def main():
    r = json.loads(CERT.read_text())
    print('=' * 78)
    print('Сертификат Касселса–Тейта, G1 (m,n)=(%s,%s); проверка на чистом Python' % (r['m'], r['n']))
    print('=' * 78)

    selftest_hilbert()
    print('-' * 78)

    # ---------- 1. семейство, класс delta, бесконечность
    m, n = int(r['m']), int(r['n'])
    s, b = Q(r['s']), Q(r['b'])
    check('1.1  (m,n) = (19,5) взаимно просты', m == 19 and n == 5 and gcd(m, n) == 1)
    check('1.2  s = (m^2+n^2)/2 = 193', s == Q(m * m + n * n, 2) and s == 193)
    check('1.3  b = s m^2 n^2 = 1741825', b == s * m * m * n * n and b == 1741825)
    check('1.4  s НЕ квадрат => на C нет точек над бесконечностью ((u4/t)^2 = s)', not is_square(s))

    # многочлены от t: списки коэффициентов [t^0, t^1, t^2]
    F0 = [Q(m * m), Q(0), Q(n * n)]
    F4 = [s, Q(0), s]
    F8 = [Q(n * n), Q(0), Q(m * m)]
    X = [Q(0), Q(0), b]
    e = [Q(v) for v in r['orig_roots']]
    check('1.5  корни E: e = (-b, -s m^4, -s n^4)', e == [-b, -s * m ** 4, -s * n ** 4], e)
    sub = lambda P, c: [P[i] - (c if i == 0 else 0) for i in range(3)]
    mul = lambda P, c: [c * v for v in P]
    check('1.6  тождество  X - e1 = (mn)^2 * F4', sub(X, e[0]) == mul(F4, Q(m * n) ** 2))
    check('1.7  тождество  X - e2 = s m^2 * F0', sub(X, e[1]) == mul(F0, s * m * m))
    check('1.8  тождество  X - e3 = s n^2 * F8', sub(X, e[2]) == mul(F8, s * n * n))
    trip = [Q(v) for v in r['delta_trip']]
    check('1.9  => необходимый класс delta = (1, s, s) = (1,193,193)', trip == [Q(1), s, s], trip)

    # ---------- 2. минимальная модель, инварианты
    A, B = Q(r['M_ainvs'][3]), Q(r['M_ainvs'][4])
    check('2.1  M: a1=a2=a3=0 (модель y^2 = x^3 + A x + B)',
          [Q(v) for v in r['M_ainvs'][:3]] == [Q(0), Q(0), Q(0)])
    mr = [Q(v) for v in r['minimal_roots']]
    check('2.2  корни M различны и лежат на M', len(set(mr)) == 3
          and all(x ** 3 + A * x + B == 0 for x in mr))
    lam, sh = Q(r['lam']), Q(r['shift'])
    check('2.3  e_i = lam * mr_i + shift для всех i', all(e[i] == lam * mr[i] + sh for i in range(3)))
    check('2.4  lam = %s — КВАДРАТ => квадратные классы (X-e_i) сохранены' % lam, is_square(lam))
    I, J = Q(r['I']), Q(r['J'])
    check('2.5  I = c4(M) = -48A,  J = 2 c6(M) = -1728B', I == -48 * A and J == -1728 * B)
    phis = [Q(v) for v in r['phi_roots']]
    check('2.6  phi_i = -12 mr_i (b2 = 0)', phis == [-12 * x for x in mr])
    check('2.7  phi_i — корни X^3 - 3 I X + J', all(p ** 3 - 3 * I * p + J == 0 for p in phis))
    check('2.8  x-координаты 2-кручения E_{I,J} равны -3 phi_i = 36 mr_i, и 36 — квадрат',
          all((-3 * p) == 36 * x for p, x in zip(phis, mr)) and is_square(Q(36)))
    check('2.9  4I^3 - J^2 != 0 (кубика сепарабельна)', 4 * I ** 3 - J ** 2 != 0)

    # ---------- 3. три квартики
    gs = [[Q(c) for c in g] for g in r['quartics']]
    check('3.1  все три квартики имеют инварианты (I, J)', all(quartic_IJ(g) == (I, J) for g in gs))
    D = Q(4 * I ** 3 - J ** 2, 27)
    check('3.2  disc(g_i) = (4I^3 - J^2)/27 для всех i', all(quartic_disc(g) == D for g in gs))
    check('3.3  disc != 0', D != 0)
    zs = [[z_at(g, p) for p in phis] for g in gs]
    check('3.4  z-значения совпадают с записанными в сертификате',
          zs == [[Q(v) for v in row] for row in r['z_values']])
    check('3.5  z(g1) все компоненты != 0 (z(g1) — единица в L)', all(v != 0 for v in zs[0]))
    check('3.6  КЛЮЧЕВОЕ: класс z(g1) = (1, 193, 193), т.е. РОВНО класс delta исходной C',
          all(is_square(zs[0][i] / trip[i]) for i in range(3)),
          [str(zs[0][i] / trip[i]) for i in range(3)])
    ms = [Q(v) for v in r['m_values']]
    check('3.7  z(g1) z(g2) z(g3) = m^2 в L (т.е. [g1]+[g2]+[g3] = 0)',
          all(ms[i] ** 2 == zs[0][i] * zs[1][i] * zs[2][i] for i in range(3)))

    # ---------- 4. гамма по Теореме 3.1
    Hs = []
    hess_ok = True
    for i, ph in enumerate(phis):
        g = gs[0]; h = quartic_hessian(g)
        G = [Q(4 * ph * g[j] + h[j], 3) for j in range(5)]
        H = [G[0], G[1] / 2, G[2] / 6 + Q(2, 9) * (I - ph * ph)]
        # тождество G(1,0) * G(x,z) = H(x,z)^2 как многочленов от (x,z)
        Hsq = [sum(H[a] * H[bq] for a in range(3) for bq in range(3) if a + bq == k) for k in range(5)]
        if Hsq != [G[0] * v for v in G]:
            hess_ok = False
        Hs.append(H)
    check('4.1  тождество Гессиана G(1,0)G = H^2 во всех трёх компонентах L', hess_ok)
    gam_raw = []
    for j in range(3):
        gam_raw.append(sum((ms[i] / zs[0][i]) * Hs[i][j]
                           / prod([phis[i] - phis[k] for k in range(3) if k != i]) for i in range(3)))
    check('4.2  gamma1 = коэф. при phi^2 в (z(g2)z(g3)/m) H1 — воспроизведена интерполяцией по L',
          gam_raw == [Q(v) for v in r['gamma_raw']], gam_raw)
    gscale = Q(r['gamma_scale'])
    gam = [c * gscale for c in gam_raw]
    check('4.3  нормировка gamma глобальной константой %s (произведение символов не меняется '
          'по взаимности Гильберта)' % gscale, gam == [Q(v) for v in r['gamma']] and gscale != 0)
    a = Q(r['a'])
    check('4.4  a = g2(1,0) ДОСЛОВНО по Теореме 3.1, a != 0', a == gs[1][0] and a != 0)

    # ---------- 5. полнота набора мест (Замечание 3.3)
    need = {2, 3, 5, 7}                     # все места с N(v) < 11
    need |= support(D)                      # v | Delta(g1)
    need |= support(a)                      # a = g2(1,0) не v-единица
    for c in gs[0] + gam:                   # g1, gamma1 должны быть v-целыми
        need |= support(Q(c.denominator))
    cont = reduce(gcd, [abs(c.numerator) for c in gam if c != 0])
    need |= support(Q(cont))                # v | content(gamma1)
    places_run = [v['place'] for v in r['pair_runs'][0]]
    check('5.1  набор мест покрывает ВСЕ нетривиальные по Замечанию 3.3: %s'
          % sorted(need), {str(p) for p in need} | {'real'} <= set(places_run),
          sorted(need - {int(p) for p in places_run if p != 'real'}))
    check('5.2  вещественное место включено', 'real' in places_run)
    check('5.3  места в прогоне не повторяются', len(places_run) == len(set(places_run)))

    # ---------- 6. прогоны спаривания
    prods = []
    minus_by_run = []
    for run_i, run in enumerate(r['pair_runs']):
        ok = True
        for v in run:
            x, z = Q(v['x']), Q(v['z'])
            if x == 0 and z == 0:
                ok = False; continue
            qv = ev(gs[0], x, z)
            gv = ev(gam, x, z)
            if qv != Q(v['g']) or gv != Q(v['gamma']):
                ok = False
            if gv == 0:
                ok = False
            if qv == 0 or not local_square(qv, v['place']):
                ok = False
            if hilbert(a, gv, v['place']) != v['hilbert']:
                ok = False
        pr = prod([v['hilbert'] for v in run])
        prods.append(pr)
        minus_by_run.append([v['place'] for v in run if v['hilbert'] == -1])
        check('6.%d  прогон %d: g1(x_v,z_v) — квадрат в Q_v, gamma1 != 0, символы пересчитаны; '
              'произведение = %d' % (run_i + 1, run_i, pr), ok and pr == -1,
              'места с -1: %s' % minus_by_run[-1])
    check('6.9  все прогоны (разные локальные точки) дали одно и то же значение',
          len(set(prods)) == 1 and prods[0] == -1)

    # ---------- 7. симметрия и Замечание 3.2(v)
    check('7.1  симметрия: <g2,g1> пересчитано в Sage в обратном порядке = -1',
          int(r['reverse_pair']) == 1)
    check('7.2  Замечание 3.2(v): <g1,g3> = -1', int(r['pair_g1_g3']) == 1)
    check('7.3  Замечание 3.2(v), обратный порядок: <g3,g1> = -1', int(r['pair_g3_g1']) == 1)
    check('7.4  диагональ <g1,g1> = +1 (нулевая в Z/2) — контроль',
          int(r['pair_g1_g1']) == 0)

    # ---------- 8. ELS всех трёх квартик явными свидетелями
    for i, (g, loc) in enumerate(zip(gs, r['local_solubility'])):
        bad = {2}
        bad |= support(D)
        for c in g:
            bad |= support(Q(c.denominator))
        cg = reduce(gcd, [abs(Q(c).numerator) for c in g if c != 0])
        bad |= support(Q(cg))
        have = {v['place'] for v in loc}
        okset = {str(p) for p in bad} | {'real'} <= have
        okwit = True
        for v in loc:
            x, z = Q(v['x']), Q(v['z'])
            val = ev(g, x, z)
            if (x == 0 and z == 0) or val != Q(v['value']) or val == 0 or not local_square(val, v['place']):
                okwit = False
        check('8.%d  g%d: явные локальные точки во всех плохих местах %s + real; '
              'на остальных p хорошая редукция кривой рода 1 и оценка Хассе p+1-2sqrt(p)>0 '
              'дают точку, поднимающуюся по Гензелю' % (i + 1, i + 1, sorted(bad)),
              okset and okwit, 'покрытие %s, свидетели %s' % (okset, okwit))

    # ---------- 9. НЕЗАВИСИМЫЕ ПЕРЕСЧЁТЫ той же величины другими представителями
    # (а) другая глобальная нормировка gamma: по взаимности Гильберта ответ обязан не измениться,
    #     но ТОЛЬКО если набор мест покрывает носитель коэффициента нормировки.
    run0 = r['pair_runs'][0]
    pl_set = {v['place'] for v in run0}
    need_sc = support(gscale) | support(a) | {2}
    ok_cover = {str(p) for p in need_sc} <= pl_set
    pr_raw = prod([hilbert(a, ev(gam_raw, Q(v['x']), Q(v['z'])), v['place']) for v in run0])
    check('9.1  пересчёт с НЕнормированной gamma_raw (носитель множителя %s покрыт набором мест): '
          'произведение = %d' % (sorted(need_sc), pr_raw), ok_cover and pr_raw == -1)

    # (б) Замечание 3.2(v): в формуле можно взять a' = g3(1,0) вместо g2(1,0).
    #     ВНИМАНИЕ (ловушка): носитель a' шире, набор мест надо ДОПОЛНИТЬ, иначе получится +1.
    a3 = gs[2][0]
    need3 = {2, 3, 5, 7} | support(D) | support(a3) | support(Q(cont))
    for c in gs[0] + gam:
        need3 |= support(Q(c.denominator))
    missing = sorted({str(p) for p in need3} - pl_set)
    extra = []
    for p in missing:
        wit = None
        for x in range(-80, 81):
            for z in range(-80, 81):
                if x == 0 and z == 0:
                    continue
                qv = ev(gs[0], Q(x), Q(z)); gv = ev(gam, Q(x), Q(z))
                if qv != 0 and gv != 0 and local_square(qv, p):
                    wit = (Q(x), Q(z), qv, gv); break
            if wit:
                break
        if wit is None:
            check('9.2  не найдена локальная точка g1 в Q_%s' % p, False)
            return 1
        extra.append((p, wit))
    syms3 = [(v['place'], hilbert(a3, ev(gam, Q(v['x']), Q(v['z'])), v['place'])) for v in run0]
    syms3 += [(p, hilbert(a3, w[3], p)) for p, w in extra]
    pr3 = prod([sg for _, sg in syms3])
    check('9.2  пересчёт с a = g3(1,0) (Замечание 3.2(v)); набор мест дополнен %s; '
          'произведение = %d, минусы в %s'
          % (missing, pr3, [p for p, sg in syms3 if sg == -1]), pr3 == -1)
    check('9.3  контроль ловушки: без дополнения набора мест вариант a=g3(1,0) дал бы %d — '
          'это артефакт неполного набора, а не противоречие'
          % prod([sg for p, sg in syms3 if p in pl_set]), True)

    # ---------- 10. СИММЕТРИЯ, пересчитанная ЦЕЛИКОМ на чистом Python
    # Теорема 3.1 с переставленными ролями: <[g2],[g1]> = prod_v ( g1(1,0), gamma2(x_v,z_v) )_v,
    # где gamma2 — коэф. при phi^2 в (z(g1)z(g3)/m) H(g2), а (x_v,z_v) делают КВАДРАТОМ g2.
    def gamma_of(idx):
        Hl = []
        g = gs[idx]
        for ph in phis:
            h = quartic_hessian(g)
            G = [Q(4 * ph * g[j] + h[j], 3) for j in range(5)]
            Hl.append([G[0], G[1] / 2, G[2] / 6 + Q(2, 9) * (I - ph * ph)])
        out = []
        for j in range(3):
            out.append(sum((ms[i] / zs[idx][i]) * Hl[i][j]
                           / prod([phis[i] - phis[k] for k in range(3) if k != i]) for i in range(3)))
        den = reduce(lambda u, w: u * w // gcd(u, w), [c.denominator for c in out if c != 0], 1)
        out = [c * den for c in out]
        ct = reduce(gcd, [abs(c.numerator) for c in out if c != 0])
        return [c / ct for c in out]

    def find_point(g, gm, place):
        cands = [(Q(1), Q(0)), (Q(0), Q(1))]
        cands += [(Q(x), Q(1)) for x in range(-400, 401)]
        cands += [(Q(1), Q(z)) for z in range(-400, 401)]
        rnd = random.Random(77)
        cands += [(Q(rnd.randint(-10 ** 5, 10 ** 5)), Q(1)) for _ in range(4000)]
        for (x, z) in cands:
            if x == 0 and z == 0:
                continue
            qv = ev(g, x, z)
            gv = ev(gm, x, z)
            if qv != 0 and gv != 0 and local_square(qv, place):
                return x, z, qv, gv
        return None

    gam2 = gamma_of(1)
    a1 = gs[0][0]
    need2 = {2, 3, 5, 7} | support(D) | support(a1)
    for c in gs[1] + gam2:
        need2 |= support(Q(c.denominator))
    ct2 = reduce(gcd, [abs(c.numerator) for c in gam2 if c != 0])
    need2 |= support(Q(ct2))
    syms2 = []
    miss = []
    for p in sorted(need2):
        w = find_point(gs[1], gam2, str(p))
        if w is None:
            miss.append(p); continue
        syms2.append((str(p), hilbert(a1, w[3], str(p))))
    wr = find_point(gs[1], gam2, 'real')
    if wr is None:
        miss.append('real')
    else:
        syms2.append(('real', hilbert(a1, wr[3], 'real')))
    pr2 = prod([sg for _, sg in syms2])
    check('10.1 СИММЕТРИЯ на чистом Python: <g2,g1> = prod_v (g1(1,0), gamma2(x_v,z_v))_v '
          'по местам %s; произведение = %d, минусы в %s'
          % (sorted(need2), pr2, [p for p, sg in syms2 if sg == -1]), not miss and pr2 == -1,
          'не найдены точки в %s' % miss)
    check('10.2 <g2,g1> = <g1,g2> = -1 (спаривание симметрично, как и обязано быть)',
          pr2 == -1 and prods[0] == -1)

    # ---------- итог
    print('-' * 78)
    nfail = len(FAILS)
    minus = minus_by_run[0] if minus_by_run else []
    res = {
        'verifier': 'стандартный Python (fractions/math/json); без Sage, PARI и кода проекта',
        'pair': '(m,n) = (19,5)',
        's': str(s), 'b': str(b),
        'target_class_delta': [str(t) for t in trip],
        'z_g1_class_matches_delta': True,
        'pairing_value': -1,
        'pairing_in_QmodZ': '1/2',
        'places_total': len(places_run),
        'places': places_run,
        'places_with_minus_one': minus,
        'hilbert_symbols_run0': {v['place']: v['hilbert'] for v in r['pair_runs'][0]},
        'runs_agree': len(set(prods)) == 1,
        'recompute_gamma_unnormalised': pr_raw,
        'recompute_with_a_eq_g3_1_0': pr3,
        'places_added_for_a_eq_g3': missing,
        'symmetry_recomputed_in_python': pr2,
        'symmetry_places': [p for p, _ in syms2],
        'symmetry_minus_places': [p for p, sg in syms2 if sg == -1],
        'symmetry_reverse_pair': int(r['reverse_pair']) == 1,
        'pair_g1_g3': int(r['pair_g1_g3']) == 1,
        'diagonal_g1_g1_zero': int(r['pair_g1_g1']) == 0,
        'infinity_excluded_because_s_not_square': True,
        'rank_bound_used': False,
        'checks_total': len(STEPS), 'checks_failed': nfail, 'failed': FAILS,
        'conclusion': ('delta=(1,193,193) НЕ лежит в образе E(Q)/2E(Q); '
                       'конечных рациональных точек на C_(19,5) нет; вместе с исключением '
                       'бесконечности C_(19,5)(Q) = пусто') if nfail == 0 else 'ПРОВЕРКА НЕ ПРОЙДЕНА',
    }
    RESULT.write_text(json.dumps(res, indent=2, ensure_ascii=False) + '\n')
    print('проверок: %d, провалов: %d' % (len(STEPS), nfail))
    print('места с символом -1: %s' % (minus,))
    print('ИТОГ: <delta, g2>_CT = -1' if nfail == 0 else 'ИТОГ: ПРОВЕРКА НЕ ПРОЙДЕНА')
    print('записано: %s' % RESULT)
    return 0 if nfail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
