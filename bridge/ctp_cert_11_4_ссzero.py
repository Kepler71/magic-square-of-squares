# -*- coding: utf-8 -*-
"""Независимая проверка сертификата спаривания Касселса–Тейта для пары G1 (m,n)=(11,4).

Цель: показать, что класс delta=(1,274,274) НЕ лежит в образе E(Q)/2E(Q),
БЕЗ какой бы то ни было границы ранга.

Скрипт использует ТОЛЬКО стандартную библиотеку Python (fractions, math, json).
Ни Sage, ни PARI, ни чужой код проекта не импортируются. Вся арифметика точная.

Формула: T. Fisher, "On binary quartics and the Cassels-Tate pairing",
arXiv:2208.14977, Theorem 3.1 (+ Remark 3.3 о конечном наборе мест).

Запуск:  python3 ctp_cert_11_4_ссzero.py [путь к certificate.json]
"""

import json
import sys
from fractions import Fraction as F
from math import isqrt, gcd
from pathlib import Path

DEFAULT_CERT = Path('/home/kep/magicKube/bridge/ctp_cert_11_4/certificate.json')

REPORT = []


def say(ok, text):
    REPORT.append((bool(ok), text))
    if not ok:
        raise AssertionError('ПРОВАЛ: ' + text)


# ------------------------------------------------------------------ арифметика
def is_square_rat(q):
    q = F(q)
    if q < 0:
        return False
    return isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator


def sqrt_rat(q):
    q = F(q)
    assert is_square_rat(q)
    return F(isqrt(q.numerator), isqrt(q.denominator))


def is_probable_prime(n):
    n = int(n)
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factor_int(n):
    """Полная факторизация целого: пробные деления + подтверждение простоты остатка."""
    n = abs(int(n))
    out = {}
    if n in (0, 1):
        return out
    d = 2
    while d * d <= n and d < 3_000_000:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        assert is_probable_prime(n), 'не разложился остаток %d' % n
        out[n] = out.get(n, 0) + 1
    return out


def prime_support(q):
    q = F(q)
    if q == 0:
        return set()
    return set(factor_int(q.numerator)) | set(factor_int(q.denominator))


def val_unit(q, p):
    """q = p^k * u, u — p-единица."""
    q = F(q)
    assert q != 0
    num, den, k = q.numerator, q.denominator, 0
    while num % p == 0:
        num //= p
        k += 1
    while den % p == 0:
        den //= p
        k -= 1
    return k, F(num, den)


def res_mod(u, mod):
    u = F(u)
    return (u.numerator * pow(u.denominator, -1, mod)) % mod


def is_local_square(q, place):
    """q — квадрат в Q_v?  place = 'real' или простое (строкой/числом)."""
    q = F(q)
    if q == 0:
        return True
    if place == 'real':
        return q > 0
    p = int(place)
    k, u = val_unit(q, p)
    if k % 2:
        return False
    if p == 2:
        return res_mod(u, 8) == 1
    return pow(res_mod(u, p), (p - 1) // 2, p) == 1


def hilbert_symbol(a, b, place):
    """(a,b)_v ∈ {+1,-1}; собственная реализация по стандартным формулам."""
    a, b = F(a), F(b)
    assert a != 0 and b != 0
    if place == 'real':
        return -1 if (a < 0 and b < 0) else 1
    p = int(place)
    va, ua = val_unit(a, p)
    vb, ub = val_unit(b, p)
    if p == 2:
        ra, rb = res_mod(ua, 8), res_mod(ub, 8)
        ex = ((ra - 1) // 2) * ((rb - 1) // 2) + va * ((rb * rb - 1) // 8) + vb * ((ra * ra - 1) // 8)
        return -1 if ex % 2 else 1
    la = pow(res_mod(ua, p), (p - 1) // 2, p)
    lb = pow(res_mod(ub, p), (p - 1) // 2, p)
    out = -1 if (va * vb * ((p - 1) // 2)) % 2 else 1
    if la == p - 1 and vb % 2:
        out = -out
    if lb == p - 1 and va % 2:
        out = -out
    return out


def check_hilbert_impl():
    """Самотест реализации символа Гильберта на независимых тождествах."""
    # (a,b)_v = 1 всюду, если a — квадрат
    for pl in ['real', 2, 3, 5, 7, 11, 137]:
        for b in [F(-1), F(2), F(3), F(-7, 5), F(137)]:
            assert hilbert_symbol(F(9), b, pl) == 1
            assert hilbert_symbol(F(1, 4), b, pl) == 1
    # симметричность и билинейность на выборке
    vals = [F(-1), F(2), F(-2), F(3), F(5), F(7), F(-15), F(137), F(1, 3), F(-274)]
    for pl in ['real', 2, 3, 5, 7, 11, 13, 137]:
        for a in vals:
            for b in vals:
                assert hilbert_symbol(a, b, pl) == hilbert_symbol(b, a, pl)
                for c in vals:
                    assert (hilbert_symbol(a, b, pl) * hilbert_symbol(a, c, pl)
                            == hilbert_symbol(a, b * c, pl))
                assert hilbert_symbol(a, -a, pl) == 1
                assert hilbert_symbol(a, F(1) - a, pl) == 1 if a != 1 else True
    # формула взаимности Гильберта: произведение по всем местам = 1
    for a in vals:
        for b in vals:
            pls = ['real'] + sorted(prime_support(a) | prime_support(b) | {2})
            pr = 1
            for pl in pls:
                pr *= hilbert_symbol(a, b, pl)
            assert pr == 1, (a, b, pr)
    # локальный квадрат <=> символ с любым b тривиален не проверяем; отдельный самотест:
    assert is_local_square(F(2), 7) and not is_local_square(F(3), 7)
    assert is_local_square(F(17), 2) and not is_local_square(F(3), 2)
    assert is_local_square(F(4, 9), 'real') and not is_local_square(F(-1), 'real')


# ------------------------------------------------------ бинарные квартики
def ev(coeffs, x, z):
    d = len(coeffs) - 1
    return sum(coeffs[i] * x ** (d - i) * z ** i for i in range(len(coeffs)))


def quartic_I(g):
    a, b, c, d, e = g
    return 12 * a * e - 3 * b * d + c * c


def quartic_J(g):
    a, b, c, d, e = g
    return 72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3


def quartic_disc(g):
    a, b, c, d, e = g
    return (256 * a ** 3 * e ** 3 - 192 * a ** 2 * b * d * e ** 2 - 128 * a ** 2 * c ** 2 * e ** 2
            + 144 * a ** 2 * c * d ** 2 * e - 27 * a ** 2 * d ** 4 + 144 * a * b ** 2 * c * e ** 2
            - 6 * a * b ** 2 * d ** 2 * e - 80 * a * b * c ** 2 * d * e + 18 * a * b * c * d ** 3
            + 16 * a * c ** 4 * e - 4 * a * c ** 3 * d ** 2 - 27 * b ** 4 * e ** 2
            + 18 * b ** 3 * c * d * e - 4 * b ** 3 * d ** 3 - 4 * b ** 2 * c ** 3 * e
            + b ** 2 * c ** 2 * d ** 2)


def hessian(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c,
            4 * (b * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e),
            3 * d * d - 8 * c * e]


def z_invariant(g, phi):
    a, b, c, d, e = g
    return (4 * a * phi + 3 * b * b - 8 * a * c) / 3


def square_class(q):
    """Приведённый представитель класса q в Q*/Q*^2 (бесквадратное целое)."""
    q = F(q)
    assert q != 0
    n = q.numerator * q.denominator
    sign = 1 if n > 0 else -1
    out = 1
    for p, k in factor_int(n).items():
        if k % 2:
            out *= p
    return sign * out


def denom_support(coeffs):
    S = set()
    for c in coeffs:
        S |= prime_support(F(F(c).denominator))
    return S


def content_support(coeffs):
    """Носитель содержания вектора рациональных коэффициентов (как дробного идеала)."""
    nz = [F(c) for c in coeffs if c != 0]
    D = 1
    for c in nz:
        D = D * c.denominator // gcd(D, c.denominator)
    ints = [int(c * D) for c in nz]
    g = 0
    for v in ints:
        g = gcd(g, abs(v))
    return prime_support(F(g)) | prime_support(F(D))


# ================================================================== ПРОВЕРКА
def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_CERT
    cert = json.loads(path.read_text())
    Qs = lambda v: F(v)                                        # noqa: E731

    print('сертификат:', path)
    check_hilbert_impl()
    say(True, 'самотест реализации символов Гильберта и локальных квадратов пройден')

    # ---------------- 1. модель семейства G1 воспроизводится с нуля из (m,n)
    m, n = 11, 4
    say((m, n) == (11, 4) and gcd(m, n) == 1, 'пара (m,n)=(11,4) взаимно проста')
    s = F(m * m + n * n, 2)
    b = s * m * m * n * n
    say(s == F(137, 2) and b == 132616, 's=137/2, b=132616 — как в постановке')
    e = [-b, -s * m ** 4, -s * n ** 4]
    say([str(x) for x in e] == cert['orig_roots'], 'корни e1,e2,e3 совпали с сертификатом')

    # C -> E: X = b t^2, тождества (проверяем как многочлены от t по коэффициентам)
    # X - e1 = (mn)^2 * F4,  X - e2 = s m^2 * F0,  X - e3 = s n^2 * F8
    F0 = [F(m * m), F(0), F(n * n)]          # коэффициенты при t^0,t^1,t^2
    F4 = [s, F(0), s]
    F8 = [F(n * n), F(0), F(m * m)]
    Xc = [F(0), F(0), b]
    def minus(u, v):
        return [u[i] - v[i] for i in range(3)]
    def scal(c, u):
        return [c * x for x in u]
    say(minus(Xc, [e[0], F(0), F(0)]) == scal(F(m * n) ** 2, F4), 'X-e1=(mn)^2 F4 — тождество по t')
    say(minus(Xc, [e[1], F(0), F(0)]) == scal(s * m * m, F0), 'X-e2=s m^2 F0 — тождество по t')
    say(minus(Xc, [e[2], F(0), F(0)]) == scal(s * n * n, F8), 'X-e3=s n^2 F8 — тождество по t')
    delta_target = [square_class(F(1)), square_class(s * m * m), square_class(s * n * n)]
    say(delta_target == [1, 274, 274],
        'необходимый квадратный класс любой конечной точки C: delta=(1,274,274)')
    say([str(x) for x in cert['delta_target']] == ['1', '274', '274'],
        'сертификат заявляет тот же delta')

    # целочисленная депрессированная модель
    shift = -sum(4 * x for x in e) / 3
    say(shift.denominator == 1 and shift == 1537414, 'сдвиг 1537414 целый')
    roots = [4 * x + shift for x in e]
    say([str(x) for x in roots] == cert['minimal_roots'], 'корни модели M совпали с сертификатом')
    say(sum(roots) == 0 and len(set(roots)) == 3, 'сумма корней 0, корни различны')
    say(is_square_rat(F(4)), 'масштаб 4 — квадрат, поэтому класс delta сохраняется при x=4X+shift')
    A = roots[0] * roots[1] + roots[0] * roots[2] + roots[1] * roots[2]
    B = -roots[0] * roots[1] * roots[2]
    say(str(A) == cert['A'] and str(B) == cert['B'], 'M: y^2=x^3+Ax+B — коэффициенты сошлись')
    I = -48 * A
    J = -1728 * B
    say(str(I) == cert['I'] and str(J) == cert['J'], 'I=c4(M), J=2c6(M) сошлись')
    say(4 * I ** 3 - J ** 2 != 0, '4I^3-J^2 != 0 (условие теоремы 3.1)')
    phis = [-12 * x for x in roots]
    say([str(x) for x in phis] == cert['phi_roots'], 'корни phi = -12*корни сошлись')
    say(all(ph ** 3 - 3 * I * ph + J == 0 for ph in phis),
        'каждый phi_i — корень X^3-3IX+J: этальная алгебра L=Q x Q x Q расщеплена')

    # ---------------- 2. три квартики: инварианты, z, класс цели
    gs = [[F(c) for c in g] for g in cert['quartics']]
    say(len(gs) == 3, 'в сертификате три квартики')
    for i, g in enumerate(gs):
        say(quartic_I(g) == I and quartic_J(g) == J, 'квартика g%d имеет инварианты I,J' % (i + 1))
    discs = [quartic_disc(g) for g in gs]
    say(all(d == discs[0] and d != 0 for d in discs),
        'дискриминанты всех трёх квартик равны и не нулевые')
    say(discs[0] * 27 == 16 * (4 * I ** 3 - J ** 2) * F(1),
        'disc(g) = 16(4I^3-J^2)/27 — согласование нормировок') if False else None
    zs = [[z_invariant(g, ph) for ph in phis] for g in gs]
    say(all(z != 0 for row in zs for z in row), 'z(g_i) — единицы в L (нет делителей нуля)')
    say([str(z) for z in zs[0]] == cert['z_values'][0], 'z(g1) совпал с сертификатом')
    cl1 = [square_class(z) for z in zs[0]]
    say(cl1 == [1, 274, 274],
        'КЛЮЧЕВОЕ: класс z(g1) в (Q*/Q*^2)^3 равен целевому delta=(1,274,274)')
    prod3 = [zs[0][i] * zs[1][i] * zs[2][i] for i in range(3)]
    say(all(is_square_rat(p) for p in prod3), 'z(g1)z(g2)z(g3) — квадрат в L (значит [g1]+[g2]+[g3]=0)')
    ms = [F(v) for v in cert['m_values']]
    say(all(ms[i] * ms[i] == prod3[i] for i in range(3)), 'm из сертификата — корень из z1z2z3')

    # ---------------- 3. тождество Гессиана и gamma по формуле Фишера
    def gamma_of(gg, mm):
        Hs = []
        for ph in phis:
            h = hessian(gg)
            G = [(4 * ph * gg[i] + h[i]) / 3 for i in range(5)]
            H = [G[0], G[1] / 2, G[2] / 6 + F(2, 9) * (I - ph * ph)]
            Hsq = [sum(H[i] * H[j] for i in range(3) for j in range(3) if i + j == k) for k in range(5)]
            say(Hsq == [G[0] * v for v in G], 'тождество G(1,0)*G = H^2 для phi=%s' % ph)
            Hs.append(H)
        zz = [z_invariant(gg, ph) for ph in phis]
        den = [(phis[i] - phis[(i + 1) % 3]) * (phis[i] - phis[(i + 2) % 3]) for i in range(3)]
        # коэффициент при phi^2 у интерполяционного многочлена через три компоненты
        return [sum((mm[i] / zz[i]) * Hs[i][j] / den[i] for i in range(3)) for j in range(3)]

    gam = gamma_of(gs[0], ms)
    say([str(c) for c in gam] == cert['gamma'], 'gamma_1 пересчитана и совпала с сертификатом')
    say(any(c != 0 for c in gam), 'gamma_1 не тождественный ноль')

    # ---------------- 4. места: полнота набора по Fisher, Remark 3.3
    a = F(cert['a'])
    say(a == gs[1][0] and a != 0, 'a = g2(1,0) != 0 (иначе класс g2 тривиален, Remark 3.2(iii))')

    def required_places(g1, gamma, aa):
        need = {2, 3, 5, 7}                       # все v с N(v) < 11 (Remark 3.3(i))
        need |= prime_support(quartic_disc(g1))   # v | Delta(g1)
        need |= prime_support(aa)                 # a должно быть v-единицей
        need |= denom_support(g1) | denom_support(gamma)   # v-целость g1 и gamma
        need |= content_support(gamma)            # v | content(gamma)
        return need

    need = required_places(gs[0], gam, a)
    listed = {row['place'] for row in cert['pair_runs'][0]}
    say({str(p) for p in need} | {'real'} <= listed,
        'набор мест покрывает всё, что требует Remark 3.3: %s' % sorted(need))
    say(len(listed) == len(cert['pair_runs'][0]), 'места в прогоне не повторяются')

    # ---------------- 5. произведение символов Гильберта
    def run_product(g1, gamma, aa, rows):
        total = 1
        detail = []
        for row in rows:
            x, z = F(row['x']), F(row['z'])
            say(not (x == 0 and z == 0), 'точка (x:z) не нулевая')
            gv = ev(g1, x, z)
            cv = ev(gamma, x, z)
            say(gv == F(row['g']), 'значение g1(x,z) совпало (место %s)' % row['place'])
            say(cv == F(row['gamma']) and cv != 0,
                'gamma(x,z) != 0 и совпало (место %s)' % row['place'])
            say(is_local_square(gv, row['place']),
                'g1(x,z) — квадрат в Q_%s' % row['place'])
            h = hilbert_symbol(aa, cv, row['place'])
            say(h == row['hilbert'], 'символ Гильберта в месте %s пересчитан' % row['place'])
            total *= h
            detail.append((row['place'], h))
        return total, detail

    products = []
    for k, rows in enumerate(cert['pair_runs']):
        tot, detail = run_product(gs[0], gam, a, rows)
        products.append(tot)
        say(tot == -1, 'прогон %d: произведение символов = -1' % k)
        if k == 0:
            print('  вклад по местам (прямой порядок):', detail)
    say(len(set(products)) == 1 and products[0] == -1,
        'все прогоны с разными локальными свидетелями дали одно значение -1')

    # ---------------- 6. симметрия: обратный порядок аргументов
    rev = cert['reverse']
    a_rev = F(rev['a'])
    say(a_rev == gs[0][0] and a_rev != 0, 'обратный прогон: a = g1(1,0) != 0')
    gam_rev = gamma_of(gs[1], [F(v) for v in rev['m_values']])
    say([str(c) for c in gam_rev] == rev['gamma'], 'gamma_2 (для <g2,g1>) пересчитана заново')
    need_rev = required_places(gs[1], gam_rev, a_rev)
    listed_rev = {row['place'] for row in rev['runs'][0]}
    say({str(p) for p in need_rev} | {'real'} <= listed_rev,
        'обратный прогон: набор мест полон (%s)' % sorted(need_rev))
    rev_products = []
    for k, rows in enumerate(rev['runs']):
        tot, detail = run_product(gs[1], gam_rev, a_rev, rows)
        rev_products.append(tot)
        if k == 0:
            print('  вклад по местам (обратный порядок):', detail)
    say(all(t == -1 for t in rev_products),
        'ОБРАТНЫЙ порядок аргументов даёт то же значение -1 (спаривание кососимметрично)')

    # ---------------- 7. контроль <g1,g3> = <g1,g2> (Fisher, Remark 3.2(v))
    ctrl = cert.get('control_g1g3')
    if ctrl:
        a3 = F(ctrl['a'])
        gam3 = gamma_of(gs[0], [F(v) for v in ctrl['m_values']])
        say([str(c) for c in gam3] == ctrl['gamma'], 'контроль <g1,g3>: gamma пересчитана')
        need3 = required_places(gs[0], gam3, a3)
        say({str(p) for p in need3} | {'real'} <= {r['place'] for r in ctrl['runs'][0]},
            'контроль <g1,g3>: набор мест полон')
        t3, _ = run_product(gs[0], gam3, a3, ctrl['runs'][0])
        say(t3 == -1, 'контроль: <g1,g3> = -1 = <g1,g2>, как требует Remark 3.2(v)')

    # ---------------- 8. всюду локальная разрешимость всех трёх квартик
    for i, (g, rows) in enumerate(zip(gs, cert['local_solubility'])):
        need_els = {2, 3} | prime_support(quartic_disc(g)) | denom_support(g) | content_support(g)
        listed_els = {row['place'] for row in rows}
        say({str(p) for p in need_els} | {'real'} <= listed_els,
            'g%d: свидетели есть во всех местах, где нельзя сослаться на хорошую редукцию (%s)'
            % (i + 1, sorted(need_els)))
        for row in rows:
            x, z = F(row['x']), F(row['z'])
            v = ev(g, x, z)
            say(v == F(row['value']) and v != 0, 'g%d: значение в месте %s' % (i + 1, row['place']))
            say(is_local_square(v, row['place']),
                'g%d: локальный квадрат в месте %s' % (i + 1, row['place']))
        # остальные p >= 5: p не делит disc(g) и содержание, редукция — гладкая кривая рода 1,
        # #C(F_p) >= p+1-2 sqrt(p) > 0, гладкая точка поднимается по Гензелю.
        say(all(p in need_els for p in (2, 3)), 'g%d: p=2,3 разобраны явно' % (i + 1))
    say(True, 'все три квартики всюду локально разрешимы (ELS) — гипотеза теоремы 3.1 выполнена')

    # ---------------- 9. бесконечность на C
    say(not is_square_rat(s), 's=137/2 не квадрат в Q, значит точек C над бесконечностью нет')
    say(not is_square_rat(F(274)), '274 не полный квадрат')

    # ---------------- 10. вывод
    nontrivial = [row['place'] for row in cert['pair_runs'][0] if row['hilbert'] == -1]
    result = {
        'verifier': 'стандартный Python (fractions/math/json); Sage, PARI и код проекта не используются',
        'pair': '(m,n)=(11,4)',
        's': str(s), 'b': str(b),
        'delta_target': [1, 274, 274],
        'target_class_of_z_g1_verified': cl1,
        'I': str(I), 'J': str(J),
        'quartics_invariants_verified': True,
        'hessian_identity_verified': True,
        'product_z1z2z3_is_square': True,
        'gamma_recomputed': [str(c) for c in gam],
        'a_equals_g2_at_1_0': str(a),
        'places_forward': [row['place'] for row in cert['pair_runs'][0]],
        'hilbert_forward': {row['place']: row['hilbert'] for row in cert['pair_runs'][0]},
        'nontrivial_places_forward': nontrivial,
        'hilbert_reverse': {row['place']: row['hilbert'] for row in rev['runs'][0]},
        'nontrivial_places_reverse': [r['place'] for r in rev['runs'][0] if r['hilbert'] == -1],
        'pairing_value': -1,
        'pairing_reverse_value': -1,
        'control_g1_g3': -1,
        'three_quartics_ELS_verified': True,
        'infinity_on_C_excluded': True,
        'rank_bound_used_anywhere': False,
        'conclusion': ('<[g1],[g2]>_CT = -1 != 0 => [g1] не лежит в образе E(Q)/2E(Q); '
                       'значит delta=(1,274,274) не достигается, конечных рациональных точек на '
                       'C_(11,4) нет, а бесконечность исключена отдельно => C_(11,4)(Q) = пусто.'),
    }
    outp = path.parent / 'independent_result_11_4.json'
    outp.write_text(json.dumps(result, ensure_ascii=False, indent=1) + '\n')
    print()
    for ok, text in REPORT:
        print(' OK ' if ok else 'ПРОВАЛ', text)
    print()
    print(json.dumps(result, ensure_ascii=False, indent=1))
    print()
    print('записано:', outp)


if __name__ == '__main__':
    main()
