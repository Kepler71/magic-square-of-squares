# -*- coding: utf-8 -*-
u"""СЕРТИФИКАТ спаривания Касселса–Тейта для пары G1 (m,n) = (11,4).

Цель: доказать, что класс delta = (1, 274, 274) НЕ лежит в образе E(Q)/2E(Q),
БЕЗ какой бы то ни было верхней границы ранга E(Q).

Что делает скрипт:
  * читает certificate.json = {I, J, корни phi, три квартики g1,g2,g3, контроли};
  * ЗАНОВО считает инварианты I,J каждой квартики, z(g), корень m, гессиан,
    форму H, тождество G(1,0)G = H^2, форму gamma по Теореме 3.1 Фишера;
  * САМ ищет локальные точки (x_v, z_v) (в сертификате их нет);
  * САМ определяет полный набор мест по Замечанию 3.3 Фишера (факторизация);
  * считает символы Гильберта и их произведение;
  * пересчитывает спаривание в ОБРАТНОМ порядке аргументов;
  * проверяет всюду локальную разрешимость всех трёх квартик;
  * контроль ложных срабатываний: классы 2-кручения E(Q) (они ТОЧНО в образе
    E(Q)/2E(Q)) обязаны спариваться тривиально;
  * связывает результат с исходной кривой C_(11,4) рода 5.

Никакого Sage, никакого PARI, никакого чужого кода проекта: только стандартная
библиотека Python (fractions, math, random, json).

Формула: T. Fisher, "On binary quartics and the Cassels-Tate pairing",
arXiv:2208.14977, Theorem 3.1; конечный набор мест -- Remark 3.3.

Запуск:  python3 ctp_cert_11_4_indep.py [certificate.json]
"""

import json
import random
import sys
from fractions import Fraction as F
from math import gcd, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent

# ====================================================================== отчёт
CHECKS = []


def check(cond, text):
    CHECKS.append((bool(cond), text))
    if not cond:
        raise AssertionError(u'ПРОВАЛ: ' + text)
    return True


def note(text):
    print(u'    ' + text)


# ============================================================ целые и рациональные
def _is_probable_prime(n):
    if n < 2:
        return False
    for p in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % p == 0:
            return n == p
    d, r = n - 1, 0
    while d % 2 == 0:
        d //= 2
        r += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _pollard(n):
    if n % 2 == 0:
        return 2
    while True:
        c = random.randrange(1, n)
        x = y = random.randrange(0, n)
        d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n:
            return d


def factor(n):
    u"""Полная факторизация |n| (пробное деление + Полларда ро). n != 0."""
    n = abs(int(n))
    out = {}
    if n in (0, 1):
        return out
    for p in range(2, 100000):
        if p * p > n:
            break
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        v = stack.pop()
        if v == 1:
            continue
        if _is_probable_prime(v):
            out[v] = out.get(v, 0) + 1
            continue
        d = _pollard(v)
        stack += [d, v // d]
    return out


def primes_of(q):
    u"""Простые в разложении рационального q != 0 (числитель и знаменатель)."""
    q = F(q)
    if q == 0:
        return set()
    return set(factor(q.numerator)) | set(factor(q.denominator))


def rat_sqrt(q):
    u"""Точный квадратный корень из рационального q, или None."""
    q = F(q)
    if q < 0:
        return None
    a, b = isqrt(q.numerator), isqrt(q.denominator)
    if a * a == q.numerator and b * b == q.denominator:
        return F(a, b)
    return None


def squarefree_part(q):
    u"""Класс q по модулю квадратов: бесквадратное целое того же класса."""
    q = F(q)
    if q == 0:
        return 0
    sign = -1 if q < 0 else 1
    n = abs(q.numerator * q.denominator)
    out = 1
    for p, e in factor(n).items():
        if e % 2:
            out *= p
    return sign * out


def content(coeffs):
    u"""Контент набора рациональных чисел: рациональное c с v_p(c)=min_j v_p(a_j)."""
    nz = [F(c) for c in coeffs if c != 0]
    if not nz:
        return F(0)
    den = 1
    for c in nz:
        den = den * c.denominator // gcd(den, c.denominator)
    nums = [int(c * den) for c in nz]
    g = 0
    for v in nums:
        g = gcd(g, abs(v))
    return F(g, den)


# =================================================== локальная арифметика над Q_v
def val_unit(a, p):
    u"""a != 0 рациональное -> (v_p(a), единица u с a = p^v u)."""
    a = F(a)
    num, den = a.numerator, a.denominator
    v = 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v, F(num, den)


def unit_mod(u, q):
    u"""Вычет рациональной p-единицы u по модулю q (знаменатель обратим)."""
    u = F(u)
    return (u.numerator * pow(u.denominator, -1, q)) % q


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def is_local_square(a, place):
    u"""Является ли рациональное a квадратом в Q_v (place = 'inf' или простое)."""
    a = F(a)
    if a == 0:
        return True
    if place == 'inf':
        return a > 0
    p = int(place)
    v, u = val_unit(a, p)
    if v % 2:
        return False
    if p == 2:
        return unit_mod(u, 8) == 1
    return legendre(unit_mod(u, p), p) == 1


def hilbert(a, b, place):
    u"""Символ Гильберта (a,b)_v для рациональных a,b != 0. Реализация с нуля."""
    a, b = F(a), F(b)
    assert a != 0 and b != 0
    if place == 'inf':
        return -1 if (a < 0 and b < 0) else 1
    p = int(place)
    va, ua = val_unit(a, p)
    vb, ub = val_unit(b, p)
    if p == 2:
        ra, rb = unit_mod(ua, 8), unit_mod(ub, 8)
        eps_a, eps_b = ((ra - 1) // 2) % 2, ((rb - 1) // 2) % 2
        om_a, om_b = ((ra * ra - 1) // 8) % 2, ((rb * rb - 1) // 8) % 2
        expo = (eps_a * eps_b + va * om_b + vb * om_a) % 2
        return -1 if expo else 1
    res = 1
    if (va * vb * ((p - 1) // 2)) % 2:
        res = -res
    if vb % 2 and legendre(unit_mod(ua, p), p) == -1:
        res = -res
    if va % 2 and legendre(unit_mod(ub, p), p) == -1:
        res = -res
    return res


def selftest_hilbert():
    u"""Самопроверка символа Гильберта: Штейнберг, симметрия, билинейность, взаимность."""
    random.seed(20260912)
    nplaces = set()

    def allplaces(*vals):
        s = {'inf', 2}
        for v in vals:
            s |= primes_of(v)
        return s

    for _ in range(60):
        a = F(random.randint(-60, 60) or 1, random.randint(1, 30))
        b = F(random.randint(-60, 60) or 1, random.randint(1, 30))
        pl = allplaces(a, b)
        nplaces |= pl
        # симметрия
        for v in pl:
            check(hilbert(a, b, v) == hilbert(b, a, v), u'символ Гильберта симметричен')
        # взаимность: произведение по всем местам = 1
        prodv = 1
        for v in pl:
            prodv *= hilbert(a, b, v)
        check(prodv == 1, u'взаимность Гильберта для a=%s b=%s' % (a, b))
        # (a,-a)_v = 1
        for v in allplaces(a):
            check(hilbert(a, -a, v) == 1, u'(a,-a)_v = 1')
        # Штейнберг: (a,1-a)_v = 1
        if a != 1:
            for v in allplaces(a, 1 - a):
                check(hilbert(a, 1 - a, v) == 1, u'(a,1-a)_v = 1 (Штейнберг)')
        # билинейность по второму аргументу
        c = F(random.randint(-40, 40) or 1, random.randint(1, 20))
        for v in allplaces(a, b, c):
            check(hilbert(a, b * c, v) == hilbert(a, b, v) * hilbert(a, c, v),
                  u'билинейность символа Гильберта')
    # (a,b)_v = 1 <=> z^2 = a x^2 + b y^2 разрешимо в Q_v -- грубая сверка для малых p
    for p in (2, 3, 5, 7):
        for a in (-3, -2, -1, 1, 2, 3, 5, 6, 7, p, 2 * p):
            for b in (-3, -2, -1, 1, 2, 3, 5, 6, 7, p, 2 * p):
                q = p ** (7 if p == 2 else 5)
                found = False
                for x in range(q):
                    for y in range(q):
                        if x % p == 0 and y % p == 0:
                            continue
                        val = a * x * x + b * y * y
                        if val == 0:
                            found = True
                            break
                        v, u = val_unit(F(val), p)
                        # хватает точности, если q >> val
                        if v + (3 if p == 2 else 1) <= 7 and is_local_square(F(val), p):
                            found = True
                            break
                    if found:
                        break
                if found:
                    check(hilbert(F(a), F(b), p) == 1,
                          u'(%d,%d)_%d = +1 при наличии решения' % (a, b, p))
    note(u'самопроверка символа Гильберта пройдена, мест задействовано: %d' % len(nplaces))


# ======================================================= бинарные квартики
def q_I(g):
    a, b, c, d, e = g
    return 12 * a * e - 3 * b * d + c * c


def q_J(g):
    a, b, c, d, e = g
    return 72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3


def hessian(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c,
            4 * (b * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e),
            3 * d * d - 8 * c * e]


def ev(form, x, z):
    u"""Значение однородной формы (коэффициенты по убыванию степени x)."""
    n = len(form) - 1
    return sum(form[j] * x ** (n - j) * z ** j for j in range(n + 1))


def z_comps(g, phis):
    a, b, c = g[0], g[1], g[2]
    return [(4 * a * ph + 3 * b * b - 8 * a * c) / F(3) for ph in phis]


def H_form(g, ph, I):
    u"""Квадратичная форма H (Фишер, формула (4)) при вложении phi -> ph."""
    h = hessian(g)
    G = [(4 * ph * g[j] + h[j]) / F(3) for j in range(5)]
    H = [G[0], G[1] / 2, G[2] / 6 + F(2, 9) * (I - ph * ph)]
    Hsq = [H[0] * H[0],
           2 * H[0] * H[1],
           H[1] * H[1] + 2 * H[0] * H[2],
           2 * H[1] * H[2],
           H[2] * H[2]]
    ok = all(Hsq[j] == G[0] * G[j] for j in range(5))
    return H, ok


# ================================================= поиск локальных точек
def p1_reps(p, k):
    u"""Представители P^1(Z/p^k) целыми парами (x,z)."""
    q = p ** k
    for x in range(q):
        yield (x, 1)
    for w in range(p ** (k - 1)):
        yield (1, p * w)


def find_local_point(g, gamma, place, kmax=8):
    u"""(x,z) с g(x,z) квадратом в Q_v, g != 0 и gamma(x,z) != 0."""
    cand = [(1, 0), (0, 1), (1, 1), (1, -1), (2, 1), (1, 2), (-1, 1), (3, 1), (1, 3)]
    for x in range(-40, 41):
        cand.append((x, 1))
        cand.append((1, x))
    for x in range(-12, 13):
        for z in range(1, 13):
            cand.append((x, z))
    if place == 'inf':
        cand += [(1, t) for t in range(-200, 201)]
    seen = set()
    for (x, z) in cand:
        if (x, z) == (0, 0) or (x, z) in seen:
            continue
        seen.add((x, z))
        val = ev(g, F(x), F(z))
        if val == 0 or ev(gamma, F(x), F(z)) == 0:
            continue
        if is_local_square(val, place):
            return (x, z)
    if place == 'inf':
        return None
    p = int(place)
    for k in range(1, kmax + 1):
        if p ** k > 4 * 10 ** 6:
            break
        for (x, z) in p1_reps(p, k):
            if (x, z) == (0, 0):
                continue
            val = ev(g, F(x), F(z))
            if val == 0 or ev(gamma, F(x), F(z)) == 0:
                continue
            if is_local_square(val, place):
                return (x, z)
    return None


# ============================================== множество мест (Фишер, Rem. 3.3)
def places_for(g1, gamma, a, Delta):
    u"""Все места, где вклад может быть нетривиален.

    Замечание 3.3: место v даёт +1, если одновременно N(v) >= 11, формы g1 и
    gamma v-целые, v не делит Delta(g1)*content(gamma), a = g2(1,0) -- v-единица
    и v не делит 2.  Берём дополнение к этому.
    """
    S = {'inf', 2, 3, 5, 7}
    S |= primes_of(Delta)
    S |= primes_of(a)
    for c in list(g1) + list(gamma):
        if c != 0:
            S |= set(factor(F(c).denominator))
    ct = content(gamma)
    if ct != 0:
        S |= set(factor(ct.numerator))
    return S


# =========================================================== ядро: спаривание
class Fisher(object):
    def __init__(self, I, J, phis):
        self.I, self.J, self.phis = I, J, phis
        self.Delta = F(16, 27) * (4 * I ** 3 - J ** 2)

    def check_quartic(self, g, tag):
        check(q_I(g) == self.I, u'I(%s) = I' % tag)
        check(q_J(g) == self.J, u'J(%s) = J' % tag)

    def gamma(self, g1, g2, g3, signs=(1, 1, 1)):
        u"""gamma_1 из Теоремы 3.1: коэффициент при phi^2 в (z2 z3 / m) H_1."""
        ph = self.phis
        z1 = z_comps(g1, ph)
        z2 = z_comps(g2, ph)
        z3 = z_comps(g3, ph)
        check(all(v != 0 for v in z1), u'z(g1) -- единица в L')
        check(all(v != 0 for v in z2), u'z(g2) -- единица в L')
        check(all(v != 0 for v in z3), u'z(g3) -- единица в L')
        m = []
        for i in range(3):
            r = rat_sqrt(z1[i] * z2[i] * z3[i])
            check(r is not None, u'z1*z2*z3 -- квадрат в компоненте %d' % i)
            m.append(signs[i] * r)
        Hs = []
        for i in range(3):
            H, ok = H_form(g1, ph[i], self.I)
            check(ok, u'тождество G(1,0)G = H^2 в компоненте %d' % i)
            Hs.append(H)
        # (z2 z3)/m = m/z1
        for i in range(3):
            check(z2[i] * z3[i] / m[i] == m[i] / z1[i], u'z2 z3 / m = m / z1, комп. %d' % i)
        gam = []
        for j in range(3):
            tot = F(0)
            for i in range(3):
                den = F(1)
                for k in range(3):
                    if k != i:
                        den *= (ph[i] - ph[k])
                tot += (m[i] / z1[i]) * Hs[i][j] / den
            gam.append(tot)
        return gam, m, (z1, z2, z3)

    def pair(self, g1, g2, g3, tag='', signs=(1, 1, 1), verbose=True):
        u"""<[g1],[g2]>_CT по Теореме 3.1. Возвращает (значение, детали по местам)."""
        for g, t in ((g1, tag + 'g1'), (g2, tag + 'g2'), (g3, tag + 'g3')):
            self.check_quartic(g, t)
        gam, m, zs = self.gamma(g1, g2, g3, signs=signs)
        a = F(g2[0])
        check(a != 0, u'g2(1,0) != 0 (Замечание 3.2(iii))')
        S = places_for(g1, gam, a, self.Delta)
        rows = []
        total = 1
        for place in sorted(S, key=lambda v: (v == 'inf', 0 if v == 'inf' else int(v))):
            pt = find_local_point(g1, gam, place)
            check(pt is not None, u'найдена локальная точка на g1 при v = %s' % place)
            x, z = pt
            gv = ev(g1, F(x), F(z))
            gm = ev(gam, F(x), F(z))
            check(gv != 0 and is_local_square(gv, place),
                  u'g1(x_v,z_v) -- квадрат в Q_%s' % place)
            check(gm != 0, u'gamma(x_v,z_v) != 0 при v = %s' % place)
            h = hilbert(a, gm, place)
            rows.append({'place': str(place), 'x': str(x), 'z': str(z),
                         'g1': str(gv), 'gamma': str(gm), 'hilbert': h})
            total *= h
            if verbose:
                note(u'v = %-6s (x:z) = (%s:%s)   (a, gamma)_v = %+d' % (place, x, z, h))
        return total, rows, gam, m, zs


# =========================================================== ELS квартики
def els_witnesses(g, Delta, tag):
    u"""Свидетели всюду локальной разрешимости y^2 = g(x,z)."""
    S = {'inf', 2, 3} | primes_of(Delta)
    for c in g:
        if c != 0:
            S |= set(factor(F(c).denominator))
    rows = []
    for place in sorted(S, key=lambda v: (v == 'inf', 0 if v == 'inf' else int(v))):
        pt = find_local_point(g, [F(1), F(0), F(0)], place)   # gamma != 0 тут не нужна
        check(pt is not None, u'ELS: локальная точка %s при v = %s' % (tag, place))
        x, z = pt
        val = ev(g, F(x), F(z))
        check(val != 0 and is_local_square(val, place),
              u'ELS: %s -- квадрат в Q_%s' % (tag, place))
        rows.append({'place': str(place), 'x': str(x), 'z': str(z), 'value': str(val)})
    return rows, sorted(str(v) for v in S)


# =========================================================== основной проход
def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else HERE / 'certificate.json'
    if not path.exists():
        path = HERE / 'candidates.json'
    R = json.loads(path.read_text())
    print(u'\n=== СЕРТИФИКАТ Касселса-Тейта, пара %s ===' % R['pair'])
    print(u'источник данных: %s' % path)

    print(u'\n[0] самопроверка локальной арифметики')
    selftest_hilbert()

    I = F(R['I'])
    J = F(R['J'])
    phis = [F(v) for v in R['phi_roots']]
    Delta = F(16, 27) * (4 * I ** 3 - J ** 2)
    check(Delta == F(R['Delta']), u'Delta = 16(4I^3-J^2)/27 совпадает с сертификатом')
    check(Delta != 0, u'Delta != 0')
    check(len(set(phis)) == 3, u'три различных корня phi')
    for ph in phis:
        check(ph ** 3 - 3 * I * ph + J == 0, u'phi -- корень X^3-3IX+J')

    # ---------- связь с исходной задачей (независимо от сертификата) ----------
    print(u'\n[1] связь с исходной кривой C_(11,4) рода 5')
    m, n = 11, 4
    s = F(137, 2)
    b = s * m * m * n * n
    check(b == 132616, u'b = s m^2 n^2 = 132616')
    er = [-b, -s * m ** 4, -s * n ** 4]
    check([F(v) for v in R['original_roots']] == er, u'корни исходной E совпали')
    # X = b t^2:  X-e1 = (mn)^2 F4,  X-e2 = s m^2 F0,  X-e3 = s n^2 F8
    # проверяем как тождества по t^2 =: T
    T = F(7, 5)                       # произвольное T; тождества линейны по T
    for T in (F(7, 5), F(0), F(3), F(-11, 9)):
        X = b * T
        F0 = m * m + n * n * T
        F4 = s * (1 + T)
        F8 = n * n + m * m * T
        check(X - er[0] == (m * n) ** 2 * F4, u'X-e1 = (mn)^2 F4')
        check(X - er[1] == s * m * m * F0, u'X-e2 = s m^2 F0')
        check(X - er[2] == s * n * n * F8, u'X-e3 = s n^2 F8')
    check(squarefree_part((m * n) ** 2) == 1, u'класс (mn)^2 = 1')
    check(squarefree_part(s * m * m) == 274, u'класс s m^2 = 274')
    check(squarefree_part(s * n * n) == 274, u'класс s n^2 = 274')
    check(rat_sqrt(s) is None, u's = 137/2 не квадрат => на бесконечности C точек нет')
    note(u'необходимый класс любой конечной точки C: delta = (1, 274, 274)')

    # ---------- переход к минимальной модели ----------
    mr = [F(v) for v in R['minimal_roots']]
    scale = (mr[0] - mr[1]) / (er[0] - er[1])
    check((mr[0] - mr[2]) / (er[0] - er[2]) == scale, u'масштаб моделей одинаков (1,3)')
    check((mr[1] - mr[2]) / (er[1] - er[2]) == scale, u'масштаб моделей одинаков (2,3)')
    check(rat_sqrt(scale) is not None, u'масштаб %s -- квадрат, класс delta сохраняется' % scale)
    A, B = F(R['M_ainvs'][3]), F(R['M_ainvs'][4])
    check(F(R['M_ainvs'][0]) == 0 and F(R['M_ainvs'][1]) == 0 and F(R['M_ainvs'][2]) == 0,
          u'минимальная модель имеет вид y^2 = x^3 + A x + B')
    for e in mr:
        check(e ** 3 + A * e + B == 0, u'корень минимальной модели')
    check(sum(mr) == 0, u'сумма корней = 0')
    check(I == -48 * A, u'I = c4(M)')
    check(J == -1728 * B, u'J = 2 c6(M)')
    for i in range(3):
        check(phis[i] == -12 * mr[i], u'phi_i = -12 e_i (b2 = 0)')
    note(u'x_M = %s X + %s' % (scale, mr[0] - scale * er[0]))

    Fi = Fisher(I, J, phis)

    # ---------- целевая квартика ----------
    print(u'\n[2] целевая квартика g1 и её класс')
    g1 = [F(v) for v in R['g1']]
    Fi.check_quartic(g1, 'g1')
    z1 = z_comps(g1, phis)
    cls = [squarefree_part(v) for v in z1]
    check(cls == [1, 274, 274], u'z(g1) лежит в классе (1,274,274): получено %s' % (cls,))
    note(u'g1 = %s' % ([str(v) for v in g1],))
    note(u'z(g1) = %s' % ([str(v) for v in z1],))

    # ---------- перебор кандидатов g2 ----------
    print(u'\n[3] спаривание g1 с кандидатами g2 из ell2cover')
    results = []
    winner = None
    for cand in R['candidates']:
        if 'g3' not in cand:
            continue
        g2 = [F(v) for v in cand['g2']]
        g3 = [F(v) for v in cand['g3']]
        print(u'  --- кандидат #%s, g2(1,0) = %s' % (cand['i'], g2[0]))
        val, rows, gam, mval, zs = Fi.pair(g1, g2, g3, tag='c%s:' % cand['i'])
        minus = [r['place'] for r in rows if r['hilbert'] == -1]
        note(u'ПРОИЗВЕДЕНИЕ = %+d   (минусы: %s)' % (val, minus if minus else u'нет'))
        results.append({'i': cand['i'], 'g2': [str(v) for v in g2],
                        'g3': [str(v) for v in g3],
                        'gamma': [str(v) for v in gam],
                        'value': val, 'places': rows,
                        'minus_places': minus})
        if val == -1 and winner is None:
            winner = results[-1]
    check(winner is not None, u'найден g2 с НЕНУЛЕВЫМ спариванием <g1,g2>_CT = -1')

    g2 = [F(v) for v in winner['g2']]
    g3 = [F(v) for v in winner['g3']]
    print(u'\n[4] выбран кандидат #%s: <g1,g2>_CT = -1' % winner['i'])

    # ---------- независимость от выбора корня m ----------
    print(u'\n[5] независимость от выбора квадратного корня m (Замечание 3.2(vi))')
    for signs in [(1, 1, 1), (-1, 1, 1), (1, -1, 1), (1, 1, -1), (-1, -1, 1), (-1, -1, -1)]:
        v2, _, _, _, _ = Fi.pair(g1, g2, g3, signs=signs, verbose=False)
        check(v2 == winner['value'], u'знаки m = %s дают то же значение' % (signs,))
    note(u'шесть разных корней m -- одно и то же значение %+d' % winner['value'])

    # ---------- ОБРАТНЫЙ порядок аргументов ----------
    print(u'\n[6] пересчёт в ОБРАТНОМ порядке аргументов: <g2,g1>_CT')
    rev, rev_rows, _, _, _ = Fi.pair(g2, g1, g3, tag='rev:')
    check(rev == winner['value'],
          u'<g2,g1>_CT = <g1,g2>_CT = %+d (симметрия)' % winner['value'])

    print(u'\n[7] третий аргумент: <g1,g3>_CT должно совпасть (3.2(v))')
    p13, p13_rows, _, _, _ = Fi.pair(g1, g3, g2, tag='13:')
    check(p13 == winner['value'], u'<g1,g3>_CT = <g1,g2>_CT')

    # ---------- ELS всех трёх квартик ----------
    print(u'\n[8] всюду локальная разрешимость трёх квартик')
    els = []
    for g, tag in ((g1, 'g1'), (g2, 'g2'), (g3, 'g3')):
        rows, S = els_witnesses(g, Delta, tag)
        els.append({'quartic': tag, 'places': rows, 'places_checked': S})
        note(u'%s: свидетели на %s' % (tag, ', '.join(S)))
    note(u'на остальных p: g p-целая, p не делит Delta => хорошая редукция рода 1;')
    note(u'оценка Хассе p+1-2sqrt(p) > 0 даёт точку, Гензель поднимает её в Q_p.')

    # ---------- контроль ложных срабатываний ----------
    print(u'\n[9] контроль: классы 2-кручения E(Q) обязаны спариваться тривиально')
    controls = []
    for ctl in R.get('controls', []):
        gt = [F(v) for v in ctl['g']]
        Fi.check_quartic(gt, ctl['name'])
        zc = [squarefree_part(v) for v in z_comps(gt, phis)]
        for tr in ctl.get('triples', []):
            gg2 = [F(v) for v in tr['g2']]
            gg3 = [F(v) for v in tr['g3']]
            if gg2 != g2:
                continue
            v, rr, _, _, _ = Fi.pair(gt, gg2, gg3, tag='ctl:', verbose=False)
            check(v == 1, u'контроль %s (класс %s) спарился тривиально' % (ctl['name'], zc))
            controls.append({'name': ctl['name'], 'class': [str(x) for x in zc],
                             'paired_with': 'g2 #%s' % winner['i'], 'value': v})
            note(u'%s: класс %s, <.,g2> = %+d  (ожидалось +1)' % (ctl['name'], zc, v))
    check(len(controls) >= 2, u'контролей с тривиальным ответом не меньше двух')

    # ---------- итог ----------
    print(u'\n[10] итог')
    minus = winner['minus_places']
    note(u'<[g1],[g2]>_CT = %+d   (в Q/Z: 1/2)' % winner['value'])
    note(u'мест в произведении: %d, из них -1: %d (%s)'
         % (len(winner['places']), len(minus), ', '.join(minus)))
    note(u'ранг E(Q) НЕ использован ни разу.')

    res = {
        'pair': '(m,n)=(11,4)',
        'verifier': u'стандартный Python (fractions/math/random/json); '
                    u'ни Sage, ни PARI, ни кода проекта',
        's': '137/2', 'b': '132616',
        'target_squareclass': [1, 274, 274],
        'I': str(I), 'J': str(J), 'Delta': str(Delta),
        'phi_roots': [str(v) for v in phis],
        'g1': [str(v) for v in g1],
        'g2': [str(v) for v in g2],
        'g3': [str(v) for v in g3],
        'gamma': winner['gamma'],
        'a_equals_g2_at_1_0': str(g2[0]),
        'z_g1_class_verified': [1, 274, 274],
        'places': winner['places'],
        'places_with_minus_one': minus,
        'n_places': len(winner['places']),
        'pairing_value': winner['value'],
        'pairing_in_Q_mod_Z': '1/2',
        'reverse_order_pairing': rev,
        'pair_g1_g3': p13,
        'all_candidate_values': [{'i': r['i'], 'value': r['value'],
                                  'minus_places': r['minus_places']} for r in results],
        'sqrt_m_sign_independence': True,
        'els_witnesses': els,
        'torsion_controls': controls,
        'rank_bound_used': False,
        'checks_passed': len(CHECKS),
        'conclusion': (u'<[g1],[g2]>_CT = -1 != 0, поэтому класс delta=(1,274,274) не лежит '
                       u'в ядре спаривания Касселса-Тейта, а образ E(Q)/2E(Q) лежит в этом ядре; '
                       u'значит delta не приходит ни от какой рациональной точки E. '
                       u'Конечных рациональных точек на C_(11,4) нет; бесконечность исключена '
                       u'тем, что s=137/2 не квадрат. Следовательно C_(11,4)(Q) = пусто.'),
    }
    outp = HERE / 'independent_result_11_4_indep.json'
    outp.write_text(json.dumps(res, indent=2, ensure_ascii=False) + '\n')
    print(u'\nпроверок пройдено: %d' % len(CHECKS))
    print(u'результат записан в %s' % outp)
    print(u'\nВЫВОД: %s' % res['conclusion'])
    return res


if __name__ == '__main__':
    main()
