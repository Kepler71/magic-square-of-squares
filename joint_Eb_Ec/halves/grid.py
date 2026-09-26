# -*- coding: utf-8 -*-
"""Общий модуль: сетка корней r[i][j]^2 = a + i*b + j*c, половинки точек на E_b, E_c, E_{b+c}, E_{b-c},
образы половинок при 2-спуске и их локальные классы в Q_p^*/Q_p^{*2}.

Соглашения (см. NOTE.md, §1):
  линия L с тройкой корней (alpha, gamma, beta), alpha^2 = x - t, gamma^2 = x, beta^2 = x + t;
  суммы A = alpha+gamma, B = gamma+beta, C = alpha+beta;
  половинка Q = (u, v) = (A*B, A*B*C), 2Q = P = (x, +-alpha*gamma*beta);
  delta_t(Q) = (u, u - t, u + t) = (A*B, A*C, B*C)  mod квадраты.
Линии: E_b -- строки j (i меняется), E_c -- столбцы i (j меняется),
       E_{b+c} -- диагональ (-1,-1),(0,0),(1,1); E_{b-c} -- антидиагональ (-1,1),(0,0),(1,-1).
"""
from fractions import Fraction
import random

IDX = (-1, 0, 1)


def lines_of(r):
    """r: dict (i,j) -> корень. Возвращает список (имя, t-метка, (alpha, gamma, beta))."""
    out = []
    for j in IDX:
        out.append(("Eb_row%+d" % j, "b", (r[(-1, j)], r[(0, j)], r[(1, j)])))
    for i in IDX:
        out.append(("Ec_col%+d" % i, "c", (r[(i, -1)], r[(i, 0)], r[(i, 1)])))
    out.append(("Ebpc_diag", "b+c", (r[(-1, -1)], r[(0, 0)], r[(1, 1)])))
    out.append(("Ebmc_anti", "b-c", (r[(-1, 1)], r[(0, 0)], r[(1, -1)])))
    return out


def sums(al, ga, be):
    return (al + ga, ga + be, al + be)


def delta_coords(al, ga, be):
    A, B, C = sums(al, ga, be)
    return (A * B, A * C, B * C)


# ---------------------------------------------------------------- p-адика
def vp_int(n, p):
    if n == 0:
        return 10 ** 9
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def legendre(n, p):
    n %= p
    if n == 0:
        return 0
    return 1 if pow(n, (p - 1) // 2, p) == 1 else -1


def sqrt_mod_p(n, p):
    """Тонелли–Шэнкс."""
    n %= p
    if n == 0:
        return 0
    assert legendre(n, p) == 1
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
    q, s = p - 1, 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while legendre(z, p) != -1:
        z += 1
    m, c, t, rr = s, pow(z, q, p), pow(n, q, p), pow(n, (q + 1) // 2, p)
    while t != 1:
        i, tt = 0, t
        while tt != 1:
            tt = tt * tt % p
            i += 1
        bb = pow(c, 1 << (m - i - 1), p)
        m, c, t, rr = i, bb * bb % p, t * bb * bb % p, rr * bb % p
    return rr


def padic_sqrt(n, p, N):
    """Корень из целого n в Z_p по модулю p^N (p нечётно). Возвращает (x, ok).
    Допускается v_p(n) чётная; корень = p^(v/2) * sqrt(unit)."""
    if n == 0:
        return 0, True
    v = vp_int(n, p)
    if v % 2:
        return None, False
    u = n // p ** v
    if legendre(u, p) != 1:
        return None, False
    mod = p ** N
    x = sqrt_mod_p(u, p)
    # Гензель: x <- x - (x^2-u)/(2x)
    k = 1
    while k < N:
        k = min(2 * k, N)
        m = p ** k
        x = (x - (x * x - u) * pow(2 * x, -1, m)) % m
    return (x * p ** (v // 2)) % mod, True


def local_class(x, p, N):
    """x -- целое, представляющее p-адическое число по модулю p^N (x != 0 mod p^N).
    Возвращает (v mod 2, символ Лежандра единичной части) -- класс в Q_p^*/Q_p^{*2} (p нечётно)."""
    x %= p ** N
    v = vp_int(x, p)
    if v >= N - 2:
        raise ValueError("недостаточная точность")
    u = x // p ** v
    return (v % 2, legendre(u, p))


def is_padic_square(x, p, N):
    vv, s = local_class(x, p, N)
    return vv == 0 and s == 1


# ---------------------------------------------------------------- 2-адика (только проверка квадратности)
def is_2adic_square_frac(q):
    q = Fraction(q)
    if q == 0:
        return True
    num, den = q.numerator, q.denominator
    v = vp_int(abs(num), 2) - vp_int(den, 2)
    if v % 2:
        return False
    u_num = num // 2 ** vp_int(abs(num), 2)
    u_den = den // 2 ** vp_int(den, 2)
    return (u_num * u_den) % 8 == 1


# ---------------------------------------------------------------- рациональные классы
def sqfree_class(q):
    """Бесквадратный представитель класса рационального q != 0 (целое со знаком)."""
    q = Fraction(q)
    n = q.numerator * q.denominator
    sgn = -1 if n < 0 else 1
    n = abs(n)
    out = 1
    d = 2
    while d * d <= n:
        e = 0
        while n % d == 0:
            n //= d
            e += 1
        if e % 2:
            out *= d
        d += 1 if d == 2 else 2
    out *= n
    return sgn * out


def rat_local_class(q, p):
    """Класс рационального q в Q_p^*/Q_p^{*2}, p нечётно: (v mod 2, символ единичной части)."""
    q = Fraction(q)
    num, den = q.numerator, q.denominator
    vn, vd = vp_int(abs(num), p), vp_int(den, p)
    un, ud = num // p ** vn, den // p ** vd
    return ((vn - vd) % 2, legendre(un * ud, p))
