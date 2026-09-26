"""Независимая библиотека (Claude, v4, 26.09.2026): точные квадратные классы в Q_p
для рациональных чисел. Никакого кода Codex/прежних прогонов не импортирует."""
from fractions import Fraction
from math import gcd, isqrt

IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]
# восемь линий арифметической сетки f_ij = 1 + (i r + j s) z
LINES = ([[(i, j) for j in IDX] for i in IDX] +          # строки (i фиксировано)
         [[(i, j) for i in IDX] for j in IDX] +          # столбцы (j фиксировано)
         [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
T_CELLS = [(-1, 0), (1, 1), (0, -1)]   # (1-rz)(1+(r+s)z)(1-sz)
L_CELLS = [(-1, 0), (1, -1), (0, 1)]   # (1-rz)(1+(r-s)z)(1+sz)


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    k = 3
    while k * k <= n:
        if n % k == 0:
            return False
        k += 2
    return True


def vp_int(n, p):
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def vp(x, p):
    x = Fraction(x)
    return vp_int(x.numerator, p) - vp_int(x.denominator, p)


def unit_part_mod(x, p, m):
    """единичная часть x (после деления на p^v) по модулю m (m -- степень p)."""
    x = Fraction(x)
    n, d = x.numerator, x.denominator
    while n % p == 0:
        n //= p
    while d % p == 0:
        d //= p
    return (n * pow(d, -1, m)) % m


def qp_class(x, p):
    """класс x в Q_p^*/Q_p^*2 как пара (v mod 2, e):
       p нечётно: e = 0, если единичная часть -- вычет mod p, иначе 1;
       p = 2: e = единичная часть mod 8 (1,3,5,7)."""
    x = Fraction(x)
    assert x != 0
    v = vp(x, p)
    if p == 2:
        return (v % 2, unit_part_mod(x, 2, 8))
    u = unit_part_mod(x, p, p)
    leg = pow(u, (p - 1) // 2, p)
    assert leg in (1, p - 1)
    return (v % 2, 0 if leg == 1 else 1)


def is_sq_qp(x, p):
    c = qp_class(x, p)
    return c == (0, 1) if p == 2 else c == (0, 0)


def cells(r, s, z):
    z = Fraction(z)
    return {(i, j): 1 + (i * r + j * s) * z for (i, j) in CELLS}


def prod_cells(f, cl):
    out = Fraction(1)
    for c in cl:
        out *= f[c]
    return out


def eight_products(f):
    return [prod_cells(f, ln) for ln in LINES]


def eight_ok_qp(f, p):
    return all(is_sq_qp(P, p) for P in eight_products(f))


def factor_int(n):
    n = abs(n)
    out = {}
    k = 2
    while k * k <= n:
        while n % k == 0:
            out[k] = out.get(k, 0) + 1
            n //= k
        k += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out
