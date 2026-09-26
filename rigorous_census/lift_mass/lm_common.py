# lm_common.py — общие определения для проверки массового утверждения §5 REVIEW_LIFT_CLASSES_2026-09-26.
# Claude (подагент), 26.09.2026. Скрипт lift_classes/review_20260926/p7_filter_check.py НЕ открывался и НЕ использовался.
#
# Клетки: f_ij = 1 + mu_ij * z,  mu_ij = i*r + j*s,  i, j in {-1, 0, 1}  (центр нормирован к 1).
# Восемь линий АРИФМЕТИЧЕСКОЙ сетки: 3 строки (i фиксировано), 3 столбца (j фиксировано), 2 диагонали.
# T = (1 - r z)(1 + (r+s) z)(1 - s z)  = f(-1,0) f(1,1) f(0,-1)
# L = (1 - r z)(1 + (r-s) z)(1 + s z)  = f(-1,0) f(1,-1) f(0,1)
# Класс в Q_p*/Q_p*^2 (p нечётно) кодируется числом 2*e + c: e = v_p mod 2, c = 0 если единичная часть — вычет mod p,
# иначе 1. Групповая операция — XOR. Код 0 = квадрат; код 1 = единичный невычет.
from itertools import combinations, product
from math import gcd, prod

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {c: k for k, c in enumerate(CELLS)}
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)]
         + [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)]
         + [[IDX[(-1, -1)], IDX[(0, 0)], IDX[(1, 1)]], [IDX[(-1, 1)], IDX[(0, 0)], IDX[(1, -1)]]])
T_CELLS = [IDX[(-1, 0)], IDX[(1, 1)], IDX[(0, -1)]]
L_CELLS = [IDX[(-1, 0)], IDX[(1, -1)], IDX[(0, 1)]]
assert len(LINES) == 8 and all(len(set(l)) == 3 for l in LINES)


def mus(r, s):
    return [i * r + j * s for (i, j) in CELLS]


# самопроверка определений T, L через mu
def _selftest():
    r, s = 5, 11
    m = mus(r, s)
    assert [m[k] for k in T_CELLS] == [-r, r + s, -s]
    assert [m[k] for k in L_CELLS] == [-r, r - s, s]
_selftest()


def primes_upto(n):
    sv = bytearray([1]) * (n + 1); sv[0:2] = b'\x00\x00'
    for i in range(2, int(n ** 0.5) + 1):
        if sv[i]: sv[i * i::i] = bytearray(len(sv[i * i::i]))
    return [i for i in range(n + 1) if sv[i]]


def factor(n):
    n = abs(n); f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def codex_list(n):
    """D(n) = {d>0 бесквадратное: p|d => p|n, p = 1 mod 4, d = 1 mod 24} (списки Codex, lift/RESULT.md)."""
    ps = sorted(p for p in factor(n) if p % 4 == 1)
    out = []
    for k in range(len(ps) + 1):
        for comb in combinations(ps, k):
            d = prod(comb)
            if d % 24 == 1: out.append(d)
    return sorted(out)


def chi_table(p):
    """chi[x] = 0, если x — ненулевой квадратичный вычет mod p, 1 — невычет; chi[0] = None."""
    sq = set(x * x % p for x in range(1, p))
    return [None] + [0 if x in sq else 1 for x in range(1, p)]


def relaxed_image(r, s, p, chi):
    """Верхняя оценка локального образа (класс[T], класс[L]) в Q_p при z в Z_p.
    Прямо по клеткам: для каждого вычета z0 mod p клетка 1 + mu z0 != 0 mod p — единица, её класс = (0, chi);
    клетка = 0 mod p — класс неизвестен, перебираются все 4 класса группы Q_p*/Q_p*^2.
    Отбираются наборы, в которых произведения по всем восьми линиям — квадраты.
    Возвращает (множество пар кодов (T, L), число z0, давших хоть одну точку, число z0 с нулевыми клетками)."""
    m = [x % p for x in mus(r, s)]
    img = set(); good = 0; zero_z0 = 0
    for z0 in range(p):
        known = [0] * 9; unk = []
        for k in range(9):
            v = (1 + m[k] * z0) % p
            if v: known[k] = chi[v]
            else: unk.append(k)
        if unk: zero_z0 += 1
        hit = False
        for combo in product(range(4), repeat=len(unk)):
            cl = known[:]
            for k, c in zip(unk, combo): cl[k] = c
            ok = True
            for a, b, c in LINES:
                if cl[a] ^ cl[b] ^ cl[c]:
                    ok = False; break
            if ok:
                hit = True
                img.add((cl[T_CELLS[0]] ^ cl[T_CELLS[1]] ^ cl[T_CELLS[2]],
                         cl[L_CELLS[0]] ^ cl[L_CELLS[1]] ^ cl[L_CELLS[2]]))
        good += hit
    return img, good, zero_z0


def exact_class_int(num, p, shift=0):
    """Класс числа num / p^shift (num — ненулевое целое) в Q_p*/Q_p*^2, код 2e + c."""
    assert num != 0
    v = 0
    while num % p == 0: num //= p; v += 1
    c = 0 if pow(num % p, (p - 1) // 2, p) == 1 else 1
    return 2 * ((v - shift) % 2) + c
