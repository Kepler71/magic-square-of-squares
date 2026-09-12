# -*- coding: utf-8 -*-
r"""
validate_sieve_19_5.sage
========================
Независимая проверка корректности решета из search_H_19_5.sage.

Тест 1. Таблицы квадратичных вычетов: tab[p][v mod p][u mod p] должно
        совпадать с прямым вычислением kronecker(N(u,v), p) != -1.
Тест 2. Сквозная: наивный полный перебор (is_square для каждой пары)
        против решета на 0 <= u <= v <= 400 для (19,5) — множества
        найденных точек обязаны совпасть.
Тест 3. Контроль на паре (15,8), где точка t=1 известна: решето обязано
        её пропустить и точный тест обязан её найти.
Тест 4. Симметрии: N(u,v) = N(v,u) (t -> 1/t) и чётность по t (t -> -t).
"""
import numpy as np, random, sys


def make_params(m, n):
    """s = (m^2+n^2)/2 может быть полуцелым; работаем с 2s = m^2+n^2."""
    m = Integer(m); n = Integer(n)
    s2 = m ^ 2 + n ^ 2                 # = 2s
    return m, n, s2


def Nint(m, n, s2, u, v):
    """4 * v^6 * F0(u/v)F4(u/v)F8(u/v) — всегда целое, квадрат тогда же,
    когда квадратно само v^6 F0F4F8 (множитель 4 — квадрат)."""
    u = Integer(u); v = Integer(v)
    return 2 * (m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2) * s2 * (u ^ 2 + v ^ 2) * (n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2)


def build_table(m, n, s2, p):
    qr = np.zeros(p, dtype=np.uint8)
    for x in range(p):
        qr[(x * x) % p] = 1
    T = np.zeros((p, p), dtype=np.uint8)
    c = (2 * s2) % p
    for vr in range(p):
        for ur in range(p):
            val = ((m ** 2 * vr * vr + n ** 2 * ur * ur) % p) * c % p
            val = val * ((ur * ur + vr * vr) % p) % p
            val = val * ((n ** 2 * vr * vr + m ** 2 * ur * ur) % p) % p
            T[vr][ur] = qr[val]
    return T


def sieve_found(m, n, s2, HMAX, primes):
    tabs = {p: build_table(m, n, s2, p) for p in primes}
    out = []
    for v in range(1, HMAX + 1):
        L = v + 1
        alive = np.ones(L, dtype=np.uint8)
        for p in primes:
            row = tabs[p][v % p]
            np.logical_and(alive, np.tile(row, L // p + 1)[:L], out=alive)
        for u in np.nonzero(alive)[0]:
            u = int(u)
            if gcd(u, v) != 1:
                continue
            if Nint(m, n, s2, u, v).is_square():
                out.append((u, v))
    return out


def naive_found(m, n, s2, HMAX):
    out = []
    for v in range(1, HMAX + 1):
        for u in range(0, v + 1):
            if gcd(u, v) != 1:
                continue
            if Nint(m, n, s2, u, v).is_square():
                out.append((u, v))
    return out


PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

print("=== ТЕСТ 1: таблицы QR против прямого kronecker ===")
ok = True
for (m0, n0) in [(19, 5), (15, 8), (11, 4), (13, 8)]:
    m, n, s2 = make_params(m0, n0)
    for p in PRIMES:
        T = build_table(m, n, s2, p)
        for _ in range(400):
            u = random.randrange(0, 10 ** 6)
            v = random.randrange(1, 10 ** 6)
            N = Nint(m, n, s2, u, v)
            direct = (kronecker(N, p) != -1)
            tabval = bool(T[v % p][u % p])
            if direct != tabval:
                ok = False
                print("  РАСХОЖДЕНИЕ", (m0, n0), p, u, v, direct, tabval)
print("  результат:", "OK" if ok else "ПРОВАЛ")

print()
print("=== ТЕСТ 2: наивный перебор против решета, (19,5), v<=400 ===")
m, n, s2 = make_params(19, 5)
A = naive_found(m, n, s2, 400)
B = sieve_found(m, n, s2, 400, PRIMES)
print("  наивно найдено:", A)
print("  решетом найдено:", B)
print("  результат:", "OK" if A == B else "ПРОВАЛ")

print()
print("=== ТЕСТ 3: контроль (15,8), известна точка t=1, v<=60 ===")
m, n, s2 = make_params(15, 8)
A = naive_found(m, n, s2, 60)
B = sieve_found(m, n, s2, 60, PRIMES)
print("  наивно найдено:", A)
print("  решетом найдено:", B)
N11 = Nint(m, n, s2, 1, 1)
print("  N(1,1) =", N11, " квадрат?", N11.is_square(), " sqrt =", N11.sqrt())
print("  результат:", "OK" if (A == B and (1, 1) in B) else "ПРОВАЛ")

print()
print("=== ТЕСТ 4: симметрии ===")
m, n, s2 = make_params(19, 5)
ok = True
for _ in range(2000):
    u = random.randrange(0, 10 ** 7); v = random.randrange(1, 10 ** 7)
    if Nint(m, n, s2, u, v) != Nint(m, n, s2, v, u):
        ok = False
print("  N(u,v) == N(v,u)  (t -> 1/t):", "OK" if ok else "ПРОВАЛ")
R.<t> = QQ[]
F0 = m ^ 2 + n ^ 2 * t ^ 2; F4 = QQ(s2)/2 * (1 + t ^ 2); F8 = n ^ 2 + m ^ 2 * t ^ 2
print("  f(-t) == f(t)  (t -> -t):", (F0 * F4 * F8)(-t) == (F0 * F4 * F8)(t))
print("  F0(1/t)*t^2 == F8(t):", (F0.subs(t=1 / t) * t ^ 2).numerator() == F8)

print()
print("=== ТЕСТ 5: ПОЛНЫЙ конвейер решета (TILE+FANCY, 41 простое) ===")
print("    ключевой тест: конвейер обязан ПРОПУСТИТЬ настоящую точку.")
TILE_P = [int(11), int(17), int(37)]
FANCY_P = [int(q) for q in (79, 137, 113, 149, 89, 131, 167, 103, 97,
           109, 157, 127, 173, 107, 53, 3, 7, 19, 211, 233, 199, 41, 43,
           29, 59, 61, 67, 71, 73, 83, 101, 139, 151, 163, 179, 181,
           191, 197)]

def full_pipeline(m0, n0, HMAX):
    m, n, s2 = make_params(m0, n0)
    tabs_t = {p: build_table(m, n, s2, p) for p in TILE_P}
    tabs_f = {p: build_table(m, n, s2, p) for p in FANCY_P}
    M = int(TILE_P[0]*TILE_P[1]*TILE_P[2])
    tiled = {p: np.array([np.resize(tabs_t[p][r], M) for r in range(p)])
             for p in TILE_P}
    out = []
    for v in range(int(1), int(HMAX)+int(1)):
        comb = tiled[TILE_P[0]][v % TILE_P[0]] & tiled[TILE_P[1]][v % TILE_P[1]]
        comb = comb & tiled[TILE_P[2]][v % TILE_P[2]]
        if not comb.any():
            continue
        L = int(v) + int(1)
        idx = np.nonzero(np.resize(comb, L) if L > M else comb[:L])[0]
        for p in FANCY_P:
            if idx.size == 0:
                break
            idx = idx[np.nonzero(tabs_f[p][v % p][idx % p])[0]]
        for uu in idx:
            uu = int(uu)
            if gcd(uu, v) != 1:
                continue
            if Nint(m, n, s2, uu, v).is_square():
                out.append((uu, v))
    return out

for (m0, n0, HB, expect) in [(15, 8, 300, [(1, 1)]),
                             (19, 5, 3000, [])]:
    got = full_pipeline(m0, n0, HB)
    m, n, s2 = make_params(m0, n0)
    nav = naive_found(m, n, s2, HB) if HB <= 300 else None
    print("  (%d,%d), v<=%d: конвейер -> %s ; ожидалось %s ; наивно -> %s"
          % (m0, n0, HB, got, expect, nav))
    print("    результат:", "OK" if got == expect else "ПРОВАЛ")

print()
print("=== ТЕСТ 6: (19,5), v<=4000 — конвейер против 'слабого' решета ===")
print("    слабое решето = только TILE-простые + точный is_square;")
print("    совпадение показывает, что стадия FANCY ничего не теряет.")
def weak_pipeline(m0, n0, HMAX):
    m, n, s2 = make_params(m0, n0)
    tabs_t = {p: build_table(m, n, s2, p) for p in TILE_P}
    M = int(TILE_P[0]*TILE_P[1]*TILE_P[2])
    tiled = {p: np.array([np.resize(tabs_t[p][r], M) for r in range(p)])
             for p in TILE_P}
    out = []; nex = 0
    for v in range(int(1), int(HMAX)+int(1)):
        comb = tiled[TILE_P[0]][v % TILE_P[0]] & tiled[TILE_P[1]][v % TILE_P[1]]
        comb = comb & tiled[TILE_P[2]][v % TILE_P[2]]
        if not comb.any():
            continue
        L = int(v) + int(1)
        idx = np.nonzero(np.resize(comb, L) if L > M else comb[:L])[0]
        for uu in idx:
            uu = int(uu)
            if gcd(uu, v) != 1:
                continue
            nex += 1
            if Nint(m, n, s2, uu, v).is_square():
                out.append((uu, v))
    return out, nex
w, nex = weak_pipeline(19, 5, 4000)
g = full_pipeline(19, 5, 4000)
print("  слабое решето: точных проверок %d, найдено %s" % (nex, w))
print("  полный конвейер: найдено %s" % (g,))
print("  результат:", "OK" if w == g else "ПРОВАЛ")
