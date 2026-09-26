# -*- coding: utf-8 -*-
"""
Скептическая проверка Шага 3 теоремы Codex (p = 2, 3) -- независимый код (Claude, 2026-09-26).

Локальное утверждение L_p:  b, c in Q_p^*, все девять клеток f_ij = 1 + i b + j c ненулевые,
все ВОСЕМЬ произведений по линиям арифметической сетки -- квадраты в Q_p  ==>
min(v_p(b), v_p(c)) >= t_p,  t_2 = 3, t_3 = 1.
(Для b = r z, c = s z с взаимно простыми r, s: min(v(b), v(c)) = v(z).)

Метод: узел = (m, k, B0, C0): b = B/p^m, c = C/p^m, B = B0, C = C0 mod p^k.
m = 0: min(v(b), v(c)) >= 0;  m >= 1: ровно min(v(b), v(c)) = -m (не оба B, C делятся на p).
Для каждой клетки N = p^m + iB + jC mod p^k вычисляется МАСКА -- множество квадратных классов,
которые клетка может принять при каком-либо подъёме (надмножество, т.е. проверка корректна).
Узел ПРОВАЛ, если нет присваивания классов из масок, при котором все 8 линий -- квадраты.
Решатель "есть присваивание" перебирает ВСЕ 8-наборы классов нецентральных клеток
только через независимую проверку линий (см. solve_table); параметризация S,U Codex
НЕ используется, а проверяется отдельно (selftest_SU).

Класс в Q_2^*/Q_2^*^2 = F_2^3: бит0 = v mod 2, бит1 = [u = 3 mod 4], бит2 = [u = +-3 mod 8].
Класс в Q_3^*/Q_3^*^2 = F_2^2: бит0 = v mod 2, бит1 = [u = 2 mod 3].
"""
import itertools

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {ij: n for n, ij in enumerate(CELLS)}
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)]
         + [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)]
         + [[IDX[(t, t)] for t in (-1, 0, 1)], [IDX[(t, -t)] for t in (-1, 0, 1)]])
CENTER = IDX[(0, 0)]
NONC = [n for n in range(9) if n != CENTER]
GBITS = {2: 3, 3: 2}
THR = {2: 3, 3: 1}


def vp_int(n, p, cap):
    v = 0
    while v < cap and n % p == 0:      # явная граница
        n //= p
        v += 1
    return v


def exact_class_int(n, p):
    """точный класс ненулевого целого n в Q_p^*/Q_p^*^2"""
    assert n != 0
    v = 0
    for _ in range(10000):            # явная граница
        if n % p:
            break
        n //= p
        v += 1
    else:
        raise RuntimeError('valuation cap')
    if p == 2:
        u = n % 8
        return (v & 1) | ((u % 4 == 3) << 1) | ((u in (3, 5)) << 2)
    u = n % 3
    return (v & 1) | ((u == 2) << 1)


def exact_class_frac(x, p):
    return exact_class_int(x.numerator, p) ^ exact_class_int(x.denominator, p)


def mask_of_residue(N0, k, p, m):
    """маска классов клетки x = N/p^m, если известно только N = N0 mod p^k"""
    G = 1 << GBITS[p]
    q = p ** k
    N0 %= q
    if N0 == 0:
        mask = (1 << G) - 1
    else:
        w = vp_int(N0, p, k)
        u = N0 // p ** w
        if p == 2:
            prec = k - w          # известно u mod 2^prec
            if prec >= 3:
                c = (w & 1) | ((u % 4 == 3) << 1) | ((u % 8 in (3, 5)) << 2)
                mask = 1 << c
            elif prec == 2:
                c0 = (w & 1) | ((u % 4 == 3) << 1)
                mask = (1 << c0) | (1 << (c0 | 4))
            else:  # prec == 1
                mask = 0
                for eo in range(4):
                    mask |= 1 << ((w & 1) | (eo << 1))
        else:
            c = (w & 1) | ((u % 3 == 2) << 1)
            mask = 1 << c
    if m & 1:   # деление на p^m меняет бит0 при нечётном m
        nm = 0
        for c in range(G):
            if mask >> c & 1:
                nm |= 1 << (c ^ 1)
        mask = nm
    return mask


def build_solution_list(p):
    """все 8-наборы классов нецентральных клеток (центр = 1, класс 0), при которых
    все 8 линий -- квадраты; прямой перебор G^8 с отсечением по линиям."""
    G = 1 << GBITS[p]
    sols = []
    # порядок присваивания: строки сверху вниз, отсечение по каждой полностью заданной линии
    order = NONC
    assign = [None] * 9
    assign[CENTER] = 0
    lines_done_at = {}
    for L in LINES:
        last = max(order.index(n) for n in L if n != CENTER)
        lines_done_at.setdefault(last, []).append(L)

    def rec(pos):
        if pos == len(order):
            sols.append(tuple(assign))
            return
        n = order[pos]
        for c in range(G):
            assign[n] = c
            ok = True
            for L in lines_done_at.get(pos, []):
                if assign[L[0]] ^ assign[L[1]] ^ assign[L[2]]:
                    ok = False
                    break
            if ok:
                rec(pos + 1)
        assign[n] = None
    rec(0)
    return sols


def selftest_SU(p, sols):
    """проверка Шага 1 Codex: решения = { S U S+U / U 0 U / S+U U S }"""
    G = 1 << GBITS[p]
    su = set()
    for S in range(G):
        for U in range(G):
            M = {(-1, -1): S, (-1, 0): U, (-1, 1): S ^ U,
                 (0, -1): U, (0, 0): 0, (0, 1): U,
                 (1, -1): S ^ U, (1, 0): U, (1, 1): S}
            su.add(tuple(M[ij] for ij in CELLS))
    return set(sols) == su, len(sols), len(su)


def brute_all_G8(p):
    """полностью независимый перебор всех G^8 без отсечения (контроль решателя)"""
    G = 1 << GBITS[p]
    out = set()
    for tup in itertools.product(range(G), repeat=8):
        a = list(tup[:4]) + [0] + list(tup[4:])
        if all((a[L[0]] ^ a[L[1]] ^ a[L[2]]) == 0 for L in LINES):
            out.add(tuple(a))
    return out
