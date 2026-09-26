# -*- coding: utf-8 -*-
"""s2: контроль теоремы A (нет универсальных мультипликативных соотношений между образами половинок).
Берём точки X(F_p) (девять корней r_ij, r_ij^2 = a+ib+jc), считаем символы Лежандра
  16 delta-координат половинок (по 2 на каждую из 8 линий: x и x-t),
  8 «базовых» линейных форм b, c, b+-c, b+-2c, 2b+-c,
  9 корней r_ij,
  константу.
Если бы некоторое произведение этих функций было (константа)*(квадрат) в поле функций K = Q(X),
его символ был бы постоянен на точках X(F_p) (вне конечного числа кривых). Полный F_2-ранг матрицы
символов => ни одной такой связи среди проверяемых классов нет.
Санитарный контроль: добавляем третью координату x+t одной линии -- ранг не растёт (тавтология x*(x-t)*(x+t) = квадрат).
"""
import random, json, sys
sys.path.insert(0, '.')
from grid import IDX, lines_of, sums, legendre

random.seed(20260925)


def sample_point(p):
    while True:
        r00, r10, r01 = (random.randrange(1, p) for _ in range(3))
        a = r00 * r00 % p
        b = (r10 * r10 - a) % p
        c = (r01 * r01 - a) % p
        if b == 0 or c == 0:
            continue
        cells = {(i, j): (a + i * b + j * c) % p for i in IDX for j in IDX}
        if any(legendre(v, p) != 1 for v in cells.values()):
            continue
        r = {}
        for (i, j), v in cells.items():
            s = pow(v, (p + 1) // 4, p) if p % 4 == 3 else None
            if s is None:
                # общий случай -- перебор не нужен: p % 4 == 3 выбираем ниже
                raise RuntimeError
            r[(i, j)] = s if random.random() < 0.5 else (p - s)
        r[(0, 0)] = r00 if random.random() < 0.5 else p - r00
        return a, b, c, r


def rank_f2(rows):
    """rows -- список целых (битовые маски). Ранг над F_2."""
    basis = []
    for x in rows:
        for bvec in basis:
            x = min(x, x ^ bvec)
        if x:
            basis.append(x)
    return len(basis)


def row_bits(p, a, b, c, r, extra_tautology=False):
    vals = []
    for name, tl, (al, ga, be) in lines_of(r):
        A, B, C = sums(al, ga, be)
        vals += [A * B % p, A * C % p]
        if extra_tautology and name == "Eb_row+0":
            vals.append(B * C % p)
    base = [b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c]
    vals += [x % p for x in base]
    vals += [r[(i, j)] for i in IDX for j in IDX]
    bits = []
    for x in vals:
        L = legendre(x, p)
        if L == 0:
            return None
        bits.append(0 if L == 1 else 1)
    bits.append(1)  # константа
    mask = 0
    for k, bt in enumerate(bits):
        mask |= bt << k
    return mask, len(bits)


def transpose_rank(masks, ncols):
    # ранг матрицы = ранг по столбцам; считаем по строкам (эквивалентно)
    return rank_f2(masks)


out = {}
for p in (10007, 20011, 40039, 65003):  # все p = 3 mod 4 (простой квадратный корень)
    assert p % 4 == 3
    pts = []
    while len(pts) < 400:
        a, b, c, r = sample_point(p)
        rb = row_bits(p, a, b, c, r)
        if rb is None:
            continue
        pts.append((a, b, c, r, rb))
    ncols = pts[0][4][1]
    rk = rank_f2([x[4][0] for x in pts])
    rk_taut = rank_f2([row_bits(p, *x[:4], extra_tautology=True)[0] for x in pts])
    ncols_taut = row_bits(p, *pts[0][:4], extra_tautology=True)[1]
    # только E_b и E_c (12 координат) + константа
    def sub_rank(sel):
        ms = []
        for x in pts:
            m_ = x[4][0]
            mm = 0
            for k, idx in enumerate(sel):
                mm |= ((m_ >> idx) & 1) << k
            ms.append(mm)
        return rank_f2(ms)
    sel_bc = list(range(0, 12)) + [ncols - 1]
    out[p] = {"points": len(pts), "columns": ncols, "rank": rk,
              "columns_with_tautology": ncols_taut, "rank_with_tautology": rk_taut,
              "Eb+Ec_only_columns": len(sel_bc), "Eb+Ec_only_rank": sub_rank(sel_bc)}
    print(p, out[p], flush=True)

json.dump(out, open("s2_fp_independence.json", "w"), indent=1)
