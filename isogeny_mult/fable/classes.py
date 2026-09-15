# Fable 15.09. Классы изоморфизма 3–4-клеточных множителей кривой наклона: сетка {−1,0,1}², λ = i·s + j·r.
# Кривая C_T ↔ 4-подмножество сетки (3-клеточное T ↔ T ∪ {(0,0)}); класс = орбита под гомотетиями u ↦ ±u + b, ±2u + b
# (проверено символьно над ℚ(k): generic_j.py — 54 класса, все внутриклассовые твисты тривиальны).
import itertools
from sage.all import *
import sys
sys.path.insert(0, '/home/kep/magicKube/census_genus1')
from g1census import weier

GRID = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]


def lam(g, r, s):
    return g[0] * s + g[1] * r


def canon(S):
    """каноническая форма 4-подмножества сетки: минимум по орбите гомотетий, сохраняющих сетку."""
    S = frozenset(S)
    best = None
    from fractions import Fraction as Fr
    for a in (1, -1, 2, -2, Fr(1, 2), Fr(-1, 2)):
        for bi in [Fr(x, 2) for x in range(-6, 7)]:
            for bj in [Fr(x, 2) for x in range(-6, 7)]:
                img = frozenset((a * i + bi, a * j + bj) for i, j in S)
                if all(abs(i) <= 1 and abs(j) <= 1 and i.denominator == 1 and j.denominator == 1 for i, j in img):
                    img = frozenset((int(i), int(j)) for i, j in img)
                    key = tuple(sorted(img))
                    if best is None or key < best: best = key
    return best


def all_classes():
    cl = {}
    for S in itertools.combinations(GRID, 4):
        cl.setdefault(canon(S), []).append(S)
    return cl


CLASSES = all_classes()
assert len(CLASSES) == 54 and sorted(len(v) for v in CLASSES.values()) == [1] * 5 + [2] * 38 + [4] * 10 + [5]
BIG = sorted([c for c, v in CLASSES.items() if len(v) >= 4], key=lambda c: -len(CLASSES[c]))   # 11 классов: 5, 4×10


def cells_of(S, r, s):
    """T (список λ) для 4-подмножества сетки; (0,0) ↦ клетка 1 (опускается)."""
    T = sorted(lam(g, r, s) for g in S if g != (0, 0))
    return T


def curve_of(S, r, s):
    T = cells_of(S, r, s)
    E, L, c0 = weier([QQ(c) for c in T])
    return E, L, c0, T


def sym_diff_size(S1, S2):
    return len(set(S1) ^ set(S2)) - (1 if (0, 0) in (set(S1) ^ set(S2)) else 0)   # число клеток в TΔT' (0 — не клетка)


if __name__ == '__main__':
    for c in BIG:
        print(len(CLASSES[c]), CLASSES[c])
