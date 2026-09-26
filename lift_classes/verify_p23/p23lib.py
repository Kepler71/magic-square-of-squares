# -*- coding: utf-8 -*-
"""
Независимая (Claude) проверка шага 3 теоремы Codex об остаточных классах:
p = 2, 3; из квадратности ВОСЬМИ произведений арифметической сетки
f_ij = 1 + i*b + j*c  (b = r z, c = s z)  следует min(v_p(b), v_p(c)) >= 3 (p=2), >= 1 (p=3).

Шар: b = B / p^m, c = C / p^m, где B, C заданы по модулю p^K, K = k + m
(абсолютная точность b, c равна p^k).  Клетка x = 1 + i b + j c,
X = p^m x = p^m + i B + j C  (mod p^K).

Квадратный класс в Q_p^*/Q_p^*^2 кодируем парой (v mod 2, u mod M), M = 8 (p=2), 3 (p=3);
тривиален <=> (0, 1).  Класс клетки определён по шару, если X != 0 mod p^K
и точности единичной части хватает: K - v_p(X) >= 3 (p=2), >= 1 (p=3).
Все циклы конечны: границы заданы явно.
"""

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]          # индекс 0..8
IDX = {ij: n for n, ij in enumerate(CELLS)}
# восемь линий арифметической сетки: 3 строки (i фикс.), 3 столбца (j фикс.), 2 диагонали
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)] +
         [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)] +
         [[IDX[(t, t)] for t in (-1, 0, 1)], [IDX[(t, -t)] for t in (-1, 0, 1)]])
assert len(LINES) == 8
OPP = [IDX[(-i, -j)] for (i, j) in CELLS]                           # противоположная клетка

UMOD = {2: 8, 3: 3}
NEED = {2: 3, 3: 1}          # сколько цифр единичной части нужно для класса
THR = {2: 3, 3: 1}           # утверждаемый порог min v_p(b), v_p(c)


def vp_int(n, p, cap):
    """v_p(n) для n != 0, не больше cap (явная граница цикла)."""
    v = 0
    while v < cap and n % p == 0:
        n //= p
        v += 1
    return v


def cell_class(X, m, p, K):
    """Класс клетки x = X/p^m, X известно mod p^K. None, если не определён."""
    PK = p ** K
    X %= PK
    if X == 0:
        return None
    v = vp_int(X, p, K)
    if K - v < NEED[p]:
        return None
    u = (X // p ** v) % UMOD[p]
    return ((v - m) % 2, u)


def mul(c1, c2, p):
    return ((c1[0] + c2[0]) % 2, (c1[1] * c2[1]) % UMOD[p])


TRIV = {2: (0, 1), 3: (0, 1)}


def line_class(cls, line, p):
    acc = (0, 1)
    for n in line:
        if cls[n] is None:
            return None
        acc = mul(acc, cls[n], p)
    return acc


def classify_naive(cls, p):
    """EXCL: есть определённая линия с нетривиальным классом;
       PASS: все 9 классов определены и все 8 линий тривиальны;
       UND: иначе."""
    allknown = True
    for line in LINES:
        lc = line_class(cls, line, p)
        if lc is None:
            allknown = False
        elif lc != TRIV[p]:
            return 'EXCL'
    if allknown and all(c is not None for c in cls):
        return 'PASS'
    return 'UND'


def classify_pairs(cls, p):
    """Использует лишь следствие 4 центральных линий: класс клетки = класс противоположной.
       Если в паре оба класса определены и различны -> EXCL (центральная линия определена).
       Неопределённый класс заменяется классом партнёра; затем проверяются все 8 линий.
       PASS = прошёл все НЕОБХОДИМЫЕ условия при полностью заполненных классах."""
    filled = list(cls)
    for n in range(9):
        a, b = cls[n], cls[OPP[n]]
        if a is not None and b is not None and a != b:
            return 'EXCL'
        if a is None:
            filled[n] = b
    return classify_naive(filled, p)


def ball_cells(B, C, m, p, K):
    pm = p ** m
    return [cell_class(pm + i * B + j * C, m, p, K) for (i, j) in CELLS]


def minval_ball(B, C, m, p, K):
    """min(v(b), v(c)) шара; если B = C = 0 mod p^K, возвращает K - m (нижняя оценка, 'большая')."""
    PK = p ** K
    vs = []
    for Y in (B % PK, C % PK):
        vs.append(K if Y == 0 else vp_int(Y, p, K))
    return min(vs) - m
