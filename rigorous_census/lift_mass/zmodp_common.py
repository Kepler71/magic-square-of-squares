#!/usr/bin/env python3
"""zmodp_common.py -- движок для проверки массового утверждения REVIEW_LIFT_CLASSES §5.

Написано заново (Claude, 26.09.2026); p7_filter_check.py, t3c и файлы lift_mass_direct.py /
lm_common.py из этого же каталога НЕ читались.

Клетки: f_ij = 1 + (i r + j s) z,  i,j in {-1,0,1}.
Класс в Q_p^*/Q_p^*2 (p нечётно) кодируется двумя битами:
    bit0 = v_p mod 2,  bit1 = 1  <=>  единичная часть -- невычет mod p.
Групповая операция -- XOR; квадрат <=> код 0.

Восемь условий: произведения по 3 строкам (i = const), 3 столбцам (j = const) и 2 диагоналям
арифметической сетки -- квадраты.
T = (1-rz)(1+(r+s)z)(1-sz)  = клетки (-1,0),(1,1),(0,-1)
L = (1-rz)(1+(r-s)z)(1+sz)  = клетки (-1,0),(1,-1),(0,1)

local_image_Zp(r, s, p) -- ТОЧНОЕ множество векторов классов девяти клеток при z in Z_p
(все клетки != 0), для которых восемь произведений -- квадраты в Q_p. Алгоритм:
  * ветви z = a mod p^k, начиная с k = 1 (a = вычет z mod p);
  * клетка с 1 + mu*a != 0 mod p^k имеет постоянный класс на всей ветви
    (v_p < k, единичная часть известна mod p^(k-v) >= p);
  * если ровно одна клетка = 0 mod p^k, то на ветви она пробегает p^k Z_p \\ {0}
    (1+mu z = p^k (m + mu t), mu -- p-единица), т.е. ВСЕ 4 класса, остальные постоянны
    -- ветвь даёт ровно 4 вектора, отбираем по восьми условиям;
  * если нулей >= 2, спускаемся на уровень k+1. Два нуля mu1 != mu2 на уровне k
    требуют p^k | (mu1 - mu2), |mu1 - mu2| <= 2(r+s), поэтому глубина ограничена;
    явная граница KMAX (исключение при превышении).
Никакой формулы I_p здесь нет: классы каждой клетки считаются напрямую.
"""

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {c: k for k, c in enumerate(CELLS)}
LINES8 = (
    [[(i, -1), (i, 0), (i, 1)] for i in (-1, 0, 1)]
    + [[(-1, j), (0, j), (1, j)] for j in (-1, 0, 1)]
    + [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]]
)
LINES8_IDX = [tuple(IDX[c] for c in ln) for ln in LINES8]
T_IDX = tuple(IDX[c] for c in [(-1, 0), (1, 1), (0, -1)])
L_IDX = tuple(IDX[c] for c in [(-1, 0), (1, -1), (0, 1)])
KMAX = 14

_QR_CACHE = {}


def qr_table(p):
    t = _QR_CACHE.get(p)
    if t is None:
        t = [False] * p
        for x in range(1, p):
            t[x * x % p] = True
        _QR_CACHE[p] = t
    return t


def is_prime(n):
    if n < 2:
        return False
    d = 2
    while d * d <= n:
        if n % d == 0:
            return False
        d += 1
    return True


def factor(n):
    n = abs(n)
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def codex_list(n):
    """D(n) = {d>0 squarefree: p|d => p|n, p = 1 mod 4, d = 1 mod 24}, n = r-s или r+s."""
    ps = [q for q in factor(n) if q % 4 == 1]
    ds = []
    for mask in range(1 << len(ps)):
        d = 1
        for k, q in enumerate(ps):
            if mask >> k & 1:
                d *= q
        if d % 24 == 1:
            ds.append(d)
    return sorted(ds)


def unit_code(d, p):
    """Класс p-единицы d в Q_p^*/Q_p^*2."""
    assert d % p != 0
    return 0 if qr_table(p)[d % p] else 2


def cell_mus(r, s):
    return [i * r + j * s for (i, j) in CELLS]


def local_image_Zp(r, s, p, stats=None):
    """Точное множество 9-векторов классов (кортежи), z in Z_p, клетки != 0, 8 условий в Q_p."""
    assert p % 2 == 1
    qr = qr_table(p)
    mus = cell_mus(r, s)
    out = set()
    # стек ветвей (a, k): a = z mod p^k. Верхняя граница: глубина <= KMAX, ширина <= p на уровень.
    stack = [(a0, 1) for a0 in range(p)]
    nodes = 0
    while stack:
        a, k = stack.pop()
        nodes += 1
        pk = p ** k
        cls = [0] * 9
        zeros = []
        for idx in range(9):
            val = (1 + mus[idx] * a) % pk
            if val == 0:
                zeros.append(idx)
                continue
            v = 0
            while val % p == 0:
                val //= p
                v += 1
            cls[idx] = (v & 1) | (0 if qr[val % p] else 2)
        # клетки с одинаковым mu тождественно равны (бывает только при s = 2r, т.е. r/s = 1/2):
        # у них общий свободный класс; спускаемся, только если среди нулей >= 2 РАЗНЫХ mu.
        zero_mus = set(mus[idx] for idx in zeros)
        if len(zero_mus) >= 2:
            if k >= KMAX:
                raise RuntimeError("KMAX exceeded r=%d s=%d p=%d" % (r, s, p))
            for t in range(p):
                stack.append((a + t * pk, k + 1))
            continue
        if zeros:
            opts = (0, 1, 2, 3)
        else:
            opts = (None,)
        for o in opts:
            if o is not None:
                for zi in zeros:
                    cls[zi] = o
            ok = True
            for (x, y, w) in LINES8_IDX:
                if cls[x] ^ cls[y] ^ cls[w]:
                    ok = False
                    break
            if ok:
                out.add(tuple(cls))
    if stats is not None:
        stats["nodes"] = stats.get("nodes", 0) + nodes
    return out


def tl_of(vec):
    return (vec[T_IDX[0]] ^ vec[T_IDX[1]] ^ vec[T_IDX[2]],
            vec[L_IDX[0]] ^ vec[L_IDX[1]] ^ vec[L_IDX[2]])


def primes_3mod4(lo, hi):
    return [q for q in range(lo, hi + 1) if is_prime(q) and q % 4 == 3]


# ---------- независимая (от движка) точная p-адическая арифметика на рациональных числах ----------
from fractions import Fraction


def vp_int(n, p):
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v, n


def qp_class_rational(x, p):
    """Класс ненулевого рационального x в Q_p^*/Q_p^*2 (p нечётно), код как выше."""
    x = Fraction(x)
    assert x != 0
    vn, un = vp_int(x.numerator, p)
    vd, ud = vp_int(x.denominator, p)
    v = vn - vd
    u = (un * pow(ud, -1, p)) % p
    return (v & 1) | (0 if qr_table(p)[u] else 2)
