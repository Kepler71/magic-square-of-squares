#!/usr/bin/env python3
"""zres_core.py -- ядро прямой проверки (Claude, 26.09.2026, отдельная реализация; код из lift_mass/*.py
и review_20260926/p7_filter_check.py НЕ использовался и не читался).

Сетка f_ij = 1 + (i r + j s) z, i,j in {-1,0,1}.  Восемь линий: 3 строки, 3 столбца, 2 диагонали.
T = f(-1,0) f(1,1) f(0,-1),  L = f(-1,0) f(1,-1) f(0,1)  (как в lift/RESULT.md Codex).

Класс ненулевого x в Q_p* / Q_p*^2 (p нечётно) кодируется числом 0..3:
    бит 0 = v_p(x) mod 2,  бит 1 = 1, если единичная часть x/p^v -- невычет mod p.
Групповая операция -- XOR; тривиальный класс = 0.

image_Zp(p, r, s): ТОЧНОЕ множество пар ([T]_p, [L]_p) по всем z in Z_p с ненулевыми клетками,
при которых все восемь произведений -- квадраты в Q_p.  Метод -- дерево по z mod p^k:
  * клетка с (1+mu z_k) mod p^k != 0 имеет оценку v < k, её класс определён по z mod p^k;
  * если неопределённая клетка ровно одна (mu -- p-единица), то при z in z_k + p^k Z_p её значение
    пробегает ВСЁ p^k Z_p, поэтому перебираются все 4 класса (остальные клетки постоянны на узле);
  * если неопределённых клеток >= 2, узел делится на p потомков. Глубина ограничена: две клетки
    mu != mu', обе = 0 mod p^k, дают p^k | (mu - mu') (z -- единица), а 0 < |mu-mu'| <= 2(|r|+|s|).
  * клетки с ОДИНАКОВЫМ mu (сетка не инъективна: при 0<r<s это только r/s = 1/2, где r = s-r) -- одна
    величина, им даётся общий класс.
Никакой формулы для I_p не используется: классы девяти клеток и восемь произведений считаются явно.
"""
from math import gcd

GRID = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
LINES = ([[(i, j) for j in (-1, 0, 1)] for i in (-1, 0, 1)] +        # строки
         [[(i, j) for i in (-1, 0, 1)] for j in (-1, 0, 1)] +        # столбцы
         [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])   # диагонали
T_CELLS = [(-1, 0), (1, 1), (0, -1)]
L_CELLS = [(-1, 0), (1, -1), (0, 1)]
assert len(LINES) == 8 and all(len(l) == 3 for l in LINES)
KMAX = 40   # явная верхняя граница глубины дерева (реально <= log_p(1200)+1)


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [i for i in range(n + 1) if sieve[i]]


def factor(n):
    n = abs(n)
    out = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def chi_table(p):
    """chi[x] для x in 0..p-1 по критерию Эйлера (без ссылок на готовые таблицы)."""
    t = [0] * p
    for x in range(1, p):
        e = pow(x, (p - 1) // 2, p)
        assert e in (1, p - 1)
        t[x] = 1 if e == 1 else -1
    return t


def class_of_int(x, p, chi):
    """класс ненулевого целого x в Q_p (p нечётно)."""
    assert x != 0
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return (v & 1) | ((1 if chi[x % p] < 0 else 0) << 1)


def class_of_frac(num, den, p, chi):
    return class_of_int(num, p, chi) ^ class_of_int(den, p, chi)


def lines_ok(cls):
    for a, b, c in LINES:
        if cls[a] ^ cls[b] ^ cls[c]:
            return False
    return True


def tl(cls):
    cT = cls[T_CELLS[0]] ^ cls[T_CELLS[1]] ^ cls[T_CELLS[2]]
    cL = cls[L_CELLS[0]] ^ cls[L_CELLS[1]] ^ cls[L_CELLS[2]]
    return cT, cL


def image_Zp(p, r, s, chi, stats=None):
    """Точный образ ([T]_p,[L]_p) по z in Z_p (см. docstring модуля).
    Возвращает (image:set, info:dict). info['nontriv_cell'] -- встречалась ли прошедшая восемь условий
    конфигурация с нетривиальным классом хотя бы одной клетки (то есть подъём к 9 квадратам НЕ автоматический)."""
    mus = {c: c[0] * r + c[1] * s for c in GRID}
    image = set()
    nontriv_cell = False
    maxdepth = 1
    nodes = 0
    refined = 0
    stack = [(1, z0) for z0 in range(p)]
    while stack:
        k, zk = stack.pop()
        nodes += 1
        if k > maxdepth:
            maxdepth = k
        pk = p ** k
        det = {}
        und = []
        for c in GRID:
            x = (1 + mus[c] * zk) % pk
            if x == 0:
                und.append(c)
            else:
                v = 0
                while x % p == 0:
                    x //= p
                    v += 1
                det[c] = (v & 1) | ((1 if chi[x % p] < 0 else 0) << 1)
        und_mu = sorted(set(mus[c] for c in und))   # одинаковые mu (s=2r и т.п.) -- одна и та же величина
        if len(und_mu) >= 2:
            if k >= KMAX:
                raise RuntimeError("KMAX exceeded p=%d r=%d s=%d" % (p, r, s))
            # все неопределённые клетки -- с p-единичными mu (иначе клетка = 1 mod p)
            refined += 1
            for t in range(p):
                stack.append((k + 1, zk + t * pk))
            continue
        if und:
            assert und_mu[0] % p != 0
            opts = range(4)
        else:
            opts = (None,)
        for o in opts:
            for c in und:          # все клетки с этим mu равны между собой -> один класс
                det[c] = o
            if lines_ok(det):
                image.add(tl(det))
                if not nontriv_cell and any(det[c] for c in GRID):
                    nontriv_cell = True
        for c in und:
            del det[c]
    info = {"nodes": nodes, "refined": refined, "maxdepth": maxdepth, "nontriv_cell": nontriv_cell}
    if stats is not None:
        stats["nodes"] = stats.get("nodes", 0) + nodes
    return image, info


def squarefree_products(primes):
    out = [1]
    for q in primes:
        out += [d * q for d in out]
    return sorted(out)


def codex_lists(r, s):
    """D_-(r,s) (для d_T) и D_+(r,s) (для d_L) по определению теоремы (lift/RESULT.md):
    d>0 бесквадратно, p|d => p | (r-s) [соотв. r+s], p = 1 mod 4, d = 1 mod 24."""
    def D(n):
        ps = [q for q in factor(n) if q % 4 == 1]
        return [d for d in squarefree_products(ps) if d % 24 == 1]
    return D(r - s), D(r + s)
