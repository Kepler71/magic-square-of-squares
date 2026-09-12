"""
sieve_uv.py -- быстрый просеивающий поиск рациональных точек на
    H : Y^2 = F0(t)*F4(t)*F8(t),   t = u/v,
    F0 = m^2 + n^2 t^2,  F4 = ((m^2+n^2)/2)(1+t^2),  F8 = n^2 + m^2 t^2.

Целая форма (домножение на 4, Y -> Y/2):
    N(u,v) = 2*(m^2+n^2) * (u^2+v^2) * (m^2 v^2 + n^2 u^2) * (n^2 v^2 + m^2 u^2)
и H имеет точку с t = u/v  <=>  N(u,v) -- точный квадрат.

Симметрии:
    t -> -t  :  u -> -u,  N чётна по u  => достаточно u >= 0;
    t -> 1/t :  (u,v) -> (v,u),  N(u,v) = N(v,u)  => достаточно u <= v.
Поэтому перебираем только 0 <= u <= v, gcd(u,v) = 1.

Отсев: для набора модулей M (степени малых простых + простые) требуем
    N(u,v) mod M  принадлежит множеству квадратов mod M.
Это НЕОБХОДИМОЕ условие; таблицы T_M[v%M][u%M] предвычисляются.
Модули упорядочены по селективности; кандидаты сжимаются после каждого модуля.

Чистый Python/NumPy (без sage-препроцессора).
"""

import numpy as np
from math import gcd, isqrt

DEFAULT_MODULI = [64, 27, 25, 49, 121, 169, 17, 19, 23, 29, 31, 37, 41, 43,
                  47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107,
                  109, 113, 127, 131, 137, 139, 149, 151]


def N_form(u, v, m, n):
    """Целая форма степени 6; квадрат <=> t=u/v даёт точку на H."""
    m2 = m * m
    n2 = n * n
    return (2 * (m2 + n2) * (u * u + v * v)
            * (m2 * v * v + n2 * u * u)
            * (n2 * v * v + m2 * u * u))


def is_square_int(x):
    if x < 0:
        return False
    r = isqrt(x)
    return r * r == x


def build_tables(m, n, moduli=None):
    """Список (M, T, density), отсортированный по возрастанию density."""
    if moduli is None:
        moduli = DEFAULT_MODULI
    m2 = m * m
    n2 = n * n
    K = 2 * (m2 + n2)
    out = []
    for M in moduli:
        sq = np.zeros(M, dtype=bool)
        for r in range(M):
            sq[(r * r) % M] = True
        T = np.zeros((M, M), dtype=bool)
        for vr in range(M):
            for ur in range(M):
                a = (K * (ur * ur + vr * vr)) % M
                b = (m2 * vr * vr + n2 * ur * ur) % M
                c = (n2 * vr * vr + m2 * ur * ur) % M
                T[vr, ur] = sq[(a * b % M) * c % M]
        out.append((M, T, float(T.mean())))
    out.sort(key=lambda z: z[2])
    return out


def search(m, n, vmax, tables=None, report_every=None, verbose=True):
    """
    Перебор 0 <= u <= v <= vmax, gcd(u,v)=1.
    Возвращает (hits, stats), где hits -- список (u, v, N) с точным квадратом N,
    stats -- словарь со счётчиками.
    """
    if tables is None:
        tables = build_tables(m, n)
    m2, n2 = m * m, n * n
    K = 2 * (m2 + n2)

    hits = []
    n_exact = 0            # сколько раз дошли до точной проверки isqrt
    n_pairs = 0            # сколько пар (u,v) всего просмотрено (до отсева)

    M0, T0, _ = tables[0]
    rest = tables[1:]

    for v in range(1, vmax + 1):
        n_pairs += v + 1
        row = T0[v % M0]
        res = np.nonzero(row)[0]
        if res.size == 0:
            continue
        blocks = np.arange(0, v // M0 + 2, dtype=np.int64) * M0
        cand = (res.astype(np.int64)[:, None] + blocks[None, :]).ravel()
        cand = cand[cand <= v]
        for (M, T, _d) in rest:
            if cand.size == 0:
                break
            cand = cand[T[v % M][cand % M]]
        if cand.size == 0:
            continue
        for u in cand.tolist():
            if gcd(u, v) != 1:
                continue
            n_exact += 1
            N = N_form(u, v, m, n)
            if is_square_int(N):
                hits.append((u, v, N))
        if report_every and v % report_every == 0 and verbose:
            print("      v = %d  (точных проверок: %d, находок: %d)"
                  % (v, n_exact, len(hits)), flush=True)

    return hits, {"pairs_scanned": n_pairs, "exact_tests": n_exact,
                  "vmax": vmax, "moduli_used": [t[0] for t in tables]}


def lift_data(u, v, m, n):
    """
    Условия подъёма на C рода 5: F0, F4, F8 -- квадраты КАЖДЫЙ по отдельности.
    F0 = (m^2 v^2 + n^2 u^2)/v^2,  F8 = (n^2 v^2 + m^2 u^2)/v^2,
    F4 = ((m^2+n^2)/2)(u^2+v^2)/v^2  -- квадрат <=> 2(m^2+n^2)(u^2+v^2) квадрат.
    """
    m2, n2 = m * m, n * n
    A0 = m2 * v * v + n2 * u * u
    A4 = 2 * (m2 + n2) * (u * u + v * v)
    A8 = n2 * v * v + m2 * u * u
    return {"F0_num": A0, "F0_sq": is_square_int(A0),
            "F4_num": A4, "F4_sq": is_square_int(A4),
            "F8_num": A8, "F8_sq": is_square_int(A8)}


# ---------------------------------------------------------------------------
# Глубокий поиск точек на C вдоль коники "F4 -- квадрат"
# ---------------------------------------------------------------------------
def conic_param(m, n):
    """
    F4(u/v) квадрат <=> 2(m^2+n^2)(u^2+v^2) = w^2 <=> u^2+v^2 = 2(m^2+n^2) z^2.
    Для (m,n)=(11,4): u^2+v^2 = 274 z^2, базовая точка (15,7,1).
    Возвращает функцию (p,q) -> (u,v) или None, если базовой точки нет в малом поиске.
    """
    D = 2 * (m * m + n * n)
    base = None
    for z in range(1, 40):
        tgt = D * z * z
        a = 0
        while a * a * 2 <= tgt:
            b2 = tgt - a * a
            b = isqrt(b2)
            if b * b == b2 and gcd(gcd(a, b), z) == 1:
                base = (a, b, z)
                break
            a += 1
        if base:
            break
    if base is None:
        return None, None
    a0, b0, c0 = base

    def param(p, q):
        # прямая через (a0,b0,c0) с направлением (p,q,0):
        # lam = -(p^2+q^2), mu = 2*(a0*p + b0*q)
        lam = -(p * p + q * q)
        mu = 2 * (a0 * p + b0 * q)
        u = a0 * lam + p * mu
        v = b0 * lam + q * mu
        return u, v

    return base, param


def build_tables_C(m, n, base, moduli=None):
    """
    Таблицы для условия подъёма на C ВДОЛЬ КОНИКИ F4=квадрат.
    (p,q) -> (u,v) по параметризации коники; требуем, чтобы
        G0 = m^2 v^2 + n^2 u^2   и   G8 = n^2 v^2 + m^2 u^2
    ОБА были квадратами (F4 квадрат автоматически по построению коники).
    """
    if moduli is None:
        moduli = DEFAULT_MODULI
    a0, b0, _c0 = base
    m2, n2 = m * m, n * n
    out = []
    for M in moduli:
        sq = np.zeros(M, dtype=bool)
        for r in range(M):
            sq[(r * r) % M] = True
        T = np.zeros((M, M), dtype=bool)
        for qr in range(M):
            for pr in range(M):
                lam = (-(pr * pr + qr * qr)) % M
                mu = (2 * (a0 * pr + b0 * qr)) % M
                u = (a0 * lam + pr * mu) % M
                v = (b0 * lam + qr * mu) % M
                g0 = (m2 * v * v + n2 * u * u) % M
                g8 = (n2 * v * v + m2 * u * u) % M
                T[qr, pr] = sq[g0] and sq[g8]
        out.append((M, T, float(T.mean())))
    out.sort(key=lambda z: z[2])
    return out


def search_conic_C(m, n, pmax, tables=None, base=None, param=None, verbose=True,
                   report_every=None):
    """
    Глубокий поиск точек на C рода 5 вдоль коники "F4 -- квадрат".
    Перебирает (p:q), q = 1..pmax, p = -pmax..pmax, gcd(p,q)=1 (плюс (1:0)).
    Высота t = u/v растёт как ~pmax^2, т.е. покрытие много глубже прямого перебора.
    """
    if base is None or param is None:
        base, param = conic_param(m, n)
    if tables is None:
        tables = build_tables_C(m, n, base)
    hits = []
    n_exact = 0
    n_pairs = 0
    M0, T0, _ = tables[0]
    rest = tables[1:]
    width = 2 * pmax + 1
    for q in range(1, pmax + 1):
        n_pairs += width
        row = T0[q % M0]
        res = np.nonzero(row)[0]
        if res.size == 0:
            continue
        # p пробегает [-pmax, pmax]; берём представителей по модулю M0
        lo = -pmax
        start = lo + ((res.astype(np.int64) - lo) % M0)
        blocks = np.arange(0, width // M0 + 2, dtype=np.int64) * M0
        cand = (start[:, None] + blocks[None, :]).ravel()
        cand = cand[(cand >= -pmax) & (cand <= pmax)]
        for (M, T, _d) in rest:
            if cand.size == 0:
                break
            cand = cand[T[q % M][cand % M]]
        if cand.size == 0:
            continue
        for p in cand.tolist():
            if gcd(abs(p), q) != 1:
                continue
            u, v = param(p, q)
            if v == 0:
                continue
            n_exact += 1
            d = gcd(abs(u), abs(v))
            if d:
                u //= d
                v //= d
            m2, n2 = m * m, n * n
            if (is_square_int(m2 * v * v + n2 * u * u)
                    and is_square_int(n2 * v * v + m2 * u * u)
                    and is_square_int(2 * (m2 + n2) * (u * u + v * v))):
                hits.append((p, q, u, v))
        if report_every and q % report_every == 0 and verbose:
            print("      q = %d  (точных проверок: %d, находок: %d)"
                  % (q, n_exact, len(hits)), flush=True)
    return hits, {"pairs_scanned": n_pairs, "exact_tests": n_exact,
                  "pmax": pmax, "moduli_used": [t[0] for t in tables]}
