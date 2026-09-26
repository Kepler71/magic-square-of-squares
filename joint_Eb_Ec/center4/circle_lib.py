# -*- coding: utf-8 -*-
"""
circle_lib.py — «круговой» (гауссов) взгляд на общий центр.

Для центра a (целое, >0) и шага u (целое) с a-u = alpha^2, a+u = beta^2 (alpha, beta > 0)
кладём  Omega_u = alpha*beta + i*u  (гауссово целое, N(Omega_u) = a^2),
        w_u = Omega_u / a  — рациональная точка единичной окружности, Im w_u = u/a.
При квадратном центре a = g^2 имеем Omega_u = omega_u^2, omega_u = x + i y,
x = (beta+alpha)/2, y = (beta-alpha)/2 (прямоугольный треугольник с гипотенузой g).

Для простого p = 1 (mod 4), p | a, e = v_p(a), и гауссова простого pi | p:
  r_u = v_pi(Omega_u),  rb_u = v_{conj pi}(Omega_u),  r_u + rb_u = 2e;
  t_u = min(r_u, rb_u)  — p-содержание (Omega_u делится на p^t_u как гауссово целое);
  линия u «p-чистая», если t_u = 0; ориентация sigma_u = +1, если r_u = 2e, -1, если rb_u = 2e.
v_pi(w_u) = r_u - e.

Всё — точная целочисленная арифметика (python int), без плавающей точки.
"""
from math import isqrt, gcd
from functools import lru_cache
from fractions import Fraction
from sympy import factorint, Matrix


def is_sq(n):
    return n >= 0 and isqrt(n) ** 2 == n


@lru_cache(maxsize=None)
def gauss_prime(p):
    """p = 1 mod 4 -> (x, y), x^2 + y^2 = p, x odd > 0, y even > 0."""
    assert p % 4 == 1
    # корень из -1 mod p
    for g in range(2, p):
        t = pow(g, (p - 1) // 4, p)
        if t * t % p == p - 1:
            break
    # алгоритм Евклида (Эрмит–Серре / Корнаккья)
    r0, r1 = p, t
    lim = isqrt(p)
    while r1 > lim:
        r0, r1 = r1, r0 % r1
    x = r1
    y = isqrt(p - x * x)
    assert x * x + y * y == p
    if x % 2 == 0:
        x, y = y, x
    return (x, y)


def gdivides(z, q):
    """делится ли гауссово z=(x,y) на q=(c,d); вернуть частное или None"""
    x, y = z
    c, d = q
    n = c * c + d * d
    re = x * c + y * d
    im = y * c - x * d
    if re % n or im % n:
        return None
    return (re // n, im // n)


def gval(z, q):
    """v_q(z) для гауссова простого q (z != 0)"""
    assert z != (0, 0)
    k = 0
    while True:
        w = gdivides(z, q)
        if w is None:
            return k
        z = w
        k += 1


def split_primes(a):
    """простые p = 1 mod 4, делящие a, с показателями"""
    return {p: e for p, e in factorint(a).items() if p % 4 == 1}


def line_data(a, u):
    """Omega_u для центра a и шага u (a +- u — ненулевые квадраты)."""
    lo, hi = a - u, a + u
    assert lo > 0 and hi > 0 and is_sq(lo) and is_sq(hi), (a, u)
    al, be = isqrt(lo), isqrt(hi)
    Om = (al * be, u)
    assert Om[0] ** 2 + Om[1] ** 2 == a * a
    return dict(u=u, alpha=al, beta=be, Omega=Om, g=gcd(al, be))


def local_profile(a, us):
    """для каждого p | a (p = 1 mod 4): e, pi, и по каждой линии (r, rb, t, sigma, v_pi(w))."""
    prof = {}
    for p, e in split_primes(a).items():
        x, y = gauss_prime(p)
        pi, pib = (x, y), (x, -y)
        rows = {}
        for u in us:
            Om = line_data(a, u)['Omega']
            r, rb = gval(Om, pi), gval(Om, pib)
            assert r + rb == 2 * e, (a, u, p, r, rb, e)
            t = min(r, rb)
            sigma = 0 if t == e else (1 if r > rb else -1)
            rows[u] = dict(r=r, rb=rb, t=t, pure=(t == 0), sigma=sigma, vw=r - e)
        prof[p] = dict(e=e, pi=pi, lines=rows)
    return prof


def circle_rank(a, us):
    """ранг подгруппы <w_u> в S^1(Q) (по модулю mu_4) = ранг матрицы v_pi(w_u)."""
    prof = local_profile(a, us)
    ps = sorted(prof)
    if not ps:
        return 0, []
    M = [[prof[p]['lines'][u]['vw'] for p in ps] for u in us]
    return Matrix(M).rank(), M


def min_twice_ok(a, us):
    """Лемма: при линейной связи между шагами us (с коэффициентами, взаимно простыми с p)
    min_u t_u(p) достигается не менее двух раз для каждого p."""
    prof = local_profile(a, us)
    bad = []
    for p, d in prof.items():
        ts = [d['lines'][u]['t'] for u in us]
        m = min(ts)
        if ts.count(m) < 2:
            bad.append((p, ts))
    return (len(bad) == 0), bad


def balance_check(a, us, NL):
    """Теорема D (общий центр): для тройки линий us со связью sum n_u u = 0, NL = sum|n_u|,
    и класса простых P (все линии p-чистые, векторы ориентаций равны с точностью до общего знака):
       Pi_P^2 <= NL * a,   Pi_P = prod_{p in P} p^{v_p(a)}.
    Возвращает (ok, список (класс, Pi, отношение Pi^2/(NL a)))."""
    prof = local_profile(a, us)
    classes = {}
    for p, d in prof.items():
        L = d['lines']
        if not all(L[u]['pure'] for u in us):
            continue
        sig = tuple(L[u]['sigma'] for u in us)
        if sig[0] < 0:
            sig = tuple(-s for s in sig)
        classes.setdefault(sig, []).append(p)
    out = []
    ok = True
    for sig, ps in classes.items():
        Pi = 1
        for p in ps:
            Pi *= p ** prof[p]['e']
        ratio = Fraction(Pi * Pi, NL * a)
        if ratio > 1:
            ok = False
        out.append((sig, sorted(ps), Pi, ratio))
    return ok, out


def relation_Z(a, us, ns, P):
    """Явное вычисление Z из доказательства теоремы D для класса простых P (контроль):
    Z = sum n_u tau_u xi_u, где T_u = conj-или-нет(Omega_u) — единица в pi_p для всех p in P
    (pi_p выбран по ориентации первой линии), T_u = prod conj(pi_p)^{2e_p} * xi_u.
    Возвращает (Z, N(Z), Pi^2, делится ли Z на prod pi_p^{2e_p})."""
    prof = local_profile(a, us)
    # выбор гауссовых простых: для p в P берём pi_p так, чтобы первая линия имела sigma=+1
    chosen = {}
    for p in P:
        d = prof[p]
        s0 = d['lines'][us[0]]['sigma']
        x, y = d['pi']
        chosen[p] = (x, y) if s0 > 0 else (x, -y)
    Z = (0, 0)
    Pi = 1
    for p in P:
        Pi *= p ** prof[p]['e']
    for u, n in zip(us, ns):
        Om = line_data(a, u)['Omega']
        Omb = (Om[0], -Om[1])
        # T_u — тот из Om, Omb, который не делится на выбранный pi_p (для всех p в P)
        cand = []
        for T, tau in ((Om, 1), (Omb, -1)):
            if all(gval(T, chosen[p]) == 0 for p in P):
                cand.append((T, tau))
        assert len(cand) == 1, (u, cand)
        T, tau = cand[0]
        xi = T
        for p in P:
            q = chosen[p]
            qb = (q[0], -q[1])
            for _ in range(2 * prof[p]['e']):
                xi = gdivides(xi, qb)
                assert xi is not None
        Z = (Z[0] + n * tau * xi[0], Z[1] + n * tau * xi[1])
    NZ = Z[0] ** 2 + Z[1] ** 2
    div = True
    for p in P:
        q = chosen[p]
        if Z == (0, 0) or gval(Z, q) < 2 * prof[p]['e']:
            div = False
    return Z, NZ, Pi * Pi, div
