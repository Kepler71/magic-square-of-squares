"""
Fable review 2026-09-13: independent 2-isogeny Selmer computation for
E_{m,n}: y^2 = x(x+m^2)(x+n^2)   (alpha = m^2+n^2, beta = m^2 n^2)
E':      Y^2 = X(X^2 - 2 alpha X + (m^2-n^2)^2)

Quartics (Silverman-Tate III.6):
  C_d : w^2 = d z^4 + alpha z^2 + beta/d,           d | mn  (squarefree, signed)
  C'_d: w^2 = d z^4 - 2 alpha z^2 + (m^2-n^2)^2/d,  d | n^2-m^2 (squarefree, signed)
Bound: 2^rank <= |Sel_phi| * |Sel_phihat| / 4  with
  Sel_phi    = {d : C_d  everywhere locally soluble}
  Sel_phihat = {d : C'_d everywhere locally soluble}

Written from scratch without looking at criterion_proof/*.py (rule 8).
"""
from functools import lru_cache
import math

MAXDEPTH = 60
UNDETERMINED = [0]   # counter of depth exhaustion events (treated as soluble)


def vp(x, p):
    if x == 0:
        return 10**9
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def is_unit_square(u, p):
    """u is a p-adic unit (integer coprime to p); is it a square in Z_p?"""
    if p == 2:
        return u % 8 == 1
    return pow(u % p, (p - 1) // 2, p) == 1


def poly_eval(coeffs, t):
    r = 0
    for c in reversed(coeffs):
        r = r * t + c
    return r


def poly_deriv(coeffs):
    return [i * coeffs[i] for i in range(1, len(coeffs))]


def poly_shift_scale(coeffs, r, s):
    """coefficients of h(r + s*t) as polynomial in t (integer arithmetic)."""
    n = len(coeffs)
    out = [0] * n
    # h(r + s t) = sum_i c_i (r + s t)^i
    for i, c in enumerate(coeffs):
        if c == 0:
            continue
        # expand (r + s t)^i
        for j in range(i + 1):
            out[j] += c * math.comb(i, j) * (r ** (i - j)) * (s ** j)
    return out


def ball_soluble(coeffs, p, depth=0, eacc=0):
    """Is there t in Z_p with p^eacc * h(t) a square in Q_p (0 counts as a square)?
    coeffs: integer coefficients of h, low degree first.
    Returns True/False; on depth exhaustion returns True (conservative) and counts it."""
    if all(c == 0 for c in coeffs):
        return True
    e0 = min(vp(c, p) for c in coeffs if c != 0)
    pe = p ** e0
    h = [c // pe for c in coeffs]      # content 1 now
    e = eacc + e0                      # total valuation carried by this ball
    dh = poly_deriv(h)
    if depth > MAXDEPTH:
        UNDETERMINED[0] += 1
        return True
    if p == 2:
        for r in range(8):
            val = poly_eval(h, r)
            if val % 2 == 1:
                # h(t) ≡ h(r) mod 8 for all t ≡ r mod 8
                if e % 2 == 0 and val % 8 == 1:
                    return True
                continue
            # h(r) even
            if poly_eval(dh, r) % 2 == 1:
                return True     # Hensel: simple root -> exact zero -> square (0)
            # multiple root mod 2 in this class: refine t = r + 8 s
            if ball_soluble(poly_shift_scale(h, r, 8), p, depth + 1, e):
                return True
        return False
    else:
        for r in range(p):
            val = poly_eval(h, r) % p
            if val != 0:
                if e % 2 == 0 and pow(val, (p - 1) // 2, p) == 1:
                    return True
                continue
            if poly_eval(dh, r) % p != 0:
                return True     # Hensel
            if ball_soluble(poly_shift_scale(h, r, p), p, depth + 1, e):
                return True
        return False


def quartic_locally_soluble(a, b, c, p):
    """w^2 = a z^4 + b z^2 + c  soluble over Q_p (projective smooth model)?"""
    # z in Z_p
    if ball_soluble([c, 0, b, 0, a], p):
        return True
    # z = 1/t, t in p Z_p (including t=0: points at infinity, needs a square)
    g = [a, 0, b, 0, c]                     # t^4 f(1/t) = a + b t^2 + c t^4
    return ball_soluble(poly_shift_scale(g, 0, p), p)


def quartic_real_soluble(a, b, c):
    if a > 0 or c >= 0:
        return True
    # a<0, c<0: max of a u^2 + b u + c on u>=0
    if b <= 0:
        return False
    return b * b >= 4 * a * c


def prime_factors(n):
    n = abs(n)
    fs = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            fs.append(d)
            while n % d == 0:
                n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        fs.append(n)
    return fs


def squarefree_signed_divisors(primes):
    ds = [1]
    for q in primes:
        ds = ds + [d * q for d in ds]
    return [s * d for d in ds for s in (1, -1)]


def selmer_sets(m, n, extra_primes=()):
    """Exact phi- and phihat-Selmer sets (as lists of d) for E_{m,n}, gcd(m,n)=1, m != n."""
    assert math.gcd(m, n) == 1 and m != n and m > 0 and n > 0
    alpha = m * m + n * n
    beta = m * m * n * n
    D2 = (m * m - n * n) ** 2
    bad = sorted(set([2] + prime_factors(m * n) + prime_factors(n * n - m * m) + list(extra_primes)))
    S1 = []
    for d in squarefree_signed_divisors(prime_factors(m * n)):
        a, b, c = d, alpha, beta // d
        ok = quartic_real_soluble(a, b, c) and all(quartic_locally_soluble(a, b, c, p) for p in bad)
        if ok:
            S1.append(d)
    S2 = []
    for d in squarefree_signed_divisors(prime_factors(n * n - m * m)):
        a, b, c = d, -2 * alpha, D2 // d
        ok = quartic_real_soluble(a, b, c) and all(quartic_locally_soluble(a, b, c, p) for p in bad)
        if ok:
            S2.append(d)
    return S1, S2


def exact_bound(m, n):
    S1, S2 = selmer_sets(m, n)
    return int(round(math.log2(len(S1) * len(S2)))) - 2


def legendre(a, q):
    a %= q
    if a == 0:
        return 0
    return 1 if pow(a, (q - 1) // 2, q) == 1 else -1


def odd_sets(m, n):
    """S_phi, S_phihat of the odd criterion (RESULT_ODD_CRITERION §1), from the definitions."""
    assert math.gcd(m, n) == 1 and m != n
    pm = prime_factors(m * n)
    pd = prime_factors(n * n - m * m)
    q1 = [q for q in pd if q % 4 == 1]
    Sphi = [d for d in squarefree_signed_divisors(pm) if all(legendre(d, q) != -1 for q in q1)]
    podd = [p for p in pm if p != 2]
    Sphihat = []
    for d in squarefree_signed_divisors(pd):
        if d < 0:
            continue
        if any(q % 4 == 3 for q in prime_factors(d)):
            continue
        if any(legendre(d, p) == -1 for p in podd):
            continue
        Sphihat.append(d)
    return Sphi, Sphihat


def odd_bound(m, n):
    S1, S2 = odd_sets(m, n)
    return int(round(math.log2(len(S1) * len(S2)))) - 2


def is_square(x):
    if x < 0:
        return False
    r = math.isqrt(x)
    return r * r == x


def pythagorean(m, n):
    return is_square(m) and is_square(n) and is_square(m + n)


def reduce_pair(a, b):
    a, b = abs(a), abs(b)
    g = math.gcd(a, b)
    return (a // g, b // g)


def twelve_pairs(r, s):
    """The 12 pairs of RESULT_ODD_CRITERION §3 for slope k = r/s (derived independently in phase A)."""
    raw = [(r, s), (r, s - r), (s, s - r), (s, r + s), (r, r + s),                       # cubic triples
           (s - r, r + s), (s, abs(2 * r - s)), (s, 2 * r + s), (r, 2 * s - r), (r, 2 * s + r),
           (s, 2 * r), (r, 2 * s)]                                                         # quartic classes
    out = []
    for a, b in raw:
        m, n = reduce_pair(a, b)
        if m == 0 or n == 0 or m == n:
            continue
        if m > n:
            m, n = n, m
        out.append((m, n))
    return out


def slopes(smin, smax):
    for s in range(smin, smax + 1):
        for r in range(1, s):
            if math.gcd(r, s) == 1 and 2 * r != s:
                yield (r, s)
