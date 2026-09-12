#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
refute_11_4_claude.py
=====================
НЕЗАВИСИМАЯ попытка ОПРОВЕРГНУТЬ сертификат Касселса-Тейта для пары (m,n)=(11,4)
семейства G1.  Написано с нуля, чистый Python (только fractions/random).
Ничего из кода проекта, Sage, PARI не импортируется и не читается.

Что делает:
  A. Выводит delta = (1,s,s) из тождеств в Q[t] (точная полиномиальная проверка).
  B. Своя арифметика: факторизация, символ Лежандра, символ Гильберта (+ самотесты).
  C. ПОЛНОЕ вычисление локальных образов E(Q_v)/2E(Q_v) методом разбиения шаров
     (не перебор -- доказательное перечисление) и отсюда группа Селмера Sel^2.
  D. Построение бинарных квартик g для класса Селмера (коника -> параметризация ->
     нормировка инвариантов I,J) + проверка z(g) == delta.
  E. Спаривание Касселса-Тейта по Теореме 3.1 Фишера (arXiv:2208.14977).
  F. Контроли: симметрия, диагональ, кручение, независимость локального символа от
     выбора точки, полнота набора мест (перебор всех p <= PLIMIT).

РЕЗУЛЬТАТ (12.09.2026): опровергнуть НЕ УДАЛОСЬ.
    <delta, g2>_CT = -1 для 16 из 32 классов Sel^2, в том числе для класса
    (2055,137,15), использованного в сертификате проекта.
    |Sel^2| = 32, delta = (1,274,274) лежит в Sel^2, ранг матрицы CTP = 2.

Сопутствующие файлы (все в этом каталоге):
    refute_11_4_claude.log            -- основной прогон
    refute_11_4_claude_fisher34.py/.log -- калибровка на опубликованном примере 3.4
                                         Фишера (L -- кубическое поле): совпало
    refute_11_4_claude_indep.py/.log  -- независимость от m, от представителя,
                                         ловушка g3(1,0), полная матрица CTP
    refute_11_4_claude_controls.py    -- контроль на ложное срабатывание
    f8_control_claude.log, f8b_control_claude.log -- его результаты
"""

import sys, random
from fractions import Fraction as Fr
from math import gcd, isqrt

random.seed(20260912)
sys.setrecursionlimit(100000)

PLIMIT = 2000          # brute-force sweep over all primes up to this bound

def log(*a):
    print(*a); sys.stdout.flush()

def hdr(t):
    log("\n" + "=" * 78); log(t); log("=" * 78)

# ======================================================================
#  B1. factorisation, square classes
# ======================================================================
_SMALL = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]

def is_prime(n):
    if n < 2: return False
    for p in _SMALL:
        if n % p == 0: return n == p
    d, r = n - 1, 0
    while d % 2 == 0: d //= 2; r += 1
    for a in (2,3,5,7,11,13,17,19,23,29,31,37):
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1: break
        else:
            return False
    return True

def _rho(n):
    if n % 2 == 0: return 2
    while True:
        x = random.randrange(2, n); y = x; c = random.randrange(1, n); d = 1
        while d == 1:
            x = (x * x + c) % n
            y = (y * y + c) % n; y = (y * y + c) % n
            d = gcd(abs(x - y), n)
        if d != n: return d

_FCACHE = {}
def factorint(n):
    """n>0 -> dict prime:exp"""
    n = int(n)
    assert n > 0
    if n in _FCACHE: return dict(_FCACHE[n])
    orig = n
    f = {}
    for p in _SMALL:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1: continue
        if is_prime(m):
            f[m] = f.get(m, 0) + 1; continue
        d = _rho(m)
        stack.append(d); stack.append(m // d)
    _FCACHE[orig] = dict(f)
    return f

def sqfree(x):
    """squarefree integer representative of the square class of x in Q^*"""
    if isinstance(x, Fr):
        num, den = x.numerator, x.denominator
        n = num * den
    else:
        n = int(x)
    assert n != 0
    sgn = -1 if n < 0 else 1
    n = abs(n)
    r = 1
    for p, e in factorint(n).items():
        if e % 2: r *= p
    return sgn * r

def as_int_sameclass(x):
    """integer with the same square class as the rational x"""
    if isinstance(x, Fr):
        return x.numerator * x.denominator
    return int(x)

def vp(n, p):
    n = int(n)
    assert n != 0, "vp(0,p) undefined"
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v, n

def valuation(x, p):
    if isinstance(x, Fr):
        a, _ = vp(x.numerator, p); b, _ = vp(x.denominator, p)
        return a - b
    return vp(x, p)[0]

def legendre(a, p):
    a %= p
    if a == 0: return 0
    t = pow(a, (p - 1) // 2, p)
    return 1 if t == 1 else -1

def is_square_Qp(x, p):
    """x nonzero rational; is it a square in Q_p?"""
    n = as_int_sameclass(x)
    if p == 0:
        return n > 0
    if n < 0 and p == 0: return False
    v, u = vp(abs(n), p)
    if v % 2: return False
    u = u if n > 0 else -u
    if p == 2:
        return u % 8 == 1
    return legendre(u % p, p) == 1

def sqclass_p(x, p):
    """canonical label of the square class of nonzero rational x in Q_p^*/(Q_p^*)^2"""
    n = as_int_sameclass(x)
    if p == 0:
        return 1 if n > 0 else -1
    sgn = 1 if n > 0 else -1
    v, u = vp(abs(n), p)
    u *= sgn
    if p == 2:
        return (v % 2, u % 8)
    return (v % 2, legendre(u % p, p))

# ---------------- Hilbert symbol, written from scratch ----------------
def _eps2(u):     # (u-1)/2 mod 2 for odd u
    return ((u % 4) - 1) // 2
def _om2(u):      # (u^2-1)/8 mod 2 for odd u
    return 0 if (u % 8) in (1, 7) else 1

def hilbert(a, b, p):
    """Hilbert symbol (a,b)_p in {1,-1}; p==0 means the real place."""
    a = as_int_sameclass(a); b = as_int_sameclass(b)
    assert a != 0 and b != 0
    if p == 0:
        return -1 if (a < 0 and b < 0) else 1
    if p == 2:
        al, u = vp(abs(a), 2); be, v = vp(abs(b), 2)
        u = u if a > 0 else -u
        v = v if b > 0 else -v
        e = _eps2(u) * _eps2(v) + al * _om2(v) + be * _om2(u)
        return -1 if e % 2 else 1
    al, u = vp(abs(a), p); be, v = vp(abs(b), p)
    u = u if a > 0 else -u
    v = v if b > 0 else -v
    e = (al * be * ((p - 1) // 2)) % 2
    s = -1 if e else 1
    if be % 2: s *= legendre(u % p, p)
    if al % 2: s *= legendre(v % p, p)
    return s

def hilbert_selftest():
    """properties that characterise the Hilbert symbol + global reciprocity"""
    bad = 0
    places = [0, 2] + [p for p in range(3, 200) if is_prime(p)]
    for _ in range(400):
        a = random.randrange(-3000, 3000) or 7
        b = random.randrange(-3000, 3000) or 5
        # reciprocity: product over all places = 1
        S = set(places)
        for q in list(factorint(abs(a)).keys()) + list(factorint(abs(b)).keys()):
            S.add(q)
        S.add(0)
        pr = 1
        for v in S: pr *= hilbert(a, b, v)
        if pr != 1: bad += 1
    # (a,-a)_v = 1 and (a,1-a)_v = 1
    for _ in range(200):
        a = Fr(random.randrange(-500, 500) or 3, random.randrange(1, 40))
        if a == 0 or a == 1: continue
        for v in [0, 2, 3, 5, 7, 11, 13, 137]:
            if hilbert(a, -a, v) != 1: bad += 1
            if hilbert(a, 1 - a, v) != 1: bad += 1
    # bilinearity
    for _ in range(200):
        a = random.randrange(-200, 200) or 3
        b = random.randrange(-200, 200) or 5
        c = random.randrange(-200, 200) or 7
        for v in [0, 2, 3, 5, 7, 11, 137]:
            if hilbert(a * b, c, v) != hilbert(a, c, v) * hilbert(b, c, v): bad += 1
    return bad

# ======================================================================
#  Polynomials in one variable over Q (lists, index = degree)
# ======================================================================
def pmul(A, B):
    r = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B):
                r[i + j] += x * y
    return r
def padd(A, B):
    n = max(len(A), len(B)); r = [Fr(0)] * n
    for i, x in enumerate(A): r[i] += x
    for i, x in enumerate(B): r[i] += x
    return r
def psub(A, B): return padd(A, [-x for x in B])
def pscal(c, A): return [c * x for x in A]
def pzero(A): return all(x == 0 for x in A)

# ======================================================================
#  A. the family G1 and the derivation of delta
# ======================================================================
def part_A(m, n):
    hdr("A.  Семейство G1, вывод класса delta из тождеств в Q[t]  (m=%d, n=%d)" % (m, n))
    s = Fr(m * m + n * n, 2)
    b = s * m * m * n * n
    log("s = %s,  b = %s" % (s, b))
    # F0 = m^2 + n^2 t^2 ; F4 = s(1+t^2) ; F8 = n^2 + m^2 t^2   (poly in t)
    F0 = [Fr(m * m), Fr(0), Fr(n * n)]
    F4 = [s, Fr(0), s]
    F8 = [Fr(n * n), Fr(0), Fr(m * m)]
    X = [Fr(0), Fr(0), b]                 # X = b t^2
    e1, e2, e3 = -b, -s * m ** 4, -s * n ** 4
    log("e1 = %s   e2 = %s   e3 = %s" % (e1, e2, e3))
    id1 = psub(psub(X, [e1]), pscal(Fr((m * n) ** 2), F4))
    id2 = psub(psub(X, [e2]), pscal(s * m * m, F0))
    id3 = psub(psub(X, [e3]), pscal(s * n * n, F8))
    log("X - e1 == (mn)^2 * F4 :", pzero(id1))
    log("X - e2 == s m^2 * F0  :", pzero(id2))
    log("X - e3 == s n^2 * F8  :", pzero(id3))
    assert pzero(id1) and pzero(id2) and pzero(id3)
    d = (sqfree(Fr((m * n) ** 2)), sqfree(s * m * m), sqfree(s * n * n))
    log("=> delta = (sqfree((mn)^2), sqfree(s m^2), sqfree(s n^2)) = %s" % (d,))
    sq = sqfree(s)
    log("s квадрат? ", sq == 1, "  (sqfree(s) = %s)" % sq)
    log("  => точки на бесконечности C исключены  (нужно было бы s in Q^*2)")
    # V^2 = ... check E equation is the product
    return s, b, (e1, e2, e3), d

# ======================================================================
#  C. local images of E(Q_v)/2E(Q_v) and the Selmer group
# ======================================================================
class Curve:
    """E : Y^2 = (X-E1)(X-E2)(X-E3) with integral distinct roots"""
    def __init__(self, roots):
        self.e = [int(r) for r in roots]
        E1, E2, E3 = self.e
        a2 = -(E1 + E2 + E3)
        a4 = E1 * E2 + E1 * E3 + E2 * E3
        a6 = -E1 * E2 * E3
        self.a2, self.a4, self.a6 = a2, a4, a6
        b2 = 4 * a2; b4 = 2 * a4; b6 = 4 * a6
        self.b2 = b2
        c4 = b2 * b2 - 24 * b4
        c6 = -b2 ** 3 + 36 * b2 * b4 - 216 * b6
        self.c4, self.c6 = c4, c6
        self.I = Fr(c4)
        self.J = Fr(2 * c6)
        self.phi = [Fr(-(12 * r + b2)) for r in self.e]
        self.disc = 16 * ((E1 - E2) * (E1 - E3) * (E2 - E3)) ** 2
    def f(self, x):
        return (x - self.e[0]) * (x - self.e[1]) * (x - self.e[2])
    def check(self):
        ok = True
        for ph in self.phi:
            ok &= (ph ** 3 - 3 * self.I * ph + self.J == 0)
        return ok
    def torsion_classes(self):
        E1, E2, E3 = self.e
        t1 = (sqfree((E1 - E2) * (E1 - E3)), sqfree(E1 - E2), sqfree(E1 - E3))
        t2 = (sqfree(E2 - E1), sqfree((E2 - E1) * (E2 - E3)), sqfree(E2 - E3))
        t3 = (sqfree(E3 - E1), sqfree(E3 - E2), sqfree((E3 - E1) * (E3 - E2)))
        return [(1, 1, 1), t1, t2, t3]

def local_image(C, p, maxdepth=200):
    """
    COMPLETE enumeration of the image of E(Q_p) -> (Q_p^*/sq)^3, x |-> (x-e_i).
    Ball subdivision; no random search, no 'not found => empty'.
    p == 0 means the real place.
    """
    e = C.e
    if p == 0:
        img = set()
        img.add((1, 1, 1))                       # O
        pts = sorted(e)
        cands = [pts[0] - 1, pts[1] - 1, pts[2] - 1, pts[2] + 1,
                 Fr(pts[0] + pts[1], 2), Fr(pts[1] + pts[2], 2)]
        for x in cands:
            v = C.f(x)
            if v > 0:
                img.add(tuple(1 if (x - r) > 0 else -1 for r in e))
        return img

    k = 3 if p == 2 else 1
    img = set()
    img.add(tuple(sqclass_p(1, p) for _ in range(3)))   # identity (x = infinity)
    tors = C.torsion_classes()[1:]
    # start ball: v(x) >= N0 covers everything except the 'near infinity' part,
    # which only produces the identity class (proved in the report text below).
    N0 = min(valuation(r, p) for r in e) - k
    stack = [(Fr(0), N0)]
    seen = set()
    while stack:
        c, nn = stack.pop()
        if (c, nn) in seen: continue
        seen.add((c, nn))
        vs = []
        for r in e:
            d = Fr(r) - c
            vs.append(None if d == 0 else valuation(d, p))
        inside = [i for i in range(3) if vs[i] is None or vs[i] >= nn]
        if not inside:
            if all(nn >= vs[i] + k for i in range(3)):
                val = C.f(c)
                if val != 0 and is_square_Qp(val, p):
                    img.add(tuple(sqclass_p(c - r, p) for r in e))
                continue
        elif len(inside) == 1:
            i = inside[0]
            others = [j for j in range(3) if j != i]
            if all(nn >= vs[j] + k for j in others):
                img.add(tuple(sqclass_p(t, p) for t in tors[i]))
                continue
        if nn > N0 + maxdepth:
            raise RuntimeError("ball recursion too deep at p=%d" % p)
        for j in range(p):
            stack.append((c + j * Fr(p) ** nn, nn + 1))
    return img

def selmer_group(C, S, quiet=False):
    """S = list of bad primes (including 2). Returns list of delta triples."""
    gens = [-1] + S
    sub = []
    for mask in range(1 << len(gens)):
        v = 1
        for i, g in enumerate(gens):
            if mask >> i & 1: v *= g
        sub.append(v)
    places = [0] + S
    limg = {v: local_image(C, v) for v in places}
    for v in places:
        exp = 2 if v == 0 else (8 if v == 2 else 4)
        if quiet: continue
        log("   локальный образ в v=%-4s : %d элементов (ожидается %d)  %s"
            % ("R" if v == 0 else v, len(limg[v]), exp,
               "OK" if len(limg[v]) == exp else "!!! НЕСОВПАДЕНИЕ"))
    out = []
    for d1 in sub:
        for d2 in sub:
            d3 = sqfree(d1 * d2)
            ok = True
            for v in places:
                cl = (sqclass_p(d1, v), sqclass_p(d2, v), sqclass_p(d3, v))
                if cl not in limg[v]:
                    ok = False; break
            if ok: out.append((d1, d2, d3))
    return out, limg

# ======================================================================
#  D. binary quartics
# ======================================================================
def quartic_I(g):
    a, b, c, d, e = g
    return 12 * a * e - 3 * b * d + c * c
def quartic_J(g):
    a, b, c, d, e = g
    return 72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3
def quartic_disc(g):
    I = quartic_I(g); J = quartic_J(g)
    return Fr(16, 27) * (4 * I ** 3 - J ** 2)
def quartic_hessian(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c,
            4 * (b * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e),
            3 * d * d - 8 * c * e]
def qeval(g, x, z):
    a, b, c, d, e = g
    return a * x ** 4 + b * x ** 3 * z + c * x ** 2 * z ** 2 + d * x * z ** 3 + e * z ** 4

def solve_conic(A, B, Cc, box=6000, wmax=120):
    """A u^2 + B v^2 + Cc w^2 = 0 ; integer coefficients. Bounded search.
       Returns (u,v,w) integers, not all zero, or None."""
    A, B, Cc = int(A), int(B), int(Cc)
    g0 = gcd(gcd(abs(A), abs(B)), abs(Cc))
    if g0 > 1: A //= g0; B //= g0; Cc //= g0
    for (P, Q, R, perm) in ((A, B, Cc, (0, 1, 2)), (B, A, Cc, (1, 0, 2)),
                            (A, Cc, B, (0, 2, 1))):
        # solve P u^2 + Q v^2 = -R w^2
        for w in range(1, wmax + 1):
            rhs = -R * w * w
            for u in range(0, box + 1):
                t = rhs - P * u * u
                if Q == 0: continue
                if t % Q: continue
                y = t // Q
                if y < 0: continue
                r = isqrt(y)
                if r * r == y:
                    sol = [0, 0, 0]
                    sol[perm[0]] = u; sol[perm[1]] = r; sol[perm[2]] = w
                    if A * sol[0] ** 2 + B * sol[1] ** 2 + Cc * sol[2] ** 2 == 0 and any(sol):
                        return tuple(sol)
    return None

def parametrise_conic(A, B, Cc, P):
    """returns three quadratic forms (as dicts of coeffs in u^2,uv,v^2) giving the
       general rational point of A x^2 + B y^2 + Cc z^2 = 0 through P."""
    M = [A, B, Cc]
    def Q(V): return sum(M[i] * V[i] ** 2 for i in range(3))
    def Bil(U, V): return sum(M[i] * U[i] * V[i] for i in range(3))
    assert Q(P) == 0
    basis = []
    for cand in ([1,0,0],[0,1,0],[0,0,1]):
        trial = basis + [cand]
        mat = [list(P)] + trial
        if len(mat) == 3:
            det = (mat[0][0]*(mat[1][1]*mat[2][2]-mat[1][2]*mat[2][1])
                   - mat[0][1]*(mat[1][0]*mat[2][2]-mat[1][2]*mat[2][0])
                   + mat[0][2]*(mat[1][0]*mat[2][1]-mat[1][1]*mat[2][0]))
            if det != 0:
                basis = trial; break
        else:
            # independence of P and cand
            cross = [P[1]*cand[2]-P[2]*cand[1], P[2]*cand[0]-P[0]*cand[2],
                     P[0]*cand[1]-P[1]*cand[0]]
            if any(cross): basis = trial
    assert len(basis) == 2, "basis construction failed"
    Aa, Bb = basis
    # D = u*Aa + v*Bb ; Q(D) = qA u^2 + qB uv + qC v^2 ; Bil(P,D) = lA u + lB v
    qA, qB, qC = Q(Aa), 2 * Bil(Aa, Bb), Q(Bb)
    lA, lB = Bil(P, Aa), Bil(P, Bb)
    out = []
    for i in range(3):
        # comp = (qA u^2+qB uv+qC v^2)*P[i] - 2*(lA u + lB v)*(u*Aa[i]+v*Bb[i])
        c_uu = qA * P[i] - 2 * lA * Aa[i]
        c_uv = qB * P[i] - 2 * (lA * Bb[i] + lB * Aa[i])
        c_vv = qC * P[i] - 2 * lB * Bb[i]
        out.append((Fr(c_uu), Fr(c_uv), Fr(c_vv)))
    return out

def qform_mul(F, G):
    """product of two quadratic binary forms -> quartic coefficient list"""
    a0, a1, a2 = F; b0, b1, b2 = G
    return [a0 * b0,
            a0 * b1 + a1 * b0,
            a0 * b2 + a1 * b1 + a2 * b0,
            a1 * b2 + a2 * b1,
            a2 * b2]

def _quartic_from_torsor(C, d, i, j, box=6000, wmax=120):
    """Build the quartic of the torsor  x-e_l = d_l z_l^2  using the conic in
       (z_i, z_j, w).  Returns the normalised quartic or None."""
    e = C.e
    kk = 3 - i - j
    cij = e[j] - e[i]                # d_i z_i^2 - d_j z_j^2 = c_ij w^2
    sol = solve_conic(d[i], -d[j], -cij, box=box, wmax=wmax)
    if sol is None: return None, None
    par = parametrise_conic(d[i], -d[j], -cij, list(sol))
    Zi, Zj, W = par
    chk = [Fr(0)] * 5
    for t, co in ((Zi, d[i]), (Zj, -d[j]), (W, -cij)):
        q = qform_mul(t, t)
        chk = [chk[x] + co * q[x] for x in range(5)]
    assert all(x == 0 for x in chk), "parametrisation off the conic"
    cik = e[kk] - e[i]
    qi = qform_mul(Zi, Zi); qw = qform_mul(W, W)
    g = [Fr(d[kk]) * (Fr(d[i]) * qi[x] - Fr(cik) * qw[x]) for x in range(5)]
    if quartic_disc(g) == 0:  return None, None
    Ig, Jg = quartic_I(g), quartic_J(g)
    if Ig == 0 or Jg == 0: return None, None
    lam2 = Jg * C.I / (C.J * Ig)
    gn = [x / lam2 for x in g]
    if quartic_I(gn) != C.I or quartic_J(gn) != C.J:
        return None, None
    return gn, sol

def _bconv(A, B):
    r = [Fr(0)] * (len(A) + len(B) - 1)
    for i, x in enumerate(A):
        if x:
            for j, y in enumerate(B): r[i + j] += x * y
    return r

def subst(g, p, q, r, s):
    """g'(x,z) = g(p x + q z, r x + s z); index in the list = power of z."""
    A = [Fr(p), Fr(q)]; B = [Fr(r), Fr(s)]
    Ap = [[Fr(1)]]; Bp = [[Fr(1)]]
    for _ in range(4):
        Ap.append(_bconv(Ap[-1], A)); Bp.append(_bconv(Bp[-1], B))
    out = [Fr(0)] * 5
    for i in range(5):
        term = _bconv(Ap[4 - i], Bp[i])
        for k, v in enumerate(term): out[k] += g[i] * v
    return out

def make_z_unit(C, g, tries=60):
    """Fisher p.3: by a change of coordinates we may assume z(g) is a unit in L.
       Proper equivalence g -> g o gamma, gamma in SL_2(Z), keeps I, J and the class."""
    zz = z_invariant(C, g)
    if all(t != 0 for t in zz): return g
    mats = []
    for t in range(-12, 13):
        mats.append((1, 0, t, 1)); mats.append((1, t, 0, 1))
        mats.append((0, 1, -1, t)); mats.append((t, 1, -1, 0))
    for (p, q, r, s) in mats:
        if p * s - q * r not in (1, -1): continue
        g2 = subst(g, p, q, r, s)
        if quartic_I(g2) != C.I or quartic_J(g2) != C.J: continue
        zz = z_invariant(C, g2)
        if all(t != 0 for t in zz): return g2
    return None

def build_quartic(C, target, selset=None, verbose=True, box=6000, wmax=120):
    """
    Return an everywhere-locally-soluble binary quartic g with invariants (I,J)
    of C and with Fisher class  z(g) == target  (checked, not assumed).

    The elementary torsor construction singles out one root, so the resulting
    2-covering MAP differs from Fisher's canonical one by a translation by a
    rational 2-torsion point; the class then differs by a torsion class.  We do
    not guess that shift: we build candidates for target * (each torsion class)
    and keep the one whose actually computed z(g) class equals target.
    """
    tors = C.torsion_classes()
    cands = []
    for t in tors:
        d = tuple(sqfree(target[l] * t[l]) for l in range(3))
        if selset is not None and d not in selset:
            continue
        cands.append(d)
    for d in cands:
        for (i, j) in ((0, 1), (0, 2), (1, 2), (1, 0), (2, 0), (2, 1)):
            gn, sol = _quartic_from_torsor(C, [int(x) for x in d], i, j, box=box, wmax=wmax)
            if gn is None: continue
            gn = make_z_unit(C, gn)
            if gn is None: continue
            cl = tuple(sqfree(t) for t in z_invariant(C, gn))
            if cl == tuple(target):
                if verbose:
                    log("   цель=%s: торсор %s, коника (%d,%d), точка %s -> z(g) класс %s OK"
                        % (tuple(target), d, i, j, sol, cl))
                return gn
    return None

def z_invariant(C, g):
    """z(g) in L = Q^3 (component i corresponds to phi_i <-> root e_i)"""
    a, b, c, dd, e = g
    return [ (4 * a * ph + 3 * b * b - 8 * a * c) / 3 for ph in C.phi ]

def H_form(C, g):
    """H(x,z) = H2 x^2 + H1 xz + H0 z^2 with H_i in L (lists of 3 rationals).
       Also verifies the identity G(1,0)*G = H^2 over each component."""
    h = quartic_hessian(g)
    Hs = []
    for idx, ph in enumerate(C.phi):
        G = [(4 * ph * g[t] + h[t]) / 3 for t in range(5)]
        G0 = G[0]
        H = (G0, G[1] / 2, G[2] / 6 + Fr(2, 9) * (C.I - ph * ph))
        # check G0 * G == H^2
        HH = qform_mul(H, H)
        lhs = [G0 * G[t] for t in range(5)]
        assert all(lhs[t] == HH[t] for t in range(5)), "G(1,0)G = H^2 failed"
        Hs.append(H)
    # regroup: coefficient of x^2, xz, z^2 each as element of L
    return [ [Hs[i][k] for i in range(3)] for k in range(3) ]

def L_interp_gamma(C, xi):
    """xi = (xi_1,xi_2,xi_3) components of an element of L; return the coefficient
       of phi^2 when xi is written as alpha + beta*phi + gamma*phi^2."""
    ph = C.phi
    g = Fr(0)
    for i in range(3):
        den = Fr(1)
        for j in range(3):
            if j != i: den *= (ph[i] - ph[j])
        g += xi[i] / den
    return g

# ======================================================================
#  E. Fisher's Theorem 3.1
# ======================================================================
def gamma1_form(C, g1, g2, g3):
    z1 = z_invariant(C, g1); z2 = z_invariant(C, g2); z3 = z_invariant(C, g3)
    for zz in (z1, z2, z3):
        assert all(t != 0 for t in zz), "z(g) is a zero divisor in L (Fisher forbids)"
    m = []
    for i in range(3):
        pr = z1[i] * z2[i] * z3[i]
        num = pr.numerator; den = pr.denominator
        assert num > 0, "z1z2z3 not a square (negative) in component %d" % i
        r1 = isqrt(num); r2 = isqrt(den)
        assert r1 * r1 == num and r2 * r2 == den, "z1z2z3 not a square in comp %d" % i
        m.append(Fr(r1, r2))
    w = [z2[i] * z3[i] / m[i] for i in range(3)]
    H1 = H_form(C, g1)
    gam = []
    for k in range(3):
        xi = [w[i] * H1[k][i] for i in range(3)]
        gam.append(L_interp_gamma(C, xi))
    return gam, (z1, z2, z3), m

def gamma_eval(gam, x, z):
    return gam[0] * x * x + gam[1] * x * z + gam[2] * z * z

def local_points(g1, gam, p, want=6):
    """rational (x,z) with g1(x,z) a nonzero square in Q_p and gamma(x,z) != 0"""
    out = []
    cands = []
    for z in range(0, 25):
        for x in range(-25, 26):
            if x == 0 and z == 0: continue
            cands.append((Fr(x), Fr(z)))
    if p not in (0,):
        for k in range(1, 9):
            cands.append((Fr(1), Fr(p) ** k)); cands.append((Fr(p) ** k, Fr(1)))
            cands.append((Fr(1), Fr(p) ** (-k))); cands.append((Fr(p) ** (-k), Fr(1)))
            cands.append((Fr(1), Fr(p) ** k + 1))
    else:
        for k in range(1, 40):
            cands.append((Fr(1), Fr(1, 10 ** k))); cands.append((Fr(1, 10 ** k), Fr(1)))
            cands.append((Fr(k), Fr(1))); cands.append((Fr(-k), Fr(1)))
    for (x, z) in cands:
        v = qeval(g1, x, z)
        if v == 0: continue
        if not is_square_Qp(v, p): continue
        if gamma_eval(gam, x, z) == 0: continue
        out.append((x, z))
        if len(out) >= want: break
    return out

def ct_pairing(C, g1, g2, g3, verbose=True, tag=""):
    """<[g1],[g2]>_CT via Fisher Thm 3.1. Returns (value, per-place dict, info)."""
    gam, zs, m = gamma1_form(C, g1, g2, g3)
    a = g2[0]                                    # g2(1,0)
    assert a != 0, "g2(1,0) = 0 : [g2] is trivial (Remark 3.2(iii))"
    # --- place set ---
    D1 = quartic_disc(g1)
    nums = [a, D1]
    nums += [c for c in g1 if c != 0]
    nums += [c for c in gam if c != 0]
    S = set([2, 3, 5, 7])
    for t in nums:
        n = as_int_sameclass(t)
        for q in factorint(abs(n)): S.add(q)
    R33 = sorted(S)                              # Remark 3.3 place set (finite part)
    allP = sorted(S | {q for q in range(2, PLIMIT) if is_prime(q)})
    res = {}
    detail = {}
    for p in [0] + allP:
        pts = local_points(g1, gam, p, want=4)
        if not pts:
            raise RuntimeError("no local point at p=%s" % p)
        vals = set()
        for (x, z) in pts:
            vals.add(hilbert(a, gamma_eval(gam, x, z), p))
        if len(vals) != 1:
            raise RuntimeError("local symbol depends on the point at p=%s !!" % p)
        res[p] = vals.pop()
        detail[p] = pts[0]
    prod = 1
    for p in res: prod *= res[p]
    minus = [("R" if p == 0 else p) for p in res if res[p] == -1]
    if verbose:
        log("   %s a=g2(1,0)=%s" % (tag, a))
        log("   %s набор мест по Замечанию 3.3 (конечная часть): %s" % (tag, R33))
        log("   %s проверено мест всего: %d (все p < %d плюс делители)" % (tag, len(res), PLIMIT))
        log("   %s места с символом -1: %s" % (tag, minus))
        log("   %s ПРОИЗВЕДЕНИЕ = %+d" % (tag, prod))
    # restricted product over the Remark 3.3 set only, as a cross-check
    prod33 = hilbert(a, gamma_eval(gam, *detail[0]), 0)
    for p in R33:
        prod33 *= res[p]
    return prod, res, dict(gamma=gam, a=a, R33=R33, prod33=prod33, z=zs, m=m)

# ======================================================================
#  main
# ======================================================================
def main():
    m, n = 11, 4
    s, b, eorig, delta = part_A(m, n)

    hdr("B.  Собственный символ Гильберта -- самотест")
    bad = hilbert_selftest()
    log("нарушений (взаимность на 400 парах + (a,-a), (a,1-a), билинейность): %d" % bad)
    assert bad == 0

    hdr("C.  Модель E с целыми корнями, локальные образы, группа Селмера")
    # scale roots by 4 (a square) : square classes are unchanged
    roots = [int(4 * r) for r in eorig]
    C = Curve(roots)
    log("корни (после умножения на 4, квадрат => классы те же): %s" % (C.e,))
    log("I = c4 = %s" % C.I)
    log("J = 2c6 = %s" % C.J)
    log("phi_i корни X^3-3IX+J :", C.check())
    assert C.check()
    log("phi = %s" % ([str(x) for x in C.phi],))
    dd = [C.e[0]-C.e[1], C.e[0]-C.e[2], C.e[1]-C.e[2]]
    log("разности корней: %s" % dd)
    badp = set()
    for t in dd:
        for q in factorint(abs(t)): badp.add(q)
    badp.add(2)
    S = sorted(badp)
    log("плохие простые S = %s" % S)
    tors = C.torsion_classes()
    log("классы кручения (образ E(Q)[2] в (Q*/Q*^2)^3): %s" % (tors,))
    log("целевой класс delta = %s ; он среди классов кручения? %s"
        % (delta, delta in tors))
    sel, limg = selmer_group(C, S)
    log("|Sel^2(E/Q)| = %d   (dim = %d)" % (len(sel), len(sel).bit_length() - 1))
    log("delta в Sel^2 ? %s" % (tuple(delta) in [tuple(x) for x in sel]))
    assert tuple(delta) in [tuple(x) for x in sel]

    hdr("D.  Построение квартик")
    d1 = tuple(delta)
    # a second Selmer class, chosen independently of the project's certificates:
    # take every non-torsion Selmer class and try to build a quartic for it
    others = [x for x in sel if x not in tors and x != d1]
    log("нетривиальных классов Селмера кроме delta и кручения: %d" % len(others))

    selset = set(tuple(x) for x in sel)
    g1 = build_quartic(C, d1, selset)
    assert g1 is not None, "не удалось построить квартику для delta"
    log("g1 = %s" % ([str(x) for x in g1],))
    z1 = z_invariant(C, g1)
    cl1 = tuple(sqfree(t) for t in z1)
    log("z(g1) = %s -> класс %s ; совпадает с delta? %s"
        % ([str(x) for x in z1], cl1, cl1 == d1))
    assert cl1 == d1
    log("I(g1)=I, J(g1)=J :", quartic_I(g1) == C.I and quartic_J(g1) == C.J)

    hdr("E.  Спаривание Касселса-Тейта <delta, g2>_CT для всех доступных g2")
    results = {}
    quartics = {d1: g1}
    for d2 in others:
        d3 = tuple(sqfree(d1[i] * d2[i]) for i in range(3))
        if d3 not in [tuple(x) for x in sel]:
            log("   пропуск %s : произведение не в Селмере (?!)" % (d2,)); continue
        g2 = quartics.get(d2) or build_quartic(C, d2, selset, verbose=False)
        g3 = quartics.get(d3) or build_quartic(C, d3, selset, verbose=False)
        if g2 is None or g3 is None:
            log("   %s : коника не решена перебором -- НЕЗАВЕРШЕНО" % (d2,)); continue
        quartics[d2] = g2; quartics[d3] = g3
        c2 = tuple(sqfree(t) for t in z_invariant(C, g2))
        c3 = tuple(sqfree(t) for t in z_invariant(C, g3))
        assert c2 == d2 and c3 == d3, "класс квартики не совпал"
        try:
            val, res, info = ct_pairing(C, g1, g2, g3, verbose=False)
        except RuntimeError as ex:
            log("   %s : %s" % (d2, ex)); continue
        results[d2] = (val, res, info, g2, g3)
        minus = sorted([("R" if p == 0 else p) for p in res if res[p] == -1],
                       key=lambda t: (t == "R", t if t != "R" else 0))
        log("   <delta, %s> = %+d   минусы: %s   (свёрнутый набор 3.3 даёт %+d)"
            % (d2, val, minus, info['prod33']))

    hdr("F.  КОНТРОЛИ (попытки сломать результат)")

    def getq(cl):
        cl = tuple(cl)
        if cl not in quartics:
            g = build_quartic(C, cl, selset, verbose=False)
            if g is None: return None
            quartics[cl] = g
        return quartics[cl]

    def pair(clA, clB):
        """<[clA],[clB]>_CT, building all three quartics as needed"""
        clA = tuple(clA); clB = tuple(clB)
        clC = tuple(sqfree(clA[i] * clB[i]) for i in range(3))
        gA, gB, gC = getq(clA), getq(clB), getq(clC)
        if gA is None or gB is None or gC is None: return None, None
        if gB[0] == 0: return None, None
        v, res, info = ct_pairing(C, gA, gB, gC, verbose=False)
        return v, (res, info)

    # ---- F1. link "quartic <-> descent class" verified on quartics WITH a point
    log("F1. Связь z(g) <-> класс спуска (X+3phi_i), на квартиках с явной точкой:")
    nchk = 0; nbad = 0
    for cl, g in list(quartics.items()):
        pt = None
        h = quartic_hessian(g)
        for zc in range(0, 40):
            for xc in range(-40, 41):
                if xc == 0 and zc == 0: continue
                v = qeval(g, Fr(xc), Fr(zc))
                if v == 0: continue
                if sqfree(v) != 1: continue
                x0c, z0c = Fr(xc), Fr(zc)
                hvc = (h[0]*x0c**4 + h[1]*x0c**3*z0c + h[2]*x0c**2*z0c**2
                       + h[3]*x0c*z0c**3 + h[4]*z0c**4)
                Xc = Fr(3) * hvc / (4 * v)
                if all(Xc + 3*ph != 0 for ph in C.phi):
                    pt = (x0c, z0c); break
            if pt: break
        if pt is None: continue
        x0, z0 = pt
        hv = (h[0]*x0**4 + h[1]*x0**3*z0 + h[2]*x0**2*z0**2 + h[3]*x0*z0**3 + h[4]*z0**4)
        gv = qeval(g, x0, z0)
        X = Fr(3) * hv / (4 * gv)
        Y2 = X**3 - 27*C.I*X - 27*C.J
        onE = (sqfree(Y2) == 1) if Y2 != 0 else True
        zz = z_invariant(C, g)
        cls_from_X = tuple(sqfree(X + 3*ph) for ph in C.phi)
        cls_from_z = tuple(sqfree(t) for t in zz)
        ok = (cls_from_X == cls_from_z) and onE
        nchk += 1
        if not ok: nbad += 1
        log("    класс z(g)=%s : точка (%s,%s), X=3h/4g лежит на E: %s, класс X+3phi=%s  %s"
            % (cls_from_z, x0, z0, onE, cls_from_X, "OK" if ok else "!!! РАСХОЖДЕНИЕ"))
        if nchk >= 6: break
    log("    проверено квартик с точкой: %d, расхождений: %d" % (nchk, nbad))

    # ---- F2. torsion classes (definitely in im E(Q)/2E(Q)) must pair to +1
    log("F2. классы кручения обязаны спариваться тривиально с delta:")
    for t in tors[1:]:
        v, _ = pair(d1, t)
        log("    <delta, %s> = %s  %s" % (t, ("%+d" % v) if v is not None else "н/д",
            "OK" if v == 1 else ("!!! ЛОЖНОЕ СРАБАТЫВАНИЕ" if v is not None else "НЕЗАВЕРШЕНО")))

    # ---- F3. non-torsion rational point of E(Q), if one can be found
    log("F3. поиск нетривиальной рациональной точки E(Q) и спаривание её класса:")
    found = []
    for X in range(-4100000, 400000):
        val = (X-C.e[0])*(X-C.e[1])*(X-C.e[2])
        if val <= 0: continue
        r = isqrt(val)
        if r*r == val:
            cl = tuple(sqfree(X - C.e[i]) for i in range(3))
            if cl not in tors:
                found.append((X, cl))
                if len(found) >= 3: break
    if found:
        for X, cl in found:
            v, _ = pair(d1, cl)
            log("    точка X=%s, класс %s ; <delta, класс> = %s  %s"
                % (X, cl, ("%+d" % v) if v is not None else "н/д",
                   "OK" if v == 1 else "!!! СЕРТИФИКАТ ЛОЖЕН"))
    else:
        log("    целых X в [-4.1e6, 4e5] с квадратным f(X) и неторсионным классом не найдено")
        log("    ОТСУТСТВИЕ НАХОДКИ НЕ ЕСТЬ ДОКАЗАТЕЛЬСТВО ОТСУТСТВИЯ (генератор может быть высоким)")

    # ---- F4. bilinearity of <delta, . >
    log("F4. билинейность <delta, . > : множество с +1 обязано быть подгруппой индекса 2")
    vals = {tuple(k): v[0] for k, v in results.items()}
    vals[(1,1,1)] = 1
    for t in tors[1:]:
        vv, _ = pair(d1, t); vals[tuple(t)] = vv
    nbil = 0; bbil = 0
    keys = [k for k in vals if vals[k] is not None]
    for A in keys:
        for B in keys:
            Cc = tuple(sqfree(A[i]*B[i]) for i in range(3))
            if Cc in vals and vals[Cc] is not None:
                nbil += 1
                if vals[Cc] != vals[A]*vals[B]: bbil += 1
    log("    проверено троек: %d, нарушений билинейности: %d  %s"
        % (nbil, bbil, "OK" if bbil == 0 else "!!! НЕ БИЛИНЕЙНО"))
    log("    классов с +1: %d, с -1: %d (линейный функционал на Sel^2 => 16/16)"
        % (sum(1 for v in vals.values() if v == 1),
           sum(1 for v in vals.values() if v == -1)))

    # ---- F5. symmetry and diagonal
    log("F5. симметрия и диагональ:")
    shown = 0
    for d2 in list(vals.keys()):
        if vals[d2] != -1: continue
        v2, _ = pair(d2, d1)
        log("    <%s, delta> = %s  (прямое %+d)  %s" % (d2, ("%+d" % v2) if v2 is not None else "н/д",
            vals[d2], "OK" if v2 == vals[d2] else "!!! РАСХОЖДЕНИЕ"))
        vd, _ = pair(d2, d2)
        log("    диагональ <%s,%s> = %s  %s" % (d2, d2, ("%+d" % vd) if vd is not None else "н/д",
            "OK" if vd == 1 else "!!! ДИАГОНАЛЬ НЕ НУЛЕВАЯ"))
        shown += 1
        if shown >= 2: break
    vdd, _ = pair(d1, d1)
    log("    диагональ <delta,delta> = %s  %s" % (("%+d" % vdd) if vdd is not None else "н/д",
        "OK" if vdd == 1 else "!!! ДИАГОНАЛЬ НЕ НУЛЕВАЯ"))

    # ---- F6. permuted delta: is the component ORDER load-bearing?
    log("F6. переставленные варианты delta (чувствительность к порядку координат):")
    for perm in [(274,1,274), (274,274,1)]:
        if perm not in selset:
            log("    %s : не в Sel^2" % (perm,)); continue
        nz2 = 0; tot = 0
        for d2 in list(vals.keys()):
            v, _ = pair(perm, d2)
            if v is None: continue
            tot += 1
            if v == -1: nz2 += 1
        log("    <%s, .> : -1 в %d случаях из %d" % (perm, nz2, tot))

    # ---- F7. places carrying -1 must lie inside the Remark 3.3 set
    log("F7. все ли места с символом -1 лежат в наборе Замечания 3.3:")
    bad7 = 0
    for d2, (val, res, info, gg2, gg3) in results.items():
        R33 = set(info['R33'])
        for pp, sv in res.items():
            if sv == -1 and pp != 0 and pp not in R33:
                bad7 += 1
                log("    !!! место %d вне набора 3.3 даёт -1 (класс %s)" % (pp, d2))
    log("    нарушений: %d  %s" % (bad7, "OK" if bad7 == 0 else "!!! ЗАМЕЧАНИЕ 3.3 НАРУШЕНО"))
    log("    символ считался во ВСЕХ местах p < %d плюс все делители ->" % PLIMIT)
    log("    ловушка 'пропущенное место' (как с p=67 у (19,5)) исключена перебором")

    # ---- F8. false-positive control on independent curves
    log("F8. контроль на ложное срабатывание: независимые кривые, у которых Sel^2")
    log("    целиком покрыт явными рациональными точками => CTP обязан быть нулевым")
    log("    (вынесен в отдельный прогон f8_control_claude.log ради времени)")

    log("")
    log("СВОДКА")
    nz = [d2 for d2, r in results.items() if r[0] == -1]
    log("  классов g2 с <delta,g2> = -1 : %d из %d посчитанных" % (len(nz), len(results)))
    key = (2055, 137, 15)
    log("  против класса сертификата проекта %s : %s"
        % (key, ("%+d" % results[key][0]) if key in results else "не посчитано"))
    if nz:
        log("  ИТОГ: <delta, .>_CT нетривиально => delta НЕ лежит в образе E(Q)/2E(Q)")
        log("        => аффинных рациональных точек на C_(11,4) нет; бесконечность")
        log("           исключена тем, что s = 137/2 не квадрат.")
        log("        ОПРОВЕРГНУТЬ ЗАЯВЛЕНИЕ НЕ УДАЛОСЬ.")
    else:
        log("  ИТОГ: спаривание тривиально везде -- ЗАЯВЛЕНИЕ ОПРОВЕРГНУТО")
    return results


# ======================================================================
#  F8: independent control curves
# ======================================================================
def control_curves():
    """Curves with full 2-torsion where Sel^2 is exhausted by explicit rational
       points.  Then im(E(Q)/2E(Q)) = Sel^2 and the Cassels-Tate pairing must be
       identically trivial.  Any -1 produced here would expose a bug."""
    tested = 0
    for roots in ([1, 2, -1], [1, 5, -3], [3, 4, -1], [1, 10, -2], [2, 3, -5],
                  [1, 9, -7], [5, 6, -1], [1, 4, -11], [2, 7, -3], [1, 3, -5]):
        try:
            C = Curve(roots)
            if C.disc == 0: continue
            badp = set([2])
            dd = [C.e[0]-C.e[1], C.e[0]-C.e[2], C.e[1]-C.e[2]]
            for t in dd:
                for q in factorint(abs(t)): badp.add(q)
            S = sorted(badp)
            sel, _ = selmer_group(C, S, quiet=True)
            selset = set(tuple(x) for x in sel)
            img = set([(1,1,1)] + [tuple(t) for t in C.torsion_classes()])
            for X in range(-1500, 1501):
                val = C.f(X)
                if val <= 0: continue
                r = isqrt(val)
                if r*r == val:
                    img.add(tuple(sqfree(X - C.e[i]) for i in range(3)))
            for den in (2, 3, 4, 5, 6, 7, 8):
                for num in range(-1500, 1501):
                    x = Fr(num, den*den)
                    val = C.f(x)
                    if val <= 0: continue
                    if sqfree(val) == 1:
                        img.add(tuple(sqfree(x - C.e[i]) for i in range(3)))
            img = {c for c in img if c in selset}
            if len(img) != len(selset):
                log("    корни %s : |Sel^2|=%d, точками покрыто %d -- не годится в контроль"
                    % (roots, len(selset), len(img)))
                continue
            quartics = {}
            worst = 1; done = 0
            cls_list = sorted(selset)
            for A in cls_list:
                for B in cls_list:
                    Cc = tuple(sqfree(A[i]*B[i]) for i in range(3))
                    gs = []; okall = True
                    for cl in (A, B, Cc):
                        if cl not in quartics:
                            gq = build_quartic(C, cl, selset, verbose=False, box=400, wmax=20)
                            if gq is None: okall = False; break
                            quartics[cl] = gq
                        gs.append(quartics[cl])
                    if not okall or gs[1][0] == 0: continue
                    try:
                        v, _, _ = ct_pairing(C, gs[0], gs[1], gs[2], verbose=False)
                    except Exception:
                        continue
                    done += 1
                    if v == -1: worst = -1
            log("    корни %s : |Sel^2|=%d, покрыт точками, посчитано пар %d, "
                "есть ли -1: %s  %s"
                % (roots, len(selset), done, "ДА" if worst == -1 else "нет",
                   "OK" if worst == 1 else "!!! ЛОЖНОЕ СРАБАТЫВАНИЕ"))
            tested += 1
            if tested >= 3: break
        except Exception as ex:
            log("    корни %s : пропущено (%s)" % (roots, ex))
    if tested == 0:
        log("    подходящих контрольных кривых не нашлось -- КОНТРОЛЬ НЕ ВЫПОЛНЕН")


if __name__ == "__main__":
    main()
