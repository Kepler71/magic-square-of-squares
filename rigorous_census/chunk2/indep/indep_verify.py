#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Claude, 26.09.2026, второй (независимый) заход по куску 2. Код предыдущего захода (../rigdem.py) и код Fable
# (joint_sieve/demjanenko.py) НЕ импортируются. Из чужих данных берутся только: список множителей T
# (joint_sieve/dem_<r>_<s>.json) и точка P0 на минимальной модели (../rig_<r>_<s>.json) — точка лишь кандидат,
# её роль (образующая по модулю кручения) заново доказывается здесь.
#
# ЧИСТЫЙ PYTHON (дроби, целые, Decimal с направленными запасами):
#   модель E1: V² = f(X) = ∏(X + L/λ), L = ∏λ, z = X/L (3-клеточные T);
#   изоморфизм Emin → E1 (X = a x + b) — точным тождеством f(a x + b) = a³ g(x);
#   групповой закон, удвоение, сертификат Безу (своё гауссово исключение 8×8), B = ln max(K,L)/3;
#   c1 по целой матрице z ↔ x (строчные и столбцовые суммы); нижняя оценка ĥ(P0) телескопом;
#   M0; кручение (верхняя граница НОД #E(F_p), явные точки проверяются своим законом);
#   2-насыщение точно над Q (критерий 2-спуска), p-насыщение для нечётных p ≤ mmax редукцией;
#   случай A (все z(nP0+t), |n| ≤ M0), случай B (своя формула X(P+t) = −X(P)); проверка восьми клеток.
# SAGE/PARI/eclib (ПО) только для: верхней границы ранга (PARI ellrank и eclib mwrank — две реализации),
#   нижней границы ĥ на точках бесконечного порядка (height_function().min) и индекса насыщения eclib,
#   списка точек кручения (проверяется своим кодом), рациональных корней многочленов случая B,
#   и контрольных сравнений (Sage height, silverman_height_bound).
import os, sys, json, time, math
from fractions import Fraction as Fr
from decimal import Decimal, getcontext

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
getcontext().prec = 70
LN2 = Decimal(2).ln()
EPS = Decimal('1e-50')
LOGF = open(os.path.join(HERE, 'indep_verify.log'), 'a')

def log(m):
    s = time.strftime('%H:%M:%S') + ' ' + str(m)
    print(s, flush=True); LOGF.write(s + '\n'); LOGF.flush()

# ---------------- логарифмы с направленными границами ----------------
def ln_int(n):
    """(lo, hi): lo ≤ ln n ≤ hi для целого n ≥ 1"""
    n = int(n); assert n >= 1
    if n == 1:
        return Decimal(0), Decimal(0)
    b = n.bit_length()
    if b <= 160:
        v = Decimal(n).ln(); return v - EPS, v + EPS
    s = b - 160; m = n >> s               # m·2^s ≤ n < (m+1)·2^s
    return Decimal(m).ln() + s * LN2 - EPS, Decimal(m + 1).ln() + s * LN2 + EPS

def hq(q):
    q = Fr(q); return ln_int(max(abs(q.numerator), q.denominator))

def isqrt_exact(n):
    if n < 0: return None
    r = math.isqrt(n); return r if r * r == n else None

def is_sq(q):
    q = Fr(q)
    return q >= 0 and isqrt_exact(q.numerator) is not None and isqrt_exact(q.denominator) is not None

def sqrt_q(q):
    q = Fr(q); return Fr(isqrt_exact(q.numerator), isqrt_exact(q.denominator))

# ---------------- эллиптическая кривая над Q (общая форма Вейерштрасса) ----------------
class Curve:
    def __init__(self, a):
        self.a1, self.a2, self.a3, self.a4, self.a6 = [Fr(x) for x in a]
    def on(self, P):
        if P is None: return True
        x, y = P
        return y * y + self.a1 * x * y + self.a3 * y == x ** 3 + self.a2 * x * x + self.a4 * x + self.a6
    def neg(self, P):
        if P is None: return None
        x, y = P; return (x, -y - self.a1 * x - self.a3)
    def add(self, P, Q):
        if P is None: return Q
        if Q is None: return P
        x1, y1 = P; x2, y2 = Q
        if x1 == x2:
            den = 2 * y1 + self.a1 * x1 + self.a3
            if y1 + y2 + self.a1 * x2 + self.a3 == 0: return None
            lam = (3 * x1 * x1 + 2 * self.a2 * x1 + self.a4 - self.a1 * y1) / den
            nu = (-x1 ** 3 + self.a4 * x1 + 2 * self.a6 - self.a3 * y1) / den
        else:
            lam = (y2 - y1) / (x2 - x1); nu = (y1 * x2 - y2 * x1) / (x2 - x1)
        x3 = lam * lam + self.a1 * lam - self.a2 - x1 - x2
        y3 = -(lam + self.a1) * x3 - nu - self.a3
        return (x3, y3)
    def mul(self, n, P):
        if n < 0: return self.mul(-n, self.neg(P))
        R = None; Q = P
        while n:
            if n & 1: R = self.add(R, Q)
            Q = self.add(Q, Q); n >>= 1
        return R
    def binv(self):
        a1, a2, a3, a4, a6 = self.a1, self.a2, self.a3, self.a4, self.a6
        b2 = a1 * a1 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3 * a3 + 4 * a6
        b8 = a1 * a1 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3 * a3 - a4 * a4
        return b2, b4, b6, b8
    def c46(self):
        b2, b4, b6, b8 = self.binv()
        return b2 * b2 - 24 * b4, -b2 ** 3 + 36 * b2 * b4 - 216 * b6
    def disc(self):
        b2, b4, b6, b8 = self.binv()
        return -b2 * b2 * b8 - 8 * b4 ** 3 - 27 * b6 * b6 + 9 * b2 * b4 * b6

# ---------------- многочлены (списки Fraction, младший коэффициент первым) ----------------
def pnorm(p):
    p = [Fr(c) for c in p]
    while p and p[-1] == 0: p.pop()
    return p
def padd(p, q):
    n = max(len(p), len(q)); return pnorm([(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)])
def pscale(p, c): return pnorm([c * x for x in p])
def pmul(p, q):
    if not p or not q: return []
    r = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0: continue
        for j, b in enumerate(q): r[i + j] += a * b
    return pnorm(r)
def peval(p, x):
    v = Fr(0)
    for c in reversed(p): v = v * x + c
    return v
def pcompose_lin(p, a, b):
    """p(a x + b)"""
    r = []; pw = [Fr(1)]
    for c in p:
        r = padd(r, pscale(pw, c)); pw = pmul(pw, [Fr(b), Fr(a)])
    return r

def primitive_int(polys):
    den = 1
    for p in polys:
        for c in p: den = den * c.denominator // math.gcd(den, c.denominator)
    ps = [[int(c * den) for c in p] for p in polys]
    g = 0
    for p in ps:
        for c in p: g = math.gcd(g, c)
    return [[c // g for c in p] for p in ps]

# ---------------- сертификат Безу для удвоения ----------------
def solve_frac(M, rhs):
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(n)] + [Fr(rhs[i])] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None: raise ValueError('вырожденная матрица Сильвестра')
        A[col], A[piv] = A[piv], A[col]
        pv = A[col][col]
        A[col] = [x / pv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                fct = A[r][col]; A[r] = [A[r][k] - fct * A[col][k] for k in range(n + 1)]
    return [A[i][n] for i in range(n)]

def ipmul(p, q):
    r = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i + j] += a * b
    return r

def dup_cert(E):
    b2, b4, b6, b8 = E.binv()
    Fq = [-b8, -2 * b6, -b4, Fr(0), Fr(1)]
    Gq = [b6, 2 * b4, b2, Fr(4)]
    F, G = primitive_int([Fq, Gq])
    d = 4
    assert F[4] != 0 and len(F) == 5
    cols = []
    for H in (F, G):
        for i in range(d):
            v = [0] * i + H; v = v + [0] * (2 * d - len(v)); cols.append(v)
    M = [[cols[c][r] for c in range(2 * d)] for r in range(2 * d)]
    sols = {j: solve_frac(M, [1 if i == j else 0 for i in range(2 * d)]) for j in (0, 2 * d - 1)}
    R = 1
    for j in sols:
        for c in sols[j]: R = R * c.denominator // math.gcd(R, c.denominator)
    wit = []
    for j in (0, 2 * d - 1):
        v = [c * R for c in sols[j]]; assert all(c.denominator == 1 for c in v)
        v = [int(c) for c in v]; A = v[:d]; B = v[d:]
        lhs = [x + y for x, y in zip(ipmul(A, F) + [0] * 8, ipmul(B, G) + [0] * 8)]
        rhs = [0] * len(lhs); rhs[j] = R
        lhs = lhs[:max(len(lhs), 8)]
        assert all(lhs[k] == (R if k == j else 0) for k in range(len(lhs))), 'тождество Безу не выполнено'
        wit.append(dict(j=j, A=[str(c) for c in A], B=[str(c) for c in B], norm=sum(abs(c) for c in v)))
    assert R != 0
    K = max(w['norm'] for w in wit)
    Lb = max(sum(abs(c) for c in F), sum(abs(c) for c in G))
    return dict(F=[str(c) for c in F], G=[str(c) for c in G], R=str(R), K=K, L=Lb, witnesses=wit)

def dup_x(E, x):
    b2, b4, b6, b8 = E.binv()
    return (x ** 4 - b4 * x * x - 2 * b6 * x - b8) / (4 * x ** 3 + b2 * x * x + 2 * b4 * x + b6)

# ---------------- модель E_T и изоморфизм с минимальной моделью ----------------
def model3(T):
    L = 1
    for t in T: L *= t
    f = [Fr(1)]
    for t in T: f = pmul(f, [Fr(L, t), Fr(1)])
    assert len(f) == 4 and f[3] == 1
    return f, L

def iso_min_to_E1(Emin, E1, f):
    """a, b, u: X = a x + b (a = u²), V = u³ (y + (a1 x + a3)/2); проверка тождеством f(a x + b) = a³ g(x)"""
    b2, b4, b6, b8 = Emin.binv()
    g = [b6 / 4, b4 / 2, b2 / 4, Fr(1)]
    c4m, c6m = Emin.c46(); c41, c61 = E1.c46()
    cands = []
    if c4m != 0 and c6m != 0:
        cands = [c61 * c4m / (c6m * c41)]
    else:
        raise ValueError('j = 0 или 1728: не предусмотрено')
    for a in cands:
        b = (a * g[2] - f[2]) / 3
        if pcompose_lin(f, a, b) == pscale(g, a ** 3) and is_sq(a):
            return a, b, sqrt_q(a)
    raise ValueError('изоморфизм не найден')

def to_E1(P, Emin, a, b, u):
    if P is None: return None
    x, y = P; return (a * x + b, u ** 3 * (y + (Emin.a1 * x + Emin.a3) / 2))

def to_min(P, Emin, a, b, u):
    if P is None: return None
    X, V = P; x = (X - b) / a; return (x, V / u ** 3 - (Emin.a1 * x + Emin.a3) / 2)

# ---------------- арифметика по модулю простого ℓ (для E1: a1 = a3 = 0, целые коэффициенты) ----------------
def count_mod(f_int, l):
    e = (l - 1) // 2; s = 0
    for x in range(l):
        v = (((f_int[3] * x + f_int[2]) * x + f_int[1]) * x + f_int[0]) % l
        if v: s += 1 if pow(v, e, l) == 1 else -1
    return l + 1 + s

def add_mod(P, Q, a2, a4, l):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % l == 0: return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * pow(2 * y1, l - 2, l) % l
    else:
        lam = (y2 - y1) * pow(x2 - x1, l - 2, l) % l
    x3 = (lam * lam - a2 - x1 - x2) % l
    return (x3, (lam * (x1 - x3) - y1) % l)

def mul_mod(n, P, a2, a4, l):
    R = None; Q = P
    while n:
        if n & 1: R = add_mod(R, Q, a2, a4, l)
        Q = add_mod(Q, Q, a2, a4, l); n >>= 1
    return R

# ---------------- ячейки квадрата ----------------
def cells(r, s): return [s, -s, r, -r, s - r, r - s, s + r, -s - r]
def full_ok(z, S): return all(is_sq(1 + l * z) for l in S)

# ---------------- проверка одного множителя ----------------
def check_factor(r, s, T, P0min_str, Emin_ainvs, sage):
    t0 = time.time()
    S = cells(r, s)
    out = dict(T=T)
    assert len(T) == 3 and sorted(T) != sorted(-x for x in T)
    assert all(t in S for t in T)
    f, L = model3(T)
    E1 = Curve([0, f[2], 0, f[1], f[0]])
    assert E1.disc() != 0
    Emin = Curve([Fr(a) for a in Emin_ainvs])
    a, b, u = iso_min_to_E1(Emin, E1, f)
    out['iso'] = dict(a=str(a), b=str(b))
    # точка P0 (кандидат из предыдущего захода) — на кривой?
    xs, ys = P0min_str.strip('()').split(',')
    P0m = (Fr(xs.strip()), Fr(ys.strip()))
    assert Emin.on(P0m)
    P01 = to_E1(P0m, Emin, a, b, u); assert E1.on(P01)
    # контроль закона сложения и формулы удвоения
    for E, P in ((Emin, P0m), (E1, P01)):
        P2 = E.add(P, P); assert E.on(P2) and dup_x(E, P[0]) == P2[0]
        P3 = E.add(P2, P); assert E.on(P3) and E.add(P3, E.neg(P)) == P2
    assert to_E1(Emin.mul(3, P0m), Emin, a, b, u) == E1.mul(3, P01)
    # ---- кручение ----
    f_int = [int(c) for c in f]; assert all(Fr(c) == c for c in f_int)
    D1 = E1.disc(); assert D1.denominator == 1; D1 = int(D1)
    g = 0; nprimes = 0
    for p in range(3, 700):
        if any(p % q == 0 for q in range(2, int(p ** 0.5) + 1)) or D1 % p == 0: continue
        g = math.gcd(g, count_mod(f_int, p)); nprimes += 1
    tors1 = [None]
    for pt in sage['tors_E1']:
        P = (Fr(pt[0]), Fr(pt[1])); assert E1.on(P)
        assert any(E1.mul(k, P) is None for k in range(1, 17)), 'точка кручения не кручения'
        if P not in tors1: tors1.append(P)
    # 2-кручение из корней f должно входить
    roots2 = [Fr(-L, t) for t in T]
    for e in roots2:
        assert peval(f, e) == 0 and (e, Fr(0)) in tors1
    # замкнутость относительно сложения
    for P in tors1:
        for Q in tors1: assert E1.add(P, Q) in tors1
    out['tors_order'] = len(tors1); out['tors_gcd_bound'] = g; out['tors_primes'] = nprimes
    tors_proved = (len(tors1) == g)
    out['tors_proved'] = tors_proved
    # ---- ранг (ПО, две реализации) ----
    out['ellrank'] = sage['ellrank']; out['mwrank_bound'] = sage['mwrank_bound']
    rank_ok = sage['ellrank'][1] == 1 or sage['mwrank_bound'] == 1
    # ---- сертификаты и константы на двух моделях ----
    Mz = {'E1': [[1, 0], [0, L]]}
    # z = (a x + b)/L  → целая примитивная матрица
    den = a.denominator * b.denominator // math.gcd(a.denominator, b.denominator)
    m = [[int(a * den), int(b * den)], [0, L * den]]
    gg = math.gcd(math.gcd(m[0][0], m[0][1]), m[1][1]); m = [[c // gg for c in row] for row in m]
    if m[1][1] < 0: m = [[-c for c in row] for row in m]
    Mz['min'] = m
    variants = {}
    for name, E, P in (('min', Emin, P0m), ('E1', E1, P01)):
        cert = dup_cert(E)
        Bhi = ln_int(max(cert['K'], cert['L']))[1] / 3
        M = Mz[name]
        # проверка матрицы: z(P) = M(x(P)) для P0, 2P0, 3P0 (z считается на E1 как X/L)
        for k in (1, 2, 3):
            Pk = E.mul(k, P); Pk1 = E1.mul(k, P01)
            zz = (M[0][0] * Pk[0] + M[0][1]) / (M[1][0] * Pk[0] + M[1][1])
            assert zz == Pk1[0] / L
        sums = [abs(M[0][0]) + abs(M[0][1]), abs(M[1][0]) + abs(M[1][1]),
                abs(M[0][0]) + abs(M[1][0]), abs(M[0][1]) + abs(M[1][1])]
        c1hi = ln_int(max(sums))[1]
        # нижняя/верхняя оценка ĥ(P0) телескопом по этой модели
        hP = sage['hP0_sage']
        kd = 6 if hP < 8 else (5 if hP < 30 else 4)
        Q = P
        for _ in range(kd): Q = E.add(Q, Q)
        hlo_k, hhi_k = hq(Q[0])
        hK = ln_int(cert['K']); hL = ln_int(cert['L'])
        h_lo = (hlo_k - hK[1] / 3) / 4 ** kd
        h_hi = (hhi_k + hL[1] / 3) / 4 ** kd
        C = 2 * Bhi + 2 * c1hi
        variants[name] = dict(K=str(cert['K']), L=str(cert['L']), R=cert['R'], B=float(Bhi), c1=float(c1hi), C=float(C),
                              M=M, kdup=kd, hP0_lo=float(h_lo), hP0_hi=float(h_hi), _C=C, _hlo=h_lo, _hhi=h_hi, _B=Bhi, cert=cert)
    # ĥ(P0): пересечение интервалов двух моделей
    h_lo = max(variants[n]['_hlo'] for n in variants); h_hi = min(variants[n]['_hhi'] for n in variants)
    assert h_lo > 0 and h_lo <= h_hi
    out['hP0_lo'] = float(h_lo); out['hP0_hi'] = float(h_hi); out['hP0_sage'] = sage['hP0_sage']
    out['norm_ok'] = bool(float(h_lo) - 1e-12 <= sage['hP0_sage'] <= float(h_hi) + 1e-12)
    best = min(variants, key=lambda n: variants[n]['_C'])
    C = variants[best]['_C']; Bbest = variants[best]['_B']
    val = C / h_lo + Decimal('1e-40')
    M0 = int(((val + 1) / 2).to_integral_value(rounding='ROUND_FLOOR'))
    out.update(model=best, B=float(Bbest), c1=variants[best]['c1'], C=float(C), M0=M0)
    out['B_silverman_sage_min'] = sage['silverman_min']; out['B_CPS_sage_min'] = sage.get('cps_min')
    out['B_vs_silverman'] = 'строгая B больше' if float(Bbest) > sage['silverman_min'] else 'строгая B не больше'
    out['variants'] = {n: {k: v for k, v in variants[n].items() if not k.startswith('_')} for n in variants}
    # ---- насыщение ----
    lam = sage['lambda_min']
    # граница индекса по СВОЕЙ верхней оценке ĥ(P0) (телескоп), λ — нижняя граница ĥ на точках бесконечного порядка (Sage, ПО)
    mmax = int(math.isqrt(int(math.floor(float(h_hi) * (1 + 1e-9) / lam)))) if lam > 0 else None
    out['lambda_min'] = lam; out['index_bound'] = mmax
    sat = {}
    # p = 2: точный критерий 2-спуска (полное 2-кручение): Q ∈ 2E(Q) ⇔ X − e_i — квадраты для всех i
    ok2 = True
    for t in tors1:
        Q = E1.add(P01, t)
        if all(is_sq(Q[0] - e) for e in roots2): ok2 = False
    sat['2'] = ok2
    primes_odd = [p for p in range(3, (mmax or 0) + 1) if all(p % q for q in range(2, p))]
    for p in primes_odd:
        okp = True
        for t in tors1:
            Q = E1.add(P01, t); hit = None
            for l in range(5, 20000):
                if any(l % q == 0 for q in range(2, int(l ** 0.5) + 1)) or D1 % l == 0: continue
                if Q[0].denominator % l == 0 or Q[1].denominator % l == 0: continue
                N = count_mod(f_int, l)
                if N % p: continue
                Ql = (Q[0].numerator * pow(Q[0].denominator, l - 2, l) % l, Q[1].numerator * pow(Q[1].denominator, l - 2, l) % l)
                if mul_mod(N // p, Ql, f_int[2] % l, f_int[1] % l, l) is not None:
                    hit = l; break
            if hit is None: okp = False
        sat[str(p)] = okp
    out['saturation_own'] = sat; out['sat_index_eclib'] = sage['sat_index_eclib']
    sat_ok = mmax is not None and all(sat.get(str(p), False) for p in [2] + primes_odd if p <= mmax) and sage['sat_index_eclib'] == 1
    out['sat_ok'] = sat_ok
    # ---- случай A ----
    vals = set(); bad = []
    def consider(z, tag):
        if z is None: return
        vals.add(z); vals.add(-z)
        if z != 0 and full_ok(z, S): bad.append([tag, str(z)])
    for t in tors1:
        if t is not None: consider(t[0] / L, 'tors')
    nP = None; nA = 0
    for n in range(1, M0 + 1):
        nP = E1.add(nP, P01)
        for t in tors1:
            Q = E1.add(nP, t)
            if Q is None: raise AssertionError('nP0 + t = O при n ≥ 1')
            consider(Q[0] / L, 'A%d' % n); nA += 1
    out['nA'] = nA
    # ---- случай B: X(P+t) = −X(P) ----
    caseB = []; zero_poly = []
    a2 = E1.a2
    for t in tors1:
        if t is None: continue
        xt, yt = t
        q2 = [xt * xt, -2 * xt, Fr(1)]
        if yt == 0:
            Qp = padd(f, pscale(q2, -(a2 + xt)))
        else:
            W = padd(padd(f, [yt * yt]), pscale(q2, -(a2 + xt)))
            Qp = padd(pmul(W, W), pscale(f, -4 * yt * yt))
        if not Qp:
            zero_poly.append(str(t)); continue
        roots = sage['roots'](Qp)
        for X0 in roots:
            if X0 == xt: continue
            fv = peval(f, X0)
            if not is_sq(fv): continue
            yv = sqrt_q(fv)
            for yy in {yv, -yv}:
                P = (X0, yy); PT = E1.add(P, t)
                eq = PT is not None and PT[0] == -X0
                caseB.append(dict(t=str(t), X=str(X0), eq=eq))
                consider(X0 / L, 'B')
    out['caseB'] = caseB; out['caseB_zero_poly'] = zero_poly
    out['n_values'] = len(vals); out['nondegenerate'] = bad
    # ---- контроль за пределами M0: симметричные z (обе точки P_T(±z) рациональны) обязаны быть в списке ----
    miss = []; symm = []
    nP = E1.mul(M0, P01)
    for n in range(M0 + 1, 2 * M0 + 7):
        nP = E1.add(nP, P01)
        for t in tors1:
            Q = E1.add(nP, t); z = Q[0] / L
            pr = Fr(1)
            for lam_ in T: pr *= 1 - lam_ * z
            if is_sq(pr):
                symm.append([n, str(z)])
                if z not in vals: miss.append([n, str(z)])
    out['ctrl_beyond'] = dict(nmax=2 * M0 + 6, symmetric=symm, missing=miss)
    # ---- численный контроль |ĥ − h(x)| ≤ B на кратных (heights от Sage — только контроль) ----
    out['dev_sample'] = sage.get('dev_sample')
    out['dev_within_B'] = bool(sage.get('dev_sample') is not None and sage['dev_sample'][best] <= float(Bbest))
    out['closed'] = bool(rank_ok and tors_proved and sat_ok and out['norm_ok'] and not bad and not zero_poly and not miss)
    out['sec'] = round(time.time() - t0, 2)
    return out, vals

def brute_ctrl(T, vals, N):
    """все z = a/b, |a| ≤ N, 1 ≤ b ≤ N, с рациональными P_T(z) и P_T(−z): обязаны лежать в vals"""
    found = []; missing = []
    for bb in range(1, N + 1):
        for aa in range(-N, N + 1):
            if math.gcd(aa, bb) != 1: continue
            p1 = bb; p2 = bb
            for lam in T: p1 *= bb + lam * aa; p2 *= bb - lam * aa
            if p1 < 0 or p2 < 0: continue
            if isqrt_exact(p1) is not None and isqrt_exact(p2) is not None:
                z = Fr(aa, bb); found.append(str(z))
                if z not in vals: missing.append(str(z))
    return dict(bound=N, symmetric_found=found, missing=missing)

# ---------------- обращения к Sage (ПО) ----------------
def sage_part(T, Emin_ainvs, P0str, E1ainvs):
    from sage.all import EllipticCurve, QQ, pari, PolynomialRing, ZZ, RealField
    from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
    E1 = EllipticCurve(QQ, [QQ(str(x)) for x in E1ainvs])
    Emin = EllipticCurve(QQ, [QQ(str(x)) for x in Emin_ainvs])
    assert Emin.is_isomorphic(E1) and Emin.is_minimal()
    res = {}
    er = pari(E1).ellrank(); res['ellrank'] = [int(er[0]), int(er[1])]
    try:
        alarm(180); mw = Emin.mwrank_curve(); res['mwrank_bound'] = int(mw.rank_bound()); cancel_alarm()
    except (AlarmInterrupt, Exception) as ex:
        cancel_alarm(); res['mwrank_bound'] = None; res['mwrank_err'] = repr(ex)[:100]
    xs, ys = P0str.strip('()').split(',')
    P0 = Emin([QQ(xs.strip()), QQ(ys.strip())])
    res['hP0_sage'] = float(P0.height())
    res['lambda_min'] = float(Emin.height_function().min(0.0001, 20))
    sat, idx, reg = Emin.saturation([P0]); res['sat_index_eclib'] = int(idx)
    res['tors_E1'] = [[str(c) for c in t.xy()] for t in E1.torsion_points() if not t.is_zero()]
    res['silverman_min'] = float(Emin.silverman_height_bound())
    try:
        res['cps_min'] = float(Emin.CPS_height_bound())
    except Exception:
        res['cps_min'] = None
    Rx = PolynomialRing(QQ, 'x')
    def roots(pl):
        from fractions import Fraction as Fr_
        return [Fr_(int(rt.numerator()), int(rt.denominator())) for rt, _ in Rx([QQ(str(c)) for c in pl]).roots(QQ)]
    res['roots'] = roots
    # численный контроль |ĥ(Q) − h(x(Q))| на Q = nP0 + t, n ≤ 6, по обеим моделям
    iso = Emin.isomorphism_to(E1)
    dev = {'min': 0.0, 'E1': 0.0}
    tors = Emin.torsion_points()
    for n in range(1, 7):
        for t in tors:
            Q = n * P0 + t
            hh = float(Q.height())
            for name, xx in (('min', Q[0]), ('E1', iso(Q)[0])):
                hx = float(RealField(200)(max(abs(xx.numerator()), xx.denominator())).log())  # без float(целого): иначе inf
                dev[name] = max(dev[name], abs(hh - hx))
    res['dev_sample'] = dev
    return res

SLOPES = ['165/439', '422/441', '425/441', '333/442', '111/445', '329/449', '438/449', '126/451',
          '172/451', '295/452', '305/461', '236/475', '165/493', '52/499', '129/499']

if __name__ == '__main__':
    BR = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    only = sys.argv[2].split(',') if len(sys.argv) > 2 else None
    t00 = time.time()
    log('=== старт: независимая проверка, перебор-контроль |a|,b ≤ %d' % BR)
    summary = {}
    for sl in SLOPES:
        if only and sl not in only: continue
        r, s = map(int, sl.split('/'))
        fab = json.load(open('/home/kep/magicKube/joint_sieve/dem_%d_%d.json' % (r, s)))
        prev = json.load(open(os.path.join(PARENT, 'rig_%d_%d.json' % (r, s))))
        prevT = {tuple(x['T']): x for x in prev['results']}
        results = []
        for x in fab['results']:
            T = [int(v) for v in x['T']]
            pv = prevT[tuple(T)]
            f, L = model3(T)
            E1ainvs = [0, f[2], 0, f[1], f[0]]
            try:
                sg = sage_part(T, pv['Emin'], pv['P0'], E1ainvs)
                res, vals = check_factor(r, s, T, pv['P0'], pv['Emin'], sg)
                res['ctrl_brute'] = brute_ctrl(T, vals, BR)
                if res['ctrl_brute']['missing']: res['closed'] = False
                res['prev_run'] = dict(M0=pv.get('M0'), B=pv.get('B_rig'), model=pv.get('model_used'), closed=pv.get('closed'))
                res['fable'] = dict(M0=x.get('M0'), B_silverman=x.get('B'), hP0=x.get('hP0'), rank=x.get('rank'), tors=x.get('tors'))
            except Exception:
                import traceback
                res = dict(T=T, error=traceback.format_exc()[-1500:], closed=False)
            short = {k: res.get(k) for k in ('ellrank', 'mwrank_bound', 'tors_order', 'tors_gcd_bound', 'tors_proved', 'model', 'B',
                                             'B_silverman_sage_min', 'c1', 'C', 'hP0_lo', 'hP0_sage', 'norm_ok', 'M0', 'index_bound',
                                             'saturation_own', 'sat_index_eclib', 'nA', 'n_values', 'nondegenerate', 'caseB_zero_poly',
                                             'dev_within_B', 'closed', 'error', 'sec')}
            log('%s T=%s: %s' % (sl, T, json.dumps(short, ensure_ascii=False, default=str)))
            if 'ctrl_brute' in res:
                log('    перебор-контроль: симметричных %s, пропущено %s; за M0: %s; сравнение: прошлый заход %s, Fable %s' % (
                    res['ctrl_brute']['symmetric_found'], res['ctrl_brute']['missing'], json.dumps(res['ctrl_beyond']),
                    res['prev_run'], res['fable']))
            results.append(res)
        closed = [x['T'] for x in results if x.get('closed')]
        json.dump(dict(slope=sl, results=results, closed_by=closed, closed_rigorous=bool(closed)),
                  open(os.path.join(HERE, 'indep_%d_%d.json' % (r, s)), 'w'), ensure_ascii=False, indent=1, default=str)
        summary[sl] = dict(closed_by=closed, n_factors=len(results),
                           M0=[x.get('M0') for x in results], B=[x.get('B') for x in results],
                           alert=[x.get('nondegenerate') for x in results if x.get('nondegenerate')])
        log('наклон %s: закрыт строго множителями %s (%.0f с с начала)' % (sl, closed, time.time() - t00))
    json.dump(summary, open(os.path.join(HERE, 'indep_summary%s.json' % ('' if not only else '_part')), 'w'), ensure_ascii=False, indent=1)
    log('=== конец, %.0f с' % (time.time() - t00))
