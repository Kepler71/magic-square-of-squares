# Claude, 26.09.2026. Строгая перепроверка закрытий наклонов методом Демьяненко–Манина (ранг 1)
# с константой |ĥ − h(x)| из сертификата Безу для отображения удвоения (метод Codex, переписан заново).
# Код Fable (joint_sieve/demjanenko.py, jsieve.py) НЕ импортируется. Из dem_<r>_<s>.json берутся только наборы T.
#
# Нормировка: h(x) = log max(|num|,|den|) для x ∈ P¹(Q), h(∞)=0; ĥ(P) = lim 4^(−n) h(x(2ⁿP)) (совпадает с Sage P.height(),
# проверяется строгим интервалом на каждой кривой).
#
# Лемма (Безу). F,G ∈ Z[X,Z] формы степени d, A_j F + B_j G = R·X^j Z^(2d−1−j) (j = 0, 2d−1), A_j,B_j ∈ Z[X,Z] степени d−1, R ≠ 0.
#   K = max_j(‖A_j‖₁+‖B_j‖₁), L = max(‖F‖₁,‖G‖₁). Тогда для всех [m:n] ∈ P¹(Q): d·h − log K ≤ h([F:G]) ≤ d·h + log L.
# Для удвоения (d=4): h(x(2P)) − 4h(x(P)) ∈ [−log K, log L]; телескоп ⇒ ĥ(P) − h(x(P)) ∈ [−log K/3, log L/3].
# B := log max(K,L)/3 — симметричная строгая граница |ĥ − h(x)| ≤ B.
from sage.all import *
import sys, json, time, os

RIF = RealIntervalField(300)
OUT = '/home/kep/magicKube/rigorous_census/chunk1'


def say(*a):
    print(time.strftime('%H:%M:%S'), *a, flush=True)


def cells(r, s):
    base = [s, r, s - r, s + r]
    return sorted(set(base + [-c for c in base]))


def all_sq(z, S):
    z = QQ(z)
    for l in S:
        v = 1 + l * z
        if v < 0 or not v.is_square():
            return False
    return True


def Hproj(x):
    """Высота точки P¹(Q) как целое: x=None означает ∞."""
    if x is None:
        return ZZ(1)
    x = QQ(x)
    return max(abs(x.numerator()), abs(x.denominator()))


def xcoord(P):
    return None if P.is_zero() else P[0]


# ---------- 1. модель кривой E_T (своя) ----------
def model(T):
    """C_T: y² = ∏_{c∈T}(1+cz), |T| ∈ {3,4}. w = c0/(c0 z+1), X = L w: Y² = X³+a2X²+a4X+a6 (целые коэффициенты).
    Возвращает E1, L, c0. Тождество h(X(z)) = (L c0³/(c0z+1)²)² ∏(1+cz) проверяется в Q(z)."""
    T = [ZZ(c) for c in T]
    c0 = T[0]
    Rw = PolynomialRing(QQ, 'w'); w = Rw.gen()
    g = Rw.prod((c0 - c) * w + c * c0 for c in T[1:])
    if len(T) == 3:
        g = c0 * w * g
    assert g.degree() == 3
    L = g.leading_coefficient()
    g0, g1, g2 = g[0], g[1], g[2]
    a2, a4, a6 = g2, g1 * L, g0 * L**2
    assert all(v in ZZ for v in (a2, a4, a6))
    E1 = EllipticCurve(QQ, [0, a2, 0, a4, a6])
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Fz = Rz.fraction_field()
    X = Fz(L * c0 / (c0 * z + 1))
    lhs = X**3 + a2 * X**2 + a4 * X + a6
    rhs = (L * c0**3 / (c0 * z + 1)**2)**2 * Rz.prod(1 + c * z for c in T)
    assert lhs == rhs, 'тождество модели не выполнено'
    return E1, QQ(L), QQ(c0)


def mobius(E1, Emin, iso, L, c0):
    """x_min = (X − r)/u² (конвенция Sage x = u²x' + r), X = L c0/(c0 z + 1) ⇒ x_min = M·z, M целая примитивная."""
    u, r_, s_, t_ = iso.tuple()
    M = matrix(QQ, [[-r_ * c0, L * c0 - r_], [u**2 * c0, u**2]])
    den = lcm([e.denominator() for e in M.list()])
    M = M * den
    gg = gcd([ZZ(e) for e in M.list()])
    M = (M / gg).change_ring(ZZ)
    assert M.det() != 0
    a, b, c, d = M.list()
    c1_int = max(abs(a) + abs(b), abs(c) + abs(d), abs(d) + abs(b), abs(c) + abs(a))   # строки M и adj(M)
    return M, ZZ(c1_int), (u, r_, s_, t_)


def z_of_xmin(M, x):
    """z = adj(M)·x; x=None — точка O. Возвращает None, если z = ∞."""
    a, b, c, d = M.list()
    if x is None:
        num, den = QQ(d), QQ(-c)
    else:
        num, den = d * x - b, -c * x + a
    return None if den == 0 else QQ(num / den)


def xmin_of_z(M, z):
    a, b, c, d = M.list()
    den = c * z + d
    return None if den == 0 else QQ((a * z + b) / den)


# ---------- 2. сертификат Безу для удвоения ----------
def dup_forms(E):
    b2, b4, b6, b8 = [ZZ(v) for v in E.b_invariants()]
    Rt = PolynomialRing(ZZ, 't'); t = Rt.gen()
    F = t**4 - b4 * t**2 - 2 * b6 * t - b8
    G = 4 * t**3 + b2 * t**2 + 2 * b4 * t + b6
    return F, G


def bezout_cert(F, G, d):
    Rt = F.parent(); t = Rt.gen(); n = 2 * d
    cols = []
    for H in (F, G):
        for i in range(d):
            v = list((t**i * H).list()); v += [0] * (n - len(v))
            assert len(v) == n
            cols.append(v)
    S = matrix(QQ, cols).transpose()
    assert S.det() != 0
    sols = {}
    for j in (0, n - 1):
        e = vector(QQ, [1 if i == j else 0 for i in range(n)])
        sols[j] = S.solve_right(e)
    R = ZZ(lcm([v.denominator() for j in sols for v in sols[j]]))
    wit = []
    for j in (0, n - 1):
        sol = sols[j] * R
        assert all(v in ZZ for v in sol)
        A = Rt([ZZ(v) for v in sol[:d]]); B = Rt([ZZ(v) for v in sol[d:]])
        assert A * F + B * G == R * t**j
        wit.append(dict(j=j, A=[str(v) for v in A.list()], B=[str(v) for v in B.list()],
                        norm=str(sum(abs(v) for v in sol))))
    K = max(ZZ(w_['norm']) for w_ in wit)
    Lc = max(sum(abs(v) for v in F.list()), sum(abs(v) for v in G.list()))
    return dict(d=d, F=[str(v) for v in F.list()], G=[str(v) for v in G.list()], R=str(R), K=str(K), L=str(Lc), witnesses=wit)


def verify_cert(cert, E):
    """Независимая от построения проверка: формы совпадают с формулой удвоения E; тождества Безу точные; K, L пересчитаны."""
    Rt = PolynomialRing(ZZ, 't'); t = Rt.gen()
    F = Rt([ZZ(v) for v in cert['F']]); G = Rt([ZZ(v) for v in cert['G']]); d = cert['d']; R = ZZ(cert['R'])
    F0, G0 = dup_forms(E)
    assert F == F0 and G == G0 and d == 4 and R != 0
    assert max(F.degree(), G.degree()) == d and F.leading_coefficient() == 1 and F.degree() == d
    norms = []
    for w_ in cert['witnesses']:
        A = Rt([ZZ(v) for v in w_['A']]); B = Rt([ZZ(v) for v in w_['B']])
        assert A.degree() <= d - 1 and B.degree() <= d - 1
        assert A * F + B * G == R * t**w_['j']
        norms.append(sum(abs(v) for v in A.list() + B.list()))
    assert sorted(w_['j'] for w_ in cert['witnesses']) == [0, 2 * d - 1]
    assert ZZ(cert['K']) == max(norms)
    assert ZZ(cert['L']) == max(sum(abs(v) for v in F.list()), sum(abs(v) for v in G.list()))
    return True


def check_dup_formula(E, pts):
    F, G = dup_forms(E)
    for P in pts:
        if P.is_zero():
            continue
        x = P[0]
        Gx = G(x)
        Q = 2 * P
        if Gx == 0:
            assert Q.is_zero()
        else:
            assert Q[0] == F(x) / Gx
    return True


def hhat_interval(P, k, logK, logL):
    Q = P
    for _ in range(k):
        Q = 2 * Q
    hx = RIF(Hproj(xcoord(Q))).log()
    lo = (hx - logK / 3) / 4**k
    hi = (hx + logL / 3) / 4**k
    return RIF(lo.lower(), hi.upper())


# ---------- 3. кручение, ранг, образующая ----------
def torsion_gcd(E, pmax=400):
    D = E.discriminant(); g = ZZ(0); used = 0
    for p in primes(3, pmax):
        if D % p == 0:
            continue
        g = gcd(g, E.change_ring(GF(p)).cardinality()); used += 1
    return g, used


def reduce_pt(Eq, P, q):
    """Редукция точки минимальной модели по модулю q (q хорошей редукции); None — если точка уходит в O."""
    if P.is_zero():
        return Eq(0)
    x, y = P.xy()
    if x.denominator() % q == 0:
        return Eq(0)
    return Eq(GF(q)(x), GF(q)(y))


def my_saturation(Emin, P0, tors, pmax, qmax=20000):
    """Для каждого простого p ≤ pmax ищется q хорошей редукции с p | #E(F_q), при котором P0 − t ∉ pE(F_q) для всех t ∈ Tors.
    Тогда P0 ∉ pE(Q) + Tors. Возвращает список непроверенных p."""
    D = Emin.discriminant(); bad = []; info = {}
    for p in primes(2, pmax + 1):
        done = False
        for q in primes(3, qmax):
            if D % q == 0:
                continue
            Eq = Emin.change_ring(GF(q)); N = Eq.cardinality()
            if N % p:
                continue
            ok = True
            for t in tors:
                Rq = reduce_pt(Eq, P0 - t, q)
                if Rq.is_zero() or Rq.division_points(p):
                    ok = False; break
            if ok:
                info[int(p)] = int(q); done = True; break
        if not done:
            bad.append(int(p))
    return bad, info


# ---------- 4. случай B ----------
def caseB(E1, L, c0, tors1):
    a1, a2, a3, a4, a6 = E1.a_invariants(); assert a1 == 0 and a3 == 0
    Rxy = PolynomialRing(QQ, ['X', 'Y']); X, Y = Rxy.gens()
    f = X**3 + a2 * X**2 + a4 * X + a6
    Rx = PolynomialRing(QQ, 'x')
    fx = Rx([a6, a4, a2, 1])
    out = []; polys = {}
    for t in tors1:
        if t.is_zero():
            continue
        xt, yt = t.xy()
        D = (X - xt)**2
        XpD = f - 2 * yt * Y + yt**2 - (a2 + X + xt) * D        # X(P+t)·D при Y² = f
        cond = c0 * L * (X * D + XpD) - 2 * X * XpD              # D·[c0 L (X+X') − 2 X X']
        assert cond.degree(Y) <= 1
        def toRx(pp):
            return Rx({e[0]: cf for e, cf in pp.dict().items()}) if pp != 0 else Rx(0)
        A = toRx(cond.coefficient({Y: 0})); B = toRx(cond.coefficient({Y: 1}))
        assert cond == A(X) + B(X) * Y
        pol = A if B == 0 else A**2 - B**2 * fx
        if pol == 0:
            out.append(dict(t=str(t), identically_zero=True)); continue
        roots = [rt for rt, _ in pol.roots(QQ)]
        polys[str(t)] = (A, B, xt, yt)
        zs = []
        for rt in roots:
            if rt == 0:
                continue
            zs.append(QQ(L / rt - 1 / c0))
        out.append(dict(t=str(t), identically_zero=False, deg=int(pol.degree()), roots=[str(v) for v in roots], z=[str(v) for v in zs]))
    return out, polys


def caseB_control(E1, L, c0, polys, primes_=(1009, 2003)):
    """Над F_p: A(X)+B(X)Y ≡ (X−xt)²·[c0L(X+X') − 2XX'] для всех точек P (X ≠ xt), X' = x(P+t) по закону группы Sage."""
    nchk = 0
    D1 = E1.discriminant()
    for p in primes_:
        if D1 % p == 0 or L.numerator() % p == 0 or L.denominator() % p == 0 or c0 % p == 0:
            continue
        Ep = E1.change_ring(GF(p)); pts = Ep.points()
        for key, (A, B, xt, yt) in polys.items():
            if xt.denominator() % p == 0 or yt.denominator() % p == 0:
                continue
            tp = Ep(GF(p)(xt), GF(p)(yt))
            Ap = A.change_ring(GF(p)); Bp = B.change_ring(GF(p))
            for P in pts:
                if P.is_zero() or P[0] == GF(p)(xt):
                    continue
                Pp = P + tp
                if Pp.is_zero():
                    continue
                Xv, Yv = P[0], P[1]; Xp = Pp[0]
                lhs = Ap(Xv) + Bp(Xv) * Yv
                rhs = (Xv - GF(p)(xt))**2 * (GF(p)(c0 * L) * (Xv + Xp) - 2 * Xv * Xp)
                assert lhs == rhs, ('случай B: тождество нарушено', p, key)
                nchk += 1
    return nchk


# ---------- 5. один множитель ----------
def run_T(S, T, sat_pmax_cap=200, do_sallows_z=None):
    t0 = time.time()
    res = dict(T=[int(c) for c in T])
    E1, L, c0 = model(T)
    Emin = E1.minimal_model(); iso = E1.isomorphism_to(Emin); iso_inv = Emin.isomorphism_to(E1)
    res['E1'] = [str(v) for v in E1.a_invariants()]; res['Emin'] = [str(v) for v in Emin.a_invariants()]
    res['L'] = str(L); res['c0'] = str(c0); res['conductor'] = str(Emin.conductor())
    # ранг
    er = pari(Emin).ellrank()
    lo, hi = int(er[0]), int(er[1])
    res['ellrank'] = [lo, hi]
    if (lo, hi) != (1, 1):
        res['ok'] = False; res['why'] = f'ellrank {lo},{hi}'; return res
    pts = [Emin([QQ(c) for c in p]) for p in er[3]]
    tors = Emin.torsion_points()
    pts = [p for p in pts if p.order() == oo]
    assert pts
    # кручение: полное (gcd #E(F_p))
    g, used = torsion_gcd(Emin)
    res['torsion'] = dict(n=len(tors), structure=str(Emin.torsion_subgroup().invariants()), gcd_Fp=int(g), primes_used=used,
                          proven_full=(int(g) == len(tors)))
    # насыщение: eclib (как Sage), плюс своё p-насыщение до границы индекса из λ (Sage height_function().min)
    P = pts[0]
    x, y = P.xy(); dd = x.denominator().lcm(y.denominator())
    from sage.libs.eclib.all import mwrank_MordellWeil
    mw = mwrank_MordellWeil(Emin.mwrank_curve(), False)
    mw.process([(x * dd, y * dd, dd)])
    ok_sat, index, unsat = mw.saturate()
    P0 = Emin(mw.points()[0])
    assert P0.order() == oo
    res['P0'] = str(P0.xy()); res['eclib_saturation'] = dict(ok=bool(ok_sat), index=int(index), unsat=str(unsat))
    # сертификат Безу
    F, G = dup_forms(Emin)
    cert = bezout_cert(F, G, 4)
    verify_cert(cert, Emin)
    Kc, Lc = ZZ(cert['K']), ZZ(cert['L'])
    logK, logL = RIF(Kc).log(), RIF(Lc).log()
    B = RIF(max(Kc, Lc)).log() / 3
    res['cert'] = cert
    res['logK_over3'] = float(logK.upper() / 3); res['logL_over3'] = float(logL.upper() / 3); res['B_bezout'] = float(B.upper())
    sil = float(Emin.silverman_height_bound()); cps = float(Emin.CPS_height_bound())
    res['B_sage_silverman'] = sil; res['B_sage_CPS'] = cps
    # формула удвоения против закона группы Sage
    check_dup_formula(Emin, [P0, 2 * P0, 3 * P0 + tors[-1], P0 + tors[-1]] + list(tors))
    # строгий интервал ĥ(P0) и сверка нормировки с Sage
    hiv = hhat_interval(P0, 7, logK, logL)
    hs = P0.height()
    res['hhat_P0_interval'] = [float(hiv.lower()), float(hiv.upper())]; res['sage_height_P0'] = float(hs)
    res['normalization_ok'] = bool(hiv.lower() <= hs <= hiv.upper())
    assert res['normalization_ok'], ('нормировка', hiv, hs)
    # своё p-насыщение
    lam = Emin.height_function().min(0.0001, 20)
    pmax = int(floor(sqrt(float(hiv.upper()) / lam))) if lam > 0 else None
    res['lambda_sage'] = float(lam); res['index_bound'] = pmax
    if pmax is not None and pmax <= sat_pmax_cap:
        bad, info = my_saturation(Emin, P0, tors, pmax)
        res['my_saturation'] = dict(pmax=pmax, unproved_primes=bad, witnesses_q=info)
    else:
        res['my_saturation'] = dict(pmax=pmax, skipped=True)
    # Мёбиус z ↦ x_min и c1
    M, c1i, isotuple = mobius(E1, Emin, iso, L, c0)
    c1 = RIF(c1i).log()
    res['M'] = [str(v) for v in M.list()]; res['c1'] = float(c1.upper())
    # проверка Мёбиуса на точках: x_min(Q) = M·z(Q), где z(Q) из E1: z = L/X − 1/c0
    for Q in [P0, 2 * P0, P0 + tors[-1], 3 * P0] + list(tors):
        Q1 = iso_inv(Q)
        if Q1.is_zero():
            continue
        X1 = Q1[0]
        if X1 == 0:
            continue
        zz = QQ(L / X1 - 1 / c0)
        assert xmin_of_z(M, zz) == Q[0] and z_of_xmin(M, Q[0]) == zz
    # константа и M0 (строго, интервалы MPFI)
    C = 2 * B + 2 * c1
    Cup = C.upper(); hlo = hiv.lower()
    M0 = int(((RIF(Cup) / RIF(hlo) + 1) / 2).upper().floor())
    # проверка: при |n| > M0: (2|n|−1)·ĥ_lo > C_up
    assert (RIF(2 * (M0 + 1) - 1) * RIF(hlo)).lower() > Cup
    res['C'] = float(Cup); res['M0'] = M0
    res['C_onesided'] = float(((logK + logL) / 3 + 2 * c1).upper())
    # случай A: z(nP0 + t), |n| ≤ M0, t ∈ Tors.  Дополнительно (контроль): до N_big = 2·M0+10 — перепись Y-точек
    # (z, при которых ∏_T(1−λz) тоже квадрат) и численная сверка констант на кратных.
    cand = set(); sols = []; nA = 0; Ypts = []
    N_big = 2 * M0 + 10
    worst = dict(hhat_minus_hx_min=1e9, hhat_minus_hx_max=-1e9, hx_minus_hz_absmax=0.0)
    Qn_pos = Emin(0); Qn_neg = Emin(0)
    for n in range(0, N_big + 1):
        for sgn, Qn in ([(1, Qn_pos)] if n == 0 else [(1, Qn_pos), (-1, Qn_neg)]):
            for t in tors:
                Q = Qn + t
                xq = xcoord(Q)
                zz = z_of_xmin(M, xq)
                if n <= M0:
                    nA += 1
                if zz is None:
                    continue
                # контроль модели: ∏_{c∈T}(1+cz) — квадрат
                pr = prod(1 + c * zz for c in T)
                assert pr >= 0 and pr.is_square(), 'z(Q) не даёт квадрата на T'
                if n <= M0:
                    cand.add(zz)
                prm = prod(1 - c * zz for c in T)
                if zz != 0 and prm >= 0 and prm.is_square():
                    Ypts.append((sgn * n, str(t), zz))
                if n <= 12 and not Q.is_zero():
                    hx = float(RIF(Hproj(xq)).log().center()); hz = float(RIF(Hproj(zz)).log().center())
                    hh = float(Q.height()) if Q.order() == oo else 0.0
                    worst['hhat_minus_hx_min'] = min(worst['hhat_minus_hx_min'], hh - hx)
                    worst['hhat_minus_hx_max'] = max(worst['hhat_minus_hx_max'], hh - hx)
                    worst['hx_minus_hz_absmax'] = max(worst['hx_minus_hz_absmax'], abs(hx - hz))
        Qn_pos = Qn_pos + P0; Qn_neg = Qn_neg - P0
    res['nA_points'] = nA
    res['const_check_multiples_n_le_12'] = worst
    res['const_check_pass'] = bool(-float(logK.upper()) / 3 <= worst['hhat_minus_hx_min'] and worst['hhat_minus_hx_max'] <= float(logL.upper()) / 3
                                   and worst['hx_minus_hz_absmax'] <= float(c1.upper()))
    # лемма Безу на случайных x: 4h(x) − log K ≤ h(F/G) ≤ 4h(x) + log L
    set_random_seed(int(abs(hash(tuple(T)))) % 10**6)
    Fq, Gq = dup_forms(Emin); bez_ok = True
    for _ in range(300):
        xr = QQ(ZZ.random_element(-10**12, 10**12)) / ZZ.random_element(1, 10**12)
        num, den = ZZ(Fq(xr) * xr.denominator()**4), ZZ(Gq(xr) * xr.denominator()**4)
        gq = gcd(num, den); hphi = RIF(max(abs(num // gq), abs(den // gq))).log() if (num, den) != (0, 0) else None
        h4 = 4 * RIF(Hproj(xr)).log()
        if not ((h4 - logK).lower() <= hphi.upper() and hphi.lower() <= (h4 + logL).upper()):
            bez_ok = False
    res['bezout_random_x_check'] = bez_ok
    # случай B (на E1)
    tors1 = [iso_inv(t) for t in tors]
    cb, polys = caseB(E1, L, c0, tors1)
    res['caseB'] = cb
    res['caseB_identically_zero'] = any(o['identically_zero'] for o in cb)
    res['caseB_Fp_checks'] = caseB_control(E1, L, c0, polys)
    for o in cb:
        for zs in o.get('z', []):
            cand.add(QQ(zs))
    res['ncand'] = len(cand)
    res['Y_points_census'] = dict(N_big=N_big, n=[(a, b, str(c)) for a, b, c in Ypts],
                                  all_in_candidates=all((c in cand) or (-c in cand) for a, b, c in Ypts))
    assert res['Y_points_census']['all_in_candidates'], 'Y-точка вне Z_A ∪ Z_B — противоречие теореме/реализации'
    sols = sorted([zz for zz in cand if zz != 0 and (all_sq(zz, S) or all_sq(-zz, S))])
    res['nondeg_solutions'] = [str(v) for v in sols]
    if do_sallows_z is not None:
        res['control_z_found'] = (QQ(do_sallows_z) in cand) or (QQ(-do_sallows_z) in cand)
    res['ok'] = True
    res['closes'] = (not res['caseB_identically_zero']) and len(sols) == 0 and res['torsion']['proven_full'] and bool(ok_sat)
    res['time'] = float('%.1f' % (time.time() - t0))
    return res


def run_slope(r, s, Ts):
    S = cells(r, s)
    out = []
    for T in Ts:
        say(f'  {r}/{s} T={T} ...')
        try:
            res = run_T(S, T)
        except Exception as e:
            import traceback; traceback.print_exc()
            res = dict(T=T, ok=False, why='исключение: ' + repr(e))
        if res.get('ok'):
            say(f"    ранг {res['ellrank']}, tors {res['torsion']['n']} (gcd {res['torsion']['gcd_Fp']}), eclib ok={res['eclib_saturation']['ok']} idx={res['eclib_saturation']['index']}, "
                f"моё насыщ. pmax={res['my_saturation'].get('pmax')} непров.={res['my_saturation'].get('unproved_primes')}")
            say(f"    B_bez={res['B_bezout']:.2f} (Sage Silverman {res['B_sage_silverman']:.2f}, CPS {res['B_sage_CPS']:.2f}), c1={res['c1']:.2f}, "
                f"ĥ(P0)∈[{res['hhat_P0_interval'][0]:.5f},{res['hhat_P0_interval'][1]:.5f}] (Sage {res['sage_height_P0']:.5f}), C={res['C']:.2f}, M0={res['M0']}, "
                f"точек A={res['nA_points']}, канд.={res['ncand']}, B-тожд.0={res['caseB_identically_zero']}, F_p-сверок B={res['caseB_Fp_checks']}, "
                f"невырожд.={res['nondeg_solutions']}, закрывает={res['closes']}, {res['time']} с")
            say(f"    контроль констант: ĥ−h(x) ∈ [{res['const_check_multiples_n_le_12']['hhat_minus_hx_min']:.2f},{res['const_check_multiples_n_le_12']['hhat_minus_hx_max']:.2f}] "
                f"при [−{res['logK_over3']:.2f}, {res['logL_over3']:.2f}]; |h(x)−h(z)| ≤ {res['const_check_multiples_n_le_12']['hx_minus_hz_absmax']:.2f} при c1={res['c1']:.2f}; "
                f"прошло={res['const_check_pass']}; Безу на случайных x: {res['bezout_random_x_check']}; Y-точек до |n|≤{res['Y_points_census']['N_big']}: {len(res['Y_points_census']['n'])}, все в кандидатах: {res['Y_points_census']['all_in_candidates']}")
        else:
            say(f"    пропуск: {res.get('why')}")
        out.append(res)
    closed = any(o.get('ok') and o.get('closes') for o in out)
    alert = [o['nondeg_solutions'] for o in out if o.get('ok') and o['nondeg_solutions']]
    return dict(slope=f'{r}/{s}', cells=[int(c) for c in S], results=out, closed_rigorously=closed, ALERT=alert)


if __name__ == '__main__':
    say('старт', sys.argv[1:])
    if sys.argv[1] == 'sallows':
        # положительный контроль: S = {±825, ±578}, z* = 168/425² (все 4 клетки квадраты)
        S = [-825, -578, 578, 825]; zs = QQ(168) / 425**2
        assert all_sq(zs, S)
        out = []
        for T in ([-825, -578, 578], [-825, -578, 825], [-578, 578, 825], [-825, 578, 825]):
            try:
                res = run_T(S, T, do_sallows_z=zs)
            except Exception as e:
                import traceback; traceback.print_exc(); res = dict(T=T, ok=False, why=repr(e))
            say('  sallows', T, {k: res.get(k) for k in ('ellrank', 'M0', 'B_bezout', 'B_sage_silverman', 'ncand', 'nondeg_solutions', 'control_z_found', 'why')})
            out.append(res)
        json.dump(dict(control='sallows', S=S, z=str(zs), results=out), open(f'{OUT}/control_sallows.json', 'w'), ensure_ascii=False, indent=1, default=str)
    else:
        for sl in sys.argv[1:]:
            r, s = map(int, sl.split('/'))
            d = json.load(open(f'/home/kep/magicKube/joint_sieve/dem_{r}_{s}.json'))
            Ts = [o['T'] for o in d['results']]
            say(f'наклон {r}/{s}: T из dem_{r}_{s}.json: {Ts}')
            t0 = time.time()
            res = run_slope(r, s, Ts)
            res['time'] = float('%.1f' % (time.time() - t0))
            json.dump(res, open(f'{OUT}/rig_{r}_{s}.json', 'w'), ensure_ascii=False, indent=1, default=str)
            say(f'ИТОГ {r}/{s}: закрыт строго = {res["closed_rigorously"]}; ALERT = {res["ALERT"]}; {res["time"]} с')
    say('готово')
