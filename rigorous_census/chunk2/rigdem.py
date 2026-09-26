# Claude, 26.09.2026. Строгая повторная проверка закрытий наклонов методом Демьяненко–Манина (ранг 1),
# с константой |ĥ − h(x)| ≤ B, ВЫВЕДЕННОЙ из сертификата Безу для удвоения (без silverman_height_bound).
# Код Fable (joint_sieve/demjanenko.py) не открывался и не импортируется.
#
# Нормировка: ĥ(P) = lim 4^(-k) h(x(2^k P)), h(a/b) = log max(|a|,|b|) (= Sage P.height(), сверяется).
# Схема для одного множителя T (|T| = 3 или 4, T ≠ −T), E_T ранга 1:
#   (1) модель E1: V² = f(X), f монична, z = μ(X) — мёбиусово;
#   (2) минимальная модель Emin, x_min ↔ X аффинно (подгонка по точкам + точная проверка);
#   (3) удвоение x ↦ F(x)/G(x) на выбранной модели; тождества Безу A_j F + B_j G = R t^j (j = 0, 2d−1)
#       ⇒ −log K ≤ h(x(2P)) − 4h(x(P)) ≤ log L ⇒ |ĥ − h(x)| ≤ B := log max(K, L)/3 (телескоп);
#   (4) c₁ = log max(строчные и столбцовые суммы |M|), M — примитивная целая матрица z ↦ x;
#   (5) ĥ(P₀) — строгая НИЖНЯЯ оценка из того же телескопа (k шагов удвоения);
#   (6) M₀ = ⌊(C/ĥ₀ + 1)/2⌋, C = 2B + 2c₁; случай A — z(nP₀+t), 1 ≤ n ≤ M₀, и z(t);
#       случай B — рациональные корни A² − fB² (сдвиг на кручение переводит z в −z);
#   (7) каждое значение v проверяется: ±v, все восемь клеток 1 + λv.
from sage.all import *
import sys, json, time, math
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt

RIF = RealIntervalField(256)
LOGF = None

def log(msg):
    s = time.strftime('%H:%M:%S') + ' ' + str(msg)
    print(s, flush=True)
    if LOGF:
        LOGF.write(s + '\n'); LOGF.flush()

def issq(q):
    q = QQ(q)
    return q >= 0 and q.is_square()

def cells(r, s):
    return [s, -s, r, -r, s - r, r - s, s + r, -s - r]

def full_ok(v, S):
    return all(issq(1 + l * v) for l in S)

def hgt(q):
    """логарифмическая высота рационального числа как интервал"""
    q = QQ(q)
    return RIF(max(abs(q.numerator()), abs(q.denominator()))).log()

def mob(M, x):
    """мёбиусово отображение по матрице M; x = None означает ∞; возвращает None для ∞"""
    a, b, c, d = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    if x is None:
        return None if c == 0 else QQ(a) / c
    den = c * x + d
    if den == 0:
        return None
    return (a * x + b) / den

def primitive_int(M):
    den = lcm([QQ(e).denominator() for e in M.list()])
    M2 = M * den
    g = gcd([ZZ(e) for e in M2.list()])
    return matrix(ZZ, 2, 2, [ZZ(e / g) for e in M2.list()])

def c1_of(M):
    """|h(Mu) − h(u)| ≤ c1: строчные суммы M (прямое) и строчные суммы adj M = столбцовые суммы M (обратное)"""
    A = [abs(e) for e in M.list()]
    rows = [A[0] + A[1], A[2] + A[3]]
    colsum = [A[0] + A[2], A[1] + A[3]]
    m = max(rows + colsum)
    return RIF(m).log(), int(m)

# ---------- модель E_T ----------
def model(T):
    T = [ZZ(t) for t in T]
    R = PolynomialRing(QQ, 'X'); X = R.gen()
    if len(T) == 3:
        L = prod(T)
        f = R.prod(X + QQ(L) / t for t in T)      # (L y)² = f(L z)
        mu = matrix(QQ, [[1, 0], [0, L]])          # z = X / L
        kind = 3; c0 = None
    else:
        c0 = min(T); others = [t for t in T if t != c0]
        L = prod(c0 - t for t in others)
        Rw = PolynomialRing(QQ, 'w'); w = Rw.gen()
        g = Rw.prod((c0 - t) * w + t * c0 for t in others)
        f = R(L ** 2 * g(X / L))
        mu = matrix(QQ, [[-1, c0 * L], [c0, 0]])   # z = −1/c0 + L/X
        kind = 4
    assert f.is_monic() and f.degree() == 3
    # проверка модели: для 5 значений z отношение ∏(1+λz)/f(μ⁻¹(z)) — квадрат рационального
    muinv = mu.adjugate()
    for zz in [QQ(1) / 7, QQ(-3) / 11, QQ(5) / 13, QQ(2) / 17, QQ(-9) / 19]:
        Xv = mob(muinv, zz)
        lhs = prod(1 + t * zz for t in T)
        if lhs == 0 or Xv is None:
            continue
        assert issq(f(Xv) / lhs), ('модель', T, zz)
    co = f.list()
    E1 = EllipticCurve(QQ, [0, co[2], 0, co[1], co[0]])
    return dict(T=T, E1=E1, f=f, mu=mu, kind=kind, c0=c0, L=L)

def zval(Mz, P):
    """z(P) по мёбиусу от абсциссы; P = O → Mz(∞)"""
    if P.is_zero():
        return mob(Mz, None)
    return mob(Mz, P[0])

# ---------- сертификат Безу ----------
def bezout_cert(Fq, Gq):
    Rq = Fq.parent(); t = Rq.gen()
    den = lcm([c.denominator() for c in Fq.list() + Gq.list()])
    F = Fq * den; G = Gq * den
    ct = gcd([ZZ(c) for c in F.list() + G.list()])
    F = F / ct; G = G / ct
    assert all(c in ZZ for c in F.list() + G.list())
    assert gcd(F, G) == 1
    d = max(F.degree(), G.degree())
    # линейное отображение (A, B) ↦ A F + B G, deg A, B ≤ d−1; базис коэффициентов 1..t^(2d−1)
    colsv = []
    for H in (F, G):
        for i in range(d):
            v = (t ** i * H).list()
            colsv.append(v + [QQ(0)] * (2 * d - len(v)))
    Syl = matrix(QQ, colsv).transpose()
    assert Syl.nrows() == 2 * d and Syl.det() != 0
    sols = {}
    for j in (0, 2 * d - 1):
        e = vector(QQ, [1 if i == j else 0 for i in range(2 * d)])
        sols[j] = Syl.solve_right(e)
    Rres = lcm([c.denominator() for j in sols for c in sols[j]])
    wit = []
    for j in (0, 2 * d - 1):
        v = Rres * sols[j]
        assert all(c in ZZ for c in v)
        A = Rq(list(v[:d])); B = Rq(list(v[d:]))
        assert A.degree() <= d - 1 and B.degree() <= d - 1
        assert A * F + B * G == Rres * t ** j
        wit.append(dict(j=int(j), A=[str(c) for c in A.list()], B=[str(c) for c in B.list()],
                        norm=int(sum(abs(c) for c in v))))
    K = max(w['norm'] for w in wit)
    Lb = max(sum(abs(c) for c in F.list()), sum(abs(c) for c in G.list()))
    return dict(d=int(d), F=[str(c) for c in F.list()], G=[str(c) for c in G.list()],
                R=str(Rres), K=int(K), L=int(Lb), witnesses=wit)

def verify_cert(c):
    """независимая перепроверка сертификата: целочисленность, степени, тождества, нормы"""
    Rq = PolynomialRing(QQ, 't'); t = Rq.gen()
    F = Rq([QQ(x) for x in c['F']]); G = Rq([QQ(x) for x in c['G']]); d = int(c['d']); R_ = ZZ(c['R'])
    assert R_ != 0 and max(F.degree(), G.degree()) == d and gcd(F, G) == 1
    assert all(x in ZZ for x in F.list() + G.list())
    assert sorted(int(w['j']) for w in c['witnesses']) == [0, 2 * d - 1]
    norms = []
    for w in c['witnesses']:
        A = Rq([QQ(x) for x in w['A']]); B = Rq([QQ(x) for x in w['B']])
        assert all(x in ZZ for x in A.list() + B.list())
        assert A.degree() <= d - 1 and B.degree() <= d - 1
        assert A * F + B * G == R_ * t ** int(w['j'])
        norms.append(sum(abs(x) for x in A.list() + B.list()))
    assert max(norms) == ZZ(c['K'])
    assert ZZ(c['L']) == max(sum(abs(x) for x in F.list()), sum(abs(x) for x in G.list()))
    return True

def dup_map(E):
    b2, b4, b6, b8 = E.b_invariants()
    Rq = PolynomialRing(QQ, 't'); t = Rq.gen()
    F = t ** 4 - b4 * t ** 2 - 2 * b6 * t - b8
    G = 4 * t ** 3 + b2 * t ** 2 + 2 * b4 * t + b6
    return F, G

def dup_bound(E):
    F, G = dup_map(E)
    c = bezout_cert(F, G)
    verify_cert(c)
    # контроль формулы удвоения на точке (если есть)
    return c

# ---------- ранг, образующая, кручение ----------
def rank_gen_tors(Emin, timeout_mw=240):
    info = {}
    er = pari(Emin).ellrank()
    lo, hi = int(er[0]), int(er[1])
    info['ellrank'] = [lo, hi]
    if (lo, hi) != (1, 1):
        return None, info
    pts = [Emin([QQ(c) for c in p]) for p in er[3]]
    pts = [p for p in pts if p.order() == oo]
    assert pts, 'ellrank [1,1] без точки бесконечного порядка'
    sat, idx, reg = Emin.saturation(pts[:1])
    P0 = sat[0]
    assert P0.order() == oo
    info['sat_index_eclib'] = int(idx)
    # граница индекса: m² ≤ ĥ(P0)/λ, λ — нижняя оценка ĥ на ВСЕХ точках бесконечного порядка
    # (Sage height_function().min = Cremona–Siksek min_gr / lcm(тамагава-экспонент)²) — (ПО)
    lam = Emin.height_function().min(0.0001, 20)
    hP = P0.height()
    mmax = int(floor(sqrt(hP / lam)))
    info['lambda_all'] = float(lam); info['index_bound'] = mmax
    ps = list(prime_range(2, mmax + 1))
    info['sat_primes'] = [int(p) for p in ps]
    # (а) PARI ellsaturation на простых ≤ mmax
    if ps:
        V = pari(Emin).ellsaturation([list(map(pari, P0.xy()))], ps[-1])
        Pp = Emin([QQ(c) for c in V[0]])
        info['pari_sat_same'] = bool((Pp - P0).order() != oo or (Pp + P0).order() != oo)
    else:
        info['pari_sat_same'] = True
    # (б) своя p-насыщенность редукцией: для каждого p ≤ mmax и t ∈ Tors найти ℓ, где P0 + t ∉ pE(F_ℓ)
    tors_ = Emin.torsion_points()
    own = {}
    D = Emin.discriminant()
    for p in ps:
        okp = True
        for tt in tors_:
            R_ = P0 + tt
            den = R_[0].denominator() * R_[1].denominator()
            hit = None
            for ell in prime_range(3, 5000):
                if D % ell == 0 or den % ell == 0:
                    continue
                El = Emin.change_ring(GF(ell))
                if El.cardinality() % p != 0:
                    continue
                Rl = El(R_[0], R_[1])
                if not Rl.division_points(p):
                    hit = int(ell); break
            if hit is None:
                okp = False
        own[str(p)] = okp
    info['own_psat'] = own
    info['own_psat_all'] = all(own.values())
    # вторая реализация верхней границы ранга: mwrank (eclib)
    try:
        alarm(timeout_mw)
        mw = Emin.mwrank_curve()
        info['mwrank_rank'] = int(mw.rank()); info['mwrank_certain'] = bool(mw.certain())
        cancel_alarm()
    except (AlarmInterrupt, Exception) as ex:
        cancel_alarm()
        info['mwrank_err'] = repr(ex)[:120]
    tors = Emin.torsion_points()
    info['tors_order'] = len(tors)
    info['tors_pari'] = int(pari(Emin).elltors()[0])
    # независимая верхняя граница кручения: НОД #E(F_p) по хорошим p, 3 ≤ p < 600
    Nd = ZZ(0)
    D = Emin.discriminant()
    for p in prime_range(3, 600):
        if D % p == 0:
            continue
        Nd = gcd(Nd, Emin.change_ring(GF(p)).cardinality())
    info['tors_gcd_Fp'] = int(Nd)
    assert info['tors_pari'] == len(tors)
    return (P0, tors), info

# ---------- основная проверка множителя ----------
def check_factor(r, s, T, kdup=6, ctrl_brute=0):
    t0 = time.time()
    S = cells(r, s)
    out = dict(T=[int(x) for x in T])
    if sorted(T) == sorted(-x for x in T):
        out['skip'] = 'T = −T'; return out
    md = model(T)
    E1, mu, f = md['E1'], md['mu'], md['f']
    if E1.discriminant() == 0:
        out['skip'] = 'вырождена'; return out
    Emin = E1.minimal_model()
    iso = E1.isomorphism_to(Emin); inv = Emin.isomorphism_to(E1)
    rg, info = rank_gen_tors(Emin)
    out.update(info)
    if rg is None:
        out['skip'] = 'ранг не [1,1]'; return out
    P0, tors = rg
    out['P0'] = str(P0.xy())
    out['Emin'] = [str(a) for a in Emin.ainvs()]
    out['E1'] = [str(a) for a in E1.ainvs()]
    # аффинная связь X = a x_min + b: подгонка по двум точкам, проверка на всех
    pts_min = [P0, 2 * P0, 3 * P0] + [P0 + tt for tt in tors if not tt.is_zero()]
    pts_min = [p for p in pts_min if not p.is_zero()]
    Xs = [(p[0], inv(p)[0]) for p in pts_min]
    (x1, X1), (x2, X2) = Xs[0], Xs[1]
    a = (X1 - X2) / (x1 - x2); b = X1 - a * x1
    assert all(Xv == a * xv + b for xv, Xv in Xs)
    assert a.is_square()
    # матрица z ↦ по x_min
    Mmin = primitive_int(mu * matrix(QQ, [[a, b], [0, 1]]))
    M1 = primitive_int(mu)
    for p in pts_min:
        assert mob(Mmin, p[0]) == mob(mu, inv(p)[0])
    # границы
    best = None
    variants = {}
    for name, EE, MM in (('min', Emin, Mmin), ('E1', E1, M1)):
        cert = dup_bound(EE)
        Bd = max(cert['K'], cert['L'])
        B = RIF(Bd).log() / 3
        BK = RIF(cert['K']).log() / 3; BL = RIF(cert['L']).log() / 3
        c1, cm = c1_of(MM)
        C = 2 * B + 2 * c1
        variants[name] = dict(K=str(cert['K']), L=str(cert['L']), R=str(cert['R']),
                              B_upper=float(B.upper()), BK=float(BK.upper()), BL=float(BL.upper()),
                              c1_upper=float(c1.upper()), M=[[int(e) for e in row] for row in MM.rows()],
                              C_upper=float(C.upper()), cert=cert)
        if best is None or C.upper() < best[3].upper():
            best = (name, EE, MM, C, B, c1)
    name, EE, MM, C, B, c1 = best
    out['model_used'] = name
    out['variants'] = variants
    # проверка сертификата удвоения на точках (формула x(2P) = F/G)
    cert = variants[name]['cert']
    Rq = PolynomialRing(QQ, 't')
    Fp = Rq([QQ(x) for x in cert['F']]); Gp = Rq([QQ(x) for x in cert['G']])
    # точки на выбранной модели
    if name == 'min':
        P0m = P0; torsm = tors; toE = lambda P: P
    else:
        P0m = inv(P0); torsm = [inv(tt) for tt in tors]; toE = lambda P: inv(P)
    for Q in [P0m, 2 * P0m, P0m + torsm[-1]]:
        if Q.is_zero():
            continue
        assert Fp(Q[0]) / Gp(Q[0]) == (2 * Q)[0]
    # строгая нижняя (и верхняя) оценка ĥ(P0) телескопом
    Q = P0m
    for _ in range(kdup):
        Q = 2 * Q
    hk = hgt(Q[0])
    errK = RIF(cert['K']).log() / (3 * 4 ** kdup)
    errL = RIF(cert['L']).log() / (3 * 4 ** kdup)
    h_lo = hk / 4 ** kdup - errK        # ĥ ≥ h_k/4^k − log K/(3·4^k)
    h_hi = hk / 4 ** kdup + errL
    h_lo = h_lo.lower(); h_hi = h_hi.upper()
    hs = P0.height()
    out['hP0_lower'] = float(h_lo); out['hP0_upper'] = float(h_hi); out['hP0_sage'] = float(hs)
    out['norm_ok'] = bool(h_lo - 1e-9 <= hs <= h_hi + 1e-9)
    assert h_lo > 0
    Cup = RIF(C.upper())
    M0 = int(floor(((Cup / RIF(h_lo) + 1) / 2).upper()))
    out['B_rig'] = float(B.upper()); out['c1'] = float(c1.upper()); out['C'] = float(C.upper()); out['M0'] = M0
    # сравнение с Sage
    try:
        out['B_silverman_sage'] = float(Emin.silverman_height_bound())
    except Exception as ex:
        out['B_silverman_sage'] = repr(ex)[:80]
    # контроль: |ĥ − h(x)| ≤ B на кратных (численно)
    worst = 0.0
    Q = P0m
    for n in range(1, 9):
        for tt in torsm:
            R_ = Q + tt
            if R_.is_zero():
                continue
            dlt = abs(float(R_.height()) - float(hgt(R_[0]).center()))
            worst = max(worst, dlt)
        Q = Q + P0m
    out['max_dev_sample'] = worst
    out['dev_within_B'] = bool(worst <= float(B.upper()))
    # ---------- случай A и кручение (на модели E1, z = μ(X)) ----------
    P01 = inv(P0); tors1 = [inv(tt) for tt in tors]
    vals = set()
    bad = []
    def consider(v, tag):
        if v is None:
            return
        v = QQ(v)
        vals.add(v); vals.add(-v)
        if v != 0 and (full_ok(v, S) or full_ok(-v, S)):
            bad.append([tag, str(v)])
    for tt in tors1:
        consider(zval(mu, tt), 'tors')
    Q = P01 * 0
    nA = 0
    for n in range(1, M0 + 1):
        Q = Q + P01
        for tt in tors1:
            R_ = Q + tt
            consider(zval(mu, R_), 'A n=%d' % n); nA += 1
        if n % 20 == 0:
            log('   T=%s: случай A n=%d/%d' % (out['T'], n, M0))
    out['nA'] = nA
    # ---------- случай B ----------
    a2 = E1.a2()
    Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
    fx = Rx(f)
    al, be, ga, de = [QQ(e) for e in mu.list()]
    lin = 2 * al * ga * x + al * de + be * ga
    caseB = []
    zero_poly = []
    for tt in tors1:
        if tt.is_zero():
            continue
        xt, yt = tt[0], tt[1]
        A = (fx + yt ** 2 - (a2 + x + xt) * (x - xt) ** 2) * lin + (x - xt) ** 2 * ((al * de + be * ga) * x + 2 * be * de)
        Bp = -2 * yt * lin
        Qp = A if Bp == 0 else A ** 2 - fx * Bp ** 2
        if Qp == 0:
            zero_poly.append(str(tt.xy())); continue
        for rt, _ in Qp.roots(QQ):
            if rt == xt:
                continue
            fv = fx(rt)
            if not issq(fv):
                continue
            yv = fv.sqrt()
            for yy in {yv, -yv}:
                P = E1(rt, yy)
                zP = zval(mu, P); zQ = zval(mu, P + tt)
                caseB.append(dict(t=str(tt.xy()), x=str(rt), z=str(zP), z_shift=str(zQ),
                                  eq=bool(zP is not None and zQ is not None and zQ == -zP)))
                consider(zP, 'B')
    out['caseB'] = caseB
    # ---------- контроль: кратные за пределами M0 (n ≤ 2 M0 + 6): симметричные z обязаны лежать в vals ----------
    Tl = [QQ(tq) for tq in md['T']]
    beyond = []; beyond_missing = []
    Q = P01 * 0
    for n in range(1, 2 * M0 + 7):
        Q = Q + P01
        for tt in tors1:
            v = zval(mu, Q + tt)
            if v is None:
                continue
            if issq(prod(1 - lam_ * v for lam_ in Tl)):
                beyond.append([n, str(v)])
                if v not in vals:
                    beyond_missing.append([n, str(v)])
    out['ctrl_beyond'] = dict(nmax=2 * M0 + 6, symmetric=beyond, missing=beyond_missing)
    out['caseB_zero_poly'] = zero_poly
    out['n_values'] = len(vals)
    out['nondegenerate'] = bad
    out['closed_algebraic'] = (not bad) and (not zero_poly)
    out['aux_ok'] = bool(out.get('own_psat_all') and out.get('pari_sat_same') and out['tors_gcd_Fp'] == out['tors_order']
                         and out['norm_ok'] and out['dev_within_B'] and not out['ctrl_beyond']['missing'])
    out['sec'] = round(time.time() - t0, 1)
    # ---------- контроль: перебор малых z с P_T(z), P_T(−z) рациональными — все должны быть в vals ----------
    if ctrl_brute:
        Tl = [int(tq) for tq in md['T']]
        k = len(Tl)
        found = []; missing = []
        for bb in range(1, ctrl_brute + 1):
            for aa in range(-ctrl_brute, ctrl_brute + 1):
                if math.gcd(aa, bb) != 1:
                    continue
                p1 = 1; p2 = 1
                for lam in Tl:
                    p1 *= bb + lam * aa; p2 *= bb - lam * aa
                if k == 3:
                    p1 *= bb; p2 *= bb
                if p1 < 0 or p2 < 0:
                    continue
                if math.isqrt(p1) ** 2 == p1 and math.isqrt(p2) ** 2 == p2:
                    v = QQ(aa) / bb
                    found.append(str(v))
                    if v != 0 and v not in vals:
                        missing.append(str(v))
        out['ctrl_brute'] = dict(bound=ctrl_brute, symmetric_found=found, missing=missing)
    out['closed'] = bool(out['closed_algebraic'] and out['aux_ok'] and not out.get('ctrl_brute', {}).get('missing'))
    return out

def algebra_ctrl_Fp(T, p):
    """контроль алгебры случая B: над F_p множество точек P с z(P+t) = −z(P) ⊂ корни Q mod p"""
    md = model(T); E1, mu, f = md['E1'], md['mu'], md['f']
    tors1 = E1.torsion_points()
    Ep = E1.change_ring(GF(p)); F = GF(p)
    al, be, ga, de = [F(QQ(e)) for e in mu.list()]
    a2 = F(E1.a2())
    Rx = PolynomialRing(F, 'x'); x = Rx.gen(); fx = Rx(f.change_ring(F))
    lin = 2 * al * ga * x + al * de + be * ga
    res = []
    for tt in tors1:
        if tt.is_zero():
            continue
        xt, yt = F(tt[0]), F(tt[1])
        A = (fx + yt ** 2 - (a2 + x + xt) * (x - xt) ** 2) * lin + (x - xt) ** 2 * ((al * de + be * ga) * x + 2 * be * de)
        Bp = -2 * yt * lin
        Qp = A if Bp == 0 else A ** 2 - fx * Bp ** 2
        tp = Ep(xt, yt)
        cnt = 0; ok = True
        for P in Ep.points():
            if P.is_zero():
                continue
            P2 = P + tp
            if P2.is_zero() or P[0] == xt:
                continue
            d1 = ga * P[0] + de; d2 = ga * P2[0] + de
            if d1 == 0 or d2 == 0:
                continue
            if (al * P2[0] + be) / d2 == -(al * P[0] + be) / d1:
                cnt += 1
                if Qp(P[0]) != 0:
                    ok = False
        res.append(dict(t=str(tt.xy()), matches=cnt, all_roots=ok))
    return res

if __name__ == '__main__':
    import os
    slope = sys.argv[1]
    Tlist = json.loads(sys.argv[2])
    brute = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    r, s = map(int, slope.split('/'))
    here = os.path.dirname(os.path.abspath(__file__))
    LOGF = open(os.path.join(here, 'log_%d_%d.txt' % (r, s)), 'a')
    log('наклон %s, множители %s' % (slope, Tlist))
    results = []
    for T in Tlist:
        log(' множитель T=%s' % T)
        try:
            res = check_factor(r, s, T, ctrl_brute=brute)
        except Exception as ex:
            import traceback
            res = dict(T=T, error=traceback.format_exc()[-1500:])
        short = {k: res.get(k) for k in ('ellrank', 'M0', 'B_rig', 'B_silverman_sage', 'c1', 'C', 'hP0_lower', 'hP0_sage',
                                          'norm_ok', 'nA', 'n_values', 'nondegenerate', 'caseB_zero_poly', 'closed', 'skip', 'error', 'sec')}
        log('  итог: ' + json.dumps(short, ensure_ascii=False))
        if 'ctrl_brute' in res:
            log('  контроль перебором: найдено %d симметричных, пропущено %s' % (len(res['ctrl_brute']['symmetric_found']), res['ctrl_brute']['missing']))
        results.append(res)
    closed = any(r_.get('closed') for r_ in results)
    alert = [r_['nondegenerate'] for r_ in results if r_.get('nondegenerate')]
    json.dump(dict(slope=slope, results=results, closed_rigorous=closed, alert=alert),
              open(os.path.join(here, 'rig_%d_%d.json' % (r, s)), 'w'), ensure_ascii=False, indent=1, default=str)
    log('наклон %s: закрыт строго = %s; ALERT = %s' % (slope, closed, alert))
