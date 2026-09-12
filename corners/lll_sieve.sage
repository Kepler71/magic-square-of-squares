# Поиск типа Бремнера на парах с КВАДРАТНЫМ центром: LLL-базис, перебор по канонической высоте, модульное сито.
# Идея сита и перебора по высоте — Codex (RESULT_FROM_CODEX_2026-09-11_LLL_HEIGHT_SIEVE.md, ..._LLL_SQUARE_CENTERS.md);
# реализация своя, независимая.
#
# Модель — bremner_type.sage (BaseCurve, conds, score). Для базы (c1,c2) точки E(Q) <-> f с квадратными c1, c2.
# Состояние = (вектор n, точка кручения T): P = sum n_i G_i + T, G — LLL-базис насыщенной подгруппы.
# Векторы: все n с n^T H n <= B (H — матрица канонических высот, PARI qfminim, плавающая точка), K штук по возрастанию высоты.
# Сито: для хорошего простого p строится таблица на E(F_p): точка -> «может ли f дать >= 7 квадратов»
#   (число квадратов mod p: 3 + [C0 □] + #{k : a_k + b_k f^2 — квадрат или 0 в F_p}); точки, где формула f не определена
#   mod p (O, y = 0, den = 0), разрешены. Если f(P) p-целое и формула определена в редукции, редукция f(P) = f(редукции P),
#   а квадрат в Q остаётся квадратом (или 0) в F_p, поэтому сито не отбрасывает настоящих попаданий.
#   Индексы: дискретный логарифм в E(F_p) = Z/N1 x Z/N2 (таблица перебором), векторы — numpy.
# Выжившие проверяются точно над Q. Исключительные точки (O, y = 0) — отдельно, формулы Codex (как exc_points.sage).
import sys, time, functools
print = functools.partial(print, flush=True)
args = sys.argv[1:]; sys.argv = ['x']
load('/home/kep/magicKube/corners/bremner_type.sage')
import numpy as np
Rx = PolynomialRing(QQ, 'X')

def exceptional_fs(bc):
    E = bc.E; a1, a2, a3, a4, a6 = E.ainvs(); q, c, d = bc.qq, bc.cc, bc.dd
    out = []
    if bc.den(bc.m0) != 0: out.append(('O', bc.num(bc.m0)/bc.den(bc.m0)))
    X = Rx.gen()
    for x0, _ in (X^3 + a2*X^2 + a4*X + a6).roots():
        N = 2*q*(x0 + c) - d^2/(2*q)
        Fx = -(3*x0^2 + 2*a2*x0 + a4); Fy = a1*x0 + a3
        if N != 0 or Fx == 0: out.append((f'({x0},0)', bc.f0)); continue
        m = bc.m0 - 2*q*Fy/Fx
        if bc.den(m) != 0: out.append((f'({x0},0)', bc.num(m)/bc.den(m)))
    return out

def get_gens(E, tlimit=180):
    """ Генераторы: PARI ellrank (+ mwrank с ограничением по времени, если найдено меньше верхней оценки), насыщение.
        Возвращает (список точек на E, lo, hi, источник, индекс насыщения). hi — верхняя оценка ранга (2-спуск PARI). """
    Em = E.minimal_model(); phi = Em.isomorphism_to(E)
    r = pari(Em.ainvs()).ellinit().ellrank()
    lo, hi = ZZ(r[0]), ZZ(r[1])
    pts = [Em(QQ(P[0]), QQ(P[1])) for P in r[3]]; src = 'pari'
    if len(pts) < hi:
        try:
            alarm(tlimit); G2 = Em.gens(proof=False); cancel_alarm()
            pts += list(G2); src += '+mwrank'
        except BaseException:
            cancel_alarm(); src += '+mwrank_timeout'
    ind = []
    for P in pts:
        if P.has_finite_order(): continue
        if Em.regulator_of_points(ind + [P]) > 1e-6: ind.append(P)
    idx = 1
    if ind:
        try:
            alarm(tlimit); ind, idx, _ = Em.saturation(ind); cancel_alarm()
        except BaseException:
            cancel_alarm(); src += '+sat_timeout'; idx = None
        ind, _ = Em.lll_reduce(ind)
    return [phi(P) for P in ind], lo, hi, src, idx

def height_vectors(E, G, K, hcap=2000):
    """ векторы n != 0 с n^T H n <= hcap, не более K (по одному из пары +-n), плюс противоположные и нулевой.
        Потолок по высоте обязателен: при ранге 1 и K = 50000 высота доходит до 10^9, точки с координатами в миллионы
        знаков делают точную проверку неподъёмной. """
    r = len(G)
    if r == 0: return np.zeros((0, 0), dtype=np.int64), 0.0
    H = E.height_pairing_matrix(G)
    Hp = pari(matrix(RR, H))
    detH = H.det(); vol = RR(pi)^(r/2)/RR(gamma(r/2 + 1))
    B = min(RR((2*K*sqrt(detH)/vol)^(2/r)), RR(hcap))
    res = Hp.qfminim(B, 4*K, 2)
    if ZZ(res[0])//2 <= 4*K:                       # весь эллипсоид высоты <= hcap помещается
        V = np.array([[int(res[2][j][i]) for i in range(r)] for j in range(len(res[2]))], dtype=np.int64)
        Hn = np.array([[float(H[a_, c_]) for c_ in range(r)] for a_ in range(r)])
        hts = np.einsum('ij,jk,ik->i', V, Hn, V)
        order = np.argsort(hts)[:K]; V = V[order]
        assert len(V) == 0 or np.abs(V).max() < 2**20
        return np.vstack([np.zeros((1, r), dtype=np.int64), V, -V]), float(hts[order][-1]) if len(order) else 0.0
    # qfminim при переполнении (больше m векторов) выдаёт НЕ самые короткие — подбираем B так, чтобы K <= cnt <= 4K
    lo_B, hi_B = None, None
    for _ in range(60):
        res = Hp.qfminim(B, 4*K, 2)
        cnt = ZZ(res[0])//2
        if cnt > 4*K: hi_B = B
        elif cnt < K: lo_B = B
        else: break
        B = (lo_B + hi_B)/2 if (lo_B is not None and hi_B is not None) else (B*1.6 if hi_B is None else B/1.6)
    else:
        raise RuntimeError("height bound search failed")
    V = np.array([[int(res[2][j][i]) for i in range(r)] for j in range(len(res[2]))], dtype=np.int64)
    Hn = np.array([[float(H[a, c]) for c in range(r)] for a in range(r)])
    hts = np.einsum('ij,jk,ik->i', V, Hn, V)
    order = np.argsort(hts)[:K]
    V = V[order]
    assert np.abs(V).max() < 2**20        # V @ DG в int64: |n_i| < 2^20, дискр. логарифмы < 2^12, ранг <= 8
    # нулевой вектор: чистое кручение (n = 0) — иначе точки кручения порядка > 2 не проверяются (замечание Astra)
    return np.vstack([np.zeros((1, r), dtype=np.int64), V, -V]), float(hts[order][-1]) if len(order) else 0.0

def p_integral(xs, p):
    return all(QQ(v).denominator() % p != 0 for v in xs)

def sieve_table(bc, cs, C0sq, p, G, tors):
    """ Для простого p: (N1, N2, dlog генераторов (r x 2), dlog кручения (t x 2), таблица allowed[N1, N2]) или None. """
    E = bc.E
    data = list(E.ainvs()) + [bc.qq, bc.cc, bc.dd, bc.m0] + list(bc.num.coefficients()) + list(bc.den.coefficients())
    if not p_integral(data, p) or E.discriminant().valuation(p) != 0 or bc.qq.valuation(p) != 0: return None
    F = GF(p); Ep = E.change_ring(F)
    gs = [g.element() for g in Ep.abelian_group().gens()]     # нормальная форма Смита (Ep.gens() — не всегда)
    if len(gs) == 1: gs = [gs[0], Ep(0)]
    N1, N2 = gs[0].order(), gs[1].order()
    if N1*N2 != Ep.order(): return None
    table = {}; allowed = np.zeros((N1, N2), dtype=bool)
    q, c, d, m0 = F(bc.qq), F(bc.cc), F(bc.dd), F(bc.m0)
    num = bc.num.change_ring(F); den = bc.den.change_ring(F)
    ck = [(F(a), F(be)) for (a, be) in cs.values()]
    base = 3 + C0sq
    Pi = Ep(0)
    for i in range(N1):
        Q = Pi
        for j in range(N2):
            key = (Q[0], Q[1], Q[2])
            if key in table: return None          # не прямое произведение
            table[key] = (i, j)
            ok = True
            if not Q.is_zero() and Q[1] != 0:
                x, y = Q.xy()
                m = (2*q*(x + c) - d^2/(2*q))/y + m0
                dm = den(m)
                if dm != 0:
                    f = num(m)/dm
                    cnt = base + sum(1 for (a, be) in ck if (a + be*f^2).is_square())
                    ok = cnt >= 7
            allowed[i, j] = ok
            Q = Q + gs[1]
        Pi = Pi + gs[0]
    def red(P):
        if P.is_zero(): return (0, 0)
        x, y = P.xy()
        if x.denominator() % p == 0: return (0, 0)
        R = Ep(F(x), F(y)); return table[(R[0], R[1], R[2])]
    DG = np.array([red(g) for g in G], dtype=np.int64).reshape(len(G), 2)
    DT = np.array([red(T) for T in tors], dtype=np.int64).reshape(len(tors), 2)
    return N1, N2, DG, DT, allowed

def run_base(b, h, c1, c2, f0, K, nprimes, pmax=3000, validate=False):
    C0, cs = conds(b, h); C0sq = int(C0.is_square())
    t0 = time.time()
    bc = BaseCurve(cs[c1], cs[c2], QQ(f0)); E = bc.E
    G, lo, hi, src, idx = get_gens(E)
    tors = E.torsion_points()
    V, hmax = height_vectors(E, G, K)
    nst = len(V)*len(tors) if len(G) else len(tors)
    t_gens = time.time() - t0
    # сито
    alive = np.ones((max(len(V), 1), len(tors)), dtype=bool)
    used = 0
    for p in primes(5, pmax):
        if used >= nprimes: break
        tb = sieve_table(bc, cs, C0sq, p, G, tors)
        if tb is None: continue
        N1, N2, DG, DT, allowed = tb
        base = V @ DG if len(G) else np.zeros((1, 2), dtype=np.int64)
        for t in range(len(tors)):
            i = (base[:, 0] + DT[t, 0]) % N1; j = (base[:, 1] + DT[t, 1]) % N2
            alive[:, t] &= allowed[i, j]
        used += 1
    surv = [(k, t) for k, t in zip(*np.nonzero(alive))]
    # точная проверка выживших
    hits = []; exact_checked = 0; bad = 0
    def f_of(k, t):
        P0 = sum((int(V[k, a])*G[a] for a in range(len(G))), E(0)) if len(G) else E(0)
        return bc.f_of_point(P0 + tors[t])
    for k, t in surv:
        f = f_of(k, t)
        if f is None: continue
        exact_checked += 1
        if not ((cs[c1][0] + cs[c1][1]*f^2).is_square() and (cs[c2][0] + cs[c2][1]*f^2).is_square()): bad += 1
        nsq, distinct, flags, cells = score(b, h, f, C0, cs)
        if nsq >= 7: hits.append((nsq, distinct, f))
    # исключительные точки
    exc_hits = []
    for where, f in exceptional_fs(bc):
        nsq, distinct, flags, cells = score(b, h, f, C0, cs)
        if nsq >= 7 and distinct: exc_hits.append((where, nsq, f))
    val = ''
    if validate:           # контроль сита: все состояния точно; ни одно отброшенное не должно давать >= 7
        fn = 0; tot = 0
        for k in range(len(V) if len(G) else 1):
            for t in range(len(tors)):
                if alive[k, t]: continue
                f = f_of(k, t)
                if f is None: continue
                tot += 1
                if score(b, h, f, C0, cs)[0] >= 7: fn += 1
        val = f" VALIDATE rejected-exact={tot} false-negatives={fn}"
    rank_s = f"rank[{lo},{hi}] gens={len(G)}" + ("" if len(G) == hi else " INCOMPLETE")
    line = (f"  {c1}/{c2}: {rank_s} src={src} sat={idx} tors={len(tors)} states={nst} hmax={hmax:.1f} primes={used} "
            f"surv={len(surv)} exact={exact_checked} badbase={bad} hits7={[(n_, d_) for n_, d_, _ in hits]} "
            f"exc={len(exc_hits)} t_gens={t_gens:.0f}s t={time.time()-t0:.0f}s{val}")
    return line, [x for x in hits if x[1]], exc_hits

def run_pair(b, h, K, nprimes, validate=False):
    C0, cs = conds(b, h)
    out = [f"({b},{h}) C0={C0} {'□' if C0.is_square() else 'non□'} K={K}"]
    allhits = []
    for c1, c2, f0 in [('TL','BR',h),('TL','D',h),('BR','D',h),('BL','TR',b),('BL','D',b),('TR','D',b)]:
        try:
            line, hits, exc = run_base(b, h, c1, c2, f0, K, nprimes, validate=validate)
            out.append(line); allhits += [(c1, c2) + tuple(x) for x in hits] + [(c1, c2, 'exc') + tuple(x) for x in exc]
        except Exception as ex:
            out.append(f"  {c1}/{c2}: ERROR {type(ex).__name__}: {str(ex)[:80]}")
    uniq = {}
    for hh in allhits: uniq.setdefault((hh[0], hh[1], abs(hh[-1])), hh)
    for hh in uniq.values(): out.append(f"   HIT {hh}")
    print("\n".join(out))

if __name__ == '__main__' or True:
    if args and args[0] == 'validate':
        run_pair(31, 17, 1500, 150, validate=True)
        run_pair(23, 7, 1500, 150, validate=True)
    elif len(args) >= 2:
        K = int(args[2]) if len(args) > 2 else 50000
        run_pair(int(args[0]), int(args[1]), K, 150)
