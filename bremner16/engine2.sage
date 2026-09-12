# Движок-2: условия «клетка = квадрат» как произвольные квадратичные многочлены от одного параметра.
# Покрывает два типа однопараметрических семейств:
#   F-семейства (как engine.sage): одна полная противоположная пара (P^2,Q^2) + третья клетка t^2
#                                  -> 3 автоматических квадрата (+центр, если C — квадрат);
#   G-семейства (новое): ДВЕ полные противоположные пары -> 4 автоматических квадрата ВСЕГДА.
# G: cell_i = P^2, cell_i' = Q^2, cell_k = R^2, cell_k' = S^2 требует P^2+Q^2 = R^2+S^2 = 2c;
#    полная рациональная параметризация этой квадрики (Сегре): P = mp+nq, Q = mq-np, R = mp-nq, S = mq+np.
#    При q = 1 параметр p пробегает Q, и все девять клеток — квадратичные многочлены от p.
# С точностью до D8 G-семейств ровно три: {рёберная,рёберная}, {угловая,угловая}, {рёберная,угловая}.
import itertools, functools, sys, time
print = functools.partial(print, flush=True)
import numpy as np

Rp = PolynomialRing(QQ, 'p'); pv = Rp.gen()
Rm = PolynomialRing(QQ, 'm'); mvar = Rm.gen()
Rx = PolynomialRing(QQ, 'X')

COEF = [(1, 1, 0), (1, -1, -1), (1, 0, 1), (1, -1, 1), (1, 0, 0), (1, 1, -1), (1, 0, -1), (1, 1, 1), (1, -1, 0)]
FFAM = {'F1': ((1, 7), 5), 'F2': ((1, 7), 0), 'F3': ((0, 8), 2), 'F4': ((0, 8), 1)}
GFAM = {'G1': ((1, 7), (3, 5)), 'G2': ((0, 8), (2, 6)), 'G3': ((1, 7), (0, 8))}
NAMES = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8']


def cells_F(ft, P, Q):
    (i, ip), j = FFAM[ft]
    M = matrix(QQ, [[COEF[k][1], COEF[k][2], COEF[k][0]] for k in (i, ip, j)])
    Mi = M.inverse()
    r0 = Mi*vector(QQ, [QQ(P)^2, QQ(Q)^2, 0]); r1 = Mi*vector(QQ, [0, 0, 1])
    a = r0[0] + r1[0]*pv^2; b = r0[1] + r1[1]*pv^2; c = r0[2] + r1[2]*pv^2
    cells = [Rp(cc*c + ca*a + cb*b) for cc, ca, cb in COEF]
    auto = [i, ip, j]
    return cells, auto, [k for k in range(9) if k not in auto and k != 4]


def cells_G(gt, m, n):
    (i, ip), (k, kp) = GFAM[gt]
    m, n = QQ(m), QQ(n)
    P = m*pv + n; Q = m - n*pv; R = m*pv - n; S = m + n*pv
    c = (m^2 + n^2)*(pv^2 + 1)/2
    M = matrix(QQ, [[COEF[i][1], COEF[i][2]], [COEF[k][1], COEF[k][2]]])
    Mi = M.inverse()
    rhs = vector(Rp, [Rp(P^2 - c), Rp(R^2 - c)])
    ab = Mi*rhs
    a, b = Rp(ab[0]), Rp(ab[1])
    cells = [Rp(cc*c + ca*a + cb*b) for cc, ca, cb in COEF]
    assert cells[i] == Rp(P^2) and cells[ip] == Rp(Q^2) and cells[k] == Rp(R^2) and cells[kp] == Rp(S^2)
    auto = [i, ip, k, kp]
    return cells, auto, [x for x in range(9) if x not in auto and x != 4]


def score9(cells, t):
    v = [f(QQ(t)) for f in cells]
    fl = [x > 0 and x.is_square() for x in v]
    return sum(fl), len(set(v)) == 9, all(x > 0 for x in v), v, fl


def degenerate_ts(cells):
    """ Значения параметра, при которых магический квадрат тривиален: ab(a^2-b^2)(a^2-4b^2)(4a^2-b^2)=0. """
    a = cells[0] - cells[4]; b = cells[2] - cells[4]
    out = set()
    for (u, v) in [(1, 0), (0, 1), (1, -1), (1, 1), (1, -2), (1, 2), (2, -1), (2, 1)]:
        f = u*a + v*b
        if f == 0:
            continue
        for r, _ in f.roots():
            out.add(QQ(r))
    return sorted(out)


class BaseCurve:
    """ Y1^2 = F1(p), Y2^2 = F2(p) (deg <= 2), рациональная точка p0 на первой конике -> эллиптическая кривая. """
    def __init__(self, F1, F2, p0, ratbound=None):
        if ratbound is None:
            ratbound = globals().get('RATBOUND', 600)
        p0 = QQ(p0)
        val = F1(p0)
        if val <= 0 or not val.is_square():
            raise ValueError("no point on conic")
        Y0 = val.sqrt()
        A = F1[2]; B = F1[1]
        if A == 0:
            raise ValueError("degenerate conic")
        self.num = p0*mvar^2 - 2*Y0*mvar + (A*p0 + B)
        self.den = mvar^2 - A
        self.F1, self.F2, self.p0, self.Y0 = F1, F2, p0, Y0
        self.Q = Rm(F2[2]*self.num^2 + F2[1]*self.num*self.den + F2[0]*self.den^2)
        if self.Q.degree() < 3 or self.Q.is_square():
            raise ValueError("degenerate quartic")
        pts = pari.hyperellratpoints(pari(self.Q), ratbound)
        pts = [(QQ(Pt[0]), QQ(Pt[1])) for Pt in pts if QQ(Pt[1]) != 0 and self.den(QQ(Pt[0])) != 0]
        if not pts:
            raise ValueError("no point on quartic")
        self.m0 = pts[0][0]
        g = Rm(self.Q(mvar + self.m0))
        qq = g[0].sqrt(); assert qq^2 == g[0]
        self.qq, self.cc, self.dd = qq, g[2], g[1]
        A1 = g[1]/qq; A2 = g[2] - g[1]^2/(4*qq^2); A3 = 2*qq*g[3]; A4 = -4*qq^2*g[4]
        try:
            self.E = EllipticCurve([A1, A2, A3, A4, A2*A4])
        except ArithmeticError:
            raise ValueError("singular curve")

    def t_of_point(self, Pt):
        if Pt.is_zero():
            return None
        x, y = Pt.xy()
        if y == 0:
            return None
        m = (2*self.qq*(x + self.cc) - self.dd^2/(2*self.qq))/y + self.m0
        if self.den(m) == 0:
            return None
        return self.num(m)/self.den(m)

    def exceptional_ts(self):
        E = self.E; a1, a2, a3, a4, a6 = E.ainvs(); q, c, d = self.qq, self.cc, self.dd
        out = []
        if self.den(self.m0) != 0:
            out.append(self.num(self.m0)/self.den(self.m0))
        X = Rx.gen()
        for x0, _ in (X^3 + a2*X^2 + a4*X + a6).roots():
            N = 2*q*(x0 + c) - d^2/(2*q)
            Fx = -(3*x0^2 + 2*a2*x0 + a4); Fy = a1*x0 + a3
            if N != 0 or Fx == 0:
                out.append(self.p0); continue
            m = self.m0 - 2*q*Fy/Fx
            if self.den(m) != 0:
                out.append(self.num(m)/self.den(m))
        return out


def get_gens(E, tlimit=90):
    Em = E.minimal_model(); phi = Em.isomorphism_to(E)
    r = pari(Em.ainvs()).ellinit().ellrank()
    lo, hi = ZZ(r[0]), ZZ(r[1])
    pts = [Em(QQ(Pt[0]), QQ(Pt[1])) for Pt in r[3]]; src = 'pari'
    if len(pts) < hi:
        try:
            alarm(tlimit); G2 = Em.gens(proof=False); cancel_alarm(); pts += list(G2); src += '+mw'
        except BaseException:
            cancel_alarm(); src += '+mwTO'
    ind = []
    for Pt in pts:
        if Pt.has_finite_order():
            continue
        if Em.regulator_of_points(ind + [Pt]) > 1e-6:
            ind.append(Pt)
    idx = 1
    if ind:
        try:
            alarm(tlimit); ind, idx, _ = Em.saturation(ind); cancel_alarm()
        except BaseException:
            cancel_alarm(); src += '+satTO'; idx = None
        ind, _ = Em.lll_reduce(ind)
    return [phi(Pt) for Pt in ind], lo, hi, src, idx


def height_vectors(E, G, K, hcap):
    r = len(G)
    if r == 0:
        return np.zeros((1, 0), dtype=np.int64), 0.0
    H = E.height_pairing_matrix(G)
    Hp = pari(matrix(RR, H))
    detH = H.det(); vol = RR(pi)^(r/2)/RR(gamma(r/2 + 1))
    B = min(RR((2*K*sqrt(abs(detH))/vol)^(2/r)), RR(hcap))
    res = Hp.qfminim(B, 4*K, 2)
    if ZZ(res[0])//2 > 4*K:
        lo_B, hi_B = None, B
        for _ in range(60):
            B = (lo_B + hi_B)/2 if lo_B is not None else B/1.6
            res = Hp.qfminim(B, 4*K, 2); cnt = ZZ(res[0])//2
            if cnt > 4*K: hi_B = B
            elif cnt < K: lo_B = B
            else: break
    V = np.array([[int(res[2][j][i]) for i in range(r)] for j in range(len(res[2]))], dtype=np.int64)
    if len(V) == 0:
        return np.zeros((1, r), dtype=np.int64), 0.0
    Hn = np.array([[float(H[a_, c_]) for c_ in range(r)] for a_ in range(r)])
    hts = np.einsum('ij,jk,ik->i', V, Hn, V)
    order = np.argsort(hts)[:K]; V = V[order]
    assert np.abs(V).max() < 2**20
    return np.vstack([np.zeros((1, r), dtype=np.int64), V, -V]), float(hts[order][-1])


def p_integral(xs, p):
    return all(QQ(v).denominator() % p != 0 for v in xs)


def sieve_table(bc, cells, free, base_cnt, p, G, tors):
    E = bc.E
    data = list(E.ainvs()) + [bc.qq, bc.cc, bc.dd, bc.m0] + list(bc.num.coefficients()) + list(bc.den.coefficients())
    data += [x for k in free for x in cells[k].coefficients()]
    if not p_integral(data, p) or E.discriminant().valuation(p) != 0 or bc.qq.valuation(p) != 0:
        return None
    F = GF(p); Ep = E.change_ring(F)
    gs = [g.element() for g in Ep.abelian_group().gens()]
    if len(gs) == 1:
        gs = [gs[0], Ep(0)]
    N1, N2 = gs[0].order(), gs[1].order()
    if N1*N2 != Ep.order() or N1*N2 > 40000:
        return None
    table = {}; allowed = np.zeros((N1, N2), dtype=bool)
    q, c, d, m0 = F(bc.qq), F(bc.cc), F(bc.dd), F(bc.m0)
    num = bc.num.change_ring(F); den = bc.den.change_ring(F)
    ck = [cells[k].change_ring(F) for k in free]
    Pi = Ep(0)
    for i in range(N1):
        Qp = Pi
        for jj in range(N2):
            key = (Qp[0], Qp[1], Qp[2])
            if key in table:
                return None
            table[key] = (i, jj)
            ok = True
            if not Qp.is_zero() and Qp[1] != 0:
                x, y = Qp.xy()
                mm = (2*q*(x + c) - d^2/(2*q))/y + m0
                dm = den(mm)
                if dm != 0:
                    tt = num(mm)/dm
                    ok = base_cnt + sum(1 for f in ck if f(tt).is_square()) >= 7
            allowed[i, jj] = ok
            Qp = Qp + gs[1]
        Pi = Pi + gs[0]

    def red(Pt):
        if Pt.is_zero():
            return (0, 0)
        x, y = Pt.xy()
        if x.denominator() % p == 0 or y.denominator() % p == 0:
            return (0, 0)
        R = Ep(F(x), F(y)); return table[(R[0], R[1], R[2])]
    DG = np.array([red(g) for g in G], dtype=np.int64).reshape(len(G), 2)
    DT = np.array([red(T) for T in tors], dtype=np.int64).reshape(len(tors), 2)
    return N1, N2, DG, DT, allowed


RATBOUND = 600


def run_base(cells, free, base_cnt, c1, c2, t0, K, nprimes, hcap, pmax=4000, validate=False):
    t_start = time.time()
    bc = BaseCurve(cells[c1], cells[c2], t0)
    E = bc.E
    G, lo, hi, src, idx = get_gens(E)
    tors = E.torsion_points()
    V, hmax = height_vectors(E, G, K, hcap)
    alive = np.ones((len(V), len(tors)), dtype=bool)
    used = 0
    for p in primes(5, pmax):
        if used >= nprimes:
            break
        tb = sieve_table(bc, cells, free, base_cnt, p, G, tors)
        if tb is None:
            continue
        N1, N2, DG, DT, allowed = tb
        basev = V @ DG if len(G) else np.zeros((len(V), 2), dtype=np.int64)
        for tt in range(len(tors)):
            i = (basev[:, 0] + DT[tt, 0]) % N1; jj = (basev[:, 1] + DT[tt, 1]) % N2
            alive[:, tt] &= allowed[i, jj]
        used += 1
    surv = list(zip(*np.nonzero(alive)))

    def t_of(k, tt):
        P0 = sum((int(V[k, a])*G[a] for a in range(len(G))), E(0)) if len(G) else E(0)
        return bc.t_of_point(P0 + tors[tt])
    hits = []; exact = 0; bad = 0
    for k, tt in surv:
        t = t_of(k, tt)
        if t is None:
            continue
        exact += 1
        if not (cells[c1](t).is_square() and cells[c2](t).is_square()):
            bad += 1
        n, dist, pos, v, fl = score9(cells, t)
        if n >= 7 and dist:
            hits.append((n, pos, t, v))
    for t in bc.exceptional_ts():
        n, dist, pos, v, fl = score9(cells, t)
        if n >= 7 and dist:
            hits.append((n, pos, t, v))
    val = ''
    if validate:
        fn = tot = 0
        for k in range(len(V)):
            for tt in range(len(tors)):
                if alive[k, tt]:
                    continue
                t = t_of(k, tt)
                if t is None:
                    continue
                tot += 1
                if score9(cells, t)[0] >= 7:
                    fn += 1
        val = f" VALIDATE rejected={tot} false-neg={fn}"
    line = (f"{NAMES[c1]}/{NAMES[c2]}@{t0}: rank[{lo},{hi}] g={len(G)}{'' if len(G)==hi else ' INC'} "
            f"src={src} sat={idx} tors={len(tors)} st={len(V)*len(tors)} hmax={hmax:.0f} pr={used} "
            f"surv={len(surv)} ex={exact} bad={bad} hits={len(hits)} t={time.time()-t_start:.0f}s{val}")
    return line, hits
