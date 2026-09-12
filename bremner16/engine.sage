# Общий движок для всех 16 конфигураций Бремнера (Claude, 2026-09-12).
#
# Параметризация Бремнера (Acta Arith. 99 (2001), (3)):
#     cell_k = c + A_k*a + B_k*b,  k = 0..8 построчно.
# Противоположные пары клеток: (0,8), (1,7), (2,6), (3,5); в каждой сумма = 2c.
#
# СЕМЕЙСТВО: берём одну полную противоположную пару (cell_i, cell_i') = (P^2, Q^2)  =>  c = (P^2+Q^2)/2,
# и третью клетку cell_j = t^2 (j не в паре, j != центр).  Тогда a, b линейны по T = t^2,
# значит ВСЕ девять клеток имеют вид alpha_k + beta_k * T.  Три клетки квадратны автоматически,
# центр — квадрат или нет (зависит только от (P,Q)), остальные пять клеток нетривиальны.
#
# С точностью до D8 есть ровно 4 типа таких семейств (стабилизатор пары в D8 имеет порядок 4):
#   F1 = (пара рёбер {1,7}, третья клетка ребро 5)   — это в точности corners/bremner_type.sage
#   F2 = (пара рёбер {1,7}, третья клетка угол 0)
#   F3 = (пара углов {0,8}, третья клетка угол 2)
#   F4 = (пара углов {0,8}, третья клетка ребро 1)
# Любой магический квадрат с >= 7 квадратными клетками имеет >= 2 полные противоположные пары
# (дырок всего две), поэтому попадает в одно из этих семейств => покрытие полное по постановке.
#
# БАЗОВАЯ КРИВАЯ: пара условий (cell_c1 = []^2, cell_c2 = []^2) — пересечение двух квадрик, род 1.
# Рациональная точка берётся из вырожденных значений t (t = P, t = Q, t = sqrt(C) при квадратном C),
# где квадрат имеет повторяющиеся клетки.  При КВАДРАТНОМ ЦЕНТРЕ при t = Q все пять свободных клеток
# — квадраты, поэтому доступны все 10 пар (в bremner_type.sage использовались только 6).
import itertools, functools, sys, time
print = functools.partial(print, flush=True)
import numpy as np

Rm = PolynomialRing(QQ, 'm'); mvar = Rm.gen()
Rx = PolynomialRing(QQ, 'X')

# (коэф. при c, при a, при b) для клеток 0..8 из (3) Бремнера
COEF = [(1, 1, 0), (1, -1, -1), (1, 0, 1), (1, -1, 1), (1, 0, 0), (1, 1, -1), (1, 0, -1), (1, 1, 1), (1, -1, 0)]
FAMILIES = {'F1': ((1, 7), 5), 'F2': ((1, 7), 0), 'F3': ((0, 8), 2), 'F4': ((0, 8), 1)}
NAMES = ['c0', 'c1', 'c2', 'c3', 'c4', 'c5', 'c6', 'c7', 'c8']


def family_cells(ft, P, Q):
    """ -> (cells, pair, j, free): cells[k] = (alpha, beta), клетка = alpha + beta*t^2. """
    (i, ip), j = FAMILIES[ft]
    M = matrix(QQ, [[COEF[k][1], COEF[k][2], COEF[k][0]] for k in (i, ip, j)])   # неизвестные a, b, c
    Mi = M.inverse()
    abc0 = Mi * vector(QQ, [QQ(P)^2, QQ(Q)^2, 0])
    abc1 = Mi * vector(QQ, [0, 0, 1])
    cells = []
    for k in range(9):
        cc, ca, cb = COEF[k]
        cells.append((cc*abc0[2] + ca*abc0[0] + cb*abc0[1], cc*abc1[2] + ca*abc1[0] + cb*abc1[1]))
    free = [k for k in range(9) if cells[k][1] != 0 and k != j]
    assert len(free) == 5 and cells[4][1] == 0
    return cells, (i, ip), j, free


def cellvals(cells, t):
    T = QQ(t)^2
    return [al + be*T for al, be in cells]


def score9(cells, t):
    """ (число квадратных клеток, все девять различны?, все положительны?, значения, флаги). """
    v = cellvals(cells, t)
    fl = [x > 0 and x.is_square() for x in v]
    return sum(fl), len(set(v)) == 9, all(x > 0 for x in v), v, fl


def degenerate_ts(cells, pair=None, j=None, C=None):
    """ Кандидаты на рациональную точку базовой коники: t >= 0, при которых магический квадрат тривиален.
        Квадрат (3) тривиален ровно при ab(a^2-b^2)(a^2-4b^2)(4a^2-b^2) = 0 (Бремнер, с. 290);
        a и b линейны по T = t^2, поэтому каждое условие даёт одно значение T. """
    a0, a1 = cells[0][0] - cells[4][0], cells[0][1] - cells[4][1]      # a = a0 + a1*T
    b0, b1 = cells[2][0] - cells[4][0], cells[2][1] - cells[4][1]      # b = b0 + b1*T
    ts = set()
    for (u, v) in [(1, 0), (0, 1), (1, -1), (1, 1), (1, -2), (1, 2), (2, -1), (2, 1)]:
        A, B = u*a1 + v*b1, u*a0 + v*b0                                # A*T + B = 0
        if A == 0:
            continue
        T = -B/A
        if T >= 0 and T.is_square():
            ts.add(T.sqrt())
    return sorted(ts)


class BaseCurve:
    """ Y1^2 = a1 + b1 t^2 (коника с точкой t0), Y2^2 = a2 + b2 t^2  ->  эллиптическая кривая (Коннелл). """
    def __init__(self, c1, c2, t0, ratbound=400):
        (a1, b1), (a2, b2) = c1, c2
        val = a1 + b1*QQ(t0)^2
        if val <= 0 or not val.is_square():
            raise ValueError("t0 не даёт квадрат в первом условии")
        Y10 = val.sqrt()
        self.a1, self.b1, self.a2, self.b2, self.f0, self.Y10 = a1, b1, a2, b2, QQ(t0), Y10
        den = mvar^2 - b1
        self.num = self.f0*den + 2*b1*self.f0 - 2*Y10*mvar
        self.den = den
        self.Q = Rm(a2*den^2 + b2*self.num^2)
        if self.Q.degree() < 3:
            raise ValueError("вырожденная квартика")
        pts = pari.hyperellratpoints(pari(self.Q), ratbound)
        pts = [(QQ(Pt[0]), QQ(Pt[1])) for Pt in pts if QQ(Pt[1]) != 0 and self.den(QQ(Pt[0])) != 0]
        if not pts:
            raise ValueError("нет точки на квартике")
        self.m0 = pts[0][0]
        g = Rm(self.Q(mvar + self.m0))
        qq = g[0].sqrt(); assert qq^2 == g[0]
        self.qq, self.aa, self.bb, self.cc, self.dd = qq, g[4], g[3], g[2], g[1]
        A1 = self.dd/qq; A2 = self.cc - self.dd^2/(4*qq^2); A3 = 2*qq*self.bb; A4 = -4*qq^2*self.aa
        self.E = EllipticCurve([A1, A2, A3, A4, A2*A4])

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
        """ Предельные значения t в точках, где формула не определена (O и точки порядка 2). """
        E = self.E; a1, a2, a3, a4, a6 = E.ainvs(); q, c, d = self.qq, self.cc, self.dd
        out = []
        if self.den(self.m0) != 0:
            out.append(self.num(self.m0)/self.den(self.m0))
        X = Rx.gen()
        for x0, _ in (X^3 + a2*X^2 + a4*X + a6).roots():
            N = 2*q*(x0 + c) - d^2/(2*q)
            Fx = -(3*x0^2 + 2*a2*x0 + a4); Fy = a1*x0 + a3
            if N != 0 or Fx == 0:
                out.append(self.f0); continue
            m = self.m0 - 2*q*Fy/Fx
            if self.den(m) != 0:
                out.append(self.num(m)/self.den(m))
        return out


def get_gens(E, tlimit=120):
    """ Генераторы: PARI ellrank + (при недоборе) mwrank, затем насыщение и LLL. """
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
    """ Все n != 0 с n^T H n <= hcap (не более K), плюс -n и нулевой вектор. Потолок обязателен. """
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
    data += [x for k in free for x in cells[k]]
    if not p_integral(data, p) or E.discriminant().valuation(p) != 0 or bc.qq.valuation(p) != 0:
        return None
    F = GF(p); Ep = E.change_ring(F)
    gs = [g.element() for g in Ep.abelian_group().gens()]
    if len(gs) == 1:
        gs = [gs[0], Ep(0)]
    N1, N2 = gs[0].order(), gs[1].order()
    if N1*N2 != Ep.order():
        return None
    table = {}; allowed = np.zeros((N1, N2), dtype=bool)
    q, c, d, m0 = F(bc.qq), F(bc.cc), F(bc.dd), F(bc.m0)
    num = bc.num.change_ring(F); den = bc.den.change_ring(F)
    ck = [(F(cells[k][0]), F(cells[k][1])) for k in free]
    Pi = Ep(0)
    for i in range(N1):
        QQp = Pi
        for jj in range(N2):
            key = (QQp[0], QQp[1], QQp[2])
            if key in table:
                return None
            table[key] = (i, jj)
            ok = True
            if not QQp.is_zero() and QQp[1] != 0:
                x, y = QQp.xy()
                m = (2*q*(x + c) - d^2/(2*q))/y + m0
                dm = den(m)
                if dm != 0:
                    f = num(m)/dm
                    ok = base_cnt + sum(1 for (al, be) in ck if (al + be*f^2).is_square()) >= 7
            allowed[i, jj] = ok
            QQp = QQp + gs[1]
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


def run_base(cells, free, base_cnt, c1, c2, t0, K, nprimes, hcap, pmax=3000, validate=False):
    t_start = time.time()
    bc = BaseCurve(cells[c1], cells[c2], t0)
    E = bc.E
    G, lo, hi, src, idx = get_gens(E)
    tors = E.torsion_points()
    V, hmax = height_vectors(E, G, K, hcap)
    nst = len(V)*len(tors)
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
        if not ((cells[c1][0] + cells[c1][1]*t^2).is_square() and (cells[c2][0] + cells[c2][1]*t^2).is_square()):
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
    line = (f"{NAMES[c1]}/{NAMES[c2]}@t0={t0}: rank[{lo},{hi}] g={len(G)}{'' if len(G)==hi else ' INC'} "
            f"src={src} sat={idx} tors={len(tors)} st={nst} hmax={hmax:.0f} pr={used} surv={len(surv)} "
            f"ex={exact} bad={bad} hits={len(hits)} t={time.time()-t_start:.0f}s{val}")
    return line, hits
