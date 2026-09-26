# Claude, 26.09.2026. Второй вычислительный путь для rig_*.json (тот же автор — не независимая проверка).
# По каждому закрывающему множителю из JSON:
#  1) E1 из T по элементарным симметрическим функциям; тождество замены (u,r,s,t): Emin -> E1 в Q[x,y];
#  2) b-инварианты по формулам; F, G удвоения; символьная проверка x(2P) = F/G по модулю уравнения кривой;
#  3) тождества Безу над Z, K, L, B = log max(K,L)/3 (MPFI);
#  4) P0 на Emin; вилка hhat(P0) из x(2^k P0) (только итерации x -> F/G);
#  5) c1 из матрицы z = (u^2 x + r)/L, C, M0;
#  6) перебор nP0 + t на модели E1 (кручение E1 считается заново), z = X/L, все 8 клеток;
#  7) случай B через результант Res_Y((Y-Yt)^2 - (a2+Xt)(X-Xt)^2, Y^2 - f(X)).
from sage.all import *
import sys, json, glob, time

RIF = RealIntervalField(256)
K2 = 7


def cells(r, s):
    v = [s, r, s - r, s + r]
    return [QQ(c) for c in v] + [QQ(-c) for c in v]


def sq(q):
    return q >= 0 and q.is_square()


def h(q):
    q = QQ(q); return RIF(max(abs(q.numerator()), q.denominator())).log()


def verify(r, s, o):
    S = cells(r, s)
    T = [QQ(t) for t in o['T']]
    L = prod(T)
    assert QQ(o['L']) == L
    c = [L / t for t in T]
    e1 = sum(c); e2 = c[0] * c[1] + c[0] * c[2] + c[1] * c[2]; e3 = c[0] * c[1] * c[2]
    E1 = EllipticCurve([0, e1, 0, e2, e3])
    assert [QQ(a) for a in o['E1']] == [QQ(a) for a in E1.ainvs()]
    Em = EllipticCurve([QQ(a) for a in o['Emin']])
    u, rr, ss, tt = [QQ(v) for v in o['urst']]
    R2 = PolynomialRing(QQ, ['x', 'y']); x, y = R2.gens()

    def W(E, X, Y):
        a1, a2, a3, a4, a6 = E.ainvs()
        return Y**2 + a1 * X * Y + a3 * Y - X**3 - a2 * X**2 - a4 * X - a6
    assert W(E1, u**2 * x + rr, u**3 * y + ss * u**2 * x + tt) == u**6 * W(Em, x, y)
    # 2) удвоение
    a1, a2, a3, a4, a6 = Em.ainvs()
    b2 = a1**2 + 4 * a2; b4 = 2 * a4 + a1 * a3; b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    Fc = [-b8, -2 * b6, -b4, 0, 1]; Gc = [b6, 2 * b4, b2, 4, 0]
    cert = o['cert']
    assert [QQ(v) for v in cert['F']] == Fc and [QQ(v) for v in cert['G']] == Gc
    Rx = PolynomialRing(QQ, 'X'); Ry = PolynomialRing(Rx, 'Y'); X = Rx.gen(); Y = Ry.gen()
    Wy = Y**2 + a1 * X * Y + a3 * Y - (X**3 + a2 * X**2 + a4 * X + a6)
    D = 2 * Y + a1 * X + a3
    N = 3 * X**2 + 2 * a2 * X + a4 - a1 * Y
    Fx = Rx(Fc); Gx = Rx(Gc)
    assert (D**2 - Gx) % Wy == 0                                     # G = (2y + a1 x + a3)^2 на кривой
    assert (N**2 + a1 * N * D - (a2 + 2 * X) * D**2 - Fx) % Wy == 0   # x(2P) G = F на кривой
    # 3) Безу
    Pz = PolynomialRing(ZZ, 't'); t = Pz.gen()
    F = Pz([ZZ(v) for v in Fc]); G = Pz([ZZ(v) for v in Gc]); Rr = ZZ(cert['R'])
    assert Rr != 0 and F.gcd(G) == 1
    norms = []
    for w in cert['wit']:
        A = Pz([ZZ(v) for v in w['A']]); B = Pz([ZZ(v) for v in w['B']])
        assert A.degree() <= 3 and B.degree() <= 3
        assert A * F + B * G == Rr * t**w['j']
        norms.append(sum(abs(v) for v in A.list() + B.list()))
    assert sorted(w['j'] for w in cert['wit']) == [0, 7]
    Kb = max(norms); Lb = max(sum(abs(v) for v in F.list()), sum(abs(v) for v in G.list()))
    assert Kb == ZZ(cert['K']) and Lb == ZZ(cert['L'])
    Bup = RIF(Lb).log() / 3; Blow = RIF(Kb).log() / 3
    B = max(Bup.upper(), Blow.upper())
    # 4) P0 и hhat
    P0 = Em([QQ(v) for v in o['P0']])
    xk = P0[0]
    for _ in range(K2):
        xk = Fx(xk) / Gx(xk)
    hk = h(xk)
    h0lo = ((hk - Blow) / 4**K2).lower()
    # 5) c1, M0
    al = u**2 / L; be = rr / L
    den = lcm(al.denominator(), be.denominator())
    a, b, d = ZZ(al * den), ZZ(be * den), ZZ(den)
    g = gcd([a, b, d]); a, b, d = a // g, b // g, d // g
    c1 = max(RIF(max(abs(a) + abs(b), abs(d))).log().upper(), RIF(max(abs(d) + abs(b), abs(a))).log().upper())
    C = 2 * RIF(B) + 2 * RIF(c1)
    M0 = ZZ(floor(((C / RIF(h0lo) + 1) / 2).upper()))
    assert M0 <= o['M_run'], (M0, o['M_run'])
    # 6) перебор на E1
    phi = lambda Pm: E1([u**2 * Pm[0] + rr, u**3 * Pm[1] + ss * u**2 * Pm[0] + tt])
    Q0 = phi(P0)
    tors = E1.torsion_points()
    assert len(tors) == o['tors']
    bad = []; n_c = 0
    Q = E1(0)
    for n in range(0, M0 + 1):
        if n:
            Q = Q + Q0
        for tp in tors:
            Qt = Q + tp
            if Qt.is_zero():
                continue
            n_c += 1
            z = Qt[0] / L
            if z != 0 and (all(sq(1 + l * z) for l in S) or all(sq(1 - l * z) for l in S)):
                bad.append(str(z))
    # 7) случай B через результант
    R2b = PolynomialRing(QQ, ['XX', 'YY']); XX, YY = R2b.gens()
    f = XX**3 + E1.a2() * XX**2 + E1.a4() * XX + E1.a6()
    rootsz = set(); badB = []
    for tp in tors:
        if tp.is_zero():
            continue
        Xt, Yt = tp[0], tp[1]
        Npol = (YY - Yt)**2 - (E1.a2() + Xt) * (XX - Xt)**2
        res = Npol.resultant(YY**2 - f, YY)
        resx = res.univariate_polynomial() if res.degree() > 0 else None
        assert resx is not None and resx != 0
        for rt, _ in resx.roots(QQ):
            z = rt / L
            rootsz.add(str(z))
            if z != 0 and (all(sq(1 + l * z) for l in S) or all(sq(1 - l * z) for l in S)):
                badB.append(str(z))
    sameB = rootsz == set(o['caseB_roots_z'])
    return dict(T=o['T'], M0=int(M0), M0_main=o['M0_rig'], B=float(B), B_main=o['B_rig'], c1=float(c1),
                h0lo=float(h0lo), cand=n_c, nondeg=bad + badB, caseB_same_roots=sameB,
                ok=(not bad and not badB))


if __name__ == '__main__':
    base = '/home/kep/magicKube/rigorous_census/chunk0'
    t0 = time.time()
    summary = {}
    for fn in sorted(glob.glob(f'{base}/rig_*_*.json')):
        d = json.load(open(fn))
        r, s = map(int, d['slope'].split('/'))
        vs = []
        for o in d['results']:
            if not o.get('closed_rigorous'):
                continue
            try:
                v = verify(r, s, o)
            except Exception as e:
                import traceback
                v = dict(T=o['T'], ok=False, err=repr(e), tb=traceback.format_exc()[-800:])
            vs.append(v)
        nok = sum(v['ok'] for v in vs)
        summary[d['slope']] = dict(verified_ok=nok, total=len(vs), details=vs)
        print(f"[{time.time()-t0:.0f}s] {d['slope']}: подтверждено {nok}/{len(vs)}; "
              f"M0 совпал: {sum(v.get('M0') == v.get('M0_main') for v in vs)}; "
              f"корни B совпали: {sum(bool(v.get('caseB_same_roots')) for v in vs)}", flush=True)
        for v in vs:
            if not v['ok']:
                print('   СБОЙ', v, flush=True)
    json.dump(summary, open(f'{base}/verify_chunk0.json', 'w'), ensure_ascii=False, indent=1)
    print('готово', flush=True)
