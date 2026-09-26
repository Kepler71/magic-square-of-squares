# Claude, 26.09.2026 (вторая попытка задачи чанка 0). ТРЕТИЙ путь: закрытие Демьяненко–Манина ранга 1 целиком
# на НЕминимальной модели E1: Y^2 = g(X) = (X+c1)(X+c2)(X+c3), c_i = L/lambda_i, X = L z, Y = L y, L = prod T.
# Отличия от rig_dem.py (его НЕ импортирую):
#  * нет минимальной модели и замены (u,r,s,t): z = X/L, c1 = log|L| (обе стороны);
#  * формула удвоения выведена заново: x(2P) = (g'^2 - 4(a2+2x) g) / (4g) и сверена с групповым законом;
#  * свидетели Безу — через xgcd над Q для (F, G) и для обращённых форм, общий R = НОК знаменателей;
#    лемма проверяется как ОДНОРОДНЫЕ тождества в Z[m,n]: A0 F + B0 G = R n^7, A7 F + B7 G = R m^7;
#  * C = B_up + B_low + 2 log|L| (ровно то, что даёт доказательство), M0 = floor((C/hhat_lo + 1)/2);
#  * кручение E1: замыкание деления на 2 от O + НОД #E1(F_p) (ellcard) для 3 <= p < 500;
#  * образующая на E1: PARI ellrank (зерна 1..5) + насыщение eclib на E1 + PARI ellsaturation до 10^4;
#  * случай B: PARI polresultant по Y и factor над Q (линейные множители).
# Из rig_<r>_<s>.json беру только: список T закрывающих множителей; для СВЕРКИ — P0 на Emin, urst, корни случая B.
# Запуск: env DOT_SAGE=/tmp/claude_rc0 python3 cross_E1.py
from sage.all import *
import json, glob, time, sys, multiprocessing as mp

RIF = RealIntervalField(256)
KD = 7
PSAT = 10**4
BASE = '/home/kep/magicKube/rigorous_census/chunk0'


def Hq(q):
    q = QQ(q)
    return max(abs(q.numerator()), q.denominator())


def lh(q):
    return RIF(Hq(q)).log()


def issq(q):
    q = QQ(q)
    return q >= 0 and q.is_square()


def cells8(r, s):
    return [QQ(v) for v in (s, r, s - r, s + r, -s, -r, r - s, -s - r)]


def e1_model(T):
    T = [ZZ(t) for t in T]
    L = T[0] * T[1] * T[2]
    c = [L // t for t in T]
    assert all(ci * t == L for ci, t in zip(c, T))
    a2 = c[0] + c[1] + c[2]; a4 = c[0] * c[1] + c[0] * c[2] + c[1] * c[2]; a6 = c[0] * c[1] * c[2]
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
    gz = (L * z + c[0]) * (L * z + c[1]) * (L * z + c[2])
    assert gz == L**2 * prod(1 + t * z for t in T)
    return EllipticCurve([0, a2, 0, a4, a6]), L, (a2, a4, a6)


def dup_forms(a2, a4, a6):
    Rx = PolynomialRing(ZZ, 'x'); x = Rx.gen()
    g = x**3 + a2 * x**2 + a4 * x + a6
    F = g.derivative()**2 - 4 * (a2 + 2 * x) * g
    G = 4 * g
    assert F.leading_coefficient() == 1 and F.degree() == 4 and G.degree() == 3
    return F, G


def bezout_xgcd(F, G):
    """свидетели через xgcd; проверка — однородные тождества в Z[m,n]"""
    Q = PolynomialRing(QQ, 't')
    fl = F.list() + [0] * (5 - len(F.list())); gl = G.list() + [0] * (5 - len(G.list()))
    d0, U0, V0 = xgcd(Q(fl), Q(gl))
    d7, U7, V7 = xgcd(Q(fl[::-1]), Q(gl[::-1]))
    assert d0 == 1 and d7 == 1
    for P in (U0, V0, U7, V7):
        assert P.degree() <= 3
    Rr = lcm([c.denominator() for P in (U0, V0, U7, V7) for c in P.list()])
    Mn = PolynomialRing(ZZ, ['m', 'n']); m, n = Mn.gens()

    def hom(cl, D, swap=False):
        cl = list(cl) + [0] * (D + 1 - len(cl))
        return sum(ZZ(cl[i]) * (n**i * m**(D - i) if swap else m**i * n**(D - i)) for i in range(D + 1))
    Fh = hom(fl, 4); Gh = hom(gl, 4)
    A0 = [Rr * c for c in U0.list()]; B0 = [Rr * c for c in V0.list()]
    A7 = [Rr * c for c in U7.list()]; B7 = [Rr * c for c in V7.list()]
    assert hom(A0, 3) * Fh + hom(B0, 3) * Gh == Rr * n**7
    assert hom(A7, 3, True) * Fh + hom(B7, 3, True) * Gh == Rr * m**7
    K = max(sum(abs(ZZ(c)) for c in A0 + B0), sum(abs(ZZ(c)) for c in A7 + B7))
    Lb = max(sum(abs(c) for c in fl), sum(abs(c) for c in gl))
    return dict(R=ZZ(Rr), K=ZZ(K), L=ZZ(Lb), A0=[str(c) for c in A0], B0=[str(c) for c in B0],
                A7=[str(c) for c in A7], B7=[str(c) for c in B7])


def torsion_E1(E):
    pts = [E(0)]; front = [E(0)]
    for _ in range(10):
        new = []
        for Q in front:
            for Rp in Q.division_points(2):
                if Rp not in pts:
                    pts.append(Rp); new.append(Rp)
        front = new
        if not front:
            break
    assert not front
    gp = 0; D = E.discriminant()
    for p in prime_range(3, 500):
        if D % p:
            gp = gcd(gp, E.change_ring(GF(p)).cardinality())
    odd = gp
    while odd % 2 == 0:
        odd //= 2
    return pts, int(gp), int(odd)


def caseB_roots(E, a2, g_pari, t):
    Xt, Yt = t[0], t[1]
    N = pari(f'(y - ({Yt}))^2 - ({a2} + ({Xt}))*(x - ({Xt}))^2')
    C = pari('y^2') - g_pari
    res = pari.polresultant(N, C, pari('y'))
    assert res != 0
    roots = []
    for fac in res.factor()[0]:
        if fac.poldegree() == 1:
            roots.append(QQ(-fac.polcoef(0) / fac.polcoef(1)))
    return roots, int(res.poldegree())


def run(T, S, okfun, ref=None, label=''):
    t0 = time.time()
    out = dict(label=label, T=[int(t) for t in T])
    E, L, (a2, a4, a6) = e1_model(T)
    out['E1'] = [str(v) for v in E.ainvs()]
    # удвоение: формула + сверка с групповым законом
    F, G = dup_forms(a2, a4, a6)
    cert = bezout_xgcd(F, G)
    Bup = RIF(cert['L']).log() / 3; Blow = RIF(cert['K']).log() / 3
    out['cert'] = {k: (str(v) if not isinstance(v, list) else v) for k, v in cert.items()}
    out['B_up'] = float(Bup.upper()); out['B_low'] = float(Blow.upper())
    # ранг и образующая на E1
    pe = pari(E)
    for att in range(5):
        pari.setrand(att + 1)
        rk = pe.ellrank()
        lo, hi = int(rk[0]), int(rk[1])
        if (lo, hi) != (1, 1) or len(rk[3]) > 0:
            break
    out['ellrank_E1'] = [lo, hi]
    if (lo, hi) != (1, 1):
        out['status'] = f'E1: ellrank {lo, hi}'
        return out
    pts = [E(list(map(QQ, p))) for p in rk[3]]
    pts = [p for p in pts if p.order() == oo]
    if not pts:
        out['status'] = 'E1: [1,1] без точки'
        return out
    sat, idx, _ = E.saturation([pts[0]])
    P0 = sat[0]
    assert P0.order() == oo
    V = pari.ellsaturation(pe, [pari([P0[0], P0[1]])], PSAT)
    Qs = E(list(map(QQ, V[0])))
    tors, gp, odd = torsion_E1(E)
    out['tors'] = len(tors); out['tors_gcd'] = gp; out['tors_odd'] = odd
    tset = set(tors)
    sat_ok = any((Qs - e * P0) in tset for e in (1, -1))
    out['sat_idx_eclib'] = int(idx); out['sat_pari_1e4'] = bool(sat_ok)
    out['P0_E1'] = [str(P0[0]), str(P0[1])]
    for Q in (P0, 2 * P0, 3 * P0 + tors[-1]):
        if not (2 * Q).is_zero():
            assert (2 * Q)[0] == QQ(F(Q[0])) / QQ(G(Q[0]))
    # сверка с путём rig_dem (Emin): образ P0_min равен ±P0_E1 по модулю кручения
    if ref is not None:
        u, rr, ss, tt = [QQ(v) for v in ref['urst']]
        Em = EllipticCurve([QQ(v) for v in ref['Emin']])
        Pm = Em([QQ(v) for v in ref['P0']])
        img = E([u**2 * Pm[0] + rr, u**3 * Pm[1] + ss * u**2 * Pm[0] + tt])
        out['same_P0_as_Emin_path'] = bool(any((img - e * P0) in tset for e in (1, -1)))
    # вилка hhat(P0) на E1
    xk = QQ(P0[0])
    for _ in range(KD):
        xk = QQ(F(xk)) / QQ(G(xk))
    hk = lh(xk)
    hlo = ((hk - Blow) / 4**KD).lower(); hhi = ((hk + Bup) / 4**KD).upper()
    out['hP0_lo'] = float(hlo); out['hP0_hi'] = float(hhi)
    out['hP0_sage'] = float(P0.height())
    out['sage_in_bracket'] = bool(hlo <= P0.height() <= hhi)
    if ref is not None:
        out['brackets_overlap'] = bool(max(hlo, ref['hP0_lo']) <= min(hhi, ref['hP0_hi']) + 1e-12)
    c1 = RIF(abs(L)).log()
    C = Bup + Blow + 2 * c1
    M0 = ZZ(floor(((C / RIF(hlo) + 1) / 2).upper()))
    out['c1'] = float(c1.upper()); out['C'] = float(C.upper()); out['M0'] = int(M0)
    try:
        out['B_silverman_E1'] = float(E.silverman_height_bound())
    except Exception as e:
        out['B_silverman_E1'] = repr(e)[:80]
    # случай A на E1
    bad = []; ncand = 0; nonsq = 0; hits = []
    Q = E(0)
    for n_ in range(0, M0 + 1):
        if n_:
            Q = Q + P0
        for t in tors:
            Qt = Q + t
            if Qt.is_zero():
                continue
            ncand += 1
            z = Qt[0] / L
            if not issq(prod(1 + l * z for l in T)):
                nonsq += 1
            if z != 0 and (okfun(z) or okfun(-z)):
                bad.append(str(z)); hits.append([str(z), n_])
    out['caseA_cand'] = ncand; out['caseA_nondeg'] = bad; out['caseA_hits'] = hits; out['ctrl_T_nonsq'] = nonsq
    # случай B на E1
    Rp = PolynomialRing(QQ, 'x')
    g_pari = pari(Rp([a6, a4, a2, 1]))
    rootsz = set(); badB = []
    for t in tors:
        if t.is_zero():
            continue
        rts, deg = caseB_roots(E, a2, g_pari, t)
        for X in rts:
            z = X / L
            rootsz.add(str(z))
            if z != 0 and (okfun(z) or okfun(-z)):
                badB.append(str(z))
    out['caseB_roots_z'] = sorted(rootsz); out['caseB_nondeg'] = badB
    if ref is not None:
        out['caseB_same_as_Emin_path'] = bool(rootsz == set(ref['caseB_roots_z']))
    ok = (not bad and not badB and nonsq == 0 and odd == 1 and sat_ok and out['sage_in_bracket'])
    if ref is not None:
        ok = ok and out['same_P0_as_Emin_path'] and out['brackets_overlap']
    out['closed_E1'] = bool(ok)
    out['sec'] = round(time.time() - t0, 2)
    return out


def work(arg):
    kind, payload = arg
    try:
        if kind == 'slope':
            sl, o = payload
            r, s = map(int, sl.split('/'))
            S = cells8(r, s)
            res = run(o['T'], S, lambda z: all(issq(1 + l * z) for l in S), ref=o, label=sl)
            res['M_run_Emin_path'] = o['M_run']; res['M0_Emin_path'] = o['M0_rig']
            return res
        T, z0 = payload
        Tq = [QQ(t) for t in T]
        okf = lambda z: issq(prod(1 + l * z for l in Tq)) and issq(prod(1 - l * z for l in Tq))
        res = run(T, None, okf, label='posctrl')
        res['z0'] = str(z0)
        tg = {str(QQ(z0)), str(-QQ(z0))}
        res['found_z0'] = bool(tg & set(res.get('caseA_nondeg', []) + res.get('caseB_nondeg', [])))
        return res
    except Exception as e:
        import traceback
        return dict(label=str(kind), T=str(payload)[:100], status='ОШИБКА', err=repr(e), tb=traceback.format_exc()[-1200:])


if __name__ == '__main__':
    tasks = []
    for fn in sorted(glob.glob(f'{BASE}/rig_*_*.json')):
        d = json.load(open(fn))
        tasks += [('slope', (d['slope'], o)) for o in d['results'] if o.get('closed_rigorous')]
    pc = json.load(open(f'{BASE}/posctrl.json'))
    tasks += [('pos', (o['T'], o['z0'])) for o in pc]
    print('задач', len(tasks), flush=True)
    t0 = time.time(); out = []
    with mp.get_context('fork').Pool(3) as pool:
        for i, o in enumerate(pool.imap_unordered(work, tasks)):
            out.append(o)
            msg = f"[{i+1}/{len(tasks)} {time.time()-t0:.0f}s] {o.get('label')} T={o.get('T')} "
            if o.get('status') == 'ОШИБКА':
                msg += 'ОШИБКА ' + o['err']
            elif 'closed_E1' in o:
                msg += (f"closed_E1={o['closed_E1']} M0_E1={o['M0']} (Emin: {o.get('M0_Emin_path')}) "
                        f"B_low={o['B_low']:.2f} B_up={o['B_up']:.2f} c1={o['c1']:.2f} cand={o['caseA_cand']}")
                if o['label'] == 'posctrl':
                    msg += f" z0={o['z0']} найден={o['found_z0']}"
            else:
                msg += str(o.get('status'))
            print(msg, flush=True)
    json.dump(out, open(f'{BASE}/cross_E1.json', 'w'), ensure_ascii=False, indent=1)
    sl = [o for o in out if o.get('label') not in ('posctrl', 'pos')]
    closed = {}
    for o in sl:
        closed.setdefault(o['label'], []).append(bool(o.get('closed_E1')))
    print('наклонов закрыто на E1:', sum(any(v) for v in closed.values()), '/', len(closed), flush=True)
    for k, v in sorted(closed.items()):
        print(f'  {k}: закрывающих на E1 {sum(v)}/{len(v)}', flush=True)
    pos = [o for o in out if o.get('label') == 'posctrl']
    print('положительный контроль: найдено ±z0 в', sum(bool(o.get('found_z0')) for o in pos), 'из', len(pos), flush=True)
    print('ошибок:', sum(o.get('status') == 'ОШИБКА' for o in out), ' время %.0fs' % (time.time() - t0), flush=True)
