# Claude, 26.09.2026. Строгая перепроверка закрытий наклонов методом Демьяненко–Манина ранга 1
# (симметрия z -> -z; FABLE_JOINT_SIEVE §7.1) на 3-клеточных множителях E_T.
# Отличие от Fable: константа B (|hhat - h(x_min)|) выводится из сертификата Безу для
# отображения удвоения x -> x(2P) и телескопической суммы, а не из silverman_height_bound;
# hhat(P0) берётся не из Sage, а строгой вилкой hhat(P0) = hhat(2^k P0)/4^k (точная арифметика + MPFI).
# Код Fable (joint_sieve/demjanenko.py) не импортируется и не читался.
# Запуск: env DOT_SAGE=/tmp/claude_rc0 python3 rig_dem.py [slopes...]
from sage.all import *
from sage.libs.eclib.interface import mwrank_EllipticCurve
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt
import sys, json, time, itertools, multiprocessing as mp

RIF = RealIntervalField(256)
KDOUB = 7                       # hhat(P0) через 2^7 P0
NPRIME_TORS = 400               # простые для верхней оценки кручения


def cells(r, s):
    v = [s, r, s - r, s + r]
    return [QQ(c) for c in v] + [QQ(-c) for c in v]


def is_sq(q):
    return q >= 0 and q.is_square()


def full_ok(z, S):
    return all(is_sq(1 + l * z) for l in S)


def hgt(q):
    """логарифмическая высота рационального числа как точки P^1 (строгая вилка)"""
    q = QQ(q)
    return RIF(max(abs(q.numerator()), q.denominator())).log()


def weq(E, X, Y):
    a1, a2, a3, a4, a6 = E.ainvs()
    return Y**2 + a1 * X * Y + a3 * Y - (X**3 + a2 * X**2 + a4 * X + a6)


# ---------- модель E_T ----------
def model3(T):
    """y^2 = prod_{l in T}(1+l z)  <->  Y^2 = g(X) = prod(X + L/l), X = L z, Y = L y, L = prod T."""
    T = [QQ(t) for t in T]
    L = prod(T)
    R = PolynomialRing(QQ, 'X'); X = R.gen()
    g = R.prod(X + L / t for t in T)
    Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
    assert g(L * z) == L**2 * prod(1 + t * z for t in T)          # контроль тождества модели
    co = g.list()
    assert co[3] == 1 and all(c in ZZ for c in co)
    E1 = EllipticCurve([0, co[2], 0, co[1], co[0]])
    return E1, L, g


def urst_to(Esrc, Edst):
    """(u,r,s,t) с X = u^2 x + r, Y = u^3 y + s u^2 x + t : Esrc(x,y) -> Edst(X,Y); проверено символьно."""
    Rxy = PolynomialRing(QQ, ['x', 'y']); x, y = Rxy.gens()
    for iso in (Esrc.isomorphism_to(Edst), Edst.isomorphism_to(Esrc)):
        for (u, r, s, t) in [iso.tuple()]:
            if weq(Edst, u**2 * x + r, u**3 * y + s * u**2 * x + t) == u**6 * weq(Esrc, x, y):
                return u, r, s, t
    # попытка обратного преобразования
    for iso in (Esrc.isomorphism_to(Edst), Edst.isomorphism_to(Esrc)):
        u, r, s, t = iso.tuple()
        ui = 1 / u; ri = -r / u**2; si = -s / u; ti = (r * s - t) / u**3
        if weq(Edst, ui**2 * x + ri, ui**3 * y + si * ui**2 * x + ti) == ui**6 * weq(Esrc, x, y):
            return ui, ri, si, ti
    raise RuntimeError('нет подходящего (u,r,s,t)')


# ---------- сертификат Безу ----------
def bezout_cert(Fc, Gc, d):
    """Fc, Gc: целые коэффициенты (от младшего) форм степени d (deg <= d).
    Ищем целые A_j, B_j степени <= d-1 и R != 0: A_j F + B_j G = R t^j, j in {0, 2d-1}.
    Возвращает R, свидетели, K = max ||A_j||_1 + ||B_j||_1, L = max(||F||_1, ||G||_1)."""
    n = 2 * d
    M = matrix(QQ, n, n)
    for i in range(d):
        for k, c in enumerate(Fc):
            M[i + k, i] += c
        for k, c in enumerate(Gc):
            M[i + k, d + i] += c
    assert M.det() != 0
    sols = {}
    for j in (0, n - 1):
        e = vector(QQ, n); e[j] = 1
        sols[j] = M.solve_right(e)
    Rres = lcm([v.denominator() for j in sols for v in sols[j]])
    Pz = PolynomialRing(ZZ, 't'); t = Pz.gen()
    F = Pz(list(Fc)); G = Pz(list(Gc))
    wit = []
    for j in (0, n - 1):
        w = [ZZ(Rres * v) for v in sols[j]]
        A = Pz(w[:d]); B = Pz(w[d:])
        assert A.degree() <= d - 1 and B.degree() <= d - 1
        assert A * F + B * G == Rres * t**j                     # точная проверка тождества над Z
        wit.append(dict(j=int(j), A=[str(c) for c in w[:d]], B=[str(c) for c in w[d:]],
                        norm=ZZ(sum(abs(c) for c in w))))
    K = max(wd['norm'] for wd in wit)
    Lb = max(sum(abs(c) for c in Fc), sum(abs(c) for c in Gc))
    return dict(R=Rres, K=K, L=Lb, wit=wit)


def dup_cert(E):
    """x(2P) = F(x)/G(x) на E (формула Сильвермана III.2.3), сертификат Безу для форм степени 4."""
    b2, b4, b6, b8 = E.b_invariants()
    assert all(b in ZZ for b in (b2, b4, b6, b8))
    Fc = [ZZ(-b8), ZZ(-2 * b6), ZZ(-b4), ZZ(0), ZZ(1)]
    Gc = [ZZ(b6), ZZ(2 * b4), ZZ(b2), ZZ(4), ZZ(0)]
    c = bezout_cert(Fc, Gc, 4)
    c['F'] = [str(v) for v in Fc]; c['G'] = [str(v) for v in Gc]
    return c, Fc, Gc


# ---------- кручение: независимая проверка ----------
def torsion_certified(E):
    """E(Q)[2^oo] — точным делением на 2; нечётная часть — оценка gcd #E(F_p) (p>=3 хорошей редукции,
    кручение вкладывается) и точное деление на l. Возвращает множество точек."""
    two = [E(0)]; front = [E(0)]
    for _ in range(8):                               # по Мазуру 2-часть <= 16; явная граница
        new = []
        for Q in front:
            for Rp in Q.division_points(2):
                if Rp not in two:
                    two.append(Rp); new.append(Rp)
        front = new
        if not front:
            break
    assert not front
    g = 0
    D = E.discriminant()
    for p in prime_range(3, NPRIME_TORS):
        if D % p == 0:
            continue
        g = gcd(g, p + 1 - E.ap(p))
    odd = g
    while odd % 2 == 0:
        odd //= 2
    oddpts = [E(0)]
    for l, _ in (list(factor(odd)) if odd > 1 else []):
        front = [E(0)]
        for _ in range(6):
            new = []
            for Q in front:
                for Rp in Q.division_points(l):
                    if Rp not in oddpts:
                        oddpts.append(Rp); new.append(Rp)
            front = new
            if not front:
                break
    tors = []
    for A in two:
        for B in oddpts:
            Q = A + B
            if Q not in tors:
                tors.append(Q)
    return tors, int(g)


# ---------- основной расчёт на одном множителе ----------
def run_T(r, s, T, S=None, okfun=None):
    t_start = time.time()
    S = cells(r, s) if S is None else S
    if okfun is None:
        okfun = lambda z: full_ok(z, S)
    T = [QQ(t) for t in T]
    assert all(t in S for t in T) and len(set(T)) == 3
    out = dict(slope=f'{r}/{s}', T=[int(t) for t in T])
    E1, L, g = model3(T)
    Em = E1.minimal_model()
    u, rr, ss, tt = urst_to(Em, E1)                  # X = u^2 x + rr
    out['E1'] = [str(a) for a in E1.ainvs()]; out['Emin'] = [str(a) for a in Em.ainvs()]
    out['urst'] = [str(u), str(rr), str(ss), str(tt)]; out['L'] = str(L)
    # ранг
    pe = pari(Em)
    # ellrank рандомизирован: фиксируем зерно; при [1,1] без точки — до 5 зёрен (детерминированно)
    for att in range(5):
        pari.setrand(att + 1)
        rk = pe.ellrank()
        lo, hi = int(rk[0]), int(rk[1])
        if (lo, hi) != (1, 1) or len(rk[3]) > 0:
            break
    out['ellrank'] = [lo, hi]; out['ellrank_seed'] = att + 1
    if (lo, hi) != (1, 1):
        out['status'] = f'не ранг [1,1]: {lo, hi}'
        return out
    try:
        alarm(120)
        M = mwrank_EllipticCurve(Em.ainvs()); M.two_descent(verbose=False)
        out['mwrank_rank_bound'] = int(M.rank_bound())
        cancel_alarm()
    except AlarmInterrupt:
        out['mwrank_rank_bound'] = 'таймаут'
    # кручение
    tors = Em.torsion_points()
    tors_c, gp = torsion_certified(Em)
    assert set(tors) == set(tors_c), 'кручение Sage != независимое'
    out['tors'] = len(tors); out['tors_gcd_Fp'] = gp
    out['tors_struct'] = [int(v) for v in Em.torsion_subgroup().invariants()]
    # образующая
    pts = [Em(list(map(QQ, p))) for p in rk[3]]
    pts = [p for p in pts if p.order() == oo]
    if not pts:                                   # PARI дал нижнюю 1 без точки — метод неприменим
        out['status'] = 'ранг [1,1] без явной точки (пропуск)'
        return out
    sat, idx, reg = Em.saturation([pts[0]])
    P0 = sat[0]
    assert P0.order() == oo
    out['P_ellrank'] = [str(pts[0][0]), str(pts[0][1])]
    out['sat_index_eclib'] = int(idx)
    out['P0'] = [str(P0[0]), str(P0[1])]
    # перекрёстно: граница индекса по нижней оценке Sage HeightFunction и PARI ellsaturation
    lam = Em.height_function().min(0.0001, 20)
    ibound = floor(sqrt(P0.height() / lam)) if lam > 0 else None
    out['lambda_sage'] = float(lam); out['index_bound'] = int(ibound) if ibound is not None else None
    if ibound is not None and ibound >= 2:
        V = pari.ellsaturation(pe, [pari([P0[0], P0[1]])], ibound)
        Q = Em(list(map(QQ, V[0])))
        out['pari_sat_same'] = bool((Q - P0).order() != oo or (Q + P0).order() != oo)
    else:
        out['pari_sat_same'] = 'индекс <= 1 тривиально'
    # сертификат удвоения
    cert, Fc, Gc = dup_cert(Em)
    Bup = RIF(cert['L']).log() / 3          # hhat - h(x) <= Bup
    Blow = RIF(cert['K']).log() / 3         # h(x) - hhat <= Blow
    Brig = max(Bup.upper(), Blow.upper())
    out['cert'] = dict(R=str(cert['R']), K=str(cert['K']), L=str(cert['L']), F=cert['F'], G=cert['G'],
                       wit=cert['wit'])
    out['B_up'] = float(Bup.upper()); out['B_low'] = float(Blow.upper()); out['B_rig'] = float(Brig)
    # контроль формулы удвоения на точках
    Px = PolynomialRing(QQ, 'x'); xx = Px.gen()
    Fp = Px(Fc); Gp = Px(Gc)
    for Q in [P0, 2 * P0, 3 * P0 + tors[-1], P0 + tors[1]]:
        if Q.is_zero() or (2 * Q).is_zero():
            continue
        assert (2 * Q)[0] == Fp(Q[0]) / Gp(Q[0])
    # строгая вилка hhat(P0)
    xk = P0[0]
    for _ in range(KDOUB):
        xk = Fp(xk) / Gp(xk)                  # x(2^k P0) — только по x
    hk = hgt(xk)
    hlo = (hk - Blow) / 4**KDOUB; hhi = (hk + Bup) / 4**KDOUB
    h0 = RIF(hlo.lower(), hhi.upper())
    assert h0.lower() > 0, 'нижний конец вилки hhat(P0) <= 0: M0 не определён'   # добавлено в 3-й попытке
    h_sage = P0.height()
    out['hP0_lo'] = float(h0.lower()); out['hP0_hi'] = float(h0.upper()); out['hP0_sage'] = float(h_sage)
    out['norm_ok'] = bool(h0.lower() <= h_sage <= h0.upper())
    # c1 из матрицы z = (u^2 x + rr)/L
    al = QQ(u**2) / L; be = QQ(rr) / L
    den = lcm(al.denominator(), be.denominator())
    a, b, d = ZZ(al * den), ZZ(be * den), ZZ(den)
    gg = gcd([a, b, d]); a, b, d = a // gg, b // gg, d // gg
    c1p = RIF(max(abs(a) + abs(b), abs(d))).log()       # h(z) <= h(x) + c1p
    c1m = RIF(max(abs(d) + abs(b), abs(a))).log()       # h(x) <= h(z) + c1m
    c1 = max(c1p.upper(), c1m.upper())
    out['matrix_z_of_x'] = [[str(a), str(b)], [0, str(d)]]
    out['c1'] = float(c1)
    zof = lambda x: (a * x + b) / d
    # константы и M0 (строго)
    C = 2 * RIF(Brig) + 2 * RIF(c1)
    M0 = ZZ(floor(((C / RIF(h0.lower()) + 1) / 2).upper()))
    out['C_rig'] = float(C.upper()); out['M0_rig'] = int(M0)
    # для сравнения: B Сильвермана (Sage), как у Fable (C = 4B + 2c1, hhat из Sage)
    BS = Em.silverman_height_bound()
    CF = 4 * BS + 2 * float(c1)
    M0F = floor((CF / float(h_sage) + 1) / 2) + 1
    out['B_silverman'] = float(BS); out['M0_fable_style'] = int(M0F)
    out['Brig_gt_BS'] = bool(Brig > BS)
    Mrun = max(int(M0), int(M0F)) + 1
    out['M_run'] = Mrun
    # ---- случай A: все nP0 + t, 0 <= n <= Mrun, t in Tors (z(-Q) = z(Q) даёт n < 0) ----
    bad = []; hits = []; ncand = 0; nonsq_T = 0
    maxdev = 0.0; maxc = 0.0; ineq_fail = 0
    Q = Em(0)
    for n in range(0, Mrun + 1):
        if n > 0:
            Q = Q + P0
        hn = RIF(n)**2 * h0
        for t in tors:
            Qt = Q + t
            if Qt.is_zero():
                continue
            ncand += 1
            z = zof(Qt[0])
            if not is_sq(prod(1 + l * z for l in T)):
                nonsq_T += 1
            if z != 0 and (okfun(z) or okfun(-z)):
                bad.append(str(z)); hits.append([str(z), n, str(t)])
            # контроль неравенств
            hx = hgt(Qt[0]); hz = hgt(z)
            dv = hn - hx                              # hhat - h(x), вилка
            if dv.lower() > Bup.upper() or (-dv).lower() > Blow.upper():
                ineq_fail += 1
            if (hz - hx).lower() > c1p.upper() or (hx - hz).lower() > c1m.upper():
                ineq_fail += 1
            if n <= 12:
                maxdev = max(maxdev, float(abs(dv).upper()))
                maxc = max(maxc, float(abs(hz - hx).upper()))
    out['caseA_cand'] = ncand; out['caseA_nondeg'] = bad; out['caseA_hits'] = hits
    out['ctrl_T_square_fail'] = nonsq_T; out['ctrl_ineq_fail'] = ineq_fail
    out['ctrl_max_dev_n12'] = maxdev; out['ctrl_max_c_n12'] = maxc
    # ---- случай B: X(P+t) + X(P) = 0 на E1 (a1=a3=0), t in Tors \ {O} ----
    a2 = E1.a2()
    RX = PolynomialRing(QQ, 'X'); X = RX.gen()
    f = RX(g)
    iso_E1 = lambda Pm: E1([u**2 * Pm[0] + rr, u**3 * Pm[1] + ss * u**2 * Pm[0] + tt])
    torsE1 = [iso_E1(t) for t in tors if not t.is_zero()]

    def caseB_poly(t, Kc):
        Xt, Yt = t[0], t[1]
        c = Kc + a2 + Xt
        A = f + Yt**2 - c * (X - Xt)**2                # условие: A(X) = 2 Yt Y
        return A if Yt == 0 else A**2 - 4 * Yt**2 * f
    caseB = []; badB = []
    for t in torsE1:
        pol = caseB_poly(t, 0)
        assert pol != 0
        for rt, _ in pol.roots(QQ):
            z = rt / L
            caseB.append(str(z))
            if z != 0 and (okfun(z) or okfun(-z)):
                badB.append(str(z))
    out['caseB_roots_z'] = sorted(set(caseB)); out['caseB_nondeg'] = badB
    # контроль машинерии случая B: известная точка находится
    okc = 0; totc = 0
    for m in (1, 2, 3):
        Pm = iso_E1(m * P0)
        for t in torsE1:
            Kc = (Pm + t)[0] + Pm[0]
            pol = caseB_poly(t, Kc)
            totc += 1
            if pol(Pm[0]) == 0:
                okc += 1
    out['ctrl_caseB'] = f'{okc}/{totc}'
    # контроль случая B по модулю p: перебор точек E1(F_p)
    modp_ok = 0; modp_tot = 0
    for p in (1009, 2003):
        if E1.discriminant() % p == 0 or any(v.denominator() % p == 0 for t in torsE1 for v in t):
            continue
        Ep = E1.change_ring(GF(p))
        ptsp = Ep.points()
        for t in torsE1:
            tp = Ep([GF(p)(t[0]), GF(p)(t[1])])
            pol = caseB_poly(t, 0).change_ring(GF(p))
            for Pp in ptsp:
                if Pp.is_zero() or (Pp + tp).is_zero() or Pp == tp or Pp == -tp:
                    continue
                if (Pp + tp)[0] + Pp[0] == 0:
                    modp_tot += 1
                    if pol(Pp[0]) == 0:
                        modp_ok += 1
    out['ctrl_caseB_modp'] = f'{modp_ok}/{modp_tot}'
    ok = (not bad and not badB and ineq_fail == 0 and nonsq_T == 0 and okc == totc and modp_ok == modp_tot
          and out['norm_ok'])
    out['closed_rigorous'] = bool(ok)
    out['status'] = 'закрыт строго (ПО)' if ok else 'НЕ ЗАКРЫТ / сбой контроля'
    out['sec'] = round(time.time() - t_start, 2)
    return out


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (bool, str, float)) or o is None:
        return o
    if isinstance(o, int):
        return o
    try:
        if o in ZZ:
            return int(o)
    except Exception:
        pass
    return str(o)


def classes3(r, s):
    """все 3-подмножества S по модулю T ~ -T"""
    S = [int(c) for c in cells(r, s)]
    seen = set(); res = []
    for T in itertools.combinations(sorted(S), 3):
        key = min(tuple(sorted(T)), tuple(sorted(-t for t in T)))
        if key in seen:
            continue
        seen.add(key); res.append(list(key))
    return res


def work(arg):
    r, s, T = arg
    try:
        return run_T(r, s, T)
    except Exception as e:
        import traceback
        return dict(slope=f'{r}/{s}', T=T, status='ОШИБКА', err=repr(e), tb=traceback.format_exc()[-1500:])


SLOPES = '143/206 96/211 61/217 19/221 211/230 204/247 231/250 65/261 43/278 265/298 108/301 251/308 237/317 140/319 197/325'.split()

if __name__ == '__main__':
    import os
    base = '/home/kep/magicKube/rigorous_census/chunk0'
    slopes = sys.argv[1:] or SLOPES
    tasks = []
    for sl in slopes:
        r, s = map(int, sl.split('/'))
        fab = json.load(open(f'/home/kep/magicKube/joint_sieve/dem_{r}_{s}.json'))
        fabT = [sorted(x['T']) for x in fab['results']]
        allT = classes3(r, s)
        # сначала множители Fable, затем все остальные классы (своя перепись)
        keyf = lambda T: min(tuple(sorted(T)), tuple(sorted(-t for t in T)))
        fk = set(keyf(T) for T in fabT)
        order = fabT + ([] if os.environ.get('ONLYFABLE') else [T for T in allT if keyf(T) not in fk])
        for T in order:
            tasks.append((r, s, T))
    print('задач', len(tasks), flush=True)
    res = {}
    t0 = time.time()
    with mp.get_context('fork').Pool(3) as pool:
        for i, o in enumerate(pool.imap_unordered(work, tasks)):
            res.setdefault(o['slope'], []).append(o)
            msg = f"[{i+1}/{len(tasks)} {time.time()-t0:.0f}s] {o['slope']} T={o['T']} {o.get('status')}"
            if o.get('ellrank') == [1, 1]:
                msg += f" M0={o.get('M0_rig')} Brig={o.get('B_rig', 0):.2f} BS={o.get('B_silverman', 0):.2f} hP0={o.get('hP0_lo', 0):.4f}"
            print(msg, flush=True)
    for sl, lst in res.items():
        r, s = sl.split('/')
        json.dump(clean(dict(slope=sl, results=lst,
                       closed_rigorous=any(o.get('closed_rigorous') for o in lst))),
                  open(f'{base}/rig_{r}_{s}.json', 'w'), ensure_ascii=False, indent=1)
    print('готово', time.time() - t0, flush=True)
