# Fable 15.09.2026, после рецензии Codex (RESULT_JOINT_149_QC_FILTER). Точное перечисление значений ω_q на C(Q_q).
#
# Ключевое наблюдение: все многочлены, от валюаций (и вычетов единичных частей) которых зависит ветвление λ_q
# (алгоритм Сильвермана, копия ветви кода Bianchi) и квадратность f(x), (x ± k)/x, имеют ТОЛЬКО рациональные корни:
#   f(x): ±a, ±b, ±c;  (2Y+a1X+a3)^2 на E_{i,min} в точке φ_i(x): те же ±a,±b,±c (полное 2-кручение обеих E_i),
#   x^6·cub_2(a0/x^2): дополнительно x = 0;  клетки x ± k: ∓k (k ∈ {a,b,c,d}).
# Поэтому класс x = q^v(u0 + q^m Z_q), не содержащий корней и отделённый от них с запасом μ (μ=1 при нечётном q, μ=3 при q=2),
# имеет постоянные валюации и постоянные вычеты mod q^μ всех этих величин, т.е. ω_q и «девятиклеточность» на нём постоянны.
# Классы, содержащие корень e, дробятся до глубины, где корень один; затем перечисляются оболочки v(x−e) = v+mm, mm < M0
# (каждая — конечное число классов с постоянными валюациями), а хвост v(x−e) ≥ v+M0 насыщен: n = min(B, N/2) = N/2
# (мультипликативная) либо C < 3B (аддитивная) — значение ω там равно значению в самом корне (y = 0 для вейерштрассовых e).
# Хвосты по v (x → 0, x → ∞) насыщаются аналогично; для «nine» хвост x→0 возможен лишь при −1 ∈ Q_q^2 и k/k' ∈ Q_q^2.
# Всё это — доказуемые утверждения о конечном переборе; сертификат — распечатка порогов и проверка κ на двух глубинах.
import itertools
from sage.all import *
from omega_direct import lam_r, to_min

def _is_sq(v, q):
    return v != 0 and v.valuation() % 2 == 0 and (v / v.parent()(q)**v.valuation()).is_square()

def values_q_exact(f, q, cells, verbose=False):
    """Возвращает (V_q, V_q^nine, cert) — множества коэффициентов c (ω_q = c·log q) по всем x ∈ C(Q_q) и по x с восемью
    клетками-квадратами; cert — словарь порогов/проверок."""
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    assert a6 == 1
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    E1min = E1.minimal_model(); E2min = E2.minimal_model()
    i1 = E1.isomorphism_to(E1min); i2 = E2.isomorphism_to(E2min)
    N1 = ZZ(E1min.discriminant()).valuation(q); N2 = ZZ(E2min.discriminant()).valuation(q)
    abc = sorted(set(abs(QQ(r)) for r, _ in f.roots()))          # a, b, c
    assert len(abc) == 3
    roots = sorted(set([QQ(s*e) for e in abc for s in (1, -1)] + [QQ(s*k) for k in cells for s in (1, -1)]))
    mu = 3 if q == 2 else 1
    prec = 4*max(N1, N2) + 60
    Kq = Qp(q, prec); fq = f.change_ring(Kq)
    def cubpol(E):
        R = PolynomialRing(Kq, 'X'); X = R.gen(); a1_, a2_, a3_, a4_, a6_ = E.a_invariants()
        return 4*X**3 + (a1_**2 + 4*a2_)*X**2 + (2*a1_*a3_ + 4*a4_)*X + (a3_**2 + 4*a6_)
    cub1, cub2 = cubpol(E1min), cubpol(E2min)
    def psi3(E):
        R = PolynomialRing(Kq, 'X'); X = R.gen(); b2, b4, b6, b8 = E.b_invariants()
        return 3*X**4 + b2*X**3 + 3*b4*X**2 + 3*b6*X + b8
    psi1, psi2 = psi3(E1min), psi3(E2min)
    vals = set(); vals9 = set(); cert = {}
    def omega_at(x, y):
        x1, x2 = Kq(a6)*x**2, Kq(a0)/x**2
        r1 = lam_r(E1min, *to_min(i1, x1, Kq(a6)*y), q)
        r2 = lam_r(E2min, *to_min(i2, x2, Kq(a0)*y/x**3), q)
        return -r1 + r2 - 2*x.valuation()
    def nine_at(x, skip=None):
        return all(_is_sq((x + sg*Kq(k))/x, q) for k in cells for sg in (1, -1) if not (skip is not None and x + sg*Kq(k) == 0) and not (skip is not None and QQ(-sg*k) == skip))
    def leaf(x):
        """x — представитель класса с постоянными валюациями/вычетами: добавить значения."""
        fx = fq(x)
        if fx == 0:
            y = Kq(0)
        else:
            if not _is_sq(fx, q): return
            y = fx.sqrt()
        c = omega_at(x, y); vals.add(c)
        if nine_at(x): vals9.add(c)
    # --- порог насыщения около вейерштрассова корня e: глубина, начиная с которой обе λ_i постоянны
    # --- точные тождества факторизации (проверяются над Q): cub_1(X_1'(x)) = 4 f(x)/u1^6,  x^6 cub_2(X_2'(x)) = 4 a0^2 f(x)/u2^6
    Rx = PolynomialRing(QQ, 'x'); xq_ = Rx.gen()
    def cub_Q(E):
        a1_, a2_, a3_, a4_, a6_ = E.a_invariants()
        return lambda X: 4*X**3 + (a1_**2 + 4*a2_)*X**2 + (2*a1_*a3_ + 4*a4_)*X + (a3_**2 + 4*a6_)
    U1, R1, S1, T1 = [QQ(t) for t in i1.tuple()]; U2, R2, S2, T2 = [QQ(t) for t in i2.tuple()]
    X1Q = (QQ(a6)*xq_**2 - R1)/U1**2
    assert cub_Q(E1min)(X1Q) == 4*f/U1**6, 'тождество cub_1 нарушено'
    Fx = Rx.fraction_field()(xq_); X2Q = (QQ(a0)/Fx**2 - R2)/U2**2
    assert Fx**6 * cub_Q(E2min)(X2Q) == 4*QQ(a0)**2*f/U2**6, 'тождество cub_2 нарушено'
    vU1, vU2 = ZZ(U1).valuation(q), ZZ(U2).valuation(q)
    six = [QQ(s_*e_) for e_ in abc for s_ in (1, -1)]
    def A_and_C_at(e, which):
        """валюации A = v(3x'^2+2a2x'+a4 − a1 y') и C = v(ψ3(x')) в точной 2-точке кручения φ_i(e) на E_{i,min}; обе конечны."""
        E, iso, psi = (E1min, i1, psi1) if which == 1 else (E2min, i2, psi2)
        xi = Kq(a6)*Kq(e)**2 if which == 1 else Kq(a0)/Kq(e)**2
        xe, ye = to_min(iso, xi, Kq(0))
        a1_, a2_, a3_, a4_, a6_ = E.a_invariants()
        # y' в 2-точке кручения: 2y'+a1x'+a3 = 0
        ye = -(a1_*xe + a3_)/2
        Aval = (3*xe**2 + 2*a2_*xe + a4_ - a1_*ye).valuation(); Cval = psi(xe).valuation()
        assert Aval < prec - 5 and Cval < prec - 5, ('A или C бесконечна', q, e, which)
        return Aval, Cval, xe
    def M0_for(e):
        """порог mm: при v(x−e) ≥ v(e)+mm обе λ_i постоянны (равны значению в e). Явно:
        v(cub_i) = v(x−e) + κ_i, κ_i = Σ_{j≠e} v(e−e_j) + v(4) − 6v(u_i) [+ 2v(a0) − 6v(e) при i=2];
        нужно v(cub_i) ≥ max(N_i+1, ⌊2C_i/3⌋+2) (насыщение n = N/2, либо C < 3B), и v(x_i' − e_i') > max(A_i, C_i)
        (постоянство ветви A и значения C), где v(x_1'−e_1') = v(x−e) + v(x+e) − 2v(u1),
        v(x_2'−e_2') = v(a0) + v(x−e) + v(x+e) − 4v(e) − 2v(u2)."""
        ve = ZZ(e.numerator()).valuation(q) - ZZ(e.denominator()).valuation(q)
        others = sum((QQ(e) - ej).valuation(q) for ej in six if ej != e)   # v(e−e_j), j ≠ e  (валюации рациональных чисел)
        v2e = (2*QQ(e)).valuation(q)
        need = []
        for which, vU, N in ((1, vU1, N1), (2, vU2, N2)):
            Aval, Cval, _ = A_and_C_at(e, which)
            kappa = others + ZZ(4).valuation(q) - 6*vU + (2*ZZ(a0).valuation(q) - 6*ve if which == 2 else 0)
            need.append(max(N + 1, (2*Cval)//3 + 2) - kappa)                         # по v(cub_i)
            shift = (v2e - 2*vU) if which == 1 else (ZZ(a0).valuation(q) + v2e - 4*ve - 2*vU)
            need.append(max(Aval, Cval) + 1 - shift)                                  # по v(x_i' − e_i')
            cert.setdefault('root_thresholds', {})[f'{e}/{which}'] = dict(kappa=int(kappa), A=int(Aval), C=int(Cval), N=int(N))
        return max(need + [1]) - ve      # в единицах mm (v(x−e) = ve + mm)
    sep = max([ (e1 - e2).valuation(q) for e1 in roots for e2 in roots if e1 != e2 ] + [0])
    # --- пороги хвостов по v
    vmaxroot = max(ZZ(e.numerator()).valuation(q) - ZZ(e.denominator()).valuation(q) for e in roots)
    # x→0: φ1 → P0=(0, a6√a0) (если a0 квадрат); λ1 постоянна, когда v(x1'(x) − x1'(0)) > v(cub1(x1'(0))) (P0 не 2-кручение)
    x10 = to_min(i1, Kq(0), Kq(0))[0]; cub1_0 = cub1(x10).valuation() if cub1(x10) != 0 else 0
    u1 = i1.tuple()[0]; u2 = i2.tuple()[0]
    V0 = max(vmaxroot + 1, (cub1_0 + 2*ZZ(u1).valuation(q))//2 + 2, (ZZ(a0).valuation(q))//2 + 2)
    x2inf = to_min(i2, Kq(0), Kq(0))[0]; cub2_inf = cub2(x2inf).valuation() if cub2(x2inf) != 0 else 0
    Vinf = max(vmaxroot + 1, (cub2_inf + 2*ZZ(u2).valuation(q) + ZZ(a0).valuation(q))//2 + 2, 1) + mu + max(ZZ(k).valuation(q) for k in cells)
    Vt = max(V0, Vinf) + max(N1, N2) + 2
    cert.update(dict(N1=N1, N2=N2, mu=mu, sep=sep, Vt=Vt, roots=[str(e) for e in roots]))
    units = [u for u in range(1, q**mu) if u % q]
    additive = [E1min.c4().valuation(q) > 0 and N1 > 0, E2min.c4().valuation(q) > 0 and N2 > 0]
    cert['additive'] = additive
    RT = PolynomialRing(Kq, 'T'); T = RT.gen()
    def psi_const(v, u0, m):
        """при аддитивной редукции: v(ψ3(x_i'(x))) постоянна на классе x = q^v(u0 + q^m Z_q)? (достаточный критерий Тейлора)"""
        ok = True
        for which, psi, iso in ((1, psi1, i1), (2, psi2, i2)):
            if not additive[which-1]: continue
            xt = Kq(q)**v * (Kq(u0) + Kq(q)**m * T)
            xi = Kq(a6)*xt**2 if which == 1 else RT.fraction_field()(Kq(a0))/xt**2
            uu, rr, ss, tt = iso.tuple()
            G = psi((xi - Kq(rr))/Kq(uu)**2)
            if which == 2: G = G * xt**8          # очистить знаменатель x^8, валюация сдвигается на константу 8v
            G = RT(G) if which == 1 else RT(G.numerator()) if hasattr(G, 'numerator') else RT(G)
            v0 = G[0].valuation() if G[0] != 0 else Infinity
            if not (v0 < Infinity and all((G[j].valuation() if G[j] != 0 else Infinity) > v0 for j in range(1, G.degree()+1))):
                ok = False
        return ok
    # --- перечисление по уровням v
    def enumerate_level(v):
        Rv = [e for e in roots if ZZ(e.numerator()).valuation(q) - ZZ(e.denominator()).valuation(q) == v]
        Rv_u = {e: Kq(e) / Kq(q)**v for e in Rv}     # e' = e/q^v — единицы
        def rec(u0, m):
            # класс x = q^v (u0 + q^m Z_q); u0 — целое, единица
            inside = [e for e in Rv if (Rv_u[e] - Kq(u0)).valuation() >= m]
            if not inside:
                if all(m - (Rv_u[e] - Kq(u0)).valuation() >= mu for e in Rv) and (not any(additive) or psi_const(v, u0, m)):
                    leaf(Kq(q)**v * Kq(u0)); return
                if m > 3*max(N1, N2) + 12:
                    cert.setdefault('uncertified_classes', []).append((v, u0, m)); leaf(Kq(q)**v * Kq(u0)); return
                for d in range(q): rec(u0 + d*q**m, m+1)
                return
            if len(inside) > 1 or m < sep - v + mu + 1:
                for d in range(q): rec(u0 + d*q**m, m+1)
                return
            e = inside[0]
            M0 = max(M0_for(e), m) if abs(e) in abc else m
            cert.setdefault('near_root', {})[str(e)] = dict(v=v, m=m, M0=M0)
            for mm in range(m, M0):
                for up in units:
                    leaf(Kq(e) + Kq(q)**(v+mm) * Kq(up))
            # хвост v(x−e) ≥ v+M0: значение в самом e
            xe = Kq(e); fe = fq(xe)
            if fe == 0: y = Kq(0)
            else:
                if not _is_sq(fe, q): return
                y = fe.sqrt()
            c = omega_at(xe, y); vals.add(c)
            if nine_at(xe, skip=e): vals9.add(c)       # клетка, обращающаяся в 0 при x=e, в хвосте может быть квадратом
        for u0 in units: rec(u0, mu)
    for v in range(-Vt, Vt+1):
        enumerate_level(v)
    # --- хвост x→0 (v ≥ Vt+1): ω = ω(P0) − v(a0); nine возможно лишь при −1 ∈ Q_q^2 и k/k' ∈ Q_q^2
    A0 = Kq(a0)
    if _is_sq(A0, q):
        y = A0.sqrt(); r1 = lam_r(E1min, *to_min(i1, Kq(0), Kq(a6)*y), q); c0 = -r1 - ZZ(a0).valuation(q) + 2*ZZ(u2).valuation(q)
        vals.add(c0)
        nine0 = _is_sq(Kq(-1), q) and all(_is_sq(Kq(k)/Kq(cells[0]), q) for k in cells)
        if nine0: vals9.add(c0)
        # численная сверка хвоста: ω на уровнях Vt+1, Vt+2 (представители с единицей 1)
        chk = []
        for vv in (Vt+1, Vt+2):
            for uu in units:
                x = Kq(q)**vv * Kq(uu); fx = fq(x)
                if _is_sq(fx, q): chk.append(omega_at(x, fx.sqrt()))
        cert['tail0'] = dict(c=str(c0), nine=nine0, check=sorted(set(str(t) for t in chk)), ok=(len(chk) > 0 and all(t == c0 for t in chk)))
    # --- хвост x→∞ (v ≤ −Vt−1): ω = ω((0, a0√a6)) + v(a6); nine всегда (1 ± k/x ≡ 1 mod q^μ)
    cinf_pt = lam_r(E2min, *to_min(i2, Kq(0), Kq(a0)*Kq(a6).sqrt()), q) + ZZ(a6).valuation(q) - 2*ZZ(u1).valuation(q)
    vals.add(cinf_pt); vals9.add(cinf_pt)
    chk = []
    for vv in (Vt+1, Vt+2):
        for uu in units:
            x = Kq(q)**(-vv) * Kq(uu); fx = fq(x)
            if _is_sq(fx, q): chk.append(omega_at(x, fx.sqrt()))
    cert['tailinf'] = dict(c=str(cinf_pt), check=sorted(set(str(t) for t in chk)), ok=(len(chk) > 0 and all(t == cinf_pt for t in chk)))
    return vals, vals9, cert

def bad_primes(f):
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    S = set([2]) | set(ZZ(a0).prime_factors())
    for E in (E1, E2):
        S |= set(ZZ(E.minimal_model().discriminant()).prime_factors())
        S |= set(ZZ(E.discriminant()/E.minimal_model().discriminant()).prime_factors())
    return sorted(S)

def omega_exact(f, p, n, cells, verbose=True):
    """Ω^nine (список p-адических чисел), все теги для каждого ω, per (V_q), per9 (V_q^nine), сертификаты."""
    Kp = Qp(p, n); per = {}; per9 = {}; certs = {}
    for q in bad_primes(f):
        if q == p: continue
        S, S9, cert = values_q_exact(f, q, cells)
        per[q] = S; per9[q] = S9; certs[q] = cert
        if verbose:
            print(f'  q={q}: N=({cert["N1"]},{cert["N2"]}) Vt={cert["Vt"]} sep={cert["sep"]} V_q={sorted(S)} V_q^nine={sorted(S9)}'
                  + (f' хвост0 ok={cert["tail0"]["ok"]} nine0={cert["tail0"]["nine"]}' if 'tail0' in cert else ' хвост0: a0 не квадрат')
                  + f' хвост∞ ok={cert["tailinf"]["ok"]}', flush=True)
    qs = sorted(per); logs = {q: Kp(q).log() for q in qs}
    seen = {}; out = []; out_tags = []
    for combo in itertools.product(*[sorted(per9[q]) for q in qs]):
        w = sum(c*logs[q] for c, q in zip(combo, qs)); k = str(w)
        if k in seen: out_tags[seen[k]].append(combo); continue
        seen[k] = len(out); out.append(w); out_tags.append([combo])
    if verbose: print(f'  |Ω^nine| = {len(out)}; совпадений ω: {sum(len(t) > 1 for t in out_tags)}', flush=True)
    return out, out_tags, qs, per, per9, certs
