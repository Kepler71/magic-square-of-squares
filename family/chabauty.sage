# Эллиптический Чабо для сечения (b,h,n), класса с закруткой D и известной точкой t0 = X0^2 (обобщение s15/ell_chabauty.sage).
# Аргументы: b h n D t0 p.   ПРЕДПОЛАГАЕТ rank E1(k) = 1.
import functools, itertools, sys
print = functools.partial(print, flush=True)
b0, h0, n0, D = [ZZ(x) for x in sys.argv[1:5]]
t0_arg = QQ(sys.argv[5]); p = int(sys.argv[6])
A, C = ZZ((h0^2 + n0^2)/2), ZZ((b0^2 + n0^2)/2)
PREC, DEG = 40, 30
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen(); Frz = Rz.fraction_field()
Rt2 = PolynomialRing(QQ, 'T2'); T2 = Rt2.gen()
f = T2*(T2^2 - 1)*(A*T2 - C)*(C*T2 - A)
Fz = Rz(Frz((1 - z)^6) * Frz(f((1 + z)/(1 - z))))
c = Fz.leading_coefficient(); q4 = Rz(Fz/(c*z)); Aq, B = q4[2], q4[0]
beta = (C - A)/(C + A); assert B == beta^2
def sqf(n):
    n = ZZ(n); return sign(n) * prod(pr^(e % 2) for pr, e in n.abs().factor())
dk = sqf((C - A)*(C + A))
k = QuadraticField(dk, 'r15'); r15 = k.gen()
sb = k(beta).sqrt()
Ru = PolynomialRing(k, 'uu'); uu = Ru.gen()
cub = Ru(D * c * (uu^2 + Aq - 2*beta) * (uu + 2*sb))
a3, a2, a1, a0 = [cub[i] for i in (3, 2, 1, 0)]
E1 = EllipticCurve(k, [0, a2, 0, a1*a3, a0*a3^2])
assert a3 in QQ
assert kronecker(dk, p) == 1 and all(not E1.has_bad_reduction(P) for P in k.primes_above(p)), "p must split and be good"
A = Aq   # далее в коде A — коэффициент квартики (как в s15)
# --- генератор (образ точки t = 9/4) и насыщение
def phi(t0, W0, gam):
    z0 = (t0 - 1)/(t0 + 1); y0 = W0*(1 - z0)^3
    u0 = z0 + beta/z0; v1 = y0*(1 + gam/z0^3)/(u0 - sb)
    X_, Y_ = a3*u0, a3*D*v1
    return E1(X_, Y_) if Y_^2 == X_^3 + a2*X_^2 + a1*a3*X_ + a0*a3^2 else None
t0 = t0_arg; W0 = (f(t0)/D).sqrt()
cands = [phi(t0, s_*W0, g_) for s_ in (1, -1) for g_ in (sb^3, -sb^3)]
cands = [P for P in cands if P is not None]
print('images of t=9/4 on E1:', len(cands))
G0 = cands[0]
print("G0 height:", RR(G0.height()))
G = G0   # насыщение проверяется ЛОКАЛЬНО ниже (по Codex: достаточно l | pL), глобальные высотные границы не используются
T2tors = [E1(0)] + [E1(e, 0) for e in (E1.two_division_polynomial().roots(multiplicities=False))]
print("E1(k)[2] =", len(T2tors), "points;  x in Q for:", [Tt.is_zero() or Tt[0] in QQ for Tt in T2tors])

# --- p-адические вложения
Kp = Qp(p, PREC)
sq = Kp(dk).sqrt()
emb = [k.hom([sq], Kp), k.hom([-sq], Kp)]
Ei = [EllipticCurve(Kp, [e_(a) for a in E1.ainvs()]) for e_ in emb]
def ptmap(i, P):
    return Ei[i](0) if P.is_zero() else Ei[i](emb[i](P[0]), emb[i](P[1]))
Gi = [ptmap(i, G) for i in range(2)]
Fp = GF(p)
def red_order(i):
    Er = EllipticCurve(Fp, [Fp(a) for a in Ei[i].ainvs()])
    Pt = Gi[i]
    return Er(Fp(Pt[0]), Fp(Pt[1])).order() if Pt[0].valuation() >= 0 else 1
N = lcm(red_order(0), red_order(1))
print(f"p = {p}: orders of G mod the two primes: {red_order(0)}, {red_order(1)};  N = {N}")
# локальное насыщение (Codex, CHABAUTY_AUDIT §2): при rank 1 и gcd([Gamma:H], pL) = 1 замыкания H и Gamma в E(Q_p)^2 совпадают
def red_exponent(i):
    Er = EllipticCurve(Fp, [Fp(a) for a in Ei[i].ainvs()])
    return Er.abelian_group().exponent()
Lexp = lcm(red_exponent(0), red_exponent(1))
ells = sorted(set(ZZ(p*Lexp).prime_factors()))
loc_sat = {l: E1.saturation([G], one_prime=l)[1] for l in ells}
print(f"local saturation: L = {Lexp}, primes l | pL = {ells}, indices = {loc_sat}")
assert all(v == 1 for v in loc_sat.values()), "G not l-saturated for some l | pL"
NG = [N*Gi[i] for i in range(2)]
FG = [Ei[i].formal_group() for i in range(2)]
tpar = [-NG[i][0]/NG[i][1] for i in range(2)]
assert all(tp.valuation() >= 1 for tp in tpar)
LOG = [FG[i].log(DEG) for i in range(2)]
EXP = [LOG[i].reverse() for i in range(2)]
ell = [LOG[i](tpar[i]) for i in range(2)]
print("formal logs of NG valuations:", [x.valuation() for x in ell])
Sm.<m> = PowerSeriesRing(Kp, default_prec=DEG)
tau = [EXP[i].change_ring(Kp)(ell[i]*m) for i in range(2)]      # параметр точки m*NG

def in_E1(Pt):
    return Pt.is_zero() or Pt[0].valuation() < 0

def xseries(i, Q):
    # x(Q + [tau_i(m)]) как ряд по m; Q вне E^1
    St.<tt> = LaurentSeriesRing(Kp, default_prec=DEG)
    xR = FG[i].x(DEG).change_ring(Kp)(tt); yR = FG[i].y(DEG).change_ring(Kp)(tt)
    xQ, yQ = Q[0], Q[1]
    lam = (yR - yQ)/(xR - xQ)
    x3 = lam^2 - Ei[i].a2() - xQ - xR
    x3 = x3.power_series()
    return x3(tau[i])

Pm.<mm> = PolynomialRing(Kp)
def trunc(pol, n=DEG):
    return Pm(pol.list()[:n])
def wseries(i, Q):
    # 1/x для Q в E^1: параметр суммы = exp(log t_Q + m ell); композиция усечённых многочленов над Q_p
    St.<tt> = LaurentSeriesRing(Kp, default_prec=DEG)
    w = (1/FG[i].x(DEG).change_ring(Kp)(tt)).power_series()      # 1/x(t) = t^2 + ...
    LQ = Kp(0) if Q.is_zero() else LOG[i](-Q[0]/Q[1])
    Wp = Pm(w.polynomial().list()); Ep_ = Pm(EXP[i].change_ring(Kp).polynomial().list())
    inner = trunc(Ep_(LQ + ell[i]*mm))
    res = trunc(Wp(inner))
    return Sm(res.list(), DEG)

SHIFT = 0
def strassmann(th):
    cs = th.list()
    vals = [cq.valuation() if cq != 0 else oo for cq in cs]
    mv = min(vals)
    J = max(j for j, v in enumerate(vals) if v == mv)
    # хвост: проверяем, что последние коэффициенты заметно больше минимума (контроль точности ряда)
    # строгий хвост (Codex): для n >= DEG  v_p(c_n) >= n - v_p(n!) >= n(p-2)/(p-1); сдвиг shift для рядов, делённых на m^shift
    tail_bound = min(n_ - valuation(factorial(n_), p) for n_ in range(len(cs) + SHIFT, 3*len(cs) + 60))
    tail_ok = (tail_bound > mv) and all(v > mv for v in vals[J+1:])
    return J, mv, tail_ok

total = 0; classes = 0; cut = 0; report = []
for r in range(N):
    for Tt in T2tors:
        Q = [ptmap(i, Tt) + r*Gi[i] for i in range(2)]
        inside = [in_E1(Q[i]) for i in range(2)]
        classes += 1
        if inside[0] != inside[1]:
            cut += 1; continue                                      # отсечение: разные валюации
        if not inside[0]:
            if Fp(Q[0][0]) != Fp(Q[1][0]):
                cut += 1; continue                                  # отсечение mod p по x
            th = xseries(0, Q[0]) - xseries(1, Q[1])
            J, mv, ok = strassmann(th)
        else:
            th = wseries(0, Q[0]) - wseries(1, Q[1])
            if Q[0].is_zero():
                th = th.shift(-2)                                   # двойной нуль m=0 у w; сама m=0 — точка O
                SHIFT = 2; J, mv, ok = strassmann(th); SHIFT = 0; J += 1   # +1 за решение m=0 (точка O)
            else:
                J, mv, ok = strassmann(th)
        total += J
        report.append((r, T2tors.index(Tt), J, mv, ok))
print(f"classes: {classes}, cut by mod-{p} test: {cut}, surviving: {classes - cut}")
for rr in report:
    print(f"   r={rr[0]:3d} T#{rr[1]}: Strassmann bound {rr[2]} (min val {rr[3]}, tail ok {rr[4]})")
assert all(rr[4] for rr in report), "tail bound failed in some class"
print("all classes: rigorous tail bound OK")
print("TOTAL upper bound for #{P in E1(k): x(P) in Q} (given rank 1):", total)
# известные такие точки
known = [Tt for Tt in T2tors if Tt.is_zero() or Tt[0] in QQ] + [P for P in cands] + [P + Tt for P in cands for Tt in T2tors]
known = list(set(P for P in known if P.is_zero() or P[0] in QQ))
print("known P with x(P) in Q:", len(known), [("O" if P.is_zero() else QQ(P[0]/a3)) for P in known])

# --- корни theta в Z_p для выживших классов
def theta_poly(r, Tt):
    Q = [ptmap(i, Tt) + r*Gi[i] for i in range(2)]
    if not in_E1(Q[0]):
        th = xseries(0, Q[0]) - xseries(1, Q[1])
    else:
        th = wseries(0, Q[0]) - wseries(1, Q[1])
    return Pm(th.list()[:DEG])
print("\nroots of theta in Z_p (m = (n - r)/N), with rational guesses:")
for rr in report:
    r, ti = rr[0], rr[1]
    th = theta_poly(r, T2tors[ti])
    # нормируем, чтобы коэффициенты были целыми, и берём корни в Z_p
    v0 = min(cq.valuation() for cq in th.list() if cq != 0)
    thn = Pm([cq / Kp(p)^v0 for cq in th.list()])
    try:
        rts = [(rt, mult) for rt, mult in thn.roots() if rt.valuation() >= 0]
    except Exception as ex:
        rts = f"roots failed: {ex}"
    out = []
    if isinstance(rts, str):
        out = rts
    else:
        for rt, mult in rts:
            try:
                g = QQ(rt.rational_reconstruction()) if hasattr(rt, 'rational_reconstruction') else None
            except Exception:
                g = None
            out.append((str(rt.add_bigoh(8)), mult, g, (r + N*g) if g is not None else None))
    print(f"   r={r} T#{ti}: {out}")

# --- точность «призрачных» корней m = -1/2: алгебраическая идентификация
print("\nexactness of m = -1/2 roots:")
G3 = [(N//2)*Gi[i] for i in range(2)]
def tau_index(i):
    # 2-кручение tau с sigma_i(3G) - sigma_i(tau) в E^1
    for j, Tt in enumerate(T2tors):
        if in_E1(G3[i] - ptmap(i, Tt)): return j
    return None
ti_ = [tau_index(0), tau_index(1)]
print("   tau index per embedding:", ti_, " (N/2 =", N//2, ")")
assert N % 2 == 0
for (r, ti, J, mv, ok) in report:
    th = theta_poly(r, T2tors[ti])
    val = th(Kp(-1)/2)
    if val.valuation() < 8: continue                   # -1/2 не корень этого класса (корни найдены с точностью >= p^8)
    A = [T2tors[ti] + (r - N//2)*G - T2tors[ti_[i]] for i in range(2)]   # глобальные точки для каждого вложения
    xa = [("O" if Pt.is_zero() else Pt[0]) for Pt in A]
    same = A[0] == A[1]
    exact = same and (A[0].is_zero() or A[0][0] in QQ)
    if not same and not A[0].is_zero() and not A[1].is_zero():
        conj = k.hom([-r15])
        exact = (A[0][0] == conj(A[1][0]))
    print(f"   class r={r} T#{ti}: theta(-1/2) = O(p^{val.valuation()});  points per embedding equal: {same};  exact zero: {exact};  point: {xa[0] if same else xa}")
