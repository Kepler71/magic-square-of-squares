# ⚠️ УСТАРЕЛО 15.09 (после рецензии Codex): для режима cells (V_q^nine) здесь ошибка — x = 0 добавляется в vals9 безусловно,
# и формулы хвостов x→0, x→∞ не учитывают масштаб u изоморфизма к минимальной модели. Заменено на omega_exact.py.
# Полное V_q (без cells) годится для сверки. lam_r и to_min отсюда используются omega_exact.
# Fable 15.09.2026. Прямое вычисление множества Ω для кода Bianchi–Padurariu (в его нормировке):
#   ω_q(z) = −λ_q^{E1min}(φ1 z) + λ_q^{E2min}(φ2 z) − 2 v_q(x(z)) log_p q,   Ω = { Σ_q ω_q(z_q) : z_q ∈ C(Q_q) }.
# Код Bianchi берёт декартово произведение множеств значений λ_q на E1, на E2 и чётных v(x) (Prop. 6.5 [Bia20]) —
# для наших кривых это даёт ~10^6 элементов. Здесь ω_q считается как функция от x ∈ Q_q напрямим перебором
# классов x = q^v·u (u mod q^K), т.е. берётся только реально достижимые пары компонент. Полнота перебора
# контролируется стабилизацией при K → K+1 (численно, не доказано).
# Локальная высота λ_q: алгоритм Сильвермана в точности как в non_archimedean_local_height кода (is_minimal=True).
import itertools
from sage.all import *

def lam_r(Emin, xp, yp, q):
    """r ∈ Q с λ_q(P) = r·log q для P=(xp,yp) ∈ Emin(Q_q), Emin — минимальная модель. Копия ветвления кода Bianchi."""
    a1, a2, a3, a4, a6 = Emin.a_invariants(); b2, b4, b6, b8 = Emin.b_invariants(); c4 = Emin.c4()
    N = Emin.discriminant().valuation(q)
    A = (3*xp**2 + 2*a2*xp + a4 - a1*yp).valuation()
    B = (2*yp + a1*xp + a3).valuation()
    C = (3*xp**4 + b2*xp**3 + 3*b4*xp**2 + 3*b6*xp + b8).valuation()
    if A <= 0 or B <= 0:
        r = max(0, -xp.valuation())
    elif c4.valuation(q) == 0:
        n = min(B, QQ(N)/2); r = -n*(N-n)/N
    elif C >= 3*B:
        r = -2*B/3
    else:
        r = -C/4
    return QQ(r)

def to_min(iso, x, y):
    """перенос точки (x,y) с E на Emin по изоморфизму iso = E.isomorphism_to(Emin), tuple (u,r,s,t)."""
    u, r, s, t = iso.tuple()
    xp = (x - r)/u**2
    yp = (y - s*(x - r) - t)/u**3
    return xp, yp

def values_q(f, q, K, vmax, cells=None):
    """множество коэффициентов c (ω_q = c·log q) по всем x ∈ Q_q: классы x = q^v·u, |v| ≤ vmax, u mod q^K, плюс x=0, x=∞."""
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    E1min = E1.minimal_model(); E2min = E2.minimal_model()
    i1 = E1.isomorphism_to(E1min); i2 = E2.isomorphism_to(E2min)
    N = max(E1min.discriminant().valuation(q), E2min.discriminant().valuation(q))
    prec = 3*N + K + 2*vmax + 30
    Kq = Qp(q, prec)
    fq = f.change_ring(Kq)
    vals = set(); vals9 = set()
    def is_sq(v):
        return v != 0 and v.valuation() % 2 == 0 and (v / Kq(q)**v.valuation()).is_square()
    def add_point(x):
        fx = fq(x)
        if fx == 0:
            y = Kq(0)
        else:
            if fx.valuation() % 2: return
            if not (fx / Kq(q)**fx.valuation()).is_square(): return
            y = fx.sqrt()
        x1, y1 = Kq(a6)*x**2, Kq(a6)*y
        x2, y2 = Kq(a0)/x**2, Kq(a0)*y/x**3
        r1 = lam_r(E1min, *to_min(i1, x1, y1), q)
        r2 = lam_r(E2min, *to_min(i2, x2, y2), q)
        c = -r1 + r2 - 2*x.valuation(); vals.add(c)
        if cells is not None and all(is_sq((x + sg*Kq(k))/x) for k in cells for sg in (1, -1)):
            vals9.add(c)
    units = [u for u in range(1, q**K) if u % q]
    for v in range(-vmax, vmax+1):
        for u in units:
            add_point(Kq(q)**v * Kq(u))
    # x = 0: φ2 → O, ω = −λ1((0, a6√a0)) − v(a0) log q  (если a0 — квадрат в Q_q)
    A0 = Kq(a0)
    if A0.valuation() % 2 == 0 and (A0/Kq(q)**A0.valuation()).is_square():
        y = A0.sqrt(); r1 = lam_r(E1min, *to_min(i1, Kq(0), Kq(a6)*y), q)
        vals.add(-r1 - ZZ(a0).valuation(q)); vals9.add(-r1 - ZZ(a0).valuation(q))
    # x = ∞: φ1 → O, ω = λ2((0, a0√a6)) + v(a6) log q (если a6 — квадрат в Q_q)
    A6 = Kq(a6)
    if A6.valuation() % 2 == 0 and (A6/Kq(q)**A6.valuation()).is_square():
        y = A6.sqrt(); r2 = lam_r(E2min, *to_min(i2, Kq(0), Kq(a0)*y), q)
        vals.add(r2 + ZZ(a6).valuation(q)); vals9.add(r2 + ZZ(a6).valuation(q))
    return (vals, vals9) if cells is not None else vals


def values_q_adaptive(f, q, Nmax, vmax, cells=None):
    """Адаптивный перебор классов x = q^v·(u + q^k·Z_q) для q ≥ 5.
    Три многочлена от x: F1 = f(x); F2 = cub_1(x_1'(x)), где cub_i(X) = (2Y+a1X+a3)^2 на E_i,min, x_1' = ((a6 x^2) − r_1)/u_1^2;
    F3 = x^6·cub_2(x_2'(x)), x_2' = (a0/x^2 − r_2)/u_2^2. Ветвление λ_q (алгоритм Сильвермана) и квадратность f(x) зависят
    только от валюаций этих величин (и от v(x)). Класс — лист, если валюация каждого F_i постоянна на классе:
    G(t) = F_i(q^v(u + q^k t)) — многочлен от t, и v(G[0]) < min_{j≥1} v(G[j]); иначе дробим (k → k+1), но не глубже Nmax+2
    (у точного корня — вейерштрассова точка или узел — при k ≥ Nmax+2 ветвление уже насыщено: n = min(B, N/2) = N/2)."""
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    E1min = E1.minimal_model(); E2min = E2.minimal_model()
    i1 = E1.isomorphism_to(E1min); i2 = E2.isomorphism_to(E2min)
    prec = 3*Nmax + 6*vmax + 40
    Kq = Qp(q, prec); fq = f.change_ring(Kq)
    R = PolynomialRing(Kq, 'X'); X = R.gen()
    def cub(E):
        a1_, a2_, a3_, a4_, a6_ = E.a_invariants()
        return 4*X**3 + (a1_**2 + 4*a2_)*X**2 + (2*a1_*a3_ + 4*a4_)*X + (a3_**2 + 4*a6_)
    u1, r1, s1, t1 = i1.tuple(); u2, r2, s2, t2 = i2.tuple()
    F1 = fq
    F2 = cub(E1min)((Kq(a6)*X**2 - r1)/u1**2)
    c2 = cub(E2min); F3 = sum(c2[m]*((Kq(a0) - r2*X**2)/u2**2)**m * X**(6-2*m) for m in range(4))   # x^6·cub_2(x_2')
    Fs = [F1, F2, F3] + ([X + Kq(k) for k in cells] + [X - Kq(k) for k in cells] if cells is not None else [])
    vals = set(); vals9 = set(); leaves = [0]
    def is_sq(v):
        return v != 0 and v.valuation() % 2 == 0 and (v / Kq(q)**v.valuation()).is_square()
    def const_val(G):
        v0 = G[0].valuation()
        return all(G[m].valuation() > v0 for m in range(1, G.degree()+1))
    def visit(v, u, k):
        sub = Kq(q)**v * (Kq(u) + Kq(q)**k * X)
        if k < Nmax + 2 and not all(const_val(F(sub)) for F in Fs):
            for d in range(q):
                uu = u + d*q**k
                if uu % q: visit(v, uu, k+1)
            return
        leaves[0] += 1
        x = Kq(q)**v * Kq(u); fx = fq(x)
        if fx == 0: y = Kq(0)
        else:
            if fx.valuation() % 2: return
            if not (fx / Kq(q)**fx.valuation()).is_square(): return
            y = fx.sqrt()
        x1, x2 = Kq(a6)*x**2, Kq(a0)/x**2
        rr1 = lam_r(E1min, *to_min(i1, x1, Kq(a6)*y), q)
        rr2 = lam_r(E2min, *to_min(i2, x2, Kq(a0)*y/x**3), q)
        c = -rr1 + rr2 - 2*x.valuation(); vals.add(c)
        if cells is not None and all(is_sq((x + sg*Kq(k))/x) for k in cells for sg in (1, -1)):
            vals9.add(c)
    for v in range(-vmax, vmax+1):
        for u in range(1, q):
            visit(v, u, 1)
    A0 = Kq(a0)
    if A0.valuation() % 2 == 0 and (A0/Kq(q)**A0.valuation()).is_square():
        y = A0.sqrt(); rr1 = lam_r(E1min, *to_min(i1, Kq(0), Kq(a6)*y), q); vals.add(-rr1 - ZZ(a0).valuation(q)); vals9.add(-rr1 - ZZ(a0).valuation(q))
    A6 = Kq(a6)
    if A6.valuation() % 2 == 0 and (A6/Kq(q)**A6.valuation()).is_square():
        y = A6.sqrt(); rr2 = lam_r(E2min, *to_min(i2, Kq(0), Kq(a0)*y), q); vals.add(rr2 + ZZ(a6).valuation(q)); vals9.add(rr2 + ZZ(a6).valuation(q))
    return (vals, vals9, leaves[0]) if cells is not None else (vals, leaves[0])

def bad_primes(f):
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    S = set([2]) | set(ZZ(a0).prime_factors()) | set(ZZ(a6).prime_factors())
    for E in (E1, E2):
        S |= set(ZZ(E.minimal_model().discriminant()).prime_factors())
        S |= set(ZZ(E.discriminant()/E.minimal_model().discriminant()).prime_factors())
    return sorted(S)

def omega_direct(f, p, n, verbose=True, check_stab=True, cells=None):
    """Ω как список p-адических чисел (Qp(p,n)) плюс словарь {q: множество коэффициентов}.
    cells=[a,b,c,d]: дополнительно V_q^nine — значения, достижимые лишь точками x ∈ Q_q, у которых все восемь клеток
    (x ± k)/x, k ∈ cells, — ненулевые квадраты в Q_q; тогда Ω строится из V_q^nine (для поиска точек «девяти квадратов»)."""
    a6, a4, a2, a0 = f[6], f[4], f[2], f[0]
    E1 = EllipticCurve([0, a4, 0, a2*a6, a0*a6**2]); E2 = EllipticCurve([0, a2, 0, a0*a4, a0**2*a6])
    E1min = E1.minimal_model(); E2min = E2.minimal_model()
    Kp = Qp(p, n); per = {}; per9 = {}
    cl = cells if cells is not None else []
    for q in bad_primes(f):
        if q == p: continue
        N = max(E1min.discriminant().valuation(q), E2min.discriminant().valuation(q))
        # для мультипликативной редукции I_N: λ_q зависит от x' mod q^(N//2+1) (см. записку); +1 запас при q<=3 (аддитивная редукция)
        K = N//2 + 1 + (1 if q <= 3 else 0)
        if q == 2: K = max(K, 4)
        vmax = (ZZ(a0).valuation(q) + ZZ(a6).valuation(q) + N)//2 + 2
        if q <= 3:
            S, S9 = values_q(f, q, K, vmax, cells=cl); stab = None
            if check_stab and q**(K+1) <= 30000:
                S2, S92 = values_q(f, q, K+1, vmax+1, cells=cl); stab = (S2 == S and S92 == S9); S = S | S2; S9 = S9 | S92
        else:
            S, S9, nl = values_q_adaptive(f, q, N, vmax, cells=cl); stab = f'адаптивно, листьев {nl}'
            if check_stab and q**(K+1) <= 30000:
                S2, S92 = values_q(f, q, K+1, vmax+1, cells=cl); stab += f', прямой перебор K+1 совпал: {S2 == S and S92 == S9}'; S = S | S2; S9 = S9 | S92
        per[q] = S; per9[q] = S9
        if verbose: print(f'  q={q}: N={N}, K={K}, |V_q|={len(S)}, стабильно при K+1: {stab}, V_q={sorted(S)}' + (f', V_q^nine={sorted(S9)}' if cells is not None else ''), flush=True)
    qs = sorted(per)
    logs = {q: Kp(q).log() for q in qs}
    src = per9 if cells is not None else per
    Om = []; tags = []
    for combo in itertools.product(*[sorted(src[q]) for q in qs]):
        Om.append(sum(c*logs[q] for c, q in zip(combo, qs))); tags.append(combo)
    # убрать совпадения (как f7 в коде Bianchi); тег — список всех комбинаций, давших это ω
    seen = {}; out = []; out_tags = []
    for w, tg in zip(Om, tags):
        k = str(w)
        if k in seen: out_tags[seen[k]].append(tg); continue
        seen[k] = len(out); out.append(w); out_tags.append([tg])
    if verbose: print(f'  |Ω_direct{"^nine" if cells is not None else ""}| = {len(out)} (произведение {len(Om)})', flush=True)
    omega_direct.last_tags = out_tags; omega_direct.last_qs = qs; omega_direct.last_per9 = per9
    return out, per
