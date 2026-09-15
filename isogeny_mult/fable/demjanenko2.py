# Fable 15.09. Обобщённый метод Демьяненко–Манина для кривой наклона: класс изоморфных множителей E_T ≅ E (кратность m),
# ранг E(ℚ) = r ∈ {1, 2}. Морфизмы f_T = φ_T ∘ π_T : X_k → E; Q_T = f_T(z) ∈ E(ℚ) при рациональном решении z.
#   Нормы:  |ĥ(Q_T) − h(z)| ≤ C_T = SAFETY·B + c₁(T)   (Силверман на общей минимальной модели + мёбиусова связь x_min ↔ z).
#   Спаривания: для пары (T,T′) с |TΔT′| ≤ 4 (третья кривая рода ≤ 1) x(Q_T ± Q_T′) = u_±(R), R = P_{TΔT′}(z);
#     сумма полюсов u_+ — кручение (S = O: проверено численно на всех парах 204/247, 19/60; символьно — check_S_symbolic),
#     поэтому 4⟨Q_T,Q_T′⟩ = h(u_+(R)) − h(u_+(−R)) ± 2B ограничено; |Λ| = |h(u_+(R)) − h(u_−(R))| оценивается ЧИСЛЕННО
#     по выборке точек R (константа ЧИСЛЕННАЯ, с запасом SAFETY_P; для рода 0 — по выборке параметра коники).
#   r = 1:  Q₁ = n₁G + t₁, Q₂ = n₂G + t₂; |n₁n₂|ĥ(G) ≤ c₁₂ и n_i²ĥ(G) ≥ h − C_i ⇒ |n₁| ≤ M₀ — перебор (случай A).
#   r = 2:  хаб T₁ и соседи T₂, T₃ (|T₁ΔT_j| ≤ 4): Q₂, Q₃ почти ортогональны Q₁ ⇒ ĥ(Q₂ ∓ Q₃) ≤ K(h) → 0;
#           при h ≥ H₀ (K(H₀) < λ₁ = min ĥ на E(ℚ)∖Tors): Q₂ ∓ Q₃ = t ∈ Tors — алгебраическое уравнение z₂(A) = z₃(A − t)
#           на E (случай B: конечно много A, корни многочлена); при h < H₀ — перебор точек ĥ(Q₁) < H₀ + C₁ (случай A).
from sage.all import *
import sys, json, time, itertools, random
sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable'); sys.path.insert(0, '/home/kep/magicKube/joint_sieve')
from classes import CLASSES, BIG, cells_of
from pairing_theory import Fac, setup_pair, eval_u, naive_h
from jsieve import slope_cells

SAFETY = 2.0      # запас на границу Силвермана (как в §7)
SAFETY_P = 2.0    # запас на численную константу спаривания


def all_square(z, S):
    z = QQ(z)
    return all((1 + lam * z) >= 0 and (1 + lam * z).is_square() for lam in S)


class Common:
    """общая кривая класса: минимальная модель E, изоморфизмы φ_T: E_T → E, ранг, образующие, кручение, константы."""
    def __init__(self, Ts, verbose=True):
        self.Ts = [sorted(ZZ(t) for t in T) for T in Ts]
        self.facs = [Fac(T) for T in self.Ts]
        self.E = self.facs[0].E.minimal_model()
        self.iso = []
        for f in self.facs:
            assert f.E.is_isomorphic(self.E), 'множитель не изоморфен общей кривой'
            phi = f.E.isomorphism_to(self.E); u, r, s_, t_ = phi.tuple()
            P0 = f.E(*f.XY(QQ(0), QQ(1)))                       # точка над z = 0
            assert phi(P0)[0] == (P0[0] - r) / u**2, 'соглашение Sage об (u,r,s,t) не то'
            self.iso.append(phi)
        rk = self.E.pari_curve().ellrank(); self.lo, self.hi = int(rk[0]), int(rk[1])
        pts = [self.E(list(p)) for p in rk[3]]
        if pts:
            sat, idx, reg = self.E.saturation(pts); self.gens, self.sat_index = list(sat), ZZ(idx)
        else:
            self.gens, self.sat_index = [], ZZ(1)
        self.rank = len(self.gens); self.rank_proved = (self.lo == self.hi == self.rank)
        self.tors = self.E.torsion_points()
        self.B = float(self.E.silverman_height_bound())
        # мёбиусова связь x_min ↔ z для каждого T:  x_min = (X − r)/u², X = L c₀/(c₀ z + 1)
        self.M = []; self.c1 = []
        for f, phi in zip(self.facs, self.iso):
            u, r, s_, t_ = phi.tuple(); a = 1 / u**2; b = -r / u**2      # x_min = a X + b
            M = matrix(QQ, [[b * f.c0, a * f.L * f.c0 + b], [f.c0, 1]])   # x_min = (M00 z + M01)/(M10 z + M11)
            d = lcm([x.denominator() for x in M.list()]); M = M * d
            g = gcd([ZZ(x) for x in M.list()]); M = M / g
            self.M.append(M); self.c1.append(float(log(2 * max(abs(x) for x in M.list()))))
        self.C = [SAFETY * self.B + c for c in self.c1]
        if self.rank >= 1:
            self.G = matrix(RDF, self.rank, self.rank, lambda i, j: float((self.gens[i] + self.gens[j]).height() - self.gens[i].height() - self.gens[j].height()) / 2)
            self.lam_min = min(self.G.eigenvalues()).real() * 0.999
        if verbose:
            print(f'  E = {self.E.ainvs()} (cond {self.E.conductor()}), ранг {self.lo}/{self.hi}, найдено {self.rank}, индекс насыщ. {self.sat_index}, '
                  f'кручение {len(self.tors)}, B = {self.B:.2f}, c1 = {[round(c,2) for c in self.c1]}', flush=True)

    def z_of_point(self, i, Q):
        """z для точки Q ∈ E через φ_i⁻¹: X = u² x + r, z = −1/c₀ + L/X."""
        if Q.is_zero(): return -1 / self.facs[i].c0
        u, r, s_, t_ = self.iso[i].tuple(); X = u**2 * Q[0] + r
        if X == 0: return infinity
        return -1 / self.facs[i].c0 + self.facs[i].L / X

    def xmin_of_z(self, i, z):
        M = self.M[i]; return (M[0, 0] * z + M[0, 1]) / (M[1, 0] * z + M[1, 1])


def pairing_constant(com, i, j, nmax=14, verbose=True):
    """численная оценка |⟨Q_i,Q_j⟩| ≤ c: c = (SAFETY_P·sup|Λ| + 2B)/4, Λ = h(u_+^{min}(R)) − h(u_−^{min}(R))."""
    Ti, Tj = com.Ts[i], com.Ts[j]
    P = setup_pair(Ti, Tj); D = P['D']
    psi = P['Es'].isomorphism_to(com.E); u, r, s_, t_ = psi.tuple()   # x_min = (x_s − r)/u²
    conv = lambda x: (x - r) / u**2
    lams = []
    if len(D) >= 3:
        f3 = P['f3']; E3 = f3.E
        rk = E3.pari_curve().ellrank(); gens = [E3(list(p)) for p in rk[3]]
        if gens: gens = list(E3.saturation(gens)[0])
        tors = E3.torsion_points()
        pts = []
        if len(gens) >= 1:
            for n in range(-nmax, nmax + 1):
                for t in tors: pts.append(n * gens[0] + t)
        if len(gens) >= 2:
            for n in range(-6, 7):
                for m in range(-6, 7):
                    for t in tors[:2]: pts.append(n * gens[0] + m * gens[1] + t)
        if not gens: pts = list(tors)
        for R in pts:
            if R.is_zero(): continue
            up = eval_u(P, R, +1); um = eval_u(P, R, -1)
            if up is None or um is None: continue
            up, um = conv(up), conv(um)
            if up == 0 or um == 0: continue
            lams.append(naive_h(up) - naive_h(um))
        genus = 1; npts = len(pts)
    else:
        b = abs(D[0]); random.seed(5)
        taus = [QQ(p) / q for p in range(-20, 21) for q in range(1, 21) if gcd(p, q) == 1] + \
               [QQ(random.randint(-10**e, 10**e)) / random.randint(1, 10**e) for e in range(2, 16) for _ in range(8)]
        for tau in taus:
            if len(D) == 2: zv = 2 * tau / (b * (1 + tau**2)); y3 = (1 - tau**2) / (1 + tau**2)
            else: zv = (tau**2 - 1) / D[0]; y3 = tau
            try:
                Av = P['A'].numerator()(zv) / P['A'].denominator()(zv); Bv = P['B'].numerator()(zv) / P['B'].denominator()(zv)
            except ZeroDivisionError:
                continue
            up, um = conv(Av + Bv * y3), conv(Av - Bv * y3)
            if up == 0 or um == 0: continue
            lams.append(naive_h(up) - naive_h(um))
        genus = 0; npts = len(taus)
    sup = max(abs(x) for x in lams) if lams else float('nan')
    c = (SAFETY_P * sup + 2 * com.B) / 4
    if verbose:
        print(f'  спаривание ({Ti},{Tj}): Δ = {D} (род {genus}), выборка {len(lams)}/{npts}, sup|Λ| = {sup:.2f} (численно), c = {c:.2f}', flush=True)
    return dict(i=i, j=j, D=[int(x) for x in D], genus=genus, nsample=len(lams), supLam=sup, c=c)


def case_A_points(com, Hmax):
    """все Q ∈ E(ℚ) с ĥ(Q) ≤ Hmax (по решётке образующих + кручение)."""
    r = com.rank
    if r == 0: return list(com.tors)
    bound = int(sqrt(Hmax / com.lam_min)) + 1
    out = []
    G = com.G
    if r == 1:
        for n in range(-bound, bound + 1):
            if n * n * G[0, 0] <= Hmax:
                Q = n * com.gens[0]
                for t in com.tors: out.append(Q + t)
    elif r == 2:
        for n in range(-bound, bound + 1):
            base = n * com.gens[0]
            for m in range(-bound, bound + 1):
                v = vector(RDF, [n, m])
                if v * G * v <= Hmax:
                    Q = base + m * com.gens[1]
                    for t in com.tors: out.append(Q + t)
    else:
        raise NotImplementedError
    return out


def case_B_candidates(com, i2, i3, verbose=True):
    """z с f₂(z) = ±f₃(z) + t (t ∈ Tors): уравнение z₂(A) = z₃(A − t) на E (короткая модель)."""
    Es = com.E.short_weierstrass_model(); psi = com.E.isomorphism_to(Es)
    A4, A6 = Es.a4(), Es.a6()
    # z_j(A) = M_j(x_min(A)), x_min = u²x_s + r (обратно к psi)
    u, rr, s_, t_ = psi.tuple()
    Rx = PolynomialRing(QQ, ['x', 'y']); x, y = Rx.gens()
    Rxx = PolynomialRing(QQ, 'X'); XX = Rxx.gen()
    f = x**3 + A4 * x + A6
    def mob(j, xs):   # z_j как дробь от x_s (полином/полином): z = M⁻¹(x_min), M⁻¹ = [[M11, −M01],[−M10, M00]]
        M = com.M[j]; xm = u**2 * xs + rr
        return M[1, 1] * xm - M[0, 1], -M[1, 0] * xm + M[0, 0]
    cands = set(); trivial = []
    for t in com.tors:
        ts = psi(t)
        if ts.is_zero():
            n2, d2 = mob(i2, x); n3, d3 = mob(i3, x)
            pol = n2 * d3 - n3 * d2
            polx = Rxx(pol.univariate_polynomial().list()) if pol != 0 else Rxx(0)
            if polx == 0: trivial.append('O'); continue
            roots = [rt for rt, m in polx.roots()]
            for xr in roots:
                cands.add(com.z_of_point(i2, com.E(0)) if False else None)
                fx = xr**3 + A4 * xr + A6
                if fx.is_square():
                    yr = fx.sqrt()
                    for yy in (yr, -yr):
                        A = psi.inverse()(Es(xr, yy)) if hasattr(psi, 'inverse') else (psi**-1)(Es(xr, yy))
                        cands.add(com.z_of_point(i2, A))
            continue
        xt, yt = ts[0], ts[1]
        Ppol = f + yt**2 - (x + xt) * (x - xt)**2           # x(A − t) = (Ppol + 2 yt y)/(x − xt)²
        n2, d2 = mob(i2, x)
        M = com.M[i3]
        # z₃(A−t) = (M11·xm′ − M01)/(−M10·xm′ + M00), xm′ = u² x′ + rr, x′ = (Ppol + 2yt y)/(x−xt)²
        num3 = M[1, 1] * (u**2 * (Ppol + 2 * yt * y) + rr * (x - xt)**2) - M[0, 1] * (x - xt)**2
        den3 = -M[1, 0] * (u**2 * (Ppol + 2 * yt * y) + rr * (x - xt)**2) + M[0, 0] * (x - xt)**2
        cond = n2 * den3 - num3 * d2
        Ac = cond.coefficient({y: 0}); Bc = cond.coefficient({y: 1}); assert cond == Ac + Bc * y
        conv = lambda p: Rxx(p.univariate_polynomial().list()) if p != 0 else Rxx(0)
        Ac, Bc = conv(Ac), conv(Bc)
        fX = XX**3 + A4 * XX + A6
        pol = Ac if Bc == 0 else Ac**2 - Bc**2 * fX
        if pol == 0: trivial.append(str(t)); continue
        for xr, m in pol.roots():
            fx = xr**3 + A4 * xr + A6
            if not fx.is_square(): continue
            yr = fx.sqrt()
            for yy in (yr, -yr):
                if Bc != 0 and Ac(xr) + Bc(xr) * yy != 0: continue
                A = (psi**-1)(Es(xr, yy)); cands.add(com.z_of_point(i2, A))
        # особые точки: A = t, A = −t, A = O
        for A in (t, -t, com.E(0)): cands.add(com.z_of_point(i2, A)); cands.add(com.z_of_point(i3, A))
    cands.discard(None)
    if verbose: print(f'  случай B (пара {com.Ts[i2]}, {com.Ts[i3]}): кандидатов {len(cands)}' + (f'; ТРИВИАЛЬНО при t = {trivial}' if trivial else ''), flush=True)
    return sorted(cands, key=str), trivial


def run_class(S, Ts, verbose=True):
    """S = −S — клетки; Ts — множители одного класса изоморфизма. Возвращает отчёт."""
    t0 = time.time()
    com = Common(Ts, verbose)
    rep = dict(Ts=[[int(t) for t in T] for T in com.Ts], E=list(map(int, com.E.ainvs())), cond=int(com.E.conductor()),
               rank=[com.lo, com.hi, com.rank], sat_index=int(com.sat_index), ntors=len(com.tors), B=com.B, c1=com.c1, ok=False)
    if not com.rank_proved or com.rank > 2:
        rep['why'] = f'ранг {com.lo}/{com.hi} (найдено {com.rank}) — метод не применим'; return rep
    m = len(Ts)
    # граф пар с третьей кривой рода ≤ 1
    dsz = lambda i, j: len(set(com.Ts[i]) ^ set(com.Ts[j]))
    nbrs = {i: [j for j in range(m) if j != i and dsz(i, j) <= 4] for i in range(m)}
    if com.rank == 0:
        cands = set(com.z_of_point(0, Q) for Q in com.tors); H0 = 0; pairs = []; casesB = []
    elif com.rank == 1:
        hub = max(range(m), key=lambda i: len(nbrs[i]))
        if not nbrs[hub]: rep['why'] = 'нет пары рода ≤ 1'; return rep
        j = min(nbrs[hub], key=lambda j: dsz(hub, j))
        pc = pairing_constant(com, hub, j, verbose=verbose); pairs = [pc]
        hG = float(com.gens[0].height()); c12 = pc['c']
        # |n₁n₂|ĥ(G) ≤ c12; либо n₂ = 0 ⇒ h ≤ C_j ⇒ n₁²ĥ(G) ≤ C_hub + C_j; либо |n₁| ≤ c12/ĥ(G)
        M0 = int(max(c12 / hG, sqrt((com.C[hub] + com.C[j]) / hG))) + 1
        H0 = M0 * M0 * hG
        cands = set(com.z_of_point(hub, Q) for Q in case_A_points(com, H0 + 1e-9))
        casesB = []
        if verbose: print(f'  ранг 1: ĥ(G) = {hG:.3f}, c12 = {c12:.2f}, M0 = {M0}, кандидатов A: {len(cands)}', flush=True)
    else:
        hub = max(range(m), key=lambda i: len(nbrs[i]))
        if len(nbrs[hub]) < 2: rep['why'] = 'нет хаба с двумя соседями рода ≤ 1'; return rep
        j2, j3 = sorted(nbrs[hub], key=lambda j: dsz(hub, j))[:2]
        p2 = pairing_constant(com, hub, j2, verbose=verbose); p3 = pairing_constant(com, hub, j3, verbose=verbose); pairs = [p2, p3]
        c2, c3 = p2['c'], p3['c']; C1, C2, C3 = com.C[hub], com.C[j2], com.C[j3]
        def K(h):
            if h <= max(C1, C2, C3) + 1: return float('inf')
            a2 = c2 / sqrt(h - C1); a3 = c3 / sqrt(h - C1); a = max(a2, a3)
            if h - max(C2, C3) - a * a <= 0: return float('inf')
            return (a2 + a3)**2 + ((C2 + C3 + a * a) / (2 * sqrt(h - max(C2, C3) - a * a)))**2
        lam1 = com.lam_min
        H0 = max(C1, C2, C3) + 1
        while K(H0) >= lam1: H0 *= 1.05
        H0 = float(H0)
        if verbose: print(f'  ранг 2: λ_min(Gram) = {lam1:.3f}, c = ({c2:.2f}, {c3:.2f}), H0 = {H0:.1f} (K(H0) = {K(H0):.3f})', flush=True)
        ptsA = case_A_points(com, H0 + C1)
        cands = set(com.z_of_point(hub, Q) for Q in ptsA)
        if verbose: print(f'  случай A: точек с ĥ ≤ {H0 + C1:.1f}: {len(ptsA)}', flush=True)
        zB, trivial = case_B_candidates(com, j2, j3, verbose=verbose); casesB = [dict(pair=[j2, j3], n=len(zB), trivial=trivial)]
        cands |= set(zB)
        if trivial: rep['why'] = 'тривиальное уравнение случая B'; rep['trivial'] = trivial; return rep
    cands.discard(infinity)
    sols = [z for z in cands if z != 0 and all_square(z, S)]
    rep.update(ok=True, H0=float(H0), pairs=pairs, casesB=casesB, ncand=len(cands), nondeg_solutions=[str(z) for z in sols], time=time.time() - t0)
    if verbose: print(f'  ИТОГ класса: кандидатов {len(cands)}, невырожденных решений: {sols}; {time.time()-t0:.1f} с', flush=True)
    return rep


def run_slope(r, s, class_indices=None, max_rank=2, verbose=True):
    S = slope_cells(r, s)
    out = []
    for ci in (class_indices if class_indices is not None else range(len(BIG))):
        Ts = [cells_of(Sg, r, s) for Sg in CLASSES[BIG[ci]]]
        if verbose: print(f'класс {ci} (m = {len(Ts)}): {Ts}', flush=True)
        try:
            rep = run_class(S, Ts, verbose)
        except Exception as e:
            rep = dict(ok=False, why='ошибка: ' + repr(e)[:200]); print('  ', rep['why'], flush=True)
        rep['class'] = ci; out.append(rep)
        if rep.get('ok') and not rep['nondeg_solutions']: break
    return out


if __name__ == '__main__':
    r, s = int(sys.argv[1]), int(sys.argv[2])
    cis = [int(x) for x in sys.argv[3].split(',')] if len(sys.argv) > 3 else None
    out = run_slope(r, s, cis)
    closed = any(o.get('ok') and not o['nondeg_solutions'] for o in out)
    print(f'ИТОГ {r}/{s}: ЗАКРЫТ обобщённым методом: {closed}', flush=True)
    json.dump(dict(slope=f'{r}/{s}', results=out, closed=closed), open(f'/home/kep/magicKube/isogeny_mult/fable/dem2_{r}_{s}.json', 'w'), default=str)
