# Спаривание Касселса–Тейта на Sel^2(E/k) для E: y^2 = (x-e1)(x-e2)(x-e3) с полным 2-кручением над k (Claude, 2026-09-11).
#
# Формула — эллиптический аналог Теоремы 3.3 из J. Yan, arXiv:2109.08258 (там род 2, полное 2-кручение; для эллиптических
# кривых — Cassels 1998): для 2-накрытия C_eps и функций f_T на C_eps с div f_T = 2 D_T - 2 D
#     <eps, eta> = prod_v (f_T1(P_v), beta_1)_v (f_T2(P_v), beta_2)_v,   eta: sigma -> beta~_1 T1 + beta~_2 T2.
# Координаты Куммера eta = (x - e1, x - e2) — спаривание Вейля с T1, T2; e2(T_i,T_i) = 1, e2(T1,T2) = -1 =>
#     beta_1 = eta_2, beta_2 = eta_1.
# Модель накрытия: x - e_i = eps_i u_i^2 (eps3 = eps1 eps2), в P^3 (u0:u1:u2:u3):
#     G1 = eps1 u1^2 - eps2 u2^2 - (e2-e1) u0^2,  G2 = eps1 u1^2 - eps3 u3^2 - (e3-e1) u0^2.
# Четыре конуса пучка; конус без u_j <-> инволюция u_j -> -u_j <-> (для j=0) [-1], (для j=i) [-1] o tau_{T_i}.
# Их g^1_2 — классы D (j = 0) и D_{T_i} (j = i). Касательная плоскость L_j вдоль рациональной образующей конуса j:
#     div(L_j) = 2 D_j  =>  f_{T_i} = L_i / L_0.
# Места: S (плохие) + простые из eps, eta, коэффициентов L_j + все простые нормы <= 60 + вещественные. Для остальных v
# (хорошая редукция, всё — единицы, поле вычетов > 60) найдётся P_v с единичным f(P_v), символ = 1.
import random, functools
print = functools.partial(print, flush=True)

def _int_basis_coords(k, a):
    B = k.ring_of_integers().basis()
    M = matrix(QQ, [list(k(b)) for b in B])
    return vector(QQ, list(k(a))) * M.inverse()

def reduce_mod_pN(k, a, p, N):
    """ элемент k, p-адически близкий к a (разность в p^N O_k); a без p в знаменателе после выделения p^v """
    B = k.ring_of_integers().basis()
    c = _int_basis_coords(k, a)
    D = lcm([ci.denominator() for ci in c])
    v = D.valuation(p); Dp = D // p^v
    inv = inverse_mod(Dp, p^(N + v))
    num = [ZZ(ci * D) * inv % p^(N + v) for ci in c]
    return sum(QQ(ni) / p^v * k(b) for ni, b in zip(num, B))

class LocalSqrt:
    def __init__(self, k, P):
        self.k, self.P = k, P; self.p = P.smallest_integer()
        self.e = P.ramification_index(); self.pi = k.uniformizer(P, others='positive')
        self.F = P.residue_field()
    def sqrt(self, a, N=60):
        """ s в k, s^2 = a с относительной точностью P^N (a — локальный квадрат) """
        k, P, p = self.k, self.P, self.p
        a = k(a); v = a.valuation(P); assert v % 2 == 0
        w = a / self.pi^v
        if p != 2:
            s = self.F.lift(self.F(w).sqrt())
        else:
            M = 2*self.e + 3
            s = None
            # перебор O/P^M: c0 + c1*omega, c_i mod 2^M покрывают O/P^M (P^M ⊇ 2^M O)
            B = k.ring_of_integers().basis()
            rng = range(0, 2^(M))
            for c0 in rng:
                for c1 in rng:
                    t = c0*k(B[0]) + c1*k(B[1])
                    if t.valuation(P) == 0 and (t^2 - w).valuation(P) >= M:
                        s = t; break
                if s is not None: break
            assert s is not None
        it = 0
        while (s^2 - w).valuation(P) < N and it < 40:     # Ньютон; при p = 2 точность растёт как m -> 2m - 2e
            s = (s + w/s)/2
            s = reduce_mod_pN(k, s, p, N + 10)
            it += 1
        assert (s^2 - w).valuation(P) >= N, ((s^2 - w).valuation(P), N)
        e2 = k(2).valuation(P)
        assert N > 2*e2
        # Гензель: существует точный корень alpha из w с v(s - alpha) >= N - e2; корень из a — pi^(v/2) alpha
        self.last_prec = N - e2 + v//2
        return s * self.pi^(v//2)

class CTP:
    def __init__(self, Cv, extra_norm=60):
        self.Cv = Cv; self.k = Cv.k; self.e = Cv.e; self.cert = None
        self.extra_norm = extra_norm
        self.R4 = PolynomialRing(self.k, 'u0,u1,u2,u3'); self.u = self.R4.gens()

    def cones(self, eps1, eps2):
        k = self.k; e1, e2, e3 = self.e; u0, u1, u2, u3 = self.u
        eps3 = eps1*eps2
        G1 = eps1*u1^2 - eps2*u2^2 - (e2 - e1)*u0^2
        G2 = eps1*u1^2 - eps3*u3^2 - (e3 - e1)*u0^2
        return {3: G1, 2: G2, 1: G1 - G2, 0: (e3 - e1)*G1 - (e2 - e1)*G2}, eps3

    def small_rep(self, a):
        """ малый представитель класса a в k(S,2) и s с a = rep * s^2 """
        if not hasattr(self, '_ks'):
            self._ks = self.k.selmer_space(self.Cv.S, 2)
        KSp, gens, fromV, toV = self._ks
        a = self.k(a)
        rep = self.k(fromV(toV(a)))
        s2 = a / rep
        ok = s2.is_square()
        assert ok, "class map not exact"
        return rep, s2.sqrt()

    def solve_conic(self, coeffs):
        """ точка на a X^2 + b Y^2 + c Z^2 = 0 над k (проверяется подстановкой) """
        k = self.k
        reps = [self.small_rep(c) for c in coeffs]
        a, b, c = [r_ for r_, s_ in reps]
        order = []
        pt = None
        for (i, j, l) in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
            A, B, Cc = [a, b, c][i], [a, b, c][j], [a, b, c][l]
            t = -B/A
            if t.is_square():                      # вырожденный случай: точка сразу
                st = t.sqrt(); pt = [None]*3; pt[i] = st; pt[j] = k(1); pt[l] = k(0); break
            order.append((abs(t.norm().numerator()) * abs(t.norm().denominator()), i, j, l))
        if pt is None:
            errs = []
            for cost, i, j, l in sorted(order):    # все три расстановки по возрастанию сложности
                A, B, Cc = [a, b, c][i], [a, b, c][j], [a, b, c][l]
                try:
                    try:
                        pari.allocatemem(2*10^9, 8*10^9, silent=True)
                    except Exception:
                        pass
                    t = -B/A; m = -Cc/A
                    D = lcm([QQ(c_).denominator() for c_ in list(t)]); tp = t*D^2
                    Rx_ = PolynomialRing(k, 'Xn'); Xn = Rx_.gen()
                    Lrel = k.extension(Xn^2 - tp, 'yy')
                    ok, el = m.is_norm(Lrel, element=True, proof=False)
                    assert ok, "not a norm (conic without point?)"
                    cs_ = el.list()
                    u_ = k(cs_[0]); v_ = k(cs_[1]) * D
                    pt = [None]*3; pt[i] = u_; pt[j] = v_; pt[l] = k(1)
                    break
                except (RuntimeError, MemoryError, Exception) as ex:
                    errs.append(f"{(i,j,l)}: {type(ex).__name__}: {str(ex)[:60]}")
            if pt is None: raise RuntimeError("conic unsolved: " + "; ".join(errs))
        assert a*pt[0]^2 + b*pt[1]^2 + c*pt[2]^2 == 0
        if getattr(self, 'randomize_conic', False):      # другая точка той же коники: вторая точка пересечения с прямой
            for _ in range(20):
                v = [k(random.randint(-9, 9)) for _ in range(3)]
                Qv = a*v[0]^2 + b*v[1]^2 + c*v[2]^2
                if Qv == 0: continue
                lam = -2*(a*pt[0]*v[0] + b*pt[1]*v[1] + c*pt[2]*v[2]) / Qv
                if lam == 0: continue
                pt = [pt[t] + lam*v[t] for t in range(3)]; break
            assert a*pt[0]^2 + b*pt[1]^2 + c*pt[2]^2 == 0
        # к исходным коэффициентам: coeff = rep * s^2  =>  X_orig = X / s
        return [pt[t] / reps[t][1] for t in range(3)]

    def tangent_forms(self, eps1, eps2):
        k = self.k; u = self.u
        cones, eps3 = self.cones(eps1, eps2)
        L = {}; self.conic_points = {}; self.conic_coeffs = {}
        for j, Q in cones.items():
            vars3 = [t for t in range(4) if t != j]
            coeffs = [Q.monomial_coefficient(u[t]^2) for t in vars3]
            assert all(Q.monomial_coefficient(u[a]*u[b]) == 0 for a in range(4) for b in range(a+1, 4))
            assert Q.monomial_coefficient(u[j]^2) == 0
            pt = self.solve_conic(coeffs)
            assert sum(coeffs[t]*pt[t]^2 for t in range(3)) == 0
            self.conic_points[str(j)] = [str(v) for v in pt]; self.conic_coeffs[str(j)] = [str(v) for v in coeffs]
            L[j] = sum(2*coeffs[t]*pt[t]*u[vars3[t]] for t in range(3))
        return L, eps3

    def normalize_form(self, Lj):
        """ L_j / c, где идеал коэффициентов после деления — произведение простых из фиксированного набора мест
            (S и все простые нормы <= extra_norm; они порождают группу классов, т.к. граница Минковского меньше) """
        k = self.k
        cs = [Lj.monomial_coefficient(t) for t in self.u]
        I = k.ideal([c for c in cs if c != 0])
        if not hasattr(self, '_clgens'):
            Cl = k.class_group(); self._Cl = Cl
            base = list(self.Cv.S) + [P for q in primes(2, self.extra_norm + 1) for P in k.primes_above(q) if P.norm() <= self.extra_norm]
            # BFS по классам: для каждого класса — представитель-произведение простых из base
            reps = {Cl(k.ideal(1)): k.ideal(1)}
            frontier = [k.ideal(1)]
            while frontier and len(reps) < Cl.order():
                nf_ = []
                for J in frontier:
                    for P in base:
                        JP = J*P; c_ = Cl(JP)
                        if c_ not in reps: reps[c_] = JP; nf_.append(JP)
                frontier = nf_
            assert len(reps) == Cl.order(), "places do not generate the class group"
            self._clgens = reps
        B = self._clgens[self._Cl(I)]
        J = I / B
        assert J.is_principal()
        c = J.gens_reduced()[0]
        if not hasattr(self, 'norm_consts'): self.norm_consts = []
        self.norm_consts.append(str(c))
        return Lj / c

    def places_for(self, eps_reps, eta_reps):
        """ S + простые, делящие малые представители eps, eta (их нормы невелики) + все простые нормы <= extra_norm """
        k = self.k
        Ps = list(self.Cv.S)
        def add(P):
            if not any(P == Q for Q in Ps): Ps.append(P)
        for a in list(eps_reps) + [x for eta in eta_reps for x in eta]:
            for P, e_ in k.ideal(a).factor(): add(P)
        for q in primes(2, self.extra_norm + 1):
            for P in k.primes_above(q):
                if P.norm() <= self.extra_norm: add(P)
        return Ps

    def local_point_finite(self, P, eps, tries=20000, N=60, want_prec=False):
        k = self.k; e = self.e
        ls = LocalSqrt(k, P)
        loc = LocalSq(k, P)
        for it in range(tries):
            j = random.randint(-6, 8)
            x = random.choice(list(e) + [0]) + rand_elt(k, P.smallest_integer()^4) * ls.pi^j
            vals = [(x - e[i]) / eps[i] for i in range(3)]
            if any(v_ == 0 for v_ in vals): continue
            if all(loc.is_sq(v_) for v_ in vals):
                us = []; precs = []
                for v_ in vals:
                    us.append(ls.sqrt(v_, N=N)); precs.append(ls.last_prec)
                self._last_x = x
                if want_prec: return [k(1)] + us, [Infinity] + precs
                return [k(1)] + us
        raise RuntimeError(f"no local point at {P}")

    def local_point_real(self, emb, eps, tries=20000):
        """ точный x в k с (x - e_i)/eps_i > 0 (знаки — точно, exact_sign) """
        e = self.e; k = self.k
        for it in range(tries):
            c = random.choice(list(e))
            x = c + QQ(random.randint(-10^6, 10^6)) / random.randint(1, 10^3) * random.choice([1, QQ(1)/10^3, 10^3, 10^6])
            vals = [(x - e[i]) / eps[i] for i in range(3)]
            if all(v_ != 0 and exact_sign(k, v_, emb) > 0 for v_ in vals):
                return x, vals
        raise RuntimeError("no real point")

    def real_values_certified(self, emb, vals, L):
        """ интервальные значения L_j в точке (1, sqrt(vals)) при вещественном вложении; точность растёт, пока
            интервалы L_0, L_1, L_2 не отделены от нуля. Возвращает знаки. """
        k = self.k; r = k.gen(); D = QQ(r^2); sr = 1 if emb(r) > 0 else -1
        for prec in [256, 1024, 4096, 16384, 65536]:
            RI = RealIntervalField(prec); sD = RI(D).sqrt() * sr
            def ev(a):
                a0, a1 = list(k(a)); return RI(a0) + RI(a1) * sD
            us = [RI(1)] + [ev(v_).sqrt() for v_ in vals]
            Lv = {j: sum(ev(L[j].monomial_coefficient(self.u[i])) * us[i] for i in range(4)) for j in (0, 1, 2, 3)}
            if all(not Lv[j].contains_zero() for j in (0, 1, 2, 3)):
                return {j: (1 if Lv[j] > 0 else -1) for j in (0, 1, 2, 3)}, prec
        raise RuntimeError("real sign not certified")

    def pair(self, eps, eta_list, reps=2, verbose=False, N=80):
        """ eps = (eps1, eps2) из Sel^2; eta_list — список пар (eta1, eta2). Возвращает список значений в {0,1} (1 = -1).
            Каждое значение считается reps раз с разными локальными точками — должно совпадать. """
        k = self.k
        eps1 = self.small_rep(eps[0])[0]; eps2 = self.small_rep(eps[1])[0]
        eta_reps = [(self.small_rep(h1)[0], self.small_rep(h2)[0]) for (h1, h2) in eta_list]
        L, eps3 = self.tangent_forms(eps1, eps2)
        L = {j: self.normalize_form(L[j]) for j in L}
        epsv = [eps1, eps2, eps3]
        Ps = self.places_for([eps1, eps2], eta_reps)
        self.witness_head = {'eps_reps': [str(eps1), str(eps2), str(eps3)],
                             'eta_reps': [[str(a_), str(b_)] for a_, b_ in eta_reps],
                             'L_forms': {str(j): {str(i): str(L[j].monomial_coefficient(self.u[i])) for i in range(4)} for j in L},
                             'conic_coeffs': self.conic_coeffs, 'conic_points': self.conic_points,
                             'norm_consts': list(getattr(self, 'norm_consts', [])),
                             'places': [str(P) for P in Ps]}
        self.norm_consts = []
        RF = RealField(4000)
        results = []
        for rep in range(reps):
            tot = [0]*len(eta_list)
            for P in Ps:
                e2 = k(2).valuation(P)
                Nloc = N
                while True:
                    for _try in range(50):
                        pt, precs = self.local_point_finite(P, epsv, N=Nloc, want_prec=True)
                        Lv = {j: L[j](pt) for j in L}
                        if all(Lv[j] != 0 for j in (0, 1, 2, 3)): break
                    else:
                        raise RuntimeError("local point on L=0")
                    # доказанная точность A_j значения L_j и запас до порога квадратности
                    margins = {}
                    for j in (0, 1, 2, 3):
                        Aj = min(L[j].monomial_coefficient(self.u[i]).valuation(P) + precs[i]
                                 for i in range(1, 4) if L[j].monomial_coefficient(self.u[i]) != 0)
                        margins[j] = Aj - Lv[j].valuation(P) - 2*e2
                    if all(m_ > 0 for m_ in margins.values()): break
                    Nloc *= 2
                    if Nloc > 4000: raise RuntimeError(f"precision not certified at {P}")
                f1, f2 = Lv[1]/Lv[0], Lv[2]/Lv[0]
                syms = []
                for t, (h1, h2) in enumerate(eta_reps):
                    if getattr(self, 'cassels_form', False):
                        # форма Касселса [Cas98, Lemma 7.4]: три конуса без нормировки и без перестановки,
                        # prod_{i=1..3} (L_i(P_v), eta_i)_v,  eta_3 = eta_1 eta_2
                        h3 = k(h1)*k(h2)
                        s1 = k.hilbert_symbol(Lv[1], k(h1), P) * k.hilbert_symbol(Lv[2], k(h2), P)
                        s2 = k.hilbert_symbol(Lv[3], h3, P)
                    else:
                        s1 = k.hilbert_symbol(f1, h2, P); s2 = k.hilbert_symbol(f2, h1, P)
                    tot[t] += 0 if s1*s2 == 1 else 1
                    syms.append((int(s1), int(s2)))
                if self.cert is not None:
                    self.cert.append({'place': str(P), 'gens_of_place': [str(g) for g in P.gens()],
                                      'x': str(self._last_x), 'N': int(Nloc),
                                      'u_approx': [str(v) for v in pt], 'A_i': [str(a_) for a_ in precs],
                                      'L_values': {str(j): str(Lv[j]) for j in (0, 1, 2, 3)},
                                      'v_L': {str(j): int(Lv[j].valuation(P)) for j in (0, 1, 2, 3)},
                                      'f1': str(f1), 'f2': str(f2),
                                      'margins': {str(j): str(m_) for j, m_ in margins.items()}, 'symbols': syms})
            for emb in k.real_embeddings():
                x, vals = self.local_point_real(emb, epsv)
                sg, prec = self.real_values_certified(emb, vals, L)
                f1s, f2s = sg[1]*sg[0], sg[2]*sg[0]
                syms = []
                for t, (h1, h2) in enumerate(eta_reps):
                    if getattr(self, 'cassels_form', False):
                        h3 = k(h1)*k(h2)
                        s1 = (-1 if (sg[1] < 0 and exact_sign(k, h1, emb) < 0) else 1) * (-1 if (sg[2] < 0 and exact_sign(k, h2, emb) < 0) else 1)
                        s2 = -1 if (sg[3] < 0 and exact_sign(k, h3, emb) < 0) else 1
                    else:
                        s1 = -1 if (f1s < 0 and exact_sign(k, h2, emb) < 0) else 1
                        s2 = -1 if (f2s < 0 and exact_sign(k, h1, emb) < 0) else 1
                    tot[t] += 0 if s1*s2 == 1 else 1
                    syms.append((s1, s2))
                if self.cert is not None:
                    self.cert.append({'place': f'real sgn(r)={1 if emb(k.gen()) > 0 else -1}', 'x': str(x), 'interval_prec': prec,
                                      'vals_under_sqrt': [str(v_) for v_ in vals],
                                      'signs_L': {str(j): int(v) for j, v in sg.items()}, 'symbols': syms})
            results.append([x % 2 for x in tot])
        assert all(r == results[0] for r in results), f"CTP depends on local points: {results}"
        return results[0], len(Ps)
