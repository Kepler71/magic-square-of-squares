# Спуск через Richelot-изогению phi: J = Jac(C) -> Jhat = Jac(Chat), C: y^2 = F1 F2 F3, Chat: Delta y^2 = G1 G2 G3.
# mu_F: J(Q) -> {(a1,a2,a3) in (Q*/Q*^2)^3: a1 a2 a3 = 1},  D = sum P_i  |->  (prod F_i(P_j))_i   (для дивизоров степени 0)
# mu_G: Jhat(Q) -> то же с G_i.  Точная последовательность (ker phi, ker phihat рациональны, порядка 4):
#   2^r = |Im mu_F| * |Im mu_G| / 16;  локально |Im mu_F,v| * |Im mu_G,v| = 16 (p нечётно), 64 (p = 2), 4 (v = oo).
# Верхняя граница: r <= dim Sel_F + dim Sel_G - 4.
import functools, random
print = functools.partial(print, flush=True)
R = PolynomialRing(QQ, 't'); t = R.gen()

def hom_eval(Fq, xz):
    # квадратичная форма из многочлена степени <= 2 (степень 1 => корень в бесконечности)
    x0, z0 = xz
    c0, c1, c2 = Fq[0], Fq[1], Fq[2]
    return c2*x0^2 + c1*x0*z0 + c0*z0^2

class RichelotPair:
    def __init__(self, Fs, S_extra=()):
        self.F = [R(f) for f in Fs]
        co = lambda P: [P[0], P[1], P[2]]
        self.delta = matrix(QQ, [co(f) for f in self.F]).det(); assert self.delta != 0
        F1, F2, F3 = self.F
        self.G = [F2.derivative()*F3 - F2*F3.derivative(), F3.derivative()*F1 - F3*F1.derivative(), F1.derivative()*F2 - F1*F2.derivative()]
        self.fC = prod(self.F); self.fH = self.delta * prod(self.G)
        bad = set(ZZ(2).prime_factors())
        for q in [self.fC.discriminant(), self.fH.discriminant(), self.delta] + [f.leading_coefficient() for f in self.F + self.G if f != 0]:
            q = QQ(q)
            if q != 0: bad |= set(q.numerator().prime_factors()) | set(q.denominator().prime_factors())
        self.S = sorted(bad | set(S_extra))
    # --- значения троек на точках
    def triple_point(self, Q, xz, curve):
        # точка (x:z) на кривой; значения квадр. форм; если одна = 0, заменяем произведением двух других
        vals = [hom_eval(q, xz) for q in Q]
        if vals.count(0) > 1: return None
        for i in range(3):
            if vals[i] == 0:
                vals[i] = vals[(i+1) % 3] * vals[(i+2) % 3]
        return vals
    def triple_pair(self, Q, a):
        # пара сопряжённых точек с x-координатами — корни квадратичного a (над Q), значения = нормы
        M = NumberField(a, 'z'); z = M.gen()
        vals = []
        for q in Q:
            v = q(z)
            if v == 0: return None
            vals.append(QQ(v.norm()))
        return vals
    # --- координаты квадратичных классов
    def qclass_local(self, a, p):
        a = QQ(a)
        if p == oo: return [0 if a > 0 else 1]
        v = a.valuation(p); u = a / p^v
        n = ZZ(u.numerator()*u.denominator())
        if p == 2:
            return [v % 2, 0 if n % 4 == 1 else 1, 0 if n % 8 in (1, 7) else 1]   # (-1, 2)-базис по модулю квадратов: u mod 8
        return [v % 2, 0 if kronecker(n, p) == 1 else 1]
    def vec_local(self, vals, p):
        # (a1, a2) — a3 определяется условием произведения
        return vector(GF(2), self.qclass_local(vals[0], p) + self.qclass_local(vals[1], p))
    def vec_global(self, vals):
        out = []
        for a in vals[:2]:
            a = QQ(a); out += [0 if a > 0 else 1] + [a.valuation(p) % 2 for p in self.S]
        return vector(GF(2), out)
    # --- проверка «лежит на кривой» локально
    def is_sq(self, a, p):
        a = QQ(a)
        if a == 0: return True
        if p == oo: return a > 0
        v = a.valuation(p)
        if v % 2: return False
        u = a/p^v; n = ZZ(u.numerator()*u.denominator())
        return (n % 8 == 1) if p == 2 else (kronecker(n, p) == 1)
    def is_sq_ext(self, beta, M, p):
        # beta in M квадрат в M ⊗ Q_p
        if beta == 0: return False
        if p == oo:
            return all(e(beta) > 0 for e in M.real_embeddings())
        for P in M.primes_above(p):
            v = beta.valuation(P)
            if v % 2: return False
            pi = M.uniformizer(P, others='positive'); u = beta/pi^v
            e = P.ramification_index() if p == 2 else 0; mod = P^(2*e + 1)
            G = mod.idealstar(2); lg = mod.ideallog(u)
            if any(lg[i] % 2 for i, o in enumerate(G.gens_orders()) if o % 2 == 0): return False
        return True
    def sample(self, side, p, target_total, other_dim, tries=40000, seed=7):
        # side = 'F' (кривая C, f = fC) или 'G' (Chat, f = fH)
        Q = self.F if side == 'F' else self.G
        f = self.fC if side == 'F' else self.fH
        random.seed(int(int(seed) + (0 if p == oo else int(p)) + (0 if side == "F" else 1000)))
        rroots = f.roots(QQ, multiplicities=False)
        vecs = []; base = None; rk = 0
        def push(vals):
            nonlocal base, rk
            v = self.vec_local(vals, p)
            if base is None: base = v; return
            vecs.append(v - base)
            rk = matrix(GF(2), vecs).rank()
        for it in range(tries):
            r = random.random()
            if p == oo:
                if r < 0.5:
                    x0 = QQ(random.randint(-10^6, 10^6))/random.randint(1, 10^3)
                    if not self.is_sq(f(x0), oo): continue
                    vals = self.triple_point(Q, (x0, 1), side)
                else:
                    re_ = QQ(random.randint(-10^4, 10^4))/random.randint(1, 100); im2 = QQ(random.randint(1, 10^6))/random.randint(1, 100)
                    a = (t - re_)^2 + im2
                    if f.gcd(a) != 1: continue
                    vals = self.triple_pair(Q, a)
                    if vals is None: continue
                    # пара над R: f(z) должно быть квадратом в C — всегда; но нужна рациональность дивизора — ок
            else:
                if r < 0.4:
                    k = random.randint(-4, 4); u0 = QQ(random.randint(-p^6, p^6))
                    x0 = (random.choice(rroots) + u0*p^abs(k)) if (rroots and random.random() < 0.4) else u0*p^k
                    if f(x0) == 0 or not self.is_sq(f(x0), p): continue
                    vals = self.triple_point(Q, (x0, 1), side)
                elif r < 0.45:
                    # точки у бесконечности: (1 : z0) с малым z0
                    z0 = QQ(random.randint(-p^4, p^4)) * p^random.randint(0, 4)
                    val = sum(f[i] * z0^(6 - i) for i in range(7))   # однородная f(1, z0) степени 6
                    if val == 0 or not self.is_sq(val, p): continue
                    vals = self.triple_point(Q, (1, z0), side)
                else:
                    if random.random() < 0.5:
                        b = QQ(random.randint(-p^5, p^5))*p^random.randint(-3, 3); c = QQ(random.randint(-p^5, p^5))*p^random.randint(-4, 4)
                        a = t^2 + b*t + c
                    else:
                        x1 = QQ(random.randint(-p^5, p^5))*p^random.randint(-3, 3); x2 = QQ(random.randint(-p^5, p^5))*p^random.randint(-3, 3)
                        a = (t - x1)*(t - x2)
                    if a.discriminant() == 0 or f.gcd(a) != 1: continue
                    if a.is_irreducible():
                        M = NumberField(a, 'z'); z = M.gen()
                        if not self.is_sq_ext(M(f(z)), M, p): continue
                        vals = self.triple_pair(Q, a)
                    else:
                        rts = a.roots(QQ, multiplicities=False)
                        if not all(self.is_sq(f(r0), p) for r0 in rts): continue
                        v1 = self.triple_point(Q, (rts[0], 1), side); v2 = self.triple_point(Q, (rts[1], 1), side)
                        if v1 is None or v2 is None: continue
                        vals = [v1[i]*v2[i] for i in range(3)]
            if vals is None: continue
            push(vals)
            if rk + other_dim >= target_total: break
        n = 2*len(self.qclass_local(QQ(3), p))
        M = matrix(GF(2), vecs) if vecs else matrix(GF(2), 0, n)
        return M.row_space()
    def run(self, verbose=True):
        glob_rows = {}
        # глобальная группа: (a1, a2) in Q(S,2)^2, базис: {-1} ∪ S для каждой компоненты
        gens_q = [QQ(-1)] + [QQ(p) for p in self.S]
        nG = 2*len(gens_q)
        def gvec_local(idx, p):
            # локальный вектор базисного элемента глобальной группы
            comp, k = divmod(idx, len(gens_q))
            a = gens_q[k]
            vals = [QQ(1), QQ(1)]; vals[comp] = a
            return self.vec_local(vals + [QQ(1)], p)
        condF, condG = [], []
        for p in self.S + [oo]:
            total = 4 if (p != oo and p != 2) else (6 if p == 2 else 2)
            WF = self.sample('F', p, total, 0)
            WG = self.sample('G', p, total, WF.dimension())
            # добираем F, если нужно
            if WF.dimension() + WG.dimension() < total:
                WF = self.sample('F', p, total, WG.dimension(), seed=99)
            dF, dG = WF.dimension(), WG.dimension()
            if verbose: print(f"  v = {p}: dim Im_F = {dF}, dim Im_G = {dG}, sum = {dF + dG} (expected {total})")
            assert dF + dG == total, f"local images incomplete at {p}"
            L = matrix(GF(2), [gvec_local(i, p) for i in range(nG)])
            for W, cond in ((WF, condF), (WG, condG)):
                Wc = W.basis_matrix().right_kernel_matrix() if W.dimension() > 0 else identity_matrix(GF(2), L.ncols())
                if Wc.nrows() > 0: cond.append((L * Wc.transpose()).transpose())
        selF = block_matrix([[c] for c in condF]).right_kernel()
        selG = block_matrix([[c] for c in condG]).right_kernel()
        bound = selF.dimension() + selG.dimension() - 4
        if verbose: print(f"  dim Sel_F = {selF.dimension()}, dim Sel_G = {selG.dimension()}  =>  rank <= {bound}")
        return bound
