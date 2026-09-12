# Полный 2-спуск для якобиана кривой рода 2 y^2 = f(x), deg f = 5 (модель нечётной степени), над Q.
# Cassels–Flynn гл. 6 / Stoll: mu: J(Q)/2J(Q) -> H = ker(N: L*/L*^2 -> Q*/Q*^2), L = Q[x]/(F), F — унитарная модель,
# mu(P - oo) = x(P) - theta, mu(P1 + P2 - 2oo) = a(theta) (a — унитарный многочлен от x(P1), x(P2)).
# Локально (odd degree): mu_v инъективно на J(Q_v)/2J(Q_v), |образ| = 2^(m_v - 1) * (4 если v = 2) (p-адические),
# 2^(m_R - 1)/4 для v = oo; m_v — число неприводимых множителей F над Q_v.
# rank J(Q) <= dim Sel - (m - 1).
import functools, random, itertools
print = functools.partial(print, flush=True)
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()

def monic_model(f):
    c = f.leading_coefficient()
    F = Rx(c^4 * f(x/c))               # Y = c^2 y
    assert F.is_monic()
    den = lcm([co.denominator() for co in F.coefficients()])
    # сделать коэффициенты целыми: X -> X/d, Y -> Y/d^(5/2)... проще: подстановка X = x/d, умножение на d^5 (нечётная степень -> твист на d)
    d = 1
    while not all(co in ZZ for co in Rx(d^5 * F(x/d)).coefficients()):
        d += 1
    F2 = Rx(d^5 * F(x/d))
    if d != 1:
        # d^5 F(x/d) = (Y d^{5/2})^2 — если d не квадрат, это твист на d; компенсируем: берём d квадратом
        dd = d^2
        F2 = Rx(dd^5 * F(x/dd)); assert F2.is_monic() and all(co in ZZ for co in F2.coefficients())
    return F2

class Descent:
    def __init__(self, f, extra_primes=()):
        self.F = monic_model(f)
        F = self.F
        assert F.degree() == 5 and F.is_squarefree()
        self.facs = [g for g, e in F.factor()]
        self.m = len(self.facs)
        self.fields = []
        for g in self.facs:
            K = NumberField(g, 'th%d' % len(self.fields))
            self.fields.append((K, K.gen()))
        self.rroots = F.roots(QQ, multiplicities=False)
        bad = set(ZZ(2*F.discriminant()).prime_factors()) | set(extra_primes)
        self.S = sorted(bad)
    # --- локальные координаты квадратичного класса
    def primes_above(self, p):
        out = []
        for j, (K, th) in enumerate(self.fields):
            for P in K.primes_above(p):
                e = P.ramification_index() if p == 2 else 0
                pi = K.uniformizer(P, others='positive')
                mod = P^(2*e + 1)            # p нечётно: e = 0 -> mod = P
                G = mod.idealstar(2)
                inv = G.gens_orders()
                even = [i for i, o in enumerate(inv) if o % 2 == 0]
                out.append((j, P, pi, mod, even))
        return out
    def local_coords(self, alpha, p, pl):
        # alpha — кортеж компонент (по полям); pl — список мест над p
        vec = []
        for (j, P, pi, mod, even) in pl:
            a = alpha[j]
            v = a.valuation(P)
            vec.append(v % 2)
            u = a / pi^v
            lg = mod.ideallog(u)
            vec += [lg[i] % 2 for i in even]
        return vector(GF(2), vec)
    def real_coords(self, alpha):
        vec = []
        for j, (K, th) in enumerate(self.fields):
            for emb in K.real_embeddings():
                vec.append(0 if emb(alpha[j]) > 0 else 1)
        return vector(GF(2), vec)
    def mu_x(self, x0):   # mu(P - oo) = x0 - theta
        return tuple(K(x0) - th for (K, th) in self.fields)
    def mu_a(self, a):    # a — унитарный квадратичный многочлен над Q: mu = a(theta)
        return tuple(a(th) for (K, th) in self.fields)
    # --- проверка «F — квадрат» локально
    def is_sq_Qp(self, a, p):
        a = QQ(a)
        if a == 0: return False
        v = a.valuation(p)
        if v % 2: return False
        u = a / p^v
        if p == 2: return (u.numerator()*u.denominator()) % 8 == 1
        return kronecker(ZZ(u.numerator()*u.denominator()), p) == 1
    def sq_in_quadratic(self, a, p):
        # f mod a — квадрат в Q_p[x]/(a)? a унитарный квадратичный над Q
        if a.is_irreducible():
            M = NumberField(a, 'z'); z = M.gen()
            beta = self.F(z)
            if beta == 0: return False
            for P in M.primes_above(p):
                v = beta.valuation(P)
                if v % 2: return False
                pi = M.uniformizer(P, others='positive'); u = beta/pi^v
                e = P.ramification_index() if p == 2 else 0; mod = P^(2*e+1)
                G = mod.idealstar(2); lg = mod.ideallog(u)
                if any(lg[i] % 2 for i, o in enumerate(G.gens_orders()) if o % 2 == 0): return False
            return True
        else:
            r = a.roots(QQ, multiplicities=False)
            return all(self.is_sq_Qp(self.F(r0), p) for r0 in r)
    def expected_dim(self, p):
        mp = len(self.F.change_ring(Qp(p, 60)).factor()) if p != oo else None
        return (mp - 1) + (2 if p == 2 else 0)
    def local_image(self, p, tries=60000, seed=1):
        random.seed(int(seed) + int(p))
        pl = self.primes_above(p)
        target = self.expected_dim(p)
        rows = []
        def add(vec):
            rows.append(vec)
            return matrix(GF(2), rows).rank()
        rk = 0
        for it in range(tries):
            if it % 3 == 0:
                # точка с рациональным x: u*p^k или корень + u*p^k
                k = random.randint(-4, 4)
                u0 = QQ(random.randint(-p^6, p^6))
                if random.random() < 0.4 and self.rroots:
                    x0 = random.choice(self.rroots) + u0 * p^abs(k)
                else:
                    x0 = u0 * p^k
                val = self.F(x0)
                if val == 0 or not self.is_sq_Qp(val, p): continue
                rk = add(self.local_coords(self.mu_x(x0), p, pl))
            else:
                if random.random() < 0.5:
                    x1 = QQ(random.randint(-p^5, p^5)) * p^random.randint(-3, 3)
                    x2 = QQ(random.randint(-p^5, p^5)) * p^random.randint(-3, 3)
                    if self.rroots and random.random() < 0.3: x1 = random.choice(self.rroots) + x1*p^2
                    a = (x - x1)*(x - x2)
                else:
                    b = QQ(random.randint(-p^5, p^5)) * p^random.randint(-3, 3)
                    c = QQ(random.randint(-p^5, p^5)) * p^random.randint(-4, 4)
                    a = x^2 + b*x + c
                if a.discriminant() == 0 or any(self.F(r0) == 0 for r0 in a.roots(QQbar, multiplicities=False) if r0 in QQ): continue
                if self.F.gcd(a) != 1: continue
                if not self.sq_in_quadratic(a, p): continue
                rk = add(self.local_coords(self.mu_a(a), p, pl))
            if rk == target: break
        M = matrix(GF(2), rows) if rows else matrix(GF(2), 0, len(self.local_coords(self.mu_x(QQ(1)), p, pl)))
        return M.row_space(), rk, target, pl
    def real_image(self, tries=4000, seed=1):
        # точки из интервалов, где F > 0; пары из разных интервалов; комплексно-сопряжённые пары
        random.seed(int(seed))
        mR = len(self.F.change_ring(RR).factor())
        target = (mR - 1) - 2
        rr = sorted(self.F.roots(RR, multiplicities=False))
        cuts = [-10^8] + rr + [10^8]
        ivs = [(cuts[i], cuts[i+1]) for i in range(len(cuts)-1) if self.F(QQ((cuts[i] + cuts[i+1])/2)) > 0]
        def rnd_in(iv):
            lo, hi = iv; return RR(lo + (hi - lo)*random.uniform(0.05, 0.95)).nearby_rational(max_error=RR(hi - lo)/10^6)
        rows = []; rk = 0
        n = len(self.real_coords(self.mu_x(QQ(1))))
        for it in range(tries):
            if target <= 0: break
            r = random.random()
            if r < 0.4:
                x0 = rnd_in(random.choice(ivs))
                if self.F(x0) <= 0: continue
                rows.append(self.real_coords(self.mu_x(x0)))
            elif r < 0.8 and len(ivs) >= 1:
                x1 = rnd_in(random.choice(ivs)); x2 = rnd_in(random.choice(ivs))
                if x1 == x2 or self.F(x1) <= 0 or self.F(x2) <= 0: continue
                rows.append(self.real_coords(self.mu_a((x - x1)*(x - x2))))
            else:
                re_ = QQ(random.randint(-10^4, 10^4))/random.randint(1, 100); im2 = QQ(random.randint(1, 10^6))/random.randint(1, 100)
                a_ = (x - re_)^2 + im2                    # комплексно-сопряжённая пара — всегда допустимый дивизор над R
                if self.F.gcd(a_) != 1: continue
                rows.append(self.real_coords(self.mu_a(a_)))
            rk = matrix(GF(2), rows).rank()
            if rk == target: break
        M = matrix(GF(2), rows) if rows else matrix(GF(2), 0, n)
        return M.row_space(), rk, max(target, 0)
    def selmer(self, verbose=True):
        # глобальная группа: prod K_j(S_j, 2), условие нормы
        gens = []   # (кортеж компонент)
        for j, (K, th) in enumerate(self.fields):
            Sj = sum([K.primes_above(p) for p in self.S], [])
            for g in K.selmer_generators(Sj, 2):
                comp = [Kk(1) for (Kk, _) in self.fields]; comp[j] = K(g)
                gens.append(tuple(comp))
        n = len(gens)
        if verbose: print(f"S = {self.S}, m = {self.m}, generators of prod K(S,2): {n}")
        # норма в Q*/Q*^2: координаты по {-1} ∪ S
        def qcoords(q):
            q = QQ(q); v = [0 if q > 0 else 1] + [q.valuation(p) % 2 for p in self.S]
            return v
        Nrows = []
        for g in gens:
            nq = prod(K(g[j]).norm() for j, (K, th) in enumerate(self.fields))
            Nrows.append(qcoords(nq))
        # вектор-пространство: коэффициенты e in F2^n; условие: sum e_i N_i = 0, и локальные условия
        conds = [matrix(GF(2), Nrows).transpose()]   # (len(S)+1) x n
        for p in self.S:
            W, rk, target, pl = self.local_image(p)
            if verbose: print(f"  p = {p}: local image dim {rk} (expected {target})")
            assert rk == target, f"local image at {p} incomplete"
            locs = matrix(GF(2), [self.local_coords(g, p, pl) for g in gens])   # n x d
            # условие: e*locs in W  <=>  e*locs*Q = 0, где Q — базис ортогонального дополнения W
            Wc = W.basis_matrix().right_kernel_matrix() if W.dimension() > 0 else identity_matrix(GF(2), locs.ncols())
            if Wc.nrows() > 0: conds.append((locs * Wc.transpose()).transpose())
        Wr, rkr, tr = self.real_image()
        if verbose: print(f"  v = oo: real image dim {rkr} (expected {tr})")
        assert rkr == tr
        locr = matrix(GF(2), [self.real_coords(g) for g in gens])
        Wcr = Wr.basis_matrix().right_kernel_matrix() if Wr.dimension() > 0 else identity_matrix(GF(2), locr.ncols())
        if Wcr.nrows() > 0 and locr.ncols() > 0: conds.append((locr * Wcr.transpose()).transpose())
        Cm = block_matrix([[c] for c in conds if c.nrows() > 0])
        sel = Cm.right_kernel()
        self.gens = gens
        return sel
    def rank_bound(self):
        sel = self.selmer()
        d = sel.dimension()
        return d, d - (self.m - 1)
    def in_selmer(self, alpha, sel):
        # выразить alpha через генераторы по модулю квадратов — через локальные координаты при большом наборе простых (контроль)
        return None
