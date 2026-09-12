# 2-спуск для эллиптической кривой с ПОЛНЫМ 2-кручением над квадратичным полем k (Claude, 2026-09-11).
# E: y^2 = (x - e1)(x - e2)(x - e3), e_i в O_k.
#
# (a) Полный 2-спуск: kappa: E(k)/2E(k) -> k(S,2)^2, P -> (x - e1, x - e2) (особые значения на T1, T2).
#     Локально |E(k_v)/2E(k_v)| = 4 * 2^[k_v:Q_2] (v | 2), 4 (v конечное, v ∤ 2), 2 (v вещественное, все e_i вещественны).
#     rank E(k) <= dim Sel^2 - 2.
# (b) Спуск через 2-изогению с ядром <(e_i, 0)>: сдвиг x -> x + e_i даёт y^2 = x(x^2 + a x + b), E': Y^2 = X(X^2 - 2aX + a^2 - 4b).
#     alpha: E(k) -> k*/k*^2, (x,y) -> x, (0,0) -> b;   alpha': E'(k) -> k*/k*^2, (X,Y) -> X, (0,0) -> a^2 - 4b.
#     Локально dim Im alpha_v + dim Im alpha'_v = dim E(k_v)/2E(k_v) (E[2] рационально, phi(E[2]) = E'[phi^] = <(0,0)>;
#     у E' полное 2-кручение только если b — квадрат, для формулы это не нужно),
#     rank E(k) <= dim Sel(alpha) + dim Sel(alpha') - 2.
# Локальные образы: случайные точки E(k_v) до ТОЧНОЙ теоретической размерности (иначе — ошибка, а не ответ).
# Локальные квадратичные классы: v_P(.) mod 2 и ideallog по модулю P^(2e+1) (e = ветвление над 2; для нечётных P — P).
import random, functools
print = functools.partial(print, flush=True)

class LocalSq:
    """ k_P^*/k_P^*2 для конечного P и знаки для вещественных вложений. """
    def __init__(self, k, P):
        self.k, self.P = k, P
        self.p = P.smallest_integer()
        e = P.ramification_index() if self.p == 2 else 0
        self.pi = k.uniformizer(P, others='positive')
        self.mod = P^(2*e + 1)
        G = self.mod.idealstar(2)
        self.even = [i for i, o in enumerate(G.gens_orders()) if o % 2 == 0]
        self.dimU = len(self.even)
        self.deg = P.residue_class_degree() * P.ramification_index()   # [k_P : Q_p]
        expU = 1 + (self.deg if self.p == 2 else 0)
        assert self.dimU == expU, (P, self.dimU, expU)
        self.dim = 1 + self.dimU
    def coords(self, a):
        a = self.k(a); assert a != 0
        v = a.valuation(self.P)
        u = a / self.pi^v
        lg = self.mod.ideallog(u)
        return [v % 2] + [lg[i] % 2 for i in self.even]
    def is_sq(self, a):
        a = self.k(a)
        if a == 0: return False
        return not any(self.coords(a))

def exact_sign(k, a, emb):
    """ точный знак a = a0 + a1 r (r^2 = D) при вещественном вложении emb (sgn emb(r) = sr) — сравнением в Q """
    a = k(a); r = k.gen(); D = QQ(r^2); assert k.degree() == 2 and r^2 == D
    a0, a1 = list(a)
    sr = 1 if emb(r) > 0 else -1
    if a1 == 0: return sign(a0)
    b = a1 * sr                     # a = a0 + b sqrt(D)
    if a0 == 0: return sign(b)
    if sign(a0) == sign(b): return sign(a0)
    return sign(a0) if a0^2 > b^2 * D else sign(b)

class RealSq:
    def __init__(self, k, emb):
        self.k, self.emb, self.dim = k, emb, 1
    def coords(self, a): return [0 if exact_sign(self.k, a, self.emb) > 0 else 1]
    def is_sq(self, a): return exact_sign(self.k, a, self.emb) > 0

def rand_elt(k, bound):
    w = k.ring_of_integers().basis()
    return sum(ZZ(random.randint(-bound, bound)) * b for b in w)

class Curve3:
    def __init__(self, k, roots, extra_primes=()):
        self.k = k; self.e = [k(t) for t in roots]
        e1, e2, e3 = self.e
        self.Rx = PolynomialRing(k, 'X'); X = self.Rx.gen()
        self.f = (X - e1)*(X - e2)*(X - e3)
        disc = ((e1 - e2)*(e1 - e3)*(e2 - e3))
        ps = set(ZZ(2*disc.norm()).prime_factors()) | set(extra_primes)
        self.S = sorted(sum([k.primes_above(p) for p in sorted(ps)], []), key=lambda P: (P.smallest_integer(), str(P)))
        self.places = [LocalSq(k, P) for P in self.S] + [RealSq(k, emb) for emb in k.real_embeddings()]
        self.gens = k.selmer_generators(self.S, 2)
        self.n = len(self.gens)
        self._reps = {}

    def target_full(self, pl):
        if isinstance(pl, RealSq): return 1
        return 2 + (pl.deg if pl.p == 2 else 0)

    def sample_x(self, pl, centers):
        """ случайное x в k, «локально случайное» около центров и бесконечности для места pl """
        if isinstance(pl, RealSq):
            c = random.choice(centers)
            return c + QQ(random.randint(-10^6, 10^6)) / random.randint(1, 10^3) * random.choice([1, 10^-3, 10^3, 10^6])
        j = random.randint(-8, 12)
        return random.choice(centers + [0]) + rand_elt(self.k, pl.p^6) * pl.pi^j

    # ---------- (a) полный 2-спуск ----------
    def kappa(self, x):
        e1, e2, e3 = self.e
        if x == e1: return ((e1 - e2)*(e1 - e3), e1 - e2)
        if x == e2: return (e2 - e1, (e2 - e1)*(e2 - e3))
        return (x - e1, x - e2)

    def local_image_full(self, pl, tries=200000):
        target = self.target_full(pl)
        rows = []; rk = 0; reps = []
        # 2-кручение
        for x in self.e:
            a, b = self.kappa(x); rows.append(pl.coords(a) + pl.coords(b)); reps.append((a, b))
        rk = matrix(GF(2), rows).rank()
        it = 0
        while rk < target and it < tries:
            it += 1
            x = self.sample_x(pl, self.e)
            fx = self.f(x)
            if fx == 0 or not pl.is_sq(fx): continue
            a, b = self.kappa(x); rows.append(pl.coords(a) + pl.coords(b))
            nrk = matrix(GF(2), rows).rank()
            if nrk > rk: reps.append((a, b))
            else: rows.pop()
            rk = nrk
        if rk != target: raise RuntimeError(f"full image incomplete at {pl.__dict__.get('P', 'real')}: {rk} < {target}")
        self._reps[id(pl)] = reps
        return matrix(GF(2), rows).row_space()

    def selmer_full(self, verbose=True):
        n = self.n
        conds = []
        for pl in self.places:
            W = self.local_image_full(pl)
            Lg = [pl.coords(g) for g in self.gens]              # n x d
            d = pl.dim
            # w = (w1, w2) in F2^{2n}  ->  (sum w1_i Lg_i, sum w2_i Lg_i) in F2^{2d}
            M = matrix(GF(2), 2*n, 2*d)
            for i in range(n):
                for c in range(d):
                    M[i, c] = Lg[i][c]; M[n + i, d + c] = Lg[i][c]
            Wc = W.basis_matrix().right_kernel_matrix().transpose() if W.dimension() > 0 else identity_matrix(GF(2), 2*d)
            conds.append(M * Wc)
        Mall = block_matrix(GF(2), [conds], subdivide=False) if conds else matrix(GF(2), 2*n, 0)
        Sel = Mall.left_kernel()
        return Sel

    def vec_of_pair(self, a, b):
        """ координаты пары (a, b) из k(S,2)^2 в базисе self.gens (через локальные координаты всех мест S и доп. проверку) """
        raise NotImplementedError

    # ---------- (b) спуск через 2-изогению ----------
    def isogeny_data(self, i):
        ei = self.e[i]; ej, ek = [self.e[t] for t in range(3) if t != i]
        a = 2*ei - ej - ek; b = (ei - ej)*(ei - ek)
        return a, b

    def local_image_alpha(self, pl, a, b, tries=40000, target=None):
        """ образ alpha: y^2 = x(x^2 + a x + b), (x,y) -> x, (0,0) -> b. Возвращает подпространство (наибольшее найденное). """
        Rx = self.Rx; X = Rx.gen(); g = X*(X^2 + a*X + b)
        rts = [t for t, _ in g.roots()]
        rows = []
        for t in rts:
            rows.append(pl.coords(b if t == 0 else t))
        rk = matrix(GF(2), rows).rank()
        it = 0
        while (target is None or rk < target) and it < tries:
            it += 1
            x = self.sample_x(pl, rts)
            gx = g(x)
            if gx == 0 or not pl.is_sq(gx): continue
            rows.append(pl.coords(x)); rk = matrix(GF(2), rows).rank()
        return matrix(GF(2), rows).row_space(), rk

    def selmer_from_images(self, images):
        n = self.n; conds = []
        for pl, W in images:
            Lg = matrix(GF(2), [pl.coords(g) for g in self.gens])     # n x d
            Wc = W.basis_matrix().right_kernel_matrix().transpose() if W.dimension() > 0 else identity_matrix(GF(2), pl.dim)
            conds.append(Lg * Wc)
        Mall = block_matrix(GF(2), [conds], subdivide=False)
        return Mall.left_kernel()

    def isogeny_descent(self, i, tries=200000, verbose=False):
        a, b = self.isogeny_data(i)
        a2, b2 = -2*a, a^2 - 4*b
        X = self.Rx.gen()
        curves = [(X*(X^2 + a*X + b), b), (X*(X^2 + a2*X + b2), b2)]
        imgs1, imgs2 = [], []
        for pl in self.places:
            target = self.target_full(pl)
            rows = [[], []]; rts = []
            for c, (g, bb) in enumerate(curves):
                rt = [t for t, _ in g.roots()]; rts.append(rt)
                for t in rt: rows[c].append(pl.coords(bb if t == 0 else t))
            rk = [matrix(GF(2), rows[c]).rank() for c in range(2)]
            it = 0
            while rk[0] + rk[1] < target and it < tries:     # оба образа попеременно; сумма размерностей точно известна
                it += 1; c = it % 2; g, bb = curves[c]
                x = self.sample_x(pl, rts[c])
                gx = g(x)
                if gx == 0 or not pl.is_sq(gx): continue
                rows[c].append(pl.coords(x)); rk[c] = matrix(GF(2), rows[c]).rank()
            if rk[0] + rk[1] != target:
                raise RuntimeError(f"isogeny images incomplete at {getattr(pl, 'P', 'real')}: {rk} vs {target}")
            if verbose: print(f"   place {getattr(pl, 'P', 'real')}: dims {rk[0]} + {rk[1]} = {target}")
            imgs1.append((pl, matrix(GF(2), rows[0]).row_space())); imgs2.append((pl, matrix(GF(2), rows[1]).row_space()))
        S1 = self.selmer_from_images(imgs1); S2 = self.selmer_from_images(imgs2)
        return S1.dimension(), S2.dimension()

    # ---------- (b') изогенный спуск через точные образы ----------
    # Im alpha_v (ядро <(e_i,0)>, alpha(P) = x - e_i) = проекция точного образа kappa_v на i-ю координату
    # (x - e3 = (x - e1)(x - e2) по модулю квадратов; значения на кручении согласованы с kappa).
    # Im alpha'_v = аннулятор Im alpha_v относительно локального символа Гильберта (двойственность образов 2-изогении).
    def hilb(self, a, b, pl):
        if isinstance(pl, RealSq): return self.k.hilbert_symbol(a, b, pl.emb)
        return self.k.hilbert_symbol(a, b, pl.P)

    def isogeny_descent_dual(self, i, verbose=False):
        n = self.n
        condsA, condsB = [], []
        for pl in self.places:
            if id(pl) not in self._reps: self.local_image_full(pl)
            reps = self._reps[id(pl)]
            ds = [a if i == 0 else (b if i == 1 else a*b) for (a, b) in reps]
            W = matrix(GF(2), [pl.coords(d) for d in ds]).row_space()
            Lg = matrix(GF(2), [pl.coords(g) for g in self.gens])
            Wc = W.basis_matrix().right_kernel_matrix().transpose() if W.dimension() > 0 else identity_matrix(GF(2), pl.dim)
            condsA.append(Lg * Wc)
            # c in Im alpha'_v  <=>  (c, d)_v = 1 для всех d из базиса Im alpha_v
            dB = [d for d in ds]
            H = matrix(GF(2), n, len(dB), lambda j, t: 0 if self.hilb(self.gens[j], dB[t], pl) == 1 else 1)
            condsB.append(H)
            if verbose: print(f"   place {getattr(pl, 'P', 'real')}: dim Im alpha = {W.dimension()}, dim k_v*/k_v*2 = {pl.dim}")
        SA = block_matrix(GF(2), [condsA], subdivide=False).left_kernel()
        SB = block_matrix(GF(2), [condsB], subdivide=False).left_kernel()
        return SA.dimension(), SB.dimension()

    def check_duality(self, i, tries=200000):
        """ контроль: сэмплированный образ alpha'_v должен лежать в аннуляторе Im alpha_v, а размерности — совпадать """
        a, b = self.isogeny_data(i); a2, b2 = -2*a, a^2 - 4*b
        X = self.Rx.gen(); g2 = X*(X^2 + a2*X + b2)
        out = []
        for pl in self.places:
            reps = self._reps[id(pl)]
            ds = [a_ if i == 0 else (b_ if i == 1 else a_*b_) for (a_, b_) in reps]
            dimA = matrix(GF(2), [pl.coords(d) for d in ds]).rank()
            rt = [t for t, _ in g2.roots()]
            samples = [b2 if t == 0 else t for t in rt]
            rk = matrix(GF(2), [pl.coords(c) for c in samples]).rank(); it = 0
            while rk < pl.dim - dimA and it < tries:
                it += 1
                x = self.sample_x(pl, rt); gx = g2(x)
                if gx == 0 or not pl.is_sq(gx): continue
                samples.append(x); rk = matrix(GF(2), [pl.coords(c) for c in samples]).rank()
            orth = all(self.hilb(c, d, pl) == 1 for c in samples for d in ds)
            out.append((str(getattr(pl, 'P', 'real')), dimA, rk, pl.dim - dimA, orth))
        return out
