# Полный 2-спуск над числовым полем k для E с ОДНОЙ рациональной точкой порядка 2 (Claude, 2026-09-12).
#
# E: y^2 = f(x) = x(x^2 + a x + b),  Delta = a^2 - 4b НЕ квадрат в k.
# Этальная алгебра L = k[x]/f(x) = k x K,  K = k(sqrt Delta) = k[t]/(t^2 + a t + b), theta — корень.
# H^1(k, E[2]) = ker(N_{L/k}: L*/L*^2 -> k*/k*^2)   (E[2] = ker(N: Res_{L/k} mu_2 -> mu_2)).
# mu: E(k)/2E(k) -> L*/L*^2,  P = (x,y) -> (x, x - theta);  T = (0,0) -> (f'(0), -theta) = (b, -theta).
# N_{L/k}(mu(P)) = x * N(x - theta) = x(x^2+ax+b) = y^2 — условие нормы (НЕ автоматическое, его надо накладывать).
#
# Sel^2 = { (alpha, beta) in k(S,2) x K(S_K,2) : alpha*N(beta) in k*^2,  res_v(alpha,beta) in Im mu_v для всех v in S u {oo} }.
# S = простые над 2*b*Delta (все плохие простые и простые над 2); при v ∉ S условие автоматично:
#   класс из L(S,2) неразветвлён в v, а при хорошей редукции и v∤2 Im mu_v = H^1_ur (равенство порядков).
# Локальные размерности dim_F2 E(k_v)/2E(k_v):
#   v конечное, v∤2: dim E(k_v)[2];   v|2: dim E(k_v)[2] + [k_v:Q_2];   v вещественное: dim E(R)[2] - 1.
# Локальные образы набираются случайными точками ДО ТОЧНОЙ теоретической размерности; недобор — RuntimeError.
# Итог: rank E(k) <= dim Sel^2 - dim E(k)[2] = dim Sel^2 - 1.
#
# ЗАМЕЧАНИЕ (GRH): классовые группы полей K (степень 2[k:Q]) считаются с proof=False.
import random, functools, time
print = functools.partial(print, flush=True)
proof.number_field(False)


class LocalSq:
    """ F_P^*/F_P^*2 для конечного простого P произвольного абсолютного числового поля F. """
    def __init__(self, F, P):
        self.F, self.P = F, P
        self.p = P.smallest_integer()
        e = P.ramification_index() if self.p == 2 else 0
        self.pi = F.uniformizer(P, others='positive')
        self.mod = P^(2*e + 1)
        G = self.mod.idealstar(2)
        self.even = [i for i, o in enumerate(G.gens_orders()) if o % 2 == 0]
        self.dimU = len(self.even)
        self.deg = P.residue_class_degree() * P.ramification_index()   # [F_P : Q_p]
        expU = 1 + (self.deg if self.p == 2 else 0)
        assert self.dimU == expU, (P, self.dimU, expU)
        self.dim = 1 + self.dimU
        self.name = f"{P.smallest_integer()}:{P.gens_two()[1] if P.gens_two()[1] != 0 else ''}"
    def coords(self, a):
        a = self.F(a); assert a != 0
        v = a.valuation(self.P)
        u = a / self.pi^v
        lg = self.mod.ideallog(u)
        return [ZZ(v) % 2] + [ZZ(lg[i]) % 2 for i in self.even]
    def is_sq(self, a):
        a = self.F(a)
        if a == 0: return False
        return not any(self.coords(a))


class RealSqAA:
    """ знак при вещественном вложении emb: F -> AA (точно, алгебраические вещественные). """
    def __init__(self, F, emb):
        self.F, self.emb, self.dim = F, emb, 1
        self.name = 'real'
    def coords(self, a):
        s = self.emb(self.F(a))
        assert s != 0
        return [0 if s > 0 else 1]
    def is_sq(self, a):
        a = self.F(a)
        return a != 0 and self.emb(a) > 0


class CplxSq:
    """ комплексное место: C*/C*^2 тривиально. """
    def __init__(self): self.dim = 0; self.name = 'cplx'
    def coords(self, a): return []
    def is_sq(self, a): return True


def rand_elt(F, bound):
    w = F.ring_of_integers().basis()
    return sum(ZZ(random.randint(-bound, bound)) * bb for bb in w)


class PartialDescent:
    def __init__(self, k, a, b, extra_primes=(), verbose=True):
        t0 = time.time()
        self.verbose = verbose
        self.k = k
        self.a = k(a); self.b = k(b)
        assert self.b != 0
        self.Delta = self.a^2 - 4*self.b
        assert self.Delta != 0
        assert not self.Delta.is_square(), "Delta — квадрат в k: нужен полный 2-спуск"
        Rx = PolynomialRing(k, 'X'); X = Rx.gen(); self.Rx = Rx
        self.f = X*(X^2 + self.a*X + self.b)

        # ---- множество S: простые над 2*b*Delta ----
        nrm = QQ((2*self.b*self.Delta).norm())
        ps = set(ZZ(nrm.numerator()).prime_factors()) | set(ZZ(nrm.denominator()).prime_factors()) | {2}
        ps |= set(ZZ(p) for p in extra_primes)
        self.ps = sorted(ps)
        self.S = sorted(sum([k.primes_above(p) for p in self.ps], []),
                        key=lambda P: (P.smallest_integer(), str(P)))
        self.Vk, self.gk, self.fromk, self.tok = k.selmer_space(self.S, 2)
        self.nk = self.Vk.dimension()

        # ---- K = k(sqrt d), d — уменьшенный представитель Delta в k(S,2) ----
        d = self.fromk(self.tok(self.Delta))
        sq = (self.Delta / d)
        assert sq.is_square()
        s = sq.sqrt()
        self.d = d
        T = PolynomialRing(k, 'T').gen()
        self.K = k.extension(T^2 - d, 'w')
        self.theta = (-self.a + s*self.K.gen())/2
        assert self.theta^2 + self.a*self.theta + self.b == 0
        self.Kabs = self.K.absolute_field('aa')
        self.fromKabs, self.toKabs = self.Kabs.structure()
        assert self.fromKabs.domain() is self.Kabs
        self.SK = sorted(sum([self.Kabs.primes_above(p) for p in self.ps], []),
                         key=lambda P: (P.smallest_integer(), str(P)))
        self.VK, self.gK, self.fromK, self.toK = self.Kabs.selmer_space(self.SK, 2)
        self.nK = self.VK.dimension()
        self.n = self.nk + self.nK
        if verbose:
            print(f"  k = {k.defining_polynomial()}, a = {self.a}, b = {self.b}")
            print(f"  Delta ~ d = {d};  K = k(sqrt d), abs poly {self.Kabs.defining_polynomial()}")
            print(f"  S: {len(self.S)} простых над {self.ps};  dim k(S,2) = {self.nk}, dim K(S_K,2) = {self.nK}  ({time.time()-t0:.0f}s)")

        # ---- места ----
        self.places = []
        for P in self.S:
            plk = LocalSq(k, P)
            Ws = self.K.primes_above(P)
            wabs = []
            for W in Ws:
                Wa = W.absolute_ideal(names='aa')
                assert Wa.number_field() is self.Kabs, "рассогласование абсолютного поля"
                wabs.append(Wa)
            # контроль соответствия относительных и абсолютных простых
            for W, Wa in zip(Ws, wabs):
                for _ in range(3):
                    z = self.K(rand_elt(k, 20)) + rand_elt(k, 20)*self.K.gen()
                    if z == 0: continue
                    assert z.valuation(W) == self.toKabs(z).valuation(Wa), "валюации не совпали"
            assert (sum(Wa.ramification_index()*Wa.residue_class_degree() for Wa in wabs)
                    == 2 * P.ramification_index() * P.residue_class_degree()), "сумма e*f над P не равна 2"
            self.places.append({'kind': 'fin', 'P': P, 'plk': plk,
                                'plK': [LocalSq(self.Kabs, Wa) for Wa in wabs]})
        embs_k = k.embeddings(AA)
        embs_K = self.Kabs.embeddings(AA)
        rk = k.gen()
        for sg in embs_k:
            ext = [tau for tau in embs_K if tau(self.toKabs(self.K(rk))) == sg(rk)]
            assert len(ext) in (0, 2)
            self.places.append({'kind': 'real', 'emb': sg, 'plk': RealSqAA(k, sg),
                                'plK': [RealSqAA(self.Kabs, tau) for tau in ext] if ext else [CplxSq()]})
        for pl in self.places:
            pl['dim'] = pl['plk'].dim + sum(q.dim for q in pl['plK'])
            pl['target'] = self.local_target(pl)
        if verbose:
            print("  места: " + ", ".join(f"{self.plname(pl)}[dim {pl['dim']}, tgt {pl['target']}]" for pl in self.places))

    # ---------- служебное ----------
    def plname(self, pl):
        return str(pl['P'].gens_two()) if pl['kind'] == 'fin' else ('real ' + str(pl['emb'](self.k.gen()).n(20)))

    def local_target(self, pl):
        """ dim_F2 E(k_v)/2E(k_v) """
        if pl['kind'] == 'real':
            return (2 if pl['plk'].is_sq(self.Delta) else 1) - 1
        t2 = 1 + (1 if pl['plk'].is_sq(self.Delta) else 0)     # dim E(k_v)[2]
        if pl['plk'].p == 2:
            t2 += pl['plk'].deg
        return t2

    def coords_pair(self, alpha, beta, pl):
        """ координаты (alpha, beta) in k_v*/k_v*^2 x prod_w K_w*/K_w*^2 """
        bK = self.toKabs(self.K(beta)) if beta.parent() is not self.Kabs else beta
        out = pl['plk'].coords(alpha)
        for q in pl['plK']:
            out += q.coords(bK)
        return out

    def gen_rows(self, pl):
        """ строки образов образующих k(S,2) x K(S_K,2) в локальных координатах места pl """
        rows = []
        zK = [0]*sum(q.dim for q in pl['plK'])
        for g in self.gk:
            rows.append(pl['plk'].coords(g) + zK)
        zk = [0]*pl['plk'].dim
        for g in self.gK:
            r = list(zk)
            for q in pl['plK']:
                r += q.coords(g)
            rows.append(r)
        return matrix(GF(2), rows)

    # ---------- локальные образы ----------
    def approx_roots(self, pl, m=40):
        """ приближённые (mod P^m) корни x^2+ax+b в k_P — центры для выборки.
            Нужны, когда Delta — квадрат в k_P: тогда K ⊗ k_P = k_P x k_P и классы с нечётной
            валюацией второй компоненты достигаются только при x, близких к корню. """
        if 'roots' in pl: return pl['roots']
        P = pl['P']; k = self.k; plk = pl['plk']
        out = []
        try:
            if plk.p != 2 and plk.is_sq(self.Delta):
                # корни (-a +- sqrt(Delta))/2: достаточно извлечь корень из Delta в k_P (Ньютон по y^2 = u)
                t = ZZ(self.Delta.valuation(P)); assert t % 2 == 0
                pi = plk.pi
                u = self.Delta / pi^t
                rf = P.residue_field()
                y = k(rf.lift(rf(u).sqrt()))
                I = P^(m + 2)
                for _ in range(2*m):
                    y = (y + u/y)/2
                    try:
                        y = I.reduce(y)
                    except Exception:
                        pass
                    if (y^2 - u).valuation(P) >= m: break
                if (y^2 - u).valuation(P) >= 8:
                    sD = y * pi^(t//2)
                    out = [(-self.a + sD)/2, (-self.a - sD)/2]
                    for r in out:
                        assert (r^2 + self.a*r + self.b).valuation(P) >= 8
        except Exception:
            out = []
        pl['roots'] = out
        return out

    def sample_x_fin(self, pl):
        p = pl['plk'].p
        j = random.randint(-8, 12)
        centers = [self.k(0)] + self.approx_roots(pl)
        return random.choice(centers) + rand_elt(self.k, p^6) * pl['plk'].pi^j

    def real_roots(self, pl):
        sg = pl['emb']
        A = sg(self.a); B = sg(self.b); D = sg(self.Delta)
        if D <= 0: return [AA(0)]
        sD = AA(D).sqrt()
        return sorted([AA(0), (-A - sD)/2, (-A + sD)/2])

    def sample_x_real(self, pl):
        rts = self.real_roots(pl)
        if len(rts) == 1:
            return QQ(random.randint(1, 10^6)) / random.randint(1, 10^3)
        r1, r2, r3 = rts
        if random.random() < 0.5:
            lo, hi = QQ(r1.n(80)) , QQ(r2.n(80))     # интервал между двумя меньшими корнями
            u = QQ(random.randint(1, 10^5)) / 10^5
            return lo + (hi - lo)*u
        else:
            hi = QQ(r3.n(80))
            return hi + QQ(random.randint(1, 10^7)) / random.randint(1, 10^3)

    def local_image(self, pl, tries=300000):
        target = pl['target']
        rows = [self.coords_pair(self.b, -self.theta, pl)]      # mu(T)
        M = matrix(GF(2), rows); rk = M.rank()
        assert rk <= target, f"образ кручения вне теоретической размерности в {self.plname(pl)}"
        it = 0
        while rk < target and it < tries:
            it += 1
            x = self.sample_x_real(pl) if pl['kind'] == 'real' else self.sample_x_fin(pl)
            x = self.k(x)
            if x == 0: continue
            fx = self.f(x)
            if fx == 0 or not pl['plk'].is_sq(fx): continue
            rows.append(self.coords_pair(x, x - self.theta, pl))
            M = matrix(GF(2), rows); nrk = M.rank()
            if nrk == rk: rows.pop()
            rk = nrk
        if rk != target:
            raise RuntimeError(f"локальный образ не набран в {self.plname(pl)}: {rk} < {target}")
        return matrix(GF(2), rows).row_space()

    # ---------- Selmer ----------
    def selmer(self, verbose=None):
        verbose = self.verbose if verbose is None else verbose
        conds = []
        for pl in self.places:
            t0 = time.time()
            W = self.local_image(pl)
            Lg = self.gen_rows(pl)
            d = pl['dim']
            Wc = W.basis_matrix().right_kernel_matrix().transpose() if W.dimension() > 0 else identity_matrix(GF(2), d)
            conds.append(Lg * Wc)
            if verbose:
                print(f"    {self.plname(pl)}: dim Im mu_v = {W.dimension()} / {d}  ({time.time()-t0:.0f}s)")
        # условие нормы: alpha * N(beta) in k*^2
        Nrows = []
        for g in self.gk:
            Nrows.append(list(self.tok(g)))
        for g in self.gK:
            Nrows.append(list(self.tok(self.fromKabs(g).relative_norm())))
        conds.append(matrix(GF(2), self.n, self.nk, Nrows))
        Mall = block_matrix(GF(2), [conds], subdivide=False)
        Sel = Mall.left_kernel()
        self.Sel = Sel
        return Sel

    # ---------- контроли ----------
    def torsion_vector(self):
        v = list(self.tok(self.b)) + list(self.toK(self.toKabs(-self.theta)))
        return vector(GF(2), v)

    def check_all(self):
        out = {}
        tv = self.torsion_vector()
        out['torsion_nonzero'] = bool(tv != 0)
        out['torsion_in_Sel'] = bool(tv in self.Sel)
        # представители: alpha*N(beta) должен быть квадратом в k
        ok = True
        for v in self.Sel.basis():
            al = prod([g for g, c in zip(self.gk, list(v)[:self.nk]) if c == 1], self.k(1))
            be = prod([g for g, c in zip(self.gK, list(v)[self.nk:]) if c == 1], self.Kabs(1))
            ok = ok and bool((al * self.fromKabs(be).relative_norm()).is_square())
        out['norm_condition'] = ok
        # mu(T) действительно обнуляется нормой: b * N(-theta) = b * b = b^2
        out['N_theta'] = bool(self.fromKabs(self.toKabs(self.theta)).relative_norm() == self.b)
        # независимый контроль правила mu(T): mu(T) = 1 <=> T in 2E(k) (мю инъективно на E(k)/2E(k))
        EE = EllipticCurve(self.k, [0, self.a, 0, self.b, 0])
        T_div = len(EE(0, 0).division_points(2)) > 0
        out['T_in_2E'] = bool(T_div)
        out['mu_T_consistent'] = bool(out['torsion_nonzero'] != T_div)
        return out


def partial_descent_bound(k, a, b, extra_primes=(), verbose=True):
    C = PartialDescent(k, a, b, extra_primes=extra_primes, verbose=verbose)
    Sel = C.selmer()
    chk = C.check_all()
    # torsion_nonzero = False допустимо: это значит T = (0,0) in 2E(k) (mu инъективно на E(k)/2E(k));
    # оценка rank <= dim Sel - dim E(k)[2] = dim Sel - 1 от этого не зависит.
    assert chk['torsion_in_Sel'], f"кручение вне Sel: {chk}"
    assert chk['norm_condition'] and chk['N_theta'], f"условие нормы нарушено: {chk}"
    assert chk['mu_T_consistent'], f"правило mu(T) = (b, -theta) не согласовано с делимостью T на 2: {chk}"
    dim = Sel.dimension()
    if verbose:
        print(f"  dim Sel^2 = {dim}  =>  rank E(k) <= {dim - 1}   (контроли: {chk})")
    return dim, dim - 1, C
