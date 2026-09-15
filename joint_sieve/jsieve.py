# Fable, 15.09.2026. Совместное решето Морделла–Вейля по нескольким эллиптическим множителям
# кривой наклона X_k (все восемь клеток 1+λz — квадраты) при ОБЩЕМ z и согласованных знаках корней.
#
# Множитель: T ⊂ Λ, |T| ∈ {3,4}, C_T: y² = ∏_{λ∈T}(1+λz) ≅ E_T (модель weier из census_genus1, начало O ↔ z = −1/c0, c0 = T[0]).
# Глобально: E_T(ℚ) = ⊕ ℤ G_j ⊕ T_tors (PARI ellrank, нижняя = верхняя; Sage saturation — полное насыщение).
# Класс точки — (n_1..n_r mod N_T, координаты кручения). Локально при простом ℓ хорошей редукции: образ класса в
# E_T(𝔽_ℓ)/N_T; допустимые наборы образов = образы точек X_k(𝔽_ℓ): z̄ ∈ P¹(𝔽_ℓ) со всеми клетками-квадратами (или 0),
# знаки корней ε ∈ {±1}^8 общие для всех множителей (y_T = ∏_{λ∈T} y_λ). Над z̄ = ∞ (ℓ | знаменатель z):
# возможно лишь при −1 = □ (mod ℓ) и всех λ в одном квадратичном классе; знаки тоже согласованы.
from sage.all import *
import numpy as np, itertools, json, sys, time, os
sys.path.insert(0, '/home/kep/magicKube/census_genus1')
from g1census import weier


def slope_cells(r, s):
    return [s, -s, r, -r, s - r, r - s, s + r, -s - r]


class Factor:
    def __init__(self, T, verbose=False):
        self.T = sorted(ZZ(t) for t in T)
        E, L, c0 = weier([QQ(t) for t in self.T])
        assert E.discriminant() != 0
        self.E, self.L, self.c0 = E, ZZ(L), ZZ(c0)
        a = E.a_invariants(); assert a[0] == 0 and a[2] == 0
        self.a2, self.a4, self.a6 = ZZ(a[1]), ZZ(a[3]), ZZ(a[4])
        rk = E.pari_curve().ellrank()
        self.lo, self.hi = int(rk[0]), int(rk[1])
        pts = [E(list(p)) for p in rk[3]]
        if pts:
            sat, index, reg = E.saturation(pts)   # полное насыщение (граница индекса по высотам, Sage)
            self.gens, self.sat_index, self.reg = list(sat), ZZ(index), reg
        else:
            self.gens, self.sat_index, self.reg = [], ZZ(1), 1
        self.rank = len(self.gens)
        self.rank_proved = (self.lo == self.hi == self.rank)
        tors = E.torsion_subgroup()
        self.tgens = [E(g) for g in tors.gens()]
        self.tinv = [int(g.order()) for g in self.tgens]   # порядки в порядке gens() (invariants() может идти иначе)
        assert sorted(self.tinv) == sorted(int(n) for n in tors.invariants()) and prod(self.tinv) == tors.order()
        self.tdict = {}
        for co in itertools.product(*[range(n) for n in self.tinv]):
            P = E(0)
            for c, g in zip(co, self.tgens): P = P + c * g
            self.tdict[P] = co
        self.dim = self.rank + len(self.tinv)
        self.bad = set(ZZ(E.discriminant() * self.L * self.c0).prime_divisors())
        self.name = ','.join(str(t) for t in self.T)
        self.lcache = {}
        if verbose: print(self.name, 'rank', self.lo, self.hi, 'tors', self.tinv, 'sat index', self.sat_index, flush=True)

    # точка E_T над z с заданным знаком y = ∏ y_λ (y=+1 при z=0)
    def point(self, z, y):
        z = QQ(z); y = QQ(y)
        if z + 1 / self.c0 == 0: return self.E(0)
        w = 1 / (z + 1 / self.c0)
        X = self.L * w; Y = self.L * self.c0 * w**2 * y
        return self.E(X, Y)

    def z_of(self, P):
        if P.is_zero(): return -1 / self.c0
        if P[0] == 0: return infinity
        return -1 / self.c0 + self.L / P[0]

    # разложение P = Σ n_j G_j + t; возвращает список координат длины dim
    def decompose(self, P):
        r = self.rank
        n = []
        if r > 0:
            M = matrix(RR, r, r, lambda i, j: (self.gens[i] + self.gens[j]).height() - self.gens[i].height() - self.gens[j].height())
            v = vector(RR, [(P + G).height() - P.height() - G.height() for G in self.gens])
            sol = M.solve_right(v)
            n = [int(round(x)) for x in sol]
            assert max(abs(sol[i] - n[i]) for i in range(r)) < 1e-6, (sol, n)
        R = P
        for c, G in zip(n, self.gens): R = R - c * G
        assert R in self.tdict, 'остаток не кручение'
        return n + list(self.tdict[R])

    def neg_class(self, cls):
        r = self.rank
        return [-x for x in cls[:r]] + [(-c) % m for c, m in zip(cls[r:], self.tinv)]

    # локальные данные при ℓ: инварианты группы, логарифмы образующих (x,y), фибры над z̄ ∈ Z_ℓ
    def local(self, ell, zbars, roots, inf_roots):
        # zbars: список z̄ ∈ 𝔽_ℓ; roots: dict λ -> корень ȳ_λ (0 если клетка 0); inf_roots: dict λ -> ȳ'_λ или None
        key = ell
        if key in self.lcache: return self.lcache[key]
        F = GF(ell); El = self.E.change_ring(F); G = El.abelian_group()
        # ВАЖНО: discrete_log даёт координаты в порядке G.gens(), а G.invariants() может идти в другом порядке
        inv = [int(g.order()) for g in G.gens()]
        assert sorted(inv) == sorted(int(m) for m in G.invariants())
        if len(inv) == 1: inv = [inv[0], 1]
        elif len(inv) == 0: inv = [1, 1]
        def lg(P):
            if P.is_zero(): return (0, 0)
            v = [int(x) for x in G.discrete_log(P)]
            while len(v) < 2: v.append(0)
            return (v[0] % inv[0], v[1] % inv[1])
        glog = [lg(El(P)) for P in self.gens] + [lg(El(P)) for P in self.tgens]
        L, c0 = F(self.L), F(self.c0)
        fib = {}
        for zb in zbars:
            zz = F(zb)
            if zz + 1 / c0 == 0:
                fib[zb] = ((0, 0), (0, 0)); continue
            w = 1 / (zz + 1 / c0); X = L * w
            y = F(1)
            for t in self.T: y *= roots[zb][t]
            Y = L * c0 * w**2 * y
            P = El(X, Y)
            fib[zb] = (lg(P), lg(-P))
        if inf_roots is not None:
            if len(self.T) == 3:
                P = El(0, 0)
            else:
                y = F(1)
                for t in self.T: y *= inf_roots[t]
                Y = L * c0 * y / inf_roots['z0']**2
                P = El(0, Y)
            fib['inf'] = (lg(P), lg(-P))
        res = dict(inv=inv, glog=glog, fib=fib)
        self.lcache[key] = res
        return res


def local_zset(ell, cellset):
    """z̄ ∈ 𝔽_ℓ, при которых все клетки из cellset — квадраты или 0; корни; данные над ∞."""
    F = GF(ell)
    sq = {}
    for a in range(ell):
        sq.setdefault((a * a) % ell, a)
    zbars = []; roots = {}
    for zb in range(ell):
        rt = {}
        ok = True
        for lam in cellset:
            c = (1 + lam * zb) % ell
            if c in sq: rt[lam] = F(sq[c])
            else: ok = False; break
        if ok: zbars.append(zb); roots[zb] = rt
    inf_roots = None
    if ell % 4 == 1:
        cl = set()
        for lam in cellset:
            cl.add(F(lam).is_square())
        if len(cl) == 1:
            z0 = F(cellset[0])   # λ z0 = λ·λ0 — квадрат для всех λ (все в одном классе)
            inf_roots = {'z0': z0}
            for lam in cellset:
                inf_roots[lam] = F(sq[int(F(lam) * z0)])
    return zbars, roots, inf_roots


class JointSieve:
    def __init__(self, r, s, cellset, factors, verbose=True):
        self.r, self.s = r, s
        self.cellset = list(cellset)
        self.cidx = {lam: i for i, lam in enumerate(self.cellset)}
        self.F = factors
        for f in factors:
            assert all(t in self.cidx for t in f.T)
        self.verbose = verbose
        self.zcache = {}

    def zdata(self, ell):
        if ell not in self.zcache: self.zcache[ell] = local_zset(ell, self.cellset)
        return self.zcache[ell]

    def good_primes(self, lmin, lmax, idx):
        bad = set()
        for i in idx: bad |= self.F[i].bad
        for lam in self.cellset: bad |= set(ZZ(lam).prime_divisors())   # ℓ | λ: разбор над ∞ и нули клеток требуют ℓ ∤ λ
        return [int(l) for l in prime_range(lmin, lmax) if l > 2 and l not in bad]

    # допустимые образы при ℓ для множителей idx с модулями N: массивы (key1, key2) и параметры кодирования
    def admissible(self, ell, idx, N):
        zbars, roots, inf_roots = self.zdata(ell)
        locs = [self.F[i].local(ell, zbars, roots, inf_roots) for i in idx]
        mods = [(gcd(l['inv'][0], N[i]), gcd(l['inv'][1], N[i])) for l, i in zip(locs, idx)]
        # матрицы образов: для множителя i, столбцы класса -> (x,y) mod (a1,a2)
        m = len(idx); nz = None
        adm = set()
        chis = []
        for i in idx:
            chis.append([self.cidx[t] for t in self.F[i].T])
        def enc(u, j):
            a1, a2 = mods[j]
            return (u[0] % a1) + a1 * (u[1] % a2)
        for zb in list(zbars) + (['inf'] if inf_roots is not None else []):
            if zb == 'inf':
                nzc = [self.cidx[l] for l in self.cellset]
            else:
                nzc = [self.cidx[l] for l in self.cellset if roots[zb][l] != 0]
            pos = {c: k for k, c in enumerate(nzc)}
            fibs = [locs[j]['fib'][zb] for j in range(m)]
            for eps in itertools.product((0, 1), repeat=len(nzc)):
                tup = []
                for j in range(m):
                    sgn = sum(eps[pos[c]] for c in chis[j] if c in pos) % 2
                    tup.append(enc(fibs[j][sgn], j))
                adm.add(tuple(tup))
        return locs, mods, adm

    def sieve(self, cand, idx, N, primes, cols):
        """cand: np.array (n, D) классов; idx: индексы множителей; N: модули; cols: для каждого множителя — список столбцов."""
        stats = []
        for ell in primes:
            if cand.shape[0] == 0: break
            locs, mods, adm = self.admissible(ell, idx, N)
            m = len(idx)
            # размеры кодов
            U = [a1 * a2 for a1, a2 in mods]
            # образы кандидатов
            us = []
            for j, i in enumerate(idx):
                lg = locs[j]['glog']
                Lx = np.array([g[0] for g in lg], dtype=np.int64)
                Ly = np.array([g[1] for g in lg], dtype=np.int64)
                sub = cand[:, cols[j]]
                a1, a2 = mods[j]
                x = (sub @ Lx) % a1; y = (sub @ Ly) % a2
                us.append(x + a1 * y)
            # ключи: множители группируем так, чтобы произведение размеров кодов в группе было < 2^30;
            # группы склеиваем последовательно через индексы в отсортированных допустимых ключах (без переполнения)
            groups = []; cur = []; st = 1
            for j in range(m):
                if cur and st * U[j] >= 2**30: groups.append(cur); cur = []; st = 1
                cur.append(j); st *= U[j]
            groups.append(cur)
            A = np.array(sorted(adm), dtype=np.int64).reshape(len(adm), m)
            def gkey(arrs, grp, n):
                k = np.zeros(n, dtype=np.int64); st = 1
                for j in grp: k = k + arrs[j] * st; st *= U[j]
                return k
            nc = cand.shape[0]
            keep = np.ones(nc, dtype=bool)
            idx_c = np.zeros(nc, dtype=np.int64); idx_a = np.zeros(len(adm), dtype=np.int64)
            for grp in groups:
                kc = idx_c * 2**30 + gkey(us, grp, nc)
                ka = idx_a * 2**30 + gkey([A[:, j] for j in range(m)], grp, len(adm))
                ua = np.unique(ka)
                pos = np.searchsorted(ua, kc); posc = np.minimum(pos, len(ua) - 1)
                keep &= (ua[posc] == kc)
                idx_c = posc; idx_a = np.searchsorted(ua, ka)
                assert len(ua) < 2**30
            n0 = cand.shape[0]
            cand = cand[keep]
            stats.append((ell, n0, int(cand.shape[0]), len(adm)))
            if self.verbose: print(f'  ℓ={ell}: {n0} -> {cand.shape[0]} (|adm|={len(adm)}, |Z|={len(self.zdata(ell)[0])}{"+∞" if self.zdata(ell)[2] is not None else ""})', flush=True)
        return cand, stats

    def all_classes(self, i, Ni):
        f = self.F[i]
        ranges = [range(Ni)] * f.rank + [range(n) for n in f.tinv]
        return np.array(list(itertools.product(*ranges)), dtype=np.int64).reshape(-1, f.dim)

    def run(self, order, N0, primes, lifts=(), known=None):
        """order: порядок добавления множителей; N0: dict i->N_i; lifts: список (i, p) — поднять N_i в p раз."""
        idx = []; cols = []; N = dict(N0)
        cand = None
        log = []
        for i in order:
            f = self.F[i]
            new = self.all_classes(i, N[i])
            if cand is None:
                cand = new; cols = [list(range(f.dim))]
            else:
                D = cand.shape[1]
                cand = np.hstack([np.repeat(cand, new.shape[0], axis=0), np.tile(new, (cand.shape[0], 1))])
                cols.append(list(range(D, D + f.dim)))
            idx.append(i)
            if self.verbose: print(f'+ множитель {f.name} (ранг {f.rank}, N={N[i]}): кандидатов {cand.shape[0]}', flush=True)
            cand, st = self.sieve(cand, idx, N, primes, cols)
            log.append(dict(step=('add', f.name), N=dict(N), survivors=int(cand.shape[0]), stats=st))
            if known is not None: self.check_known(cand, idx, cols, N, known)
        for (i, p) in lifts:
            j = idx.index(i); f = self.F[i]; r = f.rank
            if r == 0: N[i] *= p; continue
            ks = np.array(list(itertools.product(range(p), repeat=r)), dtype=np.int64).reshape(-1, r)
            n0 = cand.shape[0]
            big = np.repeat(cand, ks.shape[0], axis=0)
            add = np.tile(ks, (n0, 1)) * N[i]
            big[:, cols[j][:r]] += add
            N[i] *= p
            cand = big
            if self.verbose: print(f'↑ N[{f.name}] -> {N[i]}: кандидатов {cand.shape[0]}', flush=True)
            cand, st = self.sieve(cand, idx, N, primes, cols)
            log.append(dict(step=('lift', f.name, p), N=dict(N), survivors=int(cand.shape[0]), stats=st))
            if known is not None: self.check_known(cand, idx, cols, N, known)
        return cand, idx, cols, N, log

    def class_tuple(self, cls_list, idx, N):
        out = []
        for j, i in enumerate(idx):
            f = self.F[i]; c = cls_list[i]
            out += [x % N[i] for x in c[:f.rank]] + list(c[f.rank:])
        return tuple(int(x) for x in out)

    def check_known(self, cand, idx, cols, N, known):
        S = set(map(tuple, cand.tolist()))
        for name, cls_list in known:
            t = self.class_tuple(cls_list, idx, N)
            print(f'  контроль «{name}»: класс {"ЕСТЬ" if t in S else "ПОТЕРЯН"}', flush=True)
            assert t in S, 'известная точка выброшена решетом — ошибка'

    # 256 классов вырожденных точек над z=0 (и других известных полных точек) с учётом знаков
    def sign_classes(self, z, roots):
        """roots: dict λ -> рациональный корень y_λ ≥ 0 (для всех λ ∈ cellset). Возвращает список (имя, {i: класс})."""
        base = {}
        for i, f in enumerate(self.F):
            y = prod(roots[t] for t in f.T)
            base[i] = f.decompose(f.point(z, y))
        out = []
        for eps in itertools.product((0, 1), repeat=len(self.cellset)):
            cl = {}
            for i, f in enumerate(self.F):
                sg = sum(eps[self.cidx[t]] for t in f.T) % 2
                cl[i] = f.neg_class(base[i]) if sg else base[i]
            out.append((f'z={z},ε={"".join(map(str,eps))}', cl))
        return out
