# Спаривание Касселса–Тейта на Sel^2(E/k) методом Т. Фишера через бинарные квартики.
# Claude, 2026-09-12.   БЕЗ требования 2-кручения (полного или какого-либо вообще).
#
# ИСТОЧНИК. T. Fisher, "On binary quartics and the Cassels-Tate pairing", arXiv:2208.14977,
#   Res. Number Theory 8 (2022), Art. 74 (ANTS-XV).  Теорема 3.1 (стр. 3):
#
#   Пусть I, J в k, 4I^3 - J^2 != 0;  E_{I,J}: y^2 = x^3 - 27 I x - 27 J;
#   L = k[phi], phi — корень X^3 - 3 I X + J.
#   Бинарная квартика g = a x^4 + b x^3 z + c x^2 z^2 + d x z^3 + e z^4, инварианты
#       I(g) = 12ae - 3bd + c^2,   J(g) = 72ace - 27ad^2 - 27b^2 e + 9bcd - 2c^3,
#   Sel^2(E_{I,J}/k) = {всюду локально разрешимые квартики с инвариантами I, J} / (собственная эквивалентность).
#   Кубический инвариант  z(g) = (4 a phi + 3 b^2 - 8 a c)/3 в L^*;  групповой закон = умножение z(g) в L^*/L^*2.
#   Гессиан  h = (3b^2-8ac)x^4 + 4(bc-6ad)x^3 z + 2(2c^2-24ae-3bd)x^2 z^2 + 4(cd-6be)x z^3 + (3d^2-8ce)z^4;
#       G = (4 phi g + h)/3,   H = (1/12) d^2G/dx^2 + (2/9)(I - phi^2) z^2,   G(1,0)G = H^2,  z(g) = G(1,0) = H(1,0).
#
#   ТЕОРЕМА 3.1. g1, g2, g3 — всюду локально разрешимые квартики с инвариантами I, J,
#   z(g1)z(g2)z(g3) = m^2 для некоторого m в L^*, и
#       (z(g2) z(g3) / m) H1(x,z) = alpha1(x,z) + beta1(x,z) phi + gamma1(x,z) phi^2,  alpha1,beta1,gamma1 в k[x,z].
#   Для каждого места v выберем x_v, z_v в k_v с g1(x_v,z_v) — квадратом в k_v и gamma1(x_v,z_v) != 0. Тогда при g2(1,0) != 0
#       <[g1],[g2]>_CT = prod_v ( g2(1,0), gamma1(x_v, z_v) )_v .
#   (Замечание 3.2(v): [g1]+[g2]+[g3] = 0; (vi) при E(k)[2] != 0 корень m не единствен — годится любой.
#    Замечание 3.3: вклад места тривиален, если N(v) >= 11, g1 и gamma1 v-целые, v не делит Delta(g1)*content(gamma1),
#    g2(1,0) — v-единица и v не делит 2.)
#
# ЧТО НУЖНО ДОПОЛНИТЕЛЬНО (этого в статье нет): перевод элемента группы Селмера, полученного алгебраическим
# 2-спуском в виде delta в L^*/L^*2, в бинарную квартику. Замечание 3.2(i): это ровно одна коника НАД k.
# Явная конструкция (наша, выводится ниже и проверяется в коде):
#   2-накрытие C_delta задаётся в P^3 с координатами (u : t0 : t1 : t2), t = t0 + t1*Theta + t2*Theta^2 в L:
#       delta * t^2 = x u^2 - Theta u^2,  т.е.  [delta t^2]_2 = 0,  [delta t^2]_1 + u^2 = 0,  x = [delta t^2]_0 / u^2,
#   где [.]_j — коэффициент при Theta^j, Theta = образ X в L = k[X]/f(X), f — кубика правой части y^2 = f(x).
#   Проверка: N(x - Theta) = f(x), N(delta t^2/u^2) = N(delta) (N(t)/u^3)^2, а N(delta) — квадрат, значит y^2 = f(x).
#   Коника {[delta t^2]_2 = 0} в P(L) = P^2 (эквивалентно Tr_{L/k}(delta t^2 / f'(Theta)) = 0) имеет k-точку,
#   т.к. C_delta всюду локально разрешима (Хассе–Минковский). Параметризуем её: t = t(x,z) (квадратично),
#   тогда  g(x,z) = -[delta t(x,z)^2]_1  — бинарная квартика, и y^2 = g(x,z) — модель C_delta.
#   Остаётся привести инварианты к (I, J): I(g) = lam^4 I, J(g) = lam^6 J => lam^2 = J(g) I / (J I(g)) в k^*2,
#   заменяем g на g/lam^2.
#
# ТОЧНОСТЬ. Вся локальная арифметика ТОЧНАЯ: локальные точки (x_v, z_v) берутся из k (не приближения),
# проверка «g1(x_v,z_v) — квадрат в k_v» — точная (валюация + вычет / ideallog), символы Гильберта — точные.
# p-адических приближений в методе Фишера нет вообще (в отличие от метода Касселса в ctp.sage).

import random, functools, time
print = functools.partial(print, flush=True)


# ------------------------------------------------------------------ мелочи для k = QQ или числового поля

def _is_QQ(k):
    return k is QQ


def k_rand(k, bound):
    if _is_QQ(k):
        return QQ(random.randint(-bound, bound))
    w = k.ring_of_integers().basis()
    return sum(ZZ(random.randint(-bound, bound)) * bb for bb in w)


class LocSq:
    """ k_P^*/k_P^*2: точная проверка квадратичности. k = QQ (P — простое число) или числовое поле (P — идеал). """
    def __init__(self, k, P):
        self.k, self.P = k, P
        if _is_QQ(k):
            self.p = ZZ(P); self.kind = 'Q'
            self.name = f"{self.p}"
            return
        self.kind = 'nf'
        self.p = P.smallest_integer()
        self.pi = k.uniformizer(P, others='positive')
        if self.p == 2:
            e = P.ramification_index()
            self.mod = P ^ (2 * e + 1)
            G = self.mod.idealstar(2)
            self.even = [i for i, o in enumerate(G.gens_orders()) if o % 2 == 0]
        else:
            self.rf = P.residue_field()
        self.name = str(P.gens_two())

    def val(self, a):
        if self.kind == 'Q':
            return QQ(a).valuation(self.p)
        return self.k(a).valuation(self.P)

    def is_sq(self, a):
        if self.kind == 'Q':
            a = QQ(a)
            if a == 0:
                return False
            v = a.valuation(self.p)
            u = a / self.p ^ v
            if v % 2:
                return False
            if self.p == 2:
                return (ZZ(u.numerator()) * ZZ(u.denominator())) % 8 == 1
            return kronecker(ZZ(u.numerator()) * ZZ(u.denominator()), self.p) == 1
        a = self.k(a)
        if a == 0:
            return False
        v = a.valuation(self.P)
        if v % 2:
            return False
        u = a / self.pi ^ v
        if self.p == 2:
            lg = self.mod.ideallog(u)
            return not any(ZZ(lg[i]) % 2 for i in self.even)
        return self.rf(u).is_square()


class RealPlace:
    """ вещественное место: точный знак (через вложение в AA). """
    def __init__(self, k, emb=None):
        self.k, self.emb = k, emb
        self.name = 'real' if emb is None else 'real ' + str(emb(k.gen()).n(20))

    def sgn(self, a):
        if self.emb is None:
            return sign(QQ(a))
        s = self.emb(self.k(a))
        return 1 if s > 0 else (-1 if s < 0 else 0)

    def is_sq(self, a):
        return self.sgn(a) > 0


def hilb(k, a, b, pl):
    """ символ Гильберта (a,b)_v в {1,-1}. """
    if isinstance(pl, RealPlace):
        return -1 if (pl.sgn(a) < 0 and pl.sgn(b) < 0) else 1
    if pl.kind == 'Q':
        return ZZ(hilbert_symbol(QQ(a), QQ(b), pl.p))
    return ZZ(k.hilbert_symbol(k(a), k(b), pl.P))


# ------------------------------------------------------------------ этальная алгебра L = k[X]/(cubic)

class EtaleCubic:
    def __init__(self, k, cub):
        self.k = k
        self.Rx = PolynomialRing(k, 'X')
        self.cub = self.Rx(cub)
        assert self.cub.degree() == 3 and self.cub.is_monic()
        self.L = self.Rx.quotient(self.cub, 'ph')
        self.phi = self.L.gen()
        fac = list(self.cub.factor())
        assert all(m == 1 for _, m in fac), "кубика не сепарабельна"
        self.facs = [fa for fa, _ in fac]
        self.comps = []          # (поле F, корень r в F, абсолютное поле Fa, iso Fa->F, iso F->Fa)
        for idx, fa in enumerate(self.facs):
            if fa.degree() == 1:
                F = k; r = -fa[0]
                self.comps.append({'F': F, 'r': r, 'deg': 1})
            else:
                if _is_QQ(k):
                    F = NumberField(fa, 'z%d' % idx)
                    self.comps.append({'F': F, 'r': F.gen(), 'deg': fa.degree(), 'Fa': F,
                                       'toA': (lambda z: z), 'frA': (lambda z: z)})
                else:
                    F = k.extension(fa, 'z%d' % idx)
                    Fa = F.absolute_field('aa%d' % idx)
                    frA, toA = Fa.structure()      # frA: Fa -> F,  toA: F -> Fa
                    self.comps.append({'F': F, 'r': F.gen(), 'deg': fa.degree(), 'Fa': Fa,
                                       'toA': toA, 'frA': frA})
        # матрица отображения L -> prod F_i в координатах над k
        rows = []
        for j in range(3):
            v = []
            for c in self.comps:
                v += self._kcoords(c, c['r'] ^ j)
            rows.append(v)
        self.Mcomp = matrix(k, rows)     # строка j = образ Theta^j
        assert self.Mcomp.is_invertible()

    def _kcoords(self, c, el):
        if c['deg'] == 1:
            return [self.k(el)]
        return [self.k(t) for t in list(c['F'](el))]

    def elt(self, poly_or_L):
        return self.L(poly_or_L)

    def comp(self, z, i):
        """ i-я компонента элемента z из L """
        c = self.comps[i]
        return self.L(z).lift()(c['r'])

    def from_comps(self, vals):
        """ элемент L с заданными компонентами """
        target = []
        for c, v in zip(self.comps, vals):
            target += self._kcoords(c, v)
        sol = self.Mcomp.solve_left(vector(self.k, target))
        z = self.L(self.Rx(list(sol)))
        for i, v in enumerate(vals):
            assert self.comp(z, i) == self.comps[i]['F'](v), "CRT сбой"
        return z

    def inv(self, z):
        g, u, _ = self.Rx(self.L(z).lift()).xgcd(self.cub)
        assert g.degree() == 0, "элемент L не обратим"
        return self.L(u / g[0])

    def _comp_is_square(self, c, el, root=False):
        if c['deg'] == 1:
            el = self.k(el)
            if _is_QQ(self.k):
                ok = QQ(el).is_square()
                return (ok, QQ(el).sqrt()) if root else ok
            ok = el.is_square()
            return (ok, el.sqrt()) if root else ok
        a = c['toA'](c['F'](el))
        ok = a.is_square()
        if not root:
            return ok
        return (ok, c['frA'](a.sqrt())) if ok else (False, None)

    def is_unit(self, z):
        return all(self.comp(z, i) != 0 for i in range(len(self.comps)))

    def is_square(self, z):
        return all(self._comp_is_square(c, self.comp(z, i)) for i, c in enumerate(self.comps))

    def sqrt(self, z):
        vals = []
        for i, c in enumerate(self.comps):
            ok, r = self._comp_is_square(c, self.comp(z, i), root=True)
            if not ok:
                return None
            vals.append(r)
        m = self.from_comps(vals)
        assert m * m == self.L(z), "квадратный корень в L неверен"
        return m

    def coeffs(self, z):
        """ (c0, c1, c2) в базисе 1, phi, phi^2 """
        l = self.L(z).lift()
        return [self.k(l[j]) for j in range(3)]


# ------------------------------------------------------------------ инварианты бинарных квартик

def quartic_I(g):
    a, b, c, d, e = g
    return 12 * a * e - 3 * b * d + c ^ 2


def quartic_J(g):
    a, b, c, d, e = g
    return 72 * a * c * e - 27 * a * d ^ 2 - 27 * b ^ 2 * e + 9 * b * c * d - 2 * c ^ 3


def quartic_hessian(g):
    a, b, c, d, e = g
    return [3 * b ^ 2 - 8 * a * c, 4 * (b * c - 6 * a * d), 2 * (2 * c ^ 2 - 24 * a * e - 3 * b * d),
            4 * (c * d - 6 * b * e), 3 * d ^ 2 - 8 * c * e]


def quartic_eval(g, x, z):
    a, b, c, d, e = g
    return a * x ^ 4 + b * x ^ 3 * z + c * x ^ 2 * z ^ 2 + d * x * z ^ 3 + e * z ^ 4


# ------------------------------------------------------------------ основной класс

class FisherCTP:
    def __init__(self, k, I, J, verbose=False):
        self.k = k; self.I = k(I); self.J = k(J)
        assert 4 * self.I ^ 3 - self.J ^ 2 != 0
        X = PolynomialRing(k, 'X').gen()
        self.E = EtaleCubic(k, X ^ 3 - 3 * self.I * X + self.J)
        self.L = self.E.L; self.phi = self.E.phi
        self.disc = 16 * (4 * self.I ^ 3 - self.J ^ 2) / 27
        self.verbose = verbose
        self.conic_cache = {}

    # --------- z(g), H(g) ---------
    def z_inv(self, g):
        a, b, c, d, e = g
        return (4 * a * self.phi + 3 * b ^ 2 - 8 * a * c) / 3

    def H_form(self, g):
        """ H(x,z) = H0 x^2 + H1 x z + H2 z^2, коэффициенты в L """
        h = quartic_hessian(g)
        G = [(4 * self.phi * g[j] + h[j]) / 3 for j in range(5)]
        H0 = G[0]
        H1 = G[1] / 2
        H2 = G[2] / 6 + (2 * (self.I - self.phi ^ 2)) / 9
        # контроль тождества G(1,0) G(x,z) = H(x,z)^2 (сравнение как многочленов от x,z)
        Rq = PolynomialRing(self.L, 'xx,zz'); xx, zz = Rq.gens()
        Gp = sum(G[j] * xx ^ (4 - j) * zz ^ j for j in range(5))
        Hp = H0 * xx ^ 2 + H1 * xx * zz + H2 * zz ^ 2
        assert G[0] * Gp == Hp ^ 2, "нарушено тождество G(1,0)G = H^2"
        return [H0, H1, H2]

    def _make_z_unit(self, g):
        """ Фишер (стр. 3): «By a change of coordinates we may assume that z(g) is a unit in L».
            z(g) = (4 a phi + 3b^2 - 8ac)/3 может оказаться делителем нуля; тогда заменяем g на
            собственно эквивалентную (сдвиги/обращение/перекос), пока z(g) не станет единицей. """
        k = self.k
        if self.E.is_unit(self.z_inv(g)):
            return g
        moves = [[[0, 1], [1, 0]]]
        for n in (1, -1, 2, -2, 3, -3, 5, -5):
            moves.append([[1, k(n)], [0, 1]])
            moves.append([[1, 0], [k(n), 1]])
        for mt in moves:
            try:
                gn = self.gl2(g, mt)
            except Exception:
                continue
            if self.E.is_unit(self.z_inv(gn)):
                return gn
        for _ in range(400):
            mt = [[k(1) + k_rand(k, 1), k_rand(k, 3)], [k_rand(k, 3), k(1) + k_rand(k, 1)]]
            if matrix(k, mt).determinant() == 0:
                continue
            try:
                gn = self.gl2(g, mt)
            except Exception:
                continue
            if self.E.is_unit(self.z_inv(gn)):
                return gn
        raise RuntimeError("не удалось сделать z(g) единицей в L")

    def check_quartic(self, g):
        assert quartic_I(g) == self.I, f"I(g) = {quartic_I(g)} != {self.I}"
        assert quartic_J(g) == self.J, f"J(g) = {quartic_J(g)} != {self.J}"

    # --------- 2-накрытие из delta и его квартика ---------
    def covering_forms(self, delta):
        """ [delta t^2]_0, [delta t^2]_1, [delta t^2]_2 как квадратичные формы от (t0,t1,t2) над k """
        k = self.k
        # ВНИМАНИЕ: корни x-координат 2-кручения E_{I,J} равны Theta = -3 phi (проверка:
        #   (-3phi)^3 - 27 I (-3phi) - 27 J = -27(phi^3 - 3 I phi + J) = 0),
        # поэтому базис для разложения [.]_j — (1, Theta, Theta^2), Theta = -3 phi.
        # БАЗИС: берём k-базис L, СОГЛАСОВАННЫЙ с разложением L = prod F_i (через идемпотенты).
        # Тогда c2(t) = Tr_{L/k}(delta t^2 / f'(Theta)) блочно-диагональна, её коэффициенты —
        # S-единицы (delta и f'(e_i) поддержаны на S), и коника решается. В базисе (1,Theta,Theta^2)
        # коэффициенты получаются огромными и задача становится неподъёмной.
        Rt = PolynomialRing(k, 't0,t1,t2'); tg = Rt.gens()
        RtX = PolynomialRing(Rt, 'X')
        cub = RtX([Rt(c) for c in self.E.cub.list()])
        Lt = RtX.quotient(cub, 'th'); th = Lt.gen()
        dl = self.E.coeffs(delta)
        dlt = Lt(RtX([Rt(c) for c in dl]))
        BB = self._L_basis()
        tt = sum(tg[j] * Lt(RtX([Rt(c) for c in BB[j]])) for j in range(3))
        pr = (dlt * tt * tt).lift()
        d = [Rt(pr[j]) for j in range(3)]
        return [d[0], -d[1] / 3, d[2] / 9]      # коэффициенты в базисе (1, Theta, Theta^2)

    def _L_basis(self):
        """ k-базис L, согласованный с разложением на компоненты-поля; возвращает (1,phi,phi^2)-координаты """
        if hasattr(self, '_Lbas'):
            return self._Lbas
        out = []
        for i, c in enumerate(self.E.comps):
            for m in range(c['deg']):
                vals = []
                for i2, c2 in enumerate(self.E.comps):
                    vals.append(c['r'] ^ m if i2 == i else c2['F'](0))
                out.append(self.E.coeffs(self.E.from_comps(vals)))
        assert len(out) == 3
        assert matrix(self.k, out).is_invertible()
        self._Lbas = out
        return out

    def quartic_from_delta(self, delta, tag=None):
        """ бинарная квартика с инвариантами (I,J), представляющая класс delta в L^*/L^*2 """
        k = self.k
        c0, c1, c2 = self.covering_forms(delta)
        Rt = c2.parent(); tv = Rt.gens()
        M = matrix(k, 3, 3, lambda i, j: (c2.coefficient({tv[i]: 2}) if i == j
                                          else c2.coefficient({tv[i]: 1, tv[j]: 1}) / 2))
        assert M.determinant() != 0, "вырожденная коника"
        P = self.solve_conic(M, tag=tag)
        P = self.reduce_conic_point(M, P)
        tpar, Rp = self._param_from_point(M, P)
        xx, zz = Rp.gens()
        assert sum(M[i][j] * tpar[i] * tpar[j] for i in range(3) for j in range(3)) == 0, "параметризация не на конике"
        sub = {tv[i]: tpar[i] for i in range(3)}
        gpoly = -Rp(c1.subs(sub))
        g = [k(gpoly.coefficient({xx: 4 - j, zz: j})) for j in range(5)]
        assert gpoly == sum(g[j] * xx ^ (4 - j) * zz ^ j for j in range(5)), "не бинарная квартика"
        Ig, Jg = quartic_I(g), quartic_J(g)
        assert Ig != 0 or Jg != 0
        # масштабирование до инвариантов (I, J)
        if self.I != 0 and Ig != 0:
            lam2 = Jg * self.I / (self.J * Ig) if (self.J != 0 and Jg != 0) else None
        else:
            lam2 = None
        if lam2 is None:
            raise RuntimeError("вырожденный случай I=0 или J=0 — не поддержан")
        assert (lam2 ^ 2 == Ig / self.I), f"несогласованные инварианты: lam^4 != I(g)/I"
        g = [c / lam2 for c in g]
        self.check_quartic(g)
        g = self.reduce_quartic(g)
        g = self._make_z_unit(g)
        # контроль: z(g) должен лежать в том же классе L^*/L^*2, что и delta
        rat = self.z_inv(g) * self.L(delta)
        assert self.E.is_square(rat), "z(g) не в классе delta"
        return g

    # --------- приведение квартик (сохраняет инварианты) ---------
    # g |-> (det gam)^(-2) * g о gam  для любого gam в GL_2(k): I,J сохраняются, класс в Sel^2 тот же
    # (собственная эквивалентность с lam = 1/det gam, lam*det gam = 1).
    def gl2(self, g, mat):
        k = self.k
        al, ga, be, de = k(mat[0][0]), k(mat[0][1]), k(mat[1][0]), k(mat[1][1])
        dt = al * de - be * ga
        assert dt != 0
        Rp = PolynomialRing(k, 'x,z'); xx, zz = Rp.gens()
        gp = quartic_eval(g, al * xx + ga * zz, be * xx + de * zz) / dt ^ 2
        return [k(gp.coefficient({xx: 4 - j, zz: j})) for j in range(5)]

    def _den(self, c):
        if c == 0:
            return ZZ(1)
        if _is_QQ(self.k):
            return QQ(c).denominator()
        return lcm([QQ(t).denominator() for t in list(self.k(c))])

    def _msr(self, g):
        """ размер квартики в битах (точно, без переполнения плавающей точки) """
        tot = 0
        for c in g:
            if c == 0:
                continue
            co = [QQ(c)] if _is_QQ(self.k) else [QQ(t) for t in list(self.k(c))]
            for t in co:
                if t == 0:
                    continue
                tot += ZZ(t.numerator()).nbits() + ZZ(t.denominator()).nbits()
        return tot

    # --- вложения k -> R/C с произвольной точностью (по координатам, без переполнения) ---
    def _emb_list(self):
        """ список пар (тип, значение sigma(r)) для образующей r поля k; для QQ — [('R', None)] """
        if hasattr(self, '_embl'):
            return self._embl
        k = self.k
        if _is_QQ(k):
            self._embl = [('R', None)]
        else:
            D = QQ(k.gen() ^ 2)
            assert k.degree() == 2 and k.gen() ^ 2 == D, "поддержаны только квадратичные k и QQ"
            self._embl = [('R', 1), ('R', -1)] if D > 0 else [('C', 1)]
            self._D = D
        return self._embl

    def _sig(self, c, j, prec):
        """ j-е вложение элемента c в RealField/ComplexField(prec) """
        k = self.k
        if _is_QQ(k):
            return RealField(prec)(QQ(c))
        t, sg = self._emb_list()[j]
        c0, c1 = [QQ(u) for u in list(k(c))]
        if t == 'R':
            R = RealField(prec)
            return R(c0) + R(c1) * sg * R(self._D).sqrt()
        C = ComplexField(prec)
        return C(c0) + C(c1) * C(-self._D).sqrt() * C(0, 1)

    def _from_sigmas(self, T, prec, tolbits=60, abs_err=None):
        """ элемент k, приближающий заданные значения T[j] при вложениях.
            Округление c0, c1 — с ОБЩЕЙ абсолютной погрешностью, достаточной и для наименьшего |T_j|
            (иначе при сильно разных |T_j| происходит катастрофическое сокращение). """
        k = self.k
        R = RealField(prec)
        if _is_QQ(k):
            u = R(T[0])
            if u == 0:
                return QQ(0)
            return QQ(u.nearby_rational(max_error=abs(u) / 2 ^ tolbits))
        mags = [abs(R(t)) for t in T if R(t) != 0]
        if not mags:
            return k(0)
        extra = ZZ(ceil(abs(RR(max(mags) / min(mags)).log(2)))) if min(mags) > 0 else 0
        tb = tolbits + extra
        if self._emb_list()[0][0] == 'R':
            t1, t2 = R(T[0]), R(T[1])
            c0 = (t1 + t2) / 2
            c1 = (t1 - t2) / (2 * R(self._D).sqrt())
        else:
            z = ComplexField(prec)(T[0])
            c0 = z.real()
            c1 = z.imag() / R(-self._D).sqrt()
        scale = max(abs(c0), abs(c1))
        if scale == 0:
            return k(0)
        err = scale / 2 ^ tb if abs_err is None else R(abs_err)
        def rat(u):
            u = R(u)
            if abs(u) < err:
                return QQ(0)
            return QQ(u.nearby_rational(max_error=err))
        return k(rat(c0)) + k(rat(c1)) * k.gen()

    def _prec_range(self, vals, extra=300):
        """ точность, достаточная для работы со значениями разного порядка (по динамическому диапазону) """
        R = RealField(80)
        ms = []
        for v in vals:
            a = abs(R(v))
            if a != 0:
                ms.append(a)
        if not ms:
            return 200
        hi, lo = max(ms), min(ms)
        try:
            rng = ZZ(ceil(abs(RR(hi / lo).log(2)))) + ZZ(ceil(abs(RR(hi).log(2)))) + ZZ(ceil(abs(RR(lo).log(2))))
        except Exception:
            rng = 1000
        return int(min(3 * 10 ^ 4, max(200, rng + extra)))

    def _prec_for(self, g):
        vals = []
        for c in g:
            if c == 0:
                continue
            if _is_QQ(self.k):
                vals.append(QQ(c))
            else:
                vals += [QQ(t) for t in list(self.k(c)) if t != 0]
        return self._prec_range(vals)

    def _balance(self, g, iters=60):
        """ приведение: точное «депрессирование» (b -> 0 сдвигом x -> x - b/(4a), det = 1) и
            балансировка diag(s,1). Масштаб s подбирается из симметрических функций корней:
            |sigma(s)| ~ |sigma(c/a)|^(1/2), |sigma(d/a)|^(1/3), |sigma(e/a)|^(1/4)
            (это характерные величины корней депрессированной квартики). Оба преобразования сохраняют I, J. """
        k = self.k
        for it in range(iters):
            if g[0] == 0:
                break
            if g[1] != 0:
                # сдвиг: точный (b -> 0) может раздуть знаменатели, поэтому берём лучший из
                # точного и округлённых вариантов по битовой мере
                n0 = -g[1] / (4 * g[0])
                bestt = g; bmt = self._msr(g)
                for nn in (n0, self._round(n0), self._round(n0) + 1, self._round(n0) - 1):
                    try:
                        gt = self.gl2(g, [[1, nn], [0, 1]])
                    except Exception:
                        continue
                    mt = self._msr(gt)
                    if mt < bmt:
                        bestt, bmt = gt, mt
                g = bestt
            prec = self._prec_for(g)
            embs = self._emb_list()
            best = g; bm = self._msr(g)
            for (idx, p) in ((2, 2), (3, 3), (4, 4)):
                if g[idx] == 0:
                    continue
                ratio = g[idx] / g[0]
                try:
                    T = [abs(self._sig(ratio, j, prec)) ^ (QQ(1) / p) for j in range(len(embs))]
                    s = self._from_sigmas(T, prec)
                except Exception:
                    continue
                if s == 0:
                    continue
                for ss in (s, 1 / s):
                    try:
                        gn = self.gl2(g, [[ss, 0], [0, 1]])
                    except Exception:
                        continue
                    m = self._msr(gn)
                    if m < bm:
                        best, bm = gn, m
            gr = list(reversed(g))
            if self._msr(gr) < bm:
                best, bm = gr, self._msr(gr)
            if bm >= self._msr(g):
                break
            g = best
        return g

    def _round(self, c):
        """ ближайший целый элемент O_k к c """
        k = self.k
        if _is_QQ(k):
            return QQ(round(QQ(c)))
        B = k.ring_of_integers().basis()
        M = matrix(QQ, [list(k(b)) for b in B])
        co = vector(QQ, list(k(c))) * M.inverse()
        return sum(ZZ(round(t)) * k(b) for t, b in zip(co, B))

    def _covariant_reduce(self, g, iters=25, prec_cap=40000):
        """ Приведение по ковариантной точке (Юлиа / Stoll–Cremona), обобщённое на квадратичное k.
            Для каждого вложения sigma поля k: корни alpha_i многочлена sigma(g)(x,1) в C,
            ковариантная точка z0 = x0 + i y0, где
                x0 = среднее Re(alpha_i),   y0 = sqrt( sum |alpha_i - x0|^2 / 4 )
            (это минимум формы Юлиа sum |z - alpha_i|^2 / Im z).
            Подбираем c, s в k с sigma_j(c) ~ x0^(j), sigma_j(s) ~ y0^(j) и применяем
            gamma = [[s, c],[0,1]] (то есть x -> s x + c z), что сохраняет I, J. """
        k = self.k
        embs = self._emb_list()
        for it in range(iters):
            m0 = self._msr(g)
            prec = self._prec_for(g)
            Cf = ComplexField(prec)
            X0, Y0 = [], []
            ok = True
            for j in range(len(embs)):
                Rp = PolynomialRing(Cf, 'x'); xx = Rp.gen()
                try:
                    gc = sum(Cf(self._sig(g[t], j, prec)) * xx ^ (4 - t) for t in range(5))
                    rts = gc.roots(Cf, multiplicities=False)
                except Exception:
                    ok = False; break
                if len(rts) < 3:
                    ok = False; break
                x0 = sum(r.real() for r in rts) / len(rts)
                y0 = (sum(abs(r - x0) ^ 2 for r in rts) / 4).sqrt()
                if y0 == 0:
                    ok = False; break
                X0.append(x0); Y0.append(y0)
            if not ok:
                break
            gn = None; mn = m0
            R80 = RealField(80)
            ymin = min([R80(abs(y)) for y in Y0])
            for tb in (20, 60, 160):
                try:
                    c = self._from_sigmas(X0, prec, tolbits=tb, abs_err=ymin / 2 ^ tb)
                    s = self._from_sigmas(Y0, prec, tolbits=tb)
                    if s == 0:
                        continue
                    cand = self.gl2(g, [[s, c], [0, 1]])
                except Exception:
                    continue
                mc = self._msr(cand)
                if mc < mn:
                    gn, mn = cand, mc
            if gn is None:
                gr = list(reversed(g))
                if self._msr(gr) < m0:
                    g = gr; continue
                break
            g = gn
        return g

    def reduce_quartic(self, g, rounds=40):
        k = self.k
        g = self._balance(g)
        g = self._covariant_reduce(g)
        g = self._balance(g)
        if g[0] == 0 or self._msr(list(reversed(g))) < self._msr(g):
            g = self._balance(list(reversed(g)))
        best = list(g); bm = self._msr(best)
        units = []
        if not _is_QQ(k):
            U = k.unit_group()
            for u in U.gens():
                uu = k(u)
                if uu not in (1, -1):
                    units += [uu, 1 / uu]
        cands_r = [k(2), k(1) / 2, k(3), k(1) / 3, k(5), k(1) / 5, k(7), k(1) / 7] + units
        for it in range(rounds):
            improved = False
            moves = [[[0, 1], [1, 0]]]
            a, b = best[0], best[1]
            if a != 0:
                n0 = self._round(-b / (4 * a))
                moves += [[[1, n0 + t], [0, 1]] for t in (-1, 0, 1)]
            moves += [[[1, k(t)], [0, 1]] for t in (-2, -1, 1, 2)]
            moves += [[[1, 0], [k(t), 1]] for t in (-2, -1, 1, 2)]
            e, d = best[4], best[3]
            if e != 0:
                n1 = self._round(-d / (4 * e))
                moves += [[[1, 0], [n1 + t, 1]] for t in (-1, 0, 1)]
            moves += [[[r, 0], [0, 1]] for r in cands_r]
            for mt in moves:
                try:
                    gp = self.gl2(best, mt)
                except Exception:
                    continue
                m = self._msr(gp)
                if m < bm - RR(1e-9):
                    best, bm = gp, m; improved = True
            if not improved:
                break
        assert quartic_I(best) == quartic_I(g) and quartic_J(best) == quartic_J(g), "приведение сломало инварианты"
        return best

    # --------- коники над k ---------
    def diagonalize(self, M):
        """ D (список), T: T^t M T = diag(D) """
        k = self.k; n = M.nrows()
        M0 = matrix(k, M)
        M = matrix(k, M); T = identity_matrix(k, n)
        D = []
        for s in range(n):
            # ищем ненулевой диагональный элемент в блоке s..n-1
            piv = None
            for i in range(s, n):
                if M[i, i] != 0:
                    piv = i; break
            if piv is None:
                found = None
                for i in range(s, n):
                    for j in range(i + 1, n):
                        if M[i, j] != 0:
                            found = (i, j); break
                    if found:
                        break
                if found is None:
                    D += [k(0)] * (n - s); break
                i, j = found
                M.add_multiple_of_row(i, j, 1); M.add_multiple_of_column(i, j, 1)
                T.add_multiple_of_column(i, j, 1)
                piv = i
            if piv != s:
                M.swap_rows(piv, s); M.swap_columns(piv, s); T.swap_columns(piv, s)
            d = M[s, s]
            for i in range(s + 1, n):
                if M[i, s] != 0:
                    c = -M[i, s] / d
                    M.add_multiple_of_row(i, s, c); M.add_multiple_of_column(i, s, c)
                    T.add_multiple_of_column(i, s, c)
            D.append(d)
        assert T.transpose() * M0 * T == diagonal_matrix(k, D), "диагонализация неверна"
        return D, T

    def sq_reduce(self, a):
        """ a = a' * s^2 с «маленьким» a' """
        k = self.k
        if _is_QQ(k):
            a = QQ(a)
            num = ZZ(a.numerator() * a.denominator())
            sp = num.squarefree_part()
            s = (a / sp).sqrt()
            return QQ(sp), s
        a = k(a)
        co = [QQ(t) for t in list(a)]
        den = lcm([t.denominator() for t in co])
        num = gcd([ZZ(t * den) for t in co])
        c = QQ(num) / den                 # a = c * (целый примитивный)
        cn = ZZ(c.numerator() * c.denominator())
        sp = cn.squarefree_part()
        s2 = c / sp
        s = QQ(s2).sqrt()
        return a / s ^ 2, k(s)

    def _param_from_point(self, M, P):
        """ параметризация коники {v^t M v = 0} из точки P: t(x,z) = Q(w)P - 2B(P,w)w, w = x v1 + z v2 """
        k = self.k
        bas = [vector(k, [1, 0, 0]), vector(k, [0, 1, 0]), vector(k, [0, 0, 1])]
        comp = [b for b in bas if matrix(k, [P, b]).rank() == 2]
        v1 = v2 = None
        for i in range(len(comp)):
            for j in range(i + 1, len(comp)):
                if matrix(k, [P, comp[i], comp[j]]).rank() == 3:
                    v1, v2 = comp[i], comp[j]; break
            if v1 is not None:
                break
        assert v1 is not None
        Rp = PolynomialRing(k, 'x,z'); xx, zz = Rp.gens()
        w = [xx * v1[i] + zz * v2[i] for i in range(3)]
        Qw = sum(M[i][j] * w[i] * w[j] for i in range(3) for j in range(3))
        BPw = sum(M[i][j] * P[i] * w[j] for i in range(3) for j in range(3))
        return [Qw * P[i] - 2 * BPw * w[i] for i in range(3)], Rp

    def _sl2_reduce_pd_quartic(self, h):
        """ h — положительно определённая бинарная квартика над Q ([h0..h4]).
            Возвращает gamma в SL_2(Z) такую, что h о gamma приведена (ковариантная точка Юлиа
            z0 = x0 + i y0, x0 = mean Re(корни), y0 = sqrt(sum |alpha - x0|^2 / 4), приводится
            стандартным алгоритмом в фундаментальную область). """
        # точность определяется ДИНАМИЧЕСКИМ ДИАПАЗОНОМ значений коэффициентов (а не их битовой длиной:
        # коэффициенты могут быть дробями с огромными числителем и знаменателем, но умеренными значениями)
        prec = self._prec_range([QQ(c) for c in h if c != 0])
        Cf = ComplexField(prec)
        Rp = PolynomialRing(Cf, 'x'); xx = Rp.gen()
        hp = sum(Cf(QQ(h[t])) * xx ^ (4 - t) for t in range(5))
        rts = hp.roots(Cf, multiplicities=False)
        if len(rts) < 4:
            return matrix(ZZ, [[1, 0], [0, 1]])
        x0 = sum(r.real() for r in rts) / 4
        y0 = (sum(abs(r - x0) ^ 2 for r in rts) / 4).sqrt()
        if y0 <= 0:
            return matrix(ZZ, [[1, 0], [0, 1]])
        z = Cf(x0, y0)
        gm = matrix(ZZ, [[1, 0], [0, 1]])
        for _ in range(400):
            nre = z.real()
            n = ZZ(nre.round()) if hasattr(nre, 'round') else ZZ(round(RR(nre)))
            if n != 0:
                z = z - n
                gm = matrix(ZZ, [[1, -n], [0, 1]]) * gm
            if abs(z) < 1:
                z = -1 / z
                gm = matrix(ZZ, [[0, -1], [1, 0]]) * gm
            else:
                break
        return gm.inverse()

    def reduce_conic_point(self, M, P, iters=6):
        """ приводит точку P коники: строит параметризацию, сворачивает её положительно определённой
            квартикой h(x,z) = sum_i |sigma(t_i(x,z))|^2 (рациональная!), приводит h в SL_2(Z)
            и берёт P' = t(gamma(1,0)). Именно размер P определяет знаменатели итоговой квартики. """
        k = self.k
        D = QQ(k.gen() ^ 2) if not _is_QQ(k) else QQ(0)
        def size(v):
            return self._msr(list(v))
        best = vector(k, P); bs = size(best)
        for it in range(iters):
            try:
                tp, Rp = self._param_from_point(M, best)
                xx, zz = Rp.gens()
                h = [QQ(0)] * 5
                for ti in tp:
                    if _is_QQ(k):
                        parts = [ti]
                    else:
                        parts = None
                    if _is_QQ(k):
                        sq = ti * ti
                        for j in range(5):
                            h[j] += QQ(sq.coefficient({xx: 4 - j, zz: j}))
                    else:
                        # t_i = u + v r,  sum_sigma |sigma(t_i)|^2 = 2(u^2 + |D| v^2)
                        u = Rp(0); v = Rp(0)
                        for mono, co in zip(ti.monomials(), ti.coefficients()):
                            c0, c1 = [QQ(s) for s in list(k(co))]
                            u += c0 * mono; v += c1 * mono
                        sq = u * u + abs(D) * v * v
                        for j in range(5):
                            h[j] += QQ(sq.coefficient({xx: 4 - j, zz: j}))
                if all(c == 0 for c in h):
                    break
                gm = self._sl2_reduce_pd_quartic(h)
                a_, c_ = k(gm[0][0]), k(gm[1][0])
                Pn = vector(k, [ti.subs({xx: a_, zz: c_}) for ti in tp])
                if Pn == 0:
                    break
                dd = lcm([self._den(t) for t in Pn])
                Pn = vector(k, [t * dd for t in Pn])
                gg = gcd([ZZ(s) for t in Pn for s in ([QQ(t)] if _is_QQ(k) else list(k(t))) if s != 0])
                if gg > 1:
                    Pn = vector(k, [t / gg for t in Pn])
                assert (Pn * M * Pn) == 0, "приведённая точка не на конике"
                if size(Pn) < bs:
                    best, bs = Pn, size(Pn)
                else:
                    break
            except Exception:
                break
        return best

    def reduce_form(self, M):
        """ U в GL_3(k) с «малой» M' = U^t M U (LLL по положительно определённой мажоранте).
            Pos(v) = sum_j ( |sigma_j(v)|^2 + |Mhat_j sigma_j(v)|^2 ),  Mhat = M/масштаб.
            Тогда |B'(u_i,u_j)| ограничено Pos-нормами (Коши–Буняковский), т.е. M' мала.
            Без этого шага диагонализация даёт коэффициенты чудовищного размера, их факторизация
            (для приведения по модулю квадратов) и решение коники становятся неподъёмными. """
        k = self.k
        n = 3
        if _is_QQ(k):
            omegas = [QQ(1)]; deg = 1
        else:
            omegas = [k(b) for b in k.ring_of_integers().basis()]; deg = 2
        embs = self._emb_list()
        vals = [M[i][j] for i in range(n) for j in range(n) if M[i][j] != 0]
        try:
            Dg, T = self.diagonalize(M)
            if any(d == 0 for d in Dg):
                return identity_matrix(k, n)
            Ti = T.inverse()
        except Exception:
            return identity_matrix(k, n)
        allv = vals + list(Dg) + Ti.list()
        prec = self._prec_range([QQ(t) for c in allv if c != 0
                                 for t in ([QQ(c)] if _is_QQ(k) else list(k(c))) if t != 0])
        prec = int(min(prec, 8000))
        R = RealField(prec); Cf = ComplexField(prec)
        # мажоранта: B(u,v) = sum_m D_m L_m(u) L_m(v), L = T^{-1};  H(u) = sum_m |D_m| |L_m(u)|^2
        aD = [[abs(Cf(self._sig(Dg[m], j, prec))) for m in range(n)] for j in range(len(embs))]
        tiv = [[[Cf(self._sig(Ti[m][i], j, prec)) for i in range(n)] for m in range(n)]
               for j in range(len(embs))]
        omv = [[Cf(self._sig(omegas[l], j, prec)) for l in range(deg)] for j in range(len(embs))]
        dim = n * deg
        G = matrix(R, dim, dim)
        for i in range(n):
            for l in range(deg):
                for i2 in range(n):
                    for l2 in range(deg):
                        s = R(0)
                        for j in range(len(embs)):
                            for m in range(n):
                                u1 = tiv[j][m][i] * omv[j][l]
                                u2 = tiv[j][m][i2] * omv[j][l2]
                                s += aD[j][m] * (u1.conjugate() * u2).real()
                        G[i * deg + l, i2 * deg + l2] = s
        G = (G + G.transpose()) / 2
        mx = max([abs(G[a][b]) for a in range(dim) for b in range(dim)])
        if mx == 0:
            return identity_matrix(k, n)
        Gi = matrix(ZZ, dim, dim, lambda a, b: ZZ((G[a][b] / mx * 2 ^ 200).round()))
        try:
            U6 = Gi.LLL_gram()
        except Exception:
            return identity_matrix(k, n)
        cols = []
        for c in range(dim):
            v = [sum(U6[i * deg + l, c] * omegas[l] for l in range(deg)) for i in range(n)]
            if any(t != 0 for t in v):
                cols.append(vector(k, v))
        chosen = []
        for v in cols:
            if matrix(k, chosen + [v]).rank() == len(chosen) + 1:
                chosen.append(v)
            if len(chosen) == n:
                break
        if len(chosen) < n:
            return identity_matrix(k, n)
        U = matrix(k, chosen).transpose()
        assert U.is_invertible()
        return U

    def _best_tp(self, t, abc):
        """ (tp, w) с t = tp*w^2 и наименьшей |N(tp)|. Кандидаты: приведённый представитель из
            selmer_space и все +-произведения подмножеств коэффициентов {A,B,C} диагональной формы
            (класс t = -AB совпадает с одним из них, т.к. ABC = det с точностью до квадратов). """
        k = self.k
        cands = []
        try:
            r0, s0 = self.small_rep(t) if not _is_QQ(k) else self.sq_reduce(t)
            cands.append((r0, s0))
        except Exception:
            pass
        prods = [k(1)]
        for u in abc:
            prods = prods + [p * u for p in prods]
        for cand in [sg * p for p in prods for sg in (1, -1)]:
            if cand == 0:
                continue
            q = t / cand
            try:
                if (QQ(q).is_square() if _is_QQ(k) else q.is_square()):
                    cands.append((cand, QQ(q).sqrt() if _is_QQ(k) else q.sqrt()))
            except Exception:
                pass
        if not cands:
            return (t, k(1))
        def sz(c):
            nm = QQ(c) if _is_QQ(k) else QQ(k(c).norm())
            return abs(nm.numerator()) * abs(nm.denominator())
        best = min(cands, key=lambda p: sz(p[0]))
        assert best[0] * best[1] ^ 2 == t, "_best_tp: разложение неверно"
        return best

    def small_rep(self, a):
        """ a = a' * s^2 с МАЛЫМ представителем a' класса a в k*/k*^2.
            Через k.selmer_space(S,2) с S = простые, делящие (a) (плюс над 2): образующие
            группы Сельмера уже приведены, поэтому a' — произведение нескольких малых элементов.
            Без этого решения коник (is_norm) получаются астрономическими. """
        k = self.k
        a0, s0 = self.sq_reduce(a)
        if _is_QQ(k):
            return a0, s0
        if not hasattr(self, '_srcache'):
            self._srcache = {}
        # (1) сначала пробуем ИЗВЕСТНОЕ множество простых (плохие простые кривой + малые):
        # коэффициенты коники почти всегда лежат в k(S,2), а факторизовать их нормы (100+ цифр) нельзя
        BS = getattr(self, 'base_S', None)
        if BS:
            try:
                if not hasattr(self, '_bsdata'):
                    self._bsdata = k.selmer_space(list(BS), 2)
                V, gens, fromV, toV = self._bsdata
                rep = k(fromV(toV(a0)))
                q = a0 / rep
                if q != 0 and q.is_square():
                    return rep, s0 * q.sqrt()
            except Exception:
                pass
        nn = QQ(k(a0).norm())
        if (ZZ(nn.numerator()).nbits() + ZZ(nn.denominator()).nbits()) > 200:
            return a0, s0          # факторизация нормы неподъёмна — оставляем как есть
        I = k.ideal(a0)
        try:
            ps = set([ZZ(2)])
            for P, ee in I.factor():
                ps.add(ZZ(P.smallest_integer()))
            key = tuple(sorted(ps))
            if key not in self._srcache:
                S = sorted(sum([k.primes_above(p) for p in sorted(ps)], []),
                           key=lambda P: (P.smallest_integer(), str(P)))
                self._srcache[key] = k.selmer_space(S, 2)
            V, gens, fromV, toV = self._srcache[key]
            rep = k(fromV(toV(a0)))
            q = a0 / rep
            assert q.is_square(), "small_rep: класс не совпал"
            return rep, s0 * q.sqrt()
        except Exception:
            return a0, s0

    def solve_conic(self, M, tag=None, tries=40, top=6):
        """ ненулевой v с v^t M v = 0 (M симметрична 3x3, невырождена, коника разрешима).
            Стратегия: (1) перебор малых точек; (2) МНОГО дешёвых случайных замен координат
            U в GL_3(O_k) + LLL-приведение формы + диагонализация + приведение коэффициентов по
            модулю квадратов; для каждого представления и каждой из трёх расстановок переменных
            считается tp (класс -B/A по модулю квадратов) — именно |N(tp)| определяет размер поля
            k(sqrt tp) и всю стоимость bnfisnorm; (3) is_norm запускается только для нескольких
            кандидатов с наименьшей |N(tp)|, с прерываемым (cysignals) таймаутом. """
        k = self.k
        key = (tag, str(M))
        if key in self.conic_cache:
            return self.conic_cache[key]
        for bound in (1, 2, 3):
            rng = [k(i) for i in range(-bound, bound + 1)] if _is_QQ(k) else \
                  [k_rand(k, bound) for _ in range(40)] + [k(0), k(1), k(-1)]
            for _ in range(400):
                v = vector(k, [random.choice(rng) for _ in range(3)])
                if v == 0:
                    continue
                if (v * M * v) == 0:
                    self.conic_cache[key] = v
                    return v
        cands = []
        seen = set()
        tlim = time.time() + float(getattr(self, 'conic_search_time', 60))
        for attempt in range(int(tries)):
            if attempt > 0 and time.time() > tlim:
                break
            try:
                for c in self._conic_candidates(M, attempt):
                    kk = str(c[1])
                    if kk in seen:
                        continue
                    seen.add(kk)
                    cands.append(c)
            except Exception:
                continue
        cands.sort(key=lambda c: c[0])
        if getattr(self, 'conic_verbose', False):
            print(f"    кандидатов коники: {len(cands)}; лучшие |N(tp)|: {[c[0] for c in cands[:5]]}")
        from cysignals.signals import AlarmInterrupt
        errs = []
        for c in cands[:int(top)]:
            cost, tp, w, m, i, j, l, abc, T = c
            try:
                if tp.is_square() if not _is_QQ(k) else QQ(tp).is_square():
                    st = tp.sqrt() if not _is_QQ(k) else QQ(tp).sqrt()
                    pt = [None] * 3; pt[i] = st * w; pt[j] = k(1); pt[l] = k(0)
                else:
                    Xn = PolynomialRing(k, 'Xn').gen()
                    Lrel = k.extension(Xn ^ 2 - tp, 'yy')
                    alarm(int(getattr(self, 'conic_timeout', 180)))
                    try:
                        ok, el = m.is_norm(Lrel, element=True, proof=False)
                    finally:
                        cancel_alarm()
                    assert ok, "не норма"
                    cs_ = el.list()
                    pt = [None] * 3
                    pt[i] = k(cs_[0]); pt[j] = k(cs_[1]) / w; pt[l] = k(1)
                a, b, cc = abc
                assert a * pt[0] ^ 2 + b * pt[1] ^ 2 + cc * pt[2] ^ 2 == 0, "точка не на диаг. конике"
                v = T * vector(k, pt)
                assert v != 0 and (v * M * v) == 0, "точка коники неверна"
                v = self._clear(v)
                self.conic_cache[key] = v
                return v
            except (Exception, AlarmInterrupt) as ex:
                errs.append(f"{type(ex).__name__}: {str(ex)[:60]}")
                if getattr(self, 'conic_verbose', False):
                    print("      ->", errs[-1])
        raise RuntimeError("коника не решена: " + "; ".join(errs[:6]))

    def _clear(self, v):
        k = self.k
        dd = lcm([self._den(t) for t in v])
        v = vector(k, [t * dd for t in v])
        try:
            gg = gcd([ZZ(s) for t in v for s in ([QQ(t)] if _is_QQ(k) else list(k(t))) if s != 0])
            if gg > 1:
                v = vector(k, [t / gg for t in v])
        except Exception:
            pass
        return v

    def _conic_candidates(self, M0, attempt):
        """ представления коники: (стоимость, tp, w, m, i, j, l, (a,b,c), T) """
        k = self.k
        U0 = identity_matrix(k, 3)
        M = M0
        if attempt > 0:
            for _ in range(50):
                U0 = matrix(k, 3, 3, lambda i, j: k_rand(k, 2) if i != j else k(1) + k_rand(k, 1))
                if U0.is_invertible():
                    break
            M = U0.transpose() * M0 * U0
        U = self.reduce_form(M)
        U0 = U0 * U
        Mr = U.transpose() * M * U
        dn = lcm([self._den(Mr[i][j]) for i in range(3) for j in range(3)])
        Mr = Mr * dn
        D, T0 = self.diagonalize(Mr)
        T = U0 * T0
        assert all(d != 0 for d in D)
        reps = [self.small_rep(d) for d in D]
        a, b, c = [r for r, _ in reps]
        Ts = T * diagonal_matrix(k, [1 / reps[t][1] for t in range(3)])
        out = []
        for (i, j, l) in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
            A, B, Cc = [a, b, c][i], [a, b, c][j], [a, b, c][l]
            t = -B / A; m = -Cc / A
            tp, w = self._best_tp(t, [a, b, c])
            nm = QQ(tp) if _is_QQ(k) else QQ(tp.norm())
            cost = abs(nm.numerator()) * abs(nm.denominator())
            out.append((cost, tp, w, m, i, j, l, (a, b, c), Ts))
        return out

    # --------- gamma1 ---------
    def gamma1(self, g1, g2, g3):
        """ gamma1(x,z) = коэффициент при phi^2 в (z(g2)z(g3)/m) H1;  возвращает [c0,c1,c2] (x^2, xz, z^2) """
        z1, z2, z3 = self.z_inv(g1), self.z_inv(g2), self.z_inv(g3)
        pr = z1 * z2 * z3
        m = self.E.sqrt(pr)
        assert m is not None, "z(g1)z(g2)z(g3) не квадрат в L — классы не в сумме 0"
        fac = m * self.E.inv(z1)            # = z(g2)z(g3)/m
        assert fac * z1 == m
        H = self.H_form(g1)
        gam = [self.E.coeffs(fac * Hj)[2] for Hj in H]
        assert any(c != 0 for c in gam), "gamma1 = 0 (не должно быть)"
        return gam, m

    # --------- места ---------
    def _ideal_primes(self, elts):
        k = self.k
        out = []
        if _is_QQ(k):
            s = set()
            for a in elts:
                a = QQ(a)
                if a == 0:
                    continue
                s |= set(ZZ(a.numerator()).prime_factors()) | set(ZZ(a.denominator()).prime_factors())
            return sorted(s)
        I = None
        for a in elts:
            a = k(a)
            if a == 0:
                continue
            J = k.ideal(a)
            I = J if I is None else I + J        # НОД идеалов (для content)
        return I

    def places_for(self, g1, gam, a2, extra_norm=16):
        """ Точное множество мест по Замечанию 3.3 Фишера: вклад места v тривиален, если
              (i) N(v) >= 11;  (ii) g1 и gamma1 v-целые и v не делит Delta(g1)*content(gamma1);
              (iii) a2 = g2(1,0) — v-единица и v не делит 2.
            Берём дополнение: архимедовы, v | 2, v | Delta(g1), v | знаменатели g1 и gamma1,
            v | content(gamma1), v | a2, и все v с N(v) <= extra_norm (>= 10). """
        k = self.k
        assert extra_norm >= 10
        if _is_QQ(k):
            ps = set([ZZ(2)]) | set(ZZ(p) for p in primes(extra_norm + 1))
            for c in [QQ(self.disc), QQ(a2)]:
                ps |= set(ZZ(c.numerator()).prime_factors()) | set(ZZ(c.denominator()).prime_factors())
            for c in list(g1) + list(gam):
                if c != 0:
                    ps |= set(ZZ(QQ(c).denominator()).prime_factors())
            cg = gcd([ZZ(QQ(c).numerator()) for c in gam if c != 0])
            ps |= set(ZZ(cg).prime_factors())
            return [LocSq(k, p) for p in sorted(ps)] + [RealPlace(k)]
        Ps = []
        def add_ideal(I):
            if I == 0 or I == k.ideal(1):
                return
            for P, _ in I.factor():
                Ps.append(P)
        add_ideal(k.ideal(self.disc))
        add_ideal(k.ideal(a2))
        den = k.ideal(1)
        for c in list(g1) + list(gam):
            if c != 0:
                den = den * k.ideal(c).denominator()
        add_ideal(den)
        cont = None
        for c in gam:
            if c != 0:
                J = k.ideal(c)
                cont = J if cont is None else cont + J
        if cont is not None:
            add_ideal(cont.numerator())
        rats = set([ZZ(2)]) | set(ZZ(p) for p in primes(extra_norm + 1))
        for p in sorted(rats):
            Ps += list(k.primes_above(p))
        uniq = []
        for P in Ps:
            if not any(P == Q for Q in uniq):
                uniq.append(P)
        uniq.sort(key=lambda P: (P.norm(), str(P)))
        return [LocSq(k, P) for P in uniq] + [RealPlace(k, emb) for emb in k.embeddings(AA)]

    # --------- локальные точки ---------
    def local_point(self, g1, gam, pl, tries=60000):
        """ (x, z) из k^2 с g1(x,z) — квадрат в k_v (ненулевой) и gam(x,z) != 0 """
        k = self.k
        gamv = lambda x, z: gam[0] * x ^ 2 + gam[1] * x * z + gam[2] * z ^ 2
        cand = [(k(1), k(0)), (k(0), k(1)), (k(1), k(1)), (k(1), k(-1)), (k(2), k(1)), (k(1), k(2))]
        for (x, z) in cand:
            v = quartic_eval(g1, x, z)
            if v != 0 and gamv(x, z) != 0 and pl.is_sq(v):
                return x, z
        if isinstance(pl, RealPlace):
            # интервалы положительности sigma(g1)(x,1): по вещественным корням
            emb = (lambda c: RR(QQ(c))) if pl.emb is None else (lambda c: RR(pl.emb(k(c))))
            Rr = PolynomialRing(RR, 'x'); xr = Rr.gen()
            gr = sum(emb(g1[j]) * xr ^ (4 - j) for j in range(5))
            rts = sorted([r for r in gr.roots(RR, multiplicities=False)])
            pts = []
            if rts:
                pts += [rts[0] - 1, rts[0] - 10, rts[-1] + 1, rts[-1] + 10]
                pts += [(rts[i] + rts[i + 1]) / 2 for i in range(len(rts) - 1)]
            pts += [RR(0), RR(1), RR(-1), RR(10), RR(-10), RR(100), RR(-100),
                    RR(1) / 10, RR(-1) / 10, RR(1) / 1000, RR(-1) / 1000]
            for x0 in pts:
                for jj in range(60):
                    xq = QQ(x0.nearby_rational(max_denominator=10 ^ 6)) + \
                         (QQ(random.randint(-1000, 1000)) / 10 ^ 6 if jj else 0)
                    x = k(xq); z = k(1)
                    v = quartic_eval(g1, x, z)
                    if v != 0 and gamv(x, z) != 0 and pl.is_sq(v):
                        return x, z
            raise RuntimeError(f"нет вещественной локальной точки ({pl.name}); g1 = {g1}")
        # кэш локальных точек для данной пары (g1, место)
        if not hasattr(self, '_lpcache'):
            self._lpcache = {}
        key = (str(g1), pl.name)
        pts = self._lpcache.get(key, [])
        for (x, z) in pts:
            if gamv(x, z) != 0:
                return x, z
        # (1) быстрый случайный поиск
        centers = self._approx_roots(g1, pl)
        for it in range(min(tries, 4000)):
            j = random.randint(-3, 8)
            cen = random.choice(centers + [k(0)])
            x = cen + (k_rand(k, pl.p ^ 4) * pl.pi ^ j if pl.kind == 'nf'
                       else k(random.randint(-pl.p ^ 4, pl.p ^ 4)) * QQ(pl.p) ^ j)
            z = k(1)
            v = quartic_eval(g1, x, z)
            if v != 0 and pl.is_sq(v):
                pts.append((x, z)); self._lpcache[key] = pts
                if gamv(x, z) != 0:
                    return x, z
                if len(pts) > 8:
                    break
        for (x, z) in pts:
            if gamv(x, z) != 0:
                return x, z
        # (2) систематический поиск по дереву классов вычетов.
        # Для целой квартики G и класса x in c + P^n O: G(x) = G(c) + O(P^n), поэтому при
        #     n - v(G(c)) >= 2 v(2) + 1
        # квадратичный класс G(x) постоянен на классе — класс либо принимается, либо мёртв.
        # Ветвление идёт НЕ по всем вычетам (поле вычетов бывает порядка 10^13 и больше),
        # а только по КОРНЯМ приведённого многочлена: для G_c(T) = G(c + pi^n T), m = min валюация
        # его коэффициентов, Ghat = G_c/pi^m, ghat = Ghat mod P — валюация v(G_c(T)) превышает m
        # лишь при ghat(T) = 0, значит живые подклассы отвечают корням ghat (их не больше 4).
        # «Общий» случай (валюация не растёт) ловится случайными пробами внутри класса.
        pi = pl.pi if pl.kind == 'nf' else QQ(pl.p)
        e2 = (pl.val(k(2)) if pl.kind == 'nf' else QQ(2).valuation(pl.p))
        prec = 2 * e2 + 1
        def integral(g):
            D = lcm([self._den(c) for c in g])
            return [c * D ^ 2 for c in g]
        RT = PolynomialRing(k, 'T'); TT = RT.gen()
        rf = pl.P.residue_field() if pl.kind == 'nf' else GF(pl.p)
        RU = PolynomialRing(rf, 'u'); uu = RU.gen()
        def redroots(poly):
            try:
                cs = [k(poly[t]) for t in range(poly.degree() + 1)]
                vs = [pl.val(c) for c in cs if c != 0]
                if not vs:
                    return None
                m = min(vs)
                hat = [c / pi ^ m for c in cs]
                ph = sum((rf(hat[t]) if pl.kind == 'nf' else rf(QQ(hat[t]))) * uu ^ t
                         for t in range(len(hat)))
                if ph == 0:
                    return None
                return [k(rf.lift(r)) if pl.kind == 'nf' else k(ZZ(r))
                        for r in ph.roots(multiplicities=False)]
            except Exception:
                return None
        charts = [(integral(g1), False), (integral(list(reversed(g1))), True)]
        for gg, rev in charts:
            gpoly = sum(gg[t] * TT ^ (4 - t) for t in range(5))
            live = [(k(0), 0)]
            for depth in range(0, 60):
                new = []
                for (c, n) in live:
                    for c2 in [c] + [c + pi ^ n * k_rand(k, 10 ^ 5) for _ in range(25)]:
                        val = quartic_eval(gg, c2, k(1))
                        x, z = (k(1), c2) if rev else (c2, k(1))
                        if val != 0 and pl.is_sq(val):
                            pts.append((x, z)); self._lpcache[key] = pts
                            if gamv(x, z) != 0:
                                return x, z
                    val = quartic_eval(gg, c, k(1))
                    if val != 0 and n - pl.val(val) >= prec:
                        continue                        # класс мёртв
                    rts = redroots(gpoly(c + pi ^ n * TT))
                    if not rts:
                        continue
                    for r in rts:
                        new.append((c + pi ^ n * r, n + 1))
                live = new
                if not live:
                    break
        for (x, z) in pts:
            if gamv(x, z) != 0:
                return x, z
        raise RuntimeError(f"нет локальной точки при {pl.name}; g1 = {g1}")

    RES_CAP = 4096

    def _residues(self, pl):
        """ полный набор представителей O/P; None, если поле вычетов слишком велико
            (тогда систематический перебор не имеет смысла — при N(P) большом случайная точка
            находится с вероятностью ~1/2 за попытку) """
        if not hasattr(self, '_rescache'):
            self._rescache = {}
        if pl.name in self._rescache:
            return self._rescache[pl.name]
        k = self.k
        if pl.kind == 'Q':
            out = [k(t) for t in range(pl.p)] if pl.p <= self.RES_CAP else None
        else:
            q = ZZ(pl.P.norm())
            if q > self.RES_CAP:
                out = None
            else:
                rf = pl.P.residue_field()
                out = [k(rf.lift(t)) for t in rf]
        self._rescache[pl.name] = out
        return out

    def _approx_roots(self, g, pl, m=30):
        """ приближённые (mod P^m) корни g(x,1) в k_v — центры для выборки """
        k = self.k
        key = ('rt', pl.name, str(g))
        if not hasattr(self, '_rtcache'):
            self._rtcache = {}
        if key in self._rtcache:
            return self._rtcache[key]
        out = []
        try:
            Rx = PolynomialRing(k, 'x'); x = Rx.gen()
            gp = sum(g[j] * x ^ (4 - j) for j in range(5))
            dgp = gp.derivative()
            if pl.kind == 'Q':
                p = pl.p
                for r in range(p):
                    if QQ(gp(r)).valuation(p) > 0 and QQ(dgp(r)).valuation(p) == 0:
                        y = QQ(r)
                        for _ in range(m):
                            y = y - gp(y) / dgp(y)
                            y = QQ(ZZ(y.numerator() * inverse_mod(y.denominator(), p ^ (m + 4))) % p ^ (m + 4))
                        out.append(k(y))
            else:
                P = pl.P; rf = P.residue_field()
                gr = sum(rf(g[j]) * rf['t'].gen() ^ (4 - j) for j in range(5)) if all(
                    k(c).valuation(P) >= 0 for c in g) else None
                if gr is not None and gr != 0:
                    for r, _ in gr.roots():
                        r0 = k(rf.lift(r))
                        if dgp(r0).valuation(P) != 0:
                            continue
                        y = r0
                        I = P ^ (m + 4)
                        for _ in range(3 * m):
                            y = y - gp(y) / dgp(y)
                            try:
                                y = I.reduce(y)
                            except Exception:
                                pass
                            if gp(y).valuation(P) >= m:
                                break
                        if gp(y).valuation(P) >= 8:
                            out.append(y)
        except Exception:
            out = []
        self._rtcache[key] = out
        return out

    # --------- спаривание ---------
    def pair(self, g1, g2, g3, reps=1, verbose=None, extra_norm=16):
        """ <[g1],[g2]>_CT в {0,1} (1 = -1) """
        verbose = self.verbose if verbose is None else verbose
        k = self.k
        self.check_quartic(g1); self.check_quartic(g2); self.check_quartic(g3)
        a2 = g2[0]
        assert a2 != 0, "g2(1,0) = 0 — класс тривиален, спаривание 0"
        gam, m = self.gamma1(g1, g2, g3)
        # нормировка gamma1 глобальной константой (свободна по формуле произведения)
        den = lcm([QQ(t).denominator() for c in gam for t in ([c] if _is_QQ(k) else list(k(c)))])
        gam = [c * den for c in gam]
        num = gcd([ZZ(t) for c in gam for t in ([QQ(c)] if _is_QQ(k) else list(k(c))) if t != 0])
        if num > 1:
            gam = [c / num for c in gam]
        pls = self.places_for(g1, gam, a2, extra_norm=extra_norm)
        results = []
        for rep in range(reps):
            if rep > 0:
                self._lpcache = {}          # другие локальные точки: результат обязан совпасть
            tot = 0; detail = []
            for pl in pls:
                x, z = self.local_point(g1, gam, pl)
                gv = gam[0] * x ^ 2 + gam[1] * x * z + gam[2] * z ^ 2
                s = hilb(k, a2, gv, pl)
                if s == -1:
                    tot += 1
                detail.append((pl.name, int(s)))
            results.append(tot % 2)
            if verbose:
                print(f"    места: {len(pls)}, нетривиальные: {[n for n, s in detail if s == -1]}")
        assert all(r == results[0] for r in results), f"CTP зависит от локальных точек: {results}"
        return results[0], len(pls)


# ------------------------------------------------------------------ от кривой к (I, J) и обратно

def IJ_of_curve(E):
    """ E ~ E_{I,J}: y^2 = x^3 - 27 I x - 27 J """
    I = E.c4(); J = 2 * E.c6()
    return I, J


def phi_of_root(E, e):
    """ корень e кубики y^2 = f(x) (модель a1 = a3 = 0) -> корень phi многочлена X^3 - 3IX + J.
        X = 36 x + 3 b2 — стандартный переход к y^2 = x^3 - 27 c4 x - 54 c6 = E_{I,J};
        x-координаты 2-кручения E_{I,J} равны -3 phi, значит phi = -(12 e + b2). """
    a1, a2, a3, a4, a6 = E.ainvs()
    assert a1 == 0 and a3 == 0
    b2 = 4 * a2
    return -(12 * e + b2)


# ------------------------------------------------------------------ драйверы: Sel^2 -> матрица CTP

def _comp_order(F, phis):
    """ для каждой компоненты L — индекс соответствующего корня phi (компоненты степени 1) """
    out = []
    for c in F.E.comps:
        assert c['deg'] == 1
        out.append([F.k(p) for p in phis].index(F.k(c['r'])))
    return out


def deltas_full2(F, E, roots, Cv, Sel):
    """ полное 2-кручение: Sel^2 из ek_descent.Curve3 -> список delta в L """
    k = F.k
    phis = [phi_of_root(E, k(e)) for e in roots]
    order = _comp_order(F, phis)
    out = []
    for w in Sel.basis():
        a = prod(Cv.gens[i] ^ int(w[i]) for i in range(Cv.n))
        b = prod(Cv.gens[i] ^ int(w[Cv.n + i]) for i in range(Cv.n))
        trip = [k(a), k(b), k(a * b)]
        out.append(F.E.from_comps([trip[i] for i in order]))
    return out


def deltas_partial(F, E, PD, Sel):
    """ одна точка порядка 2: Sel^2 из ek_descent_partial.PartialDescent -> список delta в L.
        E: y^2 = x(x^2 + a x + b), компоненты L: корень phi0 = -(b2) (для e = 0) и квадратичный множитель. """
    k = F.k; a = PD.a; b = PD.b
    b2 = 4 * a
    phi0 = -(12 * k(0) + b2)
    i1 = [i for i, c in enumerate(F.E.comps) if c['deg'] == 1 and F.k(c['r']) == phi0]
    assert len(i1) == 1, f"не найдена рациональная компонента L для phi0 = {phi0}"
    i1 = i1[0]
    i2 = [i for i in range(len(F.E.comps)) if i != i1]
    assert len(i2) == 1 and F.E.comps[i2[0]]['deg'] == 2, "L не k x K"
    i2 = i2[0]
    F2 = F.E.comps[i2]['F']; r2 = F.E.comps[i2]['r']
    e2 = -(r2 + b2) / 12                      # корень x^2 + a x + b в F2
    assert e2 ^ 2 + a * e2 + b == 0, "соответствие корней нарушено"
    # K = k(w), w^2 = d;  theta = c0 + c1 w  =>  w = (theta - c0)/c1;  образ w в F2: (e2 - c0)/c1
    cth = list(PD.theta)
    c0, c1 = k(cth[0]), k(cth[1])
    assert c1 != 0
    wF = (e2 - c0) / c1
    assert wF ^ 2 == F2(PD.d), "образ sqrt(d) в F2 неверен"
    out = []
    for v in Sel.basis():
        al = prod([g for g, c in zip(PD.gk, list(v)[:PD.nk]) if c == 1], k(1))
        be = prod([g for g, c in zip(PD.gK, list(v)[PD.nk:]) if c == 1], PD.Kabs(1))
        brel = PD.fromKabs(be)
        bc = list(brel)
        bF = F2(bc[0]) + F2(bc[1]) * wF
        vals = [None, None]
        vals[i1] = k(al); vals[i2] = bF
        out.append(F.E.from_comps(vals))
    return out


def ctp_matrix(F, deltas, verbose=True, with_diag=False, reps=1, sym_checks=3, cache=None):
    """ матрица спаривания Касселса–Тейта на базисе deltas (список элементов L).
        sym_checks: сколько элементов дополнительно пересчитать в ОБРАТНОМ порядке аргументов,
        <g_j, g_i> вместо <g_i, g_j>. Симметрия матрицы сама по себе тестом НЕ является
        (верхний треугольник просто копируется), а такой пересчёт — является: формула Фишера
        несимметрична по g1 и g2 (H строится по g1, а старший коэффициент берётся у g2). """
    n = len(deltas)
    t0 = time.time()
    quart = {}
    import os as _os
    if cache and _os.path.exists(cache):
        try:
            quart = load(cache)
            if verbose:
                print(f"    квартики загружены из кэша {cache}: {len(quart)}")
        except Exception:
            quart = {}
    def _need(key, dl, tag):
        if key in quart:
            return
        quart[key] = F.quartic_from_delta(dl, tag=tag)
        if cache:
            try:
                save(quart, cache)
            except Exception:
                pass
    for i in range(n):
        _need((i,), deltas[i], f"b{i}")
        if verbose:
            print(f"    квартика {i}: {quart[(i,)]}  ({time.time()-t0:.0f}s)")
    for i in range(n):
        for j in range(i + 1, n):
            _need((i, j), deltas[i] * deltas[j], f"s{i}_{j}")
    if with_diag:
        _need((), F.L(1), "triv")
    if verbose:
        print(f"    квартик построено: {len(quart)}  ({time.time()-t0:.0f}s)")
    M = matrix(GF(2), n, n)
    for i in range(n):
        for j in range(i + 1, n):
            v, npl = F.pair(quart[(i,)], quart[(j,)], quart[(i, j)], reps=reps)
            M[i, j] = v; M[j, i] = v
            if verbose:
                print(f"    <{i},{j}> = {v}  (мест {npl}, {time.time()-t0:.0f}s)")
    if with_diag:
        for i in range(n):
            v, _ = F.pair(quart[(i,)], quart[(i,)], quart[()], reps=reps)
            M[i, i] = v
            if verbose:
                print(f"    <{i},{i}> = {v}")
    # независимая проверка симметрии: пересчёт в обратном порядке аргументов
    sym_ok = []
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    random.shuffle(pairs)
    pairs.sort(key=lambda t: 0 if M[t[0], t[1]] != 0 else 1)   # сперва НЕнулевые элементы
    for (i, j) in pairs[:int(sym_checks)]:
        vr, _ = F.pair(quart[(j,)], quart[(i,)], quart[(i, j)], reps=reps)
        ok = (GF(2)(vr) == M[i, j])
        sym_ok.append(((i, j), int(M[i, j]), int(vr), bool(ok)))
        if verbose:
            print(f"    симметрия <{j},{i}> = {vr} vs <{i},{j}> = {M[i,j]}: {'OK' if ok else 'РАСХОЖДЕНИЕ!!!'}")
        assert ok, f"спаривание несимметрично на ({i},{j})"
    return M, quart, sym_ok
