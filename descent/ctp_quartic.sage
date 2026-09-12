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
        Rt = PolynomialRing(k, 't0,t1,t2'); t0, t1, t2 = Rt.gens()
        RtX = PolynomialRing(Rt, 'X')
        cub = RtX([Rt(c) for c in self.E.cub.list()])
        Lt = RtX.quotient(cub, 'th'); th = Lt.gen()
        thx = -3 * th
        dl = self.E.coeffs(delta)
        dlt = Lt(RtX([Rt(c) for c in dl]))
        tt = t0 + t1 * thx + t2 * thx ^ 2
        pr = (dlt * tt * tt).lift()
        d = [Rt(pr[j]) for j in range(3)]
        return [d[0], -d[1] / 3, d[2] / 9]      # коэффициенты в базисе (1, Theta, Theta^2)

    def quartic_from_delta(self, delta, tag=None):
        """ бинарная квартика с инвариантами (I,J), представляющая класс delta в L^*/L^*2 """
        k = self.k
        c0, c1, c2 = self.covering_forms(delta)
        Rt = c2.parent(); tv = Rt.gens()
        M = matrix(k, 3, 3, lambda i, j: (c2.coefficient({tv[i]: 2}) if i == j
                                          else c2.coefficient({tv[i]: 1, tv[j]: 1}) / 2))
        assert M.determinant() != 0, "вырожденная коника"
        P = self.solve_conic(M, tag=tag)
        # параметризация коники: t(x,z) = Q(w) P - 2 B(P,w) w,  w = x v1 + z v2
        bas = [vector(k, [1, 0, 0]), vector(k, [0, 1, 0]), vector(k, [0, 0, 1])]
        comp = [b for b in bas if matrix(k, [P, b]).rank() == 2]
        v1, v2 = None, None
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
        tpar = [Qw * P[i] - 2 * BPw * w[i] for i in range(3)]
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

    def _from_sigmas(self, T, prec):
        """ элемент k, приближающий заданные значения T[j] при вложениях """
        k = self.k
        if _is_QQ(k):
            return QQ(RealField(prec)(T[0]).nearby_rational(max_error=abs(RealField(prec)(T[0])) / 2 ^ 40 + 2 ^ (-40)))
        R = RealField(prec)
        if self._emb_list()[0][0] == 'R':
            t1, t2 = R(T[0]), R(T[1])
            c0 = (t1 + t2) / 2
            c1 = (t1 - t2) / (2 * R(self._D).sqrt())
        else:
            z = ComplexField(prec)(T[0])
            c0 = z.real()
            c1 = z.imag() / R(-self._D).sqrt()
        def rat(u):
            u = R(u)
            if u == 0:
                return QQ(0)
            return QQ(u.nearby_rational(max_error=abs(u) / 2 ^ 60))
        return k(rat(c0)) + k(rat(c1)) * k.gen()

    def _prec_for(self, g):
        return max(64, 2 * self._msr(g) // 5 + 200)

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
                g = self.gl2(g, [[1, -g[1] / (4 * g[0])], [0, 1]])
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
            prec = min(prec_cap, max(200, m0 // 2 + 200))
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
            try:
                c = self._from_sigmas(X0, prec)
                s = self._from_sigmas(Y0, prec)
            except Exception:
                break
            if s == 0:
                break
            try:
                gn = self.gl2(g, [[s, c], [0, 1]])
            except Exception:
                break
            if self._msr(gn) >= m0:
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

    def solve_conic(self, M, tag=None):
        """ ненулевой v с v^t M v = 0 (M симметрична 3x3, невырождена, коника разрешима) """
        k = self.k
        key = (tag, str(M))
        if key in self.conic_cache:
            return self.conic_cache[key]
        # маленькие точки перебором
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
        D, T = self.diagonalize(M)
        assert all(d != 0 for d in D)
        reps = [self.small_rep(d) for d in D]
        a, b, c = [r for r, _ in reps]
        pt = None
        order = []
        for (i, j, l) in [(0, 1, 2), (0, 2, 1), (1, 2, 0)]:
            A, B = [a, b, c][i], [a, b, c][j]
            t = -B / A
            if (QQ(t).is_square() if _is_QQ(k) else t.is_square()):
                st = QQ(t).sqrt() if _is_QQ(k) else t.sqrt()
                pt = [None] * 3; pt[i] = st; pt[j] = k(1); pt[l] = k(0)
                break
            nm = QQ(t) if _is_QQ(k) else t.norm()
            order.append((abs(QQ(nm).numerator()) * abs(QQ(nm).denominator()), i, j, l))
        if pt is None:
            errs = []
            for cost, i, j, l in sorted(order):
                A, B, Cc = [a, b, c][i], [a, b, c][j], [a, b, c][l]
                try:
                    t = -B / A; m = -Cc / A
                    if _is_QQ(k):
                        con = Conic(QQ, [QQ(a), QQ(b), QQ(c)])
                        ok, p0 = con.has_rational_point(point=True)
                        assert ok
                        pt = [QQ(p0[0]), QQ(p0[1]), QQ(p0[2])]
                        break
                    dd = lcm([QQ(cc).denominator() for cc in list(t)]); tp = t * dd ^ 2
                    Xn = PolynomialRing(k, 'Xn').gen()
                    Lrel = k.extension(Xn ^ 2 - tp, 'yy')
                    ok, el = m.is_norm(Lrel, element=True, proof=False)
                    assert ok, "не норма (коника без точки?)"
                    cs_ = el.list()
                    pt = [None] * 3
                    pt[i] = k(cs_[0]); pt[j] = k(cs_[1]) * dd; pt[l] = k(1)
                    break
                except Exception as ex:
                    errs.append(f"{(i,j,l)}: {type(ex).__name__}: {str(ex)[:80]}")
            if pt is None:
                raise RuntimeError("коника не решена: " + "; ".join(errs))
        assert a * pt[0] ^ 2 + b * pt[1] ^ 2 + c * pt[2] ^ 2 == 0
        y = vector(k, [pt[t] / reps[t][1] for t in range(3)])
        v = T * y
        assert v != 0 and (v * M * v) == 0, "точка коники неверна"
        self.conic_cache[key] = v
        return v

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
        # (2) систематический поиск по дереву классов вычетов с ОТСЕЧЕНИЕМ.
        # Для целой квартики G и класса x in c + P^n O: G(x) = G(c) + O(P^n), поэтому если
        #     n - v(G(c)) >= 2 v(2) + 1,
        # то квадратичный класс G(x) постоянен на классе: либо принимаем, либо класс мёртв.
        # Мёртвые классы отбрасываются; живые ветвятся. g1 всюду локально разрешима => решение найдётся.
        res = self._residues(pl)
        pi = pl.pi if pl.kind == 'nf' else QQ(pl.p)
        e2 = (pl.val(k(2)) if pl.kind == 'nf' else QQ(2).valuation(pl.p))
        prec = 2 * e2 + 1
        # целые модели (умножение на квадрат не меняет квадратичный класс)
        def integral(g):
            D = lcm([self._den(c) for c in g])
            return [c * D ^ 2 for c in g]
        charts = [(integral(g1), False), (integral(list(reversed(g1))), True)]
        for gg, rev in charts:
            live = [(k(0), 0)]
            for depth in range(1, 60):
                new = []
                for (c, n) in live:
                    for r in res:
                        c2 = c + r * pi ^ n
                        if rev and c2 == 0 and n == 0:
                            continue
                        val = quartic_eval(gg, c2, k(1))
                        x, z = (k(1), c2) if rev else (c2, k(1))
                        if val != 0 and pl.is_sq(val):
                            pts.append((x, z)); self._lpcache[key] = pts
                            if gamv(x, z) != 0:
                                return x, z
                            new.append((c2, n + 1)); continue
                        if val == 0:
                            new.append((c2, n + 1)); continue
                        if (n + 1) - pl.val(val) >= prec:
                            continue                       # класс мёртв
                        new.append((c2, n + 1))
                live = new
                if not live:
                    break
        for (x, z) in pts:
            if gamv(x, z) != 0:
                return x, z
        raise RuntimeError(f"нет локальной точки при {pl.name}; g1 = {g1}")

    def _residues(self, pl):
        if not hasattr(self, '_rescache'):
            self._rescache = {}
        if pl.name in self._rescache:
            return self._rescache[pl.name]
        k = self.k
        if pl.kind == 'Q':
            out = [k(t) for t in range(pl.p)]
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


def ctp_matrix(F, deltas, verbose=True, with_diag=False, reps=1):
    """ матрица спаривания Касселса–Тейта на базисе deltas (список элементов L) """
    n = len(deltas)
    t0 = time.time()
    quart = {}
    for i in range(n):
        quart[(i,)] = F.quartic_from_delta(deltas[i], tag=f"b{i}")
        if verbose:
            print(f"    квартика {i}: {quart[(i,)]}  ({time.time()-t0:.0f}s)")
    for i in range(n):
        for j in range(i + 1, n):
            quart[(i, j)] = F.quartic_from_delta(deltas[i] * deltas[j], tag=f"s{i}_{j}")
    if with_diag:
        quart[()] = F.quartic_from_delta(F.L(1), tag="triv")
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
    return M, quart
