# Независимая реализация эллиптического Шаботи (Bruin) для сечений семейства.
#
# Задача: найти все P in E(k) с u(P) = x(P)/scale in P^1(Q), k = Q(sqrt(dk)) вещественное,
# rank E(k) = 1, E(k) = <G> + E(k)_tors.
#
# Метод: p расщепляется в k, хорошая редукция; два вложения sigma_1, sigma_2 : k -> Q_p.
# N = lcm порядков редукций G в E(F_p) по обоим вложениям; R = N*G лежит в ядре редукции.
# Для каждого класса (T, m0), T in tors, 0 <= m0 < N:  Q = m0*G + T, и точки класса это
#      P(s) = Q + [z(s)],  z(s) = exp_F(s * log_F(t_R)),  s in Z  (аналитически s in Z_p).
# Условие x(P) in P^1(Q) влечёт sigma_1 x(P) = sigma_2 x(P), т.е. theta(s) = 0, где
#      theta(s) = phi(sigma_1 x(P(s))) - phi(sigma_2 x(P(s))),  phi(X) = 1/(X-gamma),
# gamma in Z выбрано так, что x - gamma -- единица в обоих дисках (тогда обе ветви в Z_p[[s]]).
# Число нулей в Z_p оценивается теоремой Штрассмана; хвост: v_p(c_n) >= n*v - v_p(n!),
# v = min_i v_p(t_{R,i}) >= 1  (d_m in Z_p, т.к. phi(x(Q+[t])) регулярна и ограничена на диске;
#  exp_F(W) = sum b_j W^j/j!, b_j in Z_p (Silverman IV.6.4)).
import functools, sys, time
print = functools.partial(print, flush=True)
load('/home/kep/magicKube/family/independent_chabauty/setup.sage')


def vpfact(n, p):
    s, q = 0, p
    while q <= n:
        s += n//q; q *= p
    return s


def strassmann(coeffs, v, p, M, verbose=True):
    """coeffs = [c_0,...,c_{M-1}] в Q_p; хвост v_p(c_n) >= n*v - v_p(n!) для n >= M.
       Возвращает (число нулей в Z_p с кратностью, mu, диагностика)."""
    vals = []
    for c in coeffs:
        if c == 0:
            vals.append(c.precision_absolute())   # нижняя оценка
        else:
            vals.append(c.valuation())
    mu = min(vals)
    n0 = max(i for i, x in enumerate(vals) if x == mu)
    # проверка, что хвост заведомо больше mu
    ok = all(n*v - vpfact(n, p) > mu for n in range(M, M + 400))
    # и что оценка хвоста монотонно растёт
    return n0, mu, ok, vals


class Chab:
    def __init__(self, S, G, p, prec=80, M=30):
        self.S, self.G, self.p, self.prec, self.M = S, G, p, prec, M
        E, k = S.E, S.k
        assert kronecker(S.dk, p) == 1, "p не расщепляется"
        assert all(a.is_integral() for a in E.ainvs())
        assert E.discriminant().norm() % p != 0, "плохая редукция в p"
        self.Qp = Qp(p, prec)
        self.rt = self.Qp(S.dk).sqrt()
        Fp = GF(p)
        self.Fp = Fp
        self.rb = Fp(self.rt.residue())
        assert self.rb^2 == Fp(S.dk)
        self.tors = E.torsion_points()

    def emb(self, a, sgn):
        c = self.S.k(a).list()
        return self.Qp(c[0]) + sgn*self.Qp(c[1])*self.rt

    def red(self, a, sgn):
        z = self.emb(a, sgn)
        assert z.valuation() >= 0, "не p-целое"
        return self.Fp(z.residue())

    def Ebar(self, sgn):
        return EllipticCurve(self.Fp, [self.red(a, sgn) for a in self.S.E.ainvs()])

    def redpt(self, P, sgn):
        Eb = self.Ebar(sgn)
        if P.is_zero():
            return Eb(0)
        if self.emb(P[0], sgn).valuation() < 0:
            return Eb(0)
        return Eb([self.red(P[0], sgn), self.red(P[1], sgn)])

    def setup(self):
        self.orders = [self.redpt(self.G, s).order() for s in (1, -1)]
        self.N = lcm(self.orders)
        self.R = self.N * self.G
        E = self.S.E
        M = self.M
        fg = E.formal_group()
        self.xt = fg.x(M + 4)
        self.yt = fg.y(M + 4)
        self.lg = fg.log(M + 2)
        self.ex = self.lg.reverse()
        # формальные параметры t_R по вложениям
        self.tR = []
        for sgn in (1, -1):
            xR, yR = self.emb(self.R[0], sgn), self.emb(self.R[1], sgn)
            self.tR.append(-xR/yR)
        self.v = min(t.valuation() for t in self.tR)
        # log_F(t_R)
        self.l = []
        for i, sgn in enumerate((1, -1)):
            L = sum(self.emb(c, sgn)*self.tR[i]^j for j, c in enumerate(self.lg.list()))
            self.l.append(L)
        return self

    def xseries(self, Q):
        """точный ряд x(Q+[t]) над k (Лоран при Q = O)."""
        if Q.is_zero():
            return self.xt
        a2 = self.S.E.a2()
        xq, yq = Q[0], Q[1]
        lam = (self.yt - yq)/(self.xt - xq)
        return lam^2 - a2 - xq - self.xt

    def theta(self, Q):
        """ряд theta(s) для класса точки Q; возвращает (коэффициенты, gamma)."""
        p, M = self.p, self.M
        # выбор gamma
        bad = set()
        for sgn in (1, -1):
            P = self.redpt(Q, sgn)
            if not P.is_zero():
                bad.add(P[0])
        gamma = next(ZZ(c) for c in self.Fp if c not in bad)
        Xs = self.xseries(Q)
        Phi = 1/(Xs - gamma)            # степенной ряд над k
        Phi = Phi.power_series() if hasattr(Phi, 'power_series') else Phi
        co = [Phi[j] for j in range(M)]
        Rs = PowerSeriesRing(self.Qp, 's', default_prec=M)
        s = Rs.gen()
        out = []
        for i, sgn in enumerate((1, -1)):
            # z(s) = exp_F(l_i * s)
            zs = sum(self.Qp(self.emb(c, sgn))*(self.l[i]*s)^j for j, c in enumerate(self.ex.list()))
            zs = Rs(zs).add_bigoh(M)
            f = Rs(0)
            zp = Rs(1)
            for j in range(M):
                f += self.emb(co[j], sgn)*zp
                zp = (zp*zs).add_bigoh(M)
            out.append(f)
        th = (out[0] - out[1]).add_bigoh(M)
        return [th[j] for j in range(M)], gamma, [[o[j] for j in range(M)] for o in out], co
