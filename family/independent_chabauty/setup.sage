# Независимая реализация (сессия "проверяющий"): построение E1/k для сечения (b,h,n)
# и её точек с рациональной u-координатой. Вывод получен заново из постановки:
#
#   A = (h^2+n^2)/2, C = (b^2+n^2)/2,  beta' = (C-A)/(C+A),  e = 1+beta'
#   C_J2^{(D)}: D W^2 = t(t^2-1)(A t - C)(C t - A),  z = (t-1)/(t+1)
#       =>  D W^2 = cz * z(z^2-1)(z^2-beta'^2),  cz = -4(A+C)^2     [проверено символически]
#   u = z + beta'/z,  sb = sqrt(beta') in k = Q(sqrt(dk)), dk = sqfree((C-A)(C+A))
#       z(z^2-1)(z^2-beta'^2) = z^3 (u^2 - e^2),  (z+sb)^2 = z(u+2sb)  =>  z^3 = (z+sb)^6/(u+2sb)^3
#       =>  D V^2 = cz (u^2-e^2)(u+2sb),  V = W(u+2sb)^2/(z+sb)^3
#   E: Y^2 = (x-c0*e)(x+c0*e)(x+2*sb*c0),  x = c0*u, Y = c0*V, c0 = cz/D (умноженное на квадрат).
# Ищем точки P in E(k) с x(P) in P^1(Q)  <=>  u(P) in P^1(Q).
import functools
print = functools.partial(print, flush=True)

def sqfree(n):
    n = ZZ(n)
    return sign(n) * prod(p^(ex % 2) for p, ex in n.abs().factor())

class Section:
    def __init__(self, b, h, n, D):
        assert b^2 + h^2 == 2*n^2 and gcd([b, h, n]) == 1
        self.b, self.h, self.n = ZZ(b), ZZ(h), ZZ(n)
        self.A = A = ZZ((h^2 + n^2)/2)
        self.C = C = ZZ((b^2 + n^2)/2)
        assert gcd(A, C) == 1
        self.bet = bet = QQ(C - A)/QQ(C + A)
        self.e = 1 + bet
        self.dk = dk = sqfree((C - A)*(C + A))
        self.k = k = QuadraticField(dk, 'w')
        self.w = k.gen()
        q = ZZ(sqrt(QQ((C - A)*(C + A))/dk))
        self.sb = sb = self.w * q/(C + A)          # sqrt(beta') > 0
        assert sb^2 == bet
        self.D = D = ZZ(D)
        self.cz = cz = -4*(A + C)^2
        c0 = QQ(cz)/QQ(D)
        c0 = c0 * QQ(c0.denominator())^2           # x -> lam^2 x, модель целочисленная
        # сокращаем на квадраты, пока корни остаются целыми алгебраическими
        changed = True
        while changed:
            changed = False
            for pr in ZZ(c0).prime_divisors():
                cand = c0/pr^2
                if all((cand*r).is_integral() for r in (self.e, -self.e, -2*sb)):
                    c0 = cand; changed = True; break
        self.c0 = c0
        self.roots = [c0*self.e, -c0*self.e, -2*sb*c0]
        R = PolynomialRing(k, 'X'); X = R.gen()
        cub = R(prod(X - r for r in self.roots))
        self.E = EllipticCurve(k, [0, cub[2], 0, cub[1], cub[0]])
        self.scale = c0                            # x = scale * u

    def rhs(self, x):
        return prod(self.k(x) - r for r in self.roots)

    def u_of(self, P):
        return None if P.is_zero() else P[0]/self.scale

    def point_from_u(self, u):
        x = self.k(u)*self.scale
        v = self.rhs(x)
        if not v.is_square():
            return None
        return self.E([x, v.sqrt()])

    def t_of_u(self, u):
        """все рациональные t с z+beta'/z = u, z=(t-1)/(t+1)."""
        Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
        ts = []
        for z0 in (z^2 - QQ(u)*z + self.bet).roots(QQ, multiplicities=False):
            if z0 != 1:
                ts.append((1 + z0)/(1 - z0))
        return ts
