# Поиск конфигураций «крест + k квадратных углов» (7 = 5+2, 8 = 5+3, 9 = 5+4 квадратов) через эллиптические кривые углов.
# Сечение (b,h,n): крест B(s)=b/n, H(s)=h/n; второй крест B(t), H(t). Угол i — квадрат  <=>  Y_i^2 = Q_i(t),
#   Q_{1,2} = A(t^2+1)^2 ± beta t(t^2-1),  Q_{3,4} = C(t^2+1)^2 ± beta t(t^2-1),  A = (h^2+n^2)/2, C = (b^2+n^2)/2, beta = 2n^2.
# Кривая Y^2 = Q_1(t) рода 1 имеет точку t0 (параметр сечения, вырожденная точка). Все t с квадратным углом 1 —
# образы E_1(Q); перебираем P = sum n_j G_j + T и проверяем остальные углы.
import sys, itertools, functools
print = functools.partial(print, flush=True)
R = PolynomialRing(QQ, 'T'); T = R.gen()

def quartics(b, h, n):
    A, C, B = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2), 2*n^2
    return [A*(T^2+1)^2 + B*T*(T^2-1), A*(T^2+1)^2 - B*T*(T^2-1), C*(T^2+1)^2 + B*T*(T^2-1), C*(T^2+1)^2 - B*T*(T^2-1)]

def section_params(b, h, n):
    return [QQ(n + h)/(n + b), QQ(n - h)/(n + b)]

class QuarticCurve:
    """ Y^2 = q(t), рациональная точка t0 (q(t0) = m^2 != 0). X = t - t0: Y^2 = a X^4 + bb X^3 + c X^2 + d X + m^2.
        Формулы Коннелла (Handbook, Prop. 1.2.1) — проверяются численно при построении. """
    def __init__(self, q, t0):
        self.q, self.t0 = q, t0
        g = R(q(T + t0))
        self.m = g[0].sqrt(); assert self.m^2 == g[0] and self.m != 0
        self.a, self.bb, self.c, self.d = g[4], g[3], g[2], g[1]
        qq, a, bb, c, d = self.m, self.a, self.bb, self.c, self.d
        a1 = d/qq; a2 = c - d^2/(4*qq^2); a3 = 2*qq*bb; a4 = -4*qq^2*a; a6 = a2*a4
        self.E = EllipticCurve([a1, a2, a3, a4, a6])
    def to_E(self, X, Y):
        qq, c, d = self.m, self.c, self.d
        x = (2*qq*(Y + qq) + d*X)/X^2
        y = (4*qq^2*(Y + qq) + 2*qq*(d*X + c*X^2) - d^2*X^2/(2*qq))/X^3
        return x, y
    def from_E(self, P):
        if P.is_zero(): return None
        x, y = P.xy()
        if y == 0: return None
        qq, c, d = self.m, self.c, self.d
        a1, a3 = self.E.a1(), self.E.a3()
        # общий вид для модели с a1, a3: используем y' = y + (a1 x + a3)/2-сдвиг не нужен — берём формулу Коннелла
        X = (2*qq*(x + c) - d^2/(2*qq))/y
        if X == 0: return None                    # базовая точка t0 (второй знак Y) — вырождена
        Y = -qq + X*(X*x - d)/(2*qq)
        return X, Y

def self_test(b, h, n):
    Q = quartics(b, h, n); t0 = [s_ for s_ in section_params(b, h, n) if all(q(s_).is_square() and q(s_) != 0 for q in Q)][0]
    qc = QuarticCurve(Q[0], t0)
    import random
    ok = True
    # точки квартики: из точек E (прямое/обратное)
    G = qc.E.gens(proof=False)
    for P in [k*G[0] for k in range(1, 5)] if G else []:
        XY = qc.from_E(P)
        if XY is None: continue
        X, Y = XY
        if Y^2 != qc.q(X + t0): ok = False
        x2, y2 = qc.to_E(X, Y)
        if (x2, y2) != P.xy(): ok = False
    return ok, qc

def search_section(b, h, n, N=12, max_rank=3):
    Q = quartics(b, h, n)
    t0s = [s_ for s_ in section_params(b, h, n) if all(q(s_).is_square() and q(s_) != 0 for q in Q)]
    if not t0s: return None
    t0 = t0s[0]
    orbit = set()
    for s_ in section_params(b, h, n):
        for u in (s_, -s_, 1/s_, -1/s_): orbit.add(u)
    orbit |= {QQ(0), QQ(1), QQ(-1)}              # второй крест вырожден (клетки = центр)
    results = []
    for i in range(4):
        qc = QuarticCurve(Q[i], t0)
        E = qc.E
        try:
            r = E.rank(only_use_mwrank=False)
            G = E.gens(proof=False)
        except Exception as ex:
            results.append((i, 'rank?', None)); continue
        if len(G) > max_rank: results.append((i, len(G), 'skip')); continue
        tors = E.torsion_points()
        found = {}
        for ns in itertools.product(range(-N, N+1), repeat=len(G)):
            P0 = sum((k*g for k, g in zip(ns, G)), E(0))
            for Tt in tors:
                XY = qc.from_E(P0 + Tt)
                if XY is None: continue
                t1 = XY[0] + t0
                if t1 in orbit or t1 in found: continue
                k = sum(1 for q in Q if q(t1).is_square())
                if k >= 2: found[t1] = k
        results.append((i, len(G), found))
    return t0, results

if len(sys.argv) > 1 and sys.argv[1] == 'test':
    for sec in [(17,7,13), (7,1,5), (23,7,17)]:
        ok, qc = self_test(*sec)
        print(f"self-test {sec}: maps consistent: {ok};  E1 = {qc.E.minimal_model().ainvs()}, rank {qc.E.rank()}")
