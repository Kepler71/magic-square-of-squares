# Поиск типа Бремнера: полная пара креста (b, h) (b^2 + h^2 = 2 C0), центр C0 — любой (квадрат или нет),
# клетка f^2 (f рационально) второй пары; её противоположная d^2 = 2C0 - f^2 — не обязательно квадрат.
# Квадрат:  TL  b^2  TR          TL = (h^2 + f^2)/2,   TR = (b^2 + 2h^2 - f^2)/2
#           d^2 C0   f^2         BL = (b^2 + f^2)/2,   BR = (2b^2 + h^2 - f^2)/2,   D = d^2 = b^2 + h^2 - f^2
#           BL  h^2  BR
# Число квадратов = 3 (b^2, h^2, f^2) + [C0 = □] + [D = □] + #углов.  Цель: >= 7 при различных клетках.
# Базовые кривые рода 1: пары условий с известной точкой (f = ±h для {TL, BR, D}; f = ±b для {BL, TR, D}),
# f пробегает решётку Морделла–Вейля (шаги по генераторам).
import itertools, functools, sys
print = functools.partial(print, flush=True)
Rm = PolynomialRing(QQ, 'm'); mvar = Rm.gen()

def conds(b, h):
    C0 = QQ(b^2 + h^2)/2
    # каждое условие: (alpha, beta) для alpha + beta f^2
    return C0, {'TL': (QQ(h^2)/2, QQ(1)/2), 'TR': (QQ(b^2 + 2*h^2)/2, -QQ(1)/2), 'BL': (QQ(b^2)/2, QQ(1)/2),
                'BR': (QQ(2*b^2 + h^2)/2, -QQ(1)/2), 'D': (QQ(b^2 + h^2), -QQ(1))}

class BaseCurve:
    """ Y1^2 = a1 + b1 f^2 (коника с точкой (f0, Y10)), Y2^2 = a2 + b2 f^2. Параметризация коники: f(m), затем
        Y2^2 = Q(m) (квартика, умноженная на квадрат) -> эллиптическая кривая через точку m0 (Коннелл). """
    def __init__(self, c1, c2, f0):
        (a1, b1), (a2, b2) = c1, c2
        Y10 = (a1 + b1*f0^2).sqrt(); assert Y10 in QQ and Y10 != 0
        self.a1, self.b1, self.a2, self.b2, self.f0, self.Y10 = a1, b1, a2, b2, f0, Y10
        den = mvar^2 - b1
        self.num = f0*den + 2*b1*f0 - 2*Y10*mvar          # f(m) = num/den
        self.den = den
        self.Q = Rm(a2*den^2 + b2*self.num^2)            # Y2'^2 = Q(m),  Y2' = Y2*den
        # рациональная точка на квартике: ищем малую
        pts = pari.hyperellratpoints(pari(self.Q), 200)
        pts = [(QQ(P[0]), QQ(P[1])) for P in pts if QQ(P[1]) != 0 and self.den(QQ(P[0])) != 0]
        if not pts: raise ValueError("no point on quartic")
        self.m0 = pts[0][0]
        g = Rm(self.Q(mvar + self.m0))
        qq = g[0].sqrt(); assert qq^2 == g[0]
        self.qq, self.aa, self.bb, self.cc, self.dd = qq, g[4], g[3], g[2], g[1]
        A1 = self.dd/qq; A2 = self.cc - self.dd^2/(4*qq^2); A3 = 2*qq*self.bb; A4 = -4*qq^2*self.aa; A6 = A2*A4
        self.E = EllipticCurve([A1, A2, A3, A4, A6])
    def f_of_point(self, P):
        if P.is_zero(): return None
        x, y = P.xy()
        if y == 0: return None
        X = (2*self.qq*(x + self.cc) - self.dd^2/(2*self.qq))/y
        m = X + self.m0
        if self.den(m) == 0: return None
        return self.num(m)/self.den(m)

def score(b, h, f, C0, cs):
    sq = {k: (a + be*f^2) for k, (a, be) in cs.items()}
    flags = {k: (v > 0 and v.is_square()) for k, v in sq.items()}
    cells = [sq['TL'], QQ(b^2), sq['TR'], sq['D'], C0, QQ(f^2), sq['BL'], QQ(h^2), sq['BR']]
    distinct = len(set(cells)) == 9 and all(c > 0 for c in cells)
    nsq = 3 + int(C0.is_square()) + int(flags['D']) + sum(flags[k] for k in ('TL', 'TR', 'BL', 'BR'))
    return nsq, distinct, flags, cells

def search_pair(b, h, N=6, max_rank=4, verbose=False):
    C0, cs = conds(b, h)
    bases = [('TL', 'BR', h), ('TL', 'D', h), ('BR', 'D', h), ('BL', 'TR', b), ('BL', 'D', b), ('TR', 'D', b)]
    best = []
    seen = set()
    stats = []      # по каждой базовой кривой: ('ok', rank, Nr) | ('error', сообщение) | ('skip', rank)
    for c1, c2, f0 in bases:
        try:
            bc = BaseCurve(cs[c1], cs[c2], QQ(f0))
            G = bc.E.gens(proof=False)
        except Exception as ex:
            stats.append(('error', f"{c1}/{c2}: {type(ex).__name__}: {str(ex)[:60]}"))
            if verbose: print("  base", c1, c2, "failed:", ex)
            continue
        if len(G) > max_rank:
            stats.append(('skip', f"{c1}/{c2}: rank {len(G)}")); continue   # ранг 0: проверяем только точки кручения
        tors = bc.E.torsion_points()
        Nr = {1: 12, 2: 8, 3: 5, 4: 4}.get(len(G), N) if N == 'auto' else N
        stats.append(('ok', f"{c1}/{c2}: rank {len(G)}, |n|<={Nr}"))
        for ns in itertools.product(range(-Nr, Nr+1), repeat=len(G)):
            P0 = sum((k*g for k, g in zip(ns, G)), bc.E(0))
            for Tt in tors:
                f = bc.f_of_point(P0 + Tt)
                if f is None or f in seen: continue
                seen.add(f); seen.add(-f)
                nsq, distinct, flags, cells = score(b, h, f, C0, cs)
                if distinct and nsq >= 7:
                    best.append((nsq, f, flags))
    return C0, best, stats

if len(sys.argv) > 1 and sys.argv[1] == 'control':
    C0, best, stats = search_pair(31, 17, N=6, verbose=True)
    print("stats:", stats)
    print("section (31,17), C0 =", C0, " hits with >= 7 squares:", [(n_, f_) for n_, f_, fl in best])
    print("contains f = 23/17:", any(f_ in (QQ(23)/17, -QQ(23)/17) for _, f_, _ in best))
