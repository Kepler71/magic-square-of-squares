# -*- coding: utf-8 -*-
"""
Спуск по 3-изогении для E(a,b): y^2 = x(x+M)(x+N), M=b^3(2a+b), N=a^3(a+2b).
Сдвиг x -> x + a^2 b^2 даёт  E: y^2 = x^3 + (alpha x + beta)^2,
alpha = a^2+ab+b^2, beta = a^2 b^2 (a+b)^2, T=(0,beta) — точка порядка 3 (флекс).

phi: E -> E' = E/<T>,  phihat: E' -> E.
  E(Q)/phihat(E'(Q)) --delta--> H^1(Q, mu_3) = Q*/Q*^3,      delta(P) = y - alpha x - beta
  E'(Q)/phi(E(Q))   --delta'-> H^1(Q, Z/3) = (K*/K*^3)^-,  K=Q(sqrt-3), delta'(P') = y' - l(x'), l — касательная в T'
Оценка: 3^{rank+1} = |E'(Q)/phi E(Q)| * |E(Q)/phihat E'(Q)|  <=  |Sel^phi| * |Sel^phihat|.

Локальные образы: нижняя оценка — предъявленные точки; точность доказана по двойственности Тейта:
dim im(delta_v) + dim im(delta'_v) = dim H^1(Q_v, mu_3). Если сумма найденных достигает её — оба образа точны.
Свой код (Fable), 14.09.2026.
"""
import sys, random, itertools
from sage.all import *

random.seed(int(12345))
PREC = 80

def curve_data(a, b):
    a = ZZ(a); b = ZZ(b)
    alpha = a**2 + a*b + b**2
    beta = a**2 * b**2 * (a+b)**2
    # y^2 = x^3 + alpha^2 x^2 + 2 alpha beta x + beta^2
    E = EllipticCurve(QQ, [0, alpha**2, 0, 2*alpha*beta, beta**2])
    T = E(0, beta)
    assert 3*T == E(0)
    phi = E.isogeny(T)
    Ep = phi.codomain()
    return E, Ep, phi, alpha, beta

def tangent_Tprime(Ep):
    """Точка T' в E'[phihat] (x рациональный, y = s*sqrt(-3)) и касательная y = m x + n, m,n in K=Q(sqrt-3).
    Возвращаем (x0, s, m1, n1) где y0 = s*w, m = m1*w, n = n1*w, w=sqrt(-3)."""
    a1,a2,a3,a4,a6 = Ep.ainvs(); assert a1==0 and a3==0
    rts = Ep.division_polynomial(3).roots(QQ)
    cands = []
    for x0,_ in rts:
        f0 = x0**3 + a2*x0**2 + a4*x0 + a6
        # y0^2 = f0, нужно f0 = -3 s^2
        q = -f0/3
        if q.is_square():
            cands.append((x0, q.sqrt()))
    assert len(cands) == 1, cands
    x0, s = cands[0]
    # наклон касательной: (3x0^2+2a2x0+a4)/(2 y0) = D/(2 s w) = D w/(2 s w^2) = -D w/(6 s)
    D = 3*x0**2 + 2*a2*x0 + a4
    m1 = -D/(6*s)
    # y = y0 + m (x - x0) = s w + m1 w x - m1 w x0 = m1 w x + (s - m1 x0) w
    n1 = s - m1*x0
    return x0, s, m1, n1

# ---------- классы кубов над Q_p ----------
def cube_class_Qp(u, p):
    """u — элемент Qp(p) с достаточной точностью. Возвращает вектор над GF(3): (v mod 3, [unit part])."""
    v = u.valuation()
    if v == +Infinity: raise ValueError("нуль")
    unit = u >> v  # u / p^v
    if p == 2:
        return (v % 3,)
    if p == 3:
        assert unit.precision_absolute() >= 2, "нужна точность >= 2 при p=3"
        r = ZZ(unit.residue(2))  # mod 9
        cl = {1:0, 4:1, 7:2, 8:0, 5:1, 2:2}[r]
        return (v % 3, cl)
    if p % 3 == 1:
        assert unit.precision_absolute() >= 1
        r = ZZ(unit.residue())
        Fp = GF(p); g = Fp.multiplicative_generator()
        omega = g**((p-1)//3)
        t = Fp(r)**((p-1)//3)
        cl = [omega**k for k in range(3)].index(t)
        return (v % 3, cl)
    return (v % 3,)

def cube_class_rational_Qp(d, p):
    return cube_class_Qp(Qp(p, PREC)(d), p)

# ---------- классы кубов в K_v = Q_p(sqrt-3), минус-собственная часть ----------
class KvCubes:
    """Локальная группа (K_v*/K_v*^3)^- как векторное пространство над GF(3) с явными координатами.
    p = 1 mod 3: K_v = Q_p (w -> r), вектор из cube_class_Qp (dim 2)
    p = 2 mod 3 (и p=2): K_v = Q_{p^2}, класс вычета в F_{p^2}*/кубы (dim 1)
    p = 3: K_v = Q_3(lambda), lambda^2=-3; U^1/U^4 через (s mod 9, t mod 9) (dim 3, минус-часть dim 2)"""
    def __init__(self, p):
        self.p = p
        if p % 3 == 1:
            R = Qp(p, PREC)
            self.r = R(-3).sqrt()
            self.dim = 2
        elif p == 3:
            self.dim = 3   # объемлющее U^1/U^4; минус-часть имеет dim 2
        elif p == 2:
            F = GF(4, 'z')   # z^2 + z + 1 = 0, z = zeta_3
            self.F = F
            self.z = F.gen()
            self.omega = self.z
            self.dim = 1
        else:
            F = GF(p**2, 'g')
            self.F = F
            self.w = F(-3).sqrt()
            g = F.multiplicative_generator()
            self.omega = g**((p**2-1)//3)
            self.dim = 1
    def cls(self, s, t):
        """класс элемента s + t*sqrt(-3), s,t — рациональные или p-адические (Qp(p))."""
        p = self.p
        if p % 3 == 1:
            u = s + t*self.r
            return cube_class_Qp(u, p)
        elif p == 3:
            R = Qp(3, PREC)
            s = R(s); t = R(t)
            vs = s.valuation() if s != 0 else +Infinity
            vt = t.valuation() if t != 0 else +Infinity
            v = min(2*vs, 2*vt + 1)  # v_lambda
            assert v % 3 == 0, ("валюация не кратна 3", v)
            # делим на lambda^v : lambda^2 = -3, lambda^3 = -3 lambda
            k = v // 3
            # (s + t l) / l^{3k} = (s + t l) / ((-3 l)^k) = (s+tl) / ((-3)^k l^k)
            # l^k: если k чётно, l^k = (-3)^{k/2}; если нечётно, l^k = (-3)^{(k-1)/2} l
            den = (R(-3))**k
            if k % 2 == 0:
                den *= (R(-3))**(k//2); s2, t2 = s/den, t/den
            else:
                den *= (R(-3))**((k-1)//2)
                # (s + t l)/(den * l) = (s + t l) l /(den * l^2) = (s l + t l^2)/(den*(-3)) = (-3t + s l)/(-3 den)
                s2, t2 = (-3*t)/(-3*den), s/(-3*den)
            # теперь единица: s2 должно быть единицей
            assert s2.valuation() == 0 and t2.valuation() >= 0
            s9 = ZZ(s2.residue(2)); t9 = ZZ(t2.residue(2))
            if s9 % 3 == 2:
                s9 = (-s9) % 9; t9 = (-t9) % 9
            return (s9, t9)
        elif p == 2:
            # s + t w = (s+t) + 2t zeta_3,  целый базис {1, zeta_3}
            R = Qp(2, PREC)
            A = R(s) + R(t); B = 2*R(t)
            vA = A.valuation() if A != 0 else +Infinity
            vB = B.valuation() if B != 0 else +Infinity
            v = min(vA, vB)
            assert v % 3 == 0, ("валюация не кратна 3", v)
            A2 = A >> v; B2 = B >> v
            F = self.F
            res = F(ZZ(A2.residue())) + F(ZZ(B2.residue()))*self.z
            assert res != 0
            cl = [self.omega**k for k in range(3)].index(res)
            return (cl,)
        else:
            R = Qp(p, PREC)
            s = R(s); t = R(t)
            vs = s.valuation() if s != 0 else +Infinity
            vt = t.valuation() if t != 0 else +Infinity
            v = min(vs, vt)
            assert v % 3 == 0, ("валюация не кратна 3", v)
            s2 = s >> v; t2 = t >> v
            F = self.F
            res = F(ZZ(s2.residue())) + F(ZZ(t2.residue()))*self.w
            assert res != 0
            tt = res**((p**2-1)//3)
            cl = [self.omega**k for k in range(3)].index(tt)
            return (cl,)

class Subgroup3:
    """Подгруппа в конечной абелевой группе экспоненты 3 — либо векторы над GF(3) (mul = сложение),
    либо пары mod 9 при p=3 (mul = умножение в (Z/9)[lambda])."""
    def __init__(self, p, mod9=False):
        self.p = p; self.mod9 = mod9
        self.elems = {self.one()}
    def one(self):
        if self.mod9: return (1, 0)
        return None  # определяется размерностью при первом добавлении
    def mul(self, x, y):
        if self.mod9:
            s, t = x; s2, t2 = y
            S = (s*s2 - 3*t*t2) % 9; Tt = (s*t2 + t*s2) % 9
            if S % 3 == 2: S = (-S) % 9; Tt = (-Tt) % 9
            return (S, Tt)
        return tuple((u+v) % 3 for u,v in zip(x,y))
    def add(self, x):
        if (not self.mod9) and None in self.elems:
            self.elems = {tuple([0]*len(x))}
        if x in self.elems: return False
        new = set(self.elems)
        for e in self.elems:
            new.add(self.mul(e, x)); new.add(self.mul(self.mul(e, x), x))
        self.elems = new
        return True
    def contains(self, x):
        if (not self.mod9) and None in self.elems:
            return all(c == 0 for c in x)
        return x in self.elems
    def dim(self):
        n = len(self.elems)
        d = 0
        while 3**d < n: d += 1
        assert 3**d == n
        return d

def h1_mu3_dim(p):
    if p == 3: return 2
    if p % 3 == 1: return 2
    return 1

# ---------- поиск точек ----------
def sample_xs(p, n):
    """кандидаты x: целые, а также с знаменателями p^k."""
    xs = []
    for i in range(n):
        k = int(random.choice([0,0,0,1,1,2,3]))
        num = random.randrange(int(-p**6), int(p**6)) if p < 100 else random.randrange(int(-10**8), int(10**8))
        xs.append(QQ(num) / p**k if k else QQ(num))
    return xs

def local_image_E(E, alpha, beta, p, want_dim, extra=0):
    """образ delta: E(Q_p) -> Q_p*/Q_p*^3."""
    R = Qp(p, PREC)
    G = Subgroup3(p)
    # точка T и её кратные: delta(T) = (−2β)^{-1}  (см. записку), delta(2T)=delta(-T) = -2 beta
    G.add(cube_class_Qp(R(-2*beta), p))
    tries = 0
    for x in sample_xs(p, 400 + extra):
        fx = R(x)**3 + (R(alpha)*R(x) + R(beta))**2
        if fx == 0: continue
        if not fx.is_square(): continue
        y = fx.sqrt()
        for yy in (y, -y):
            u = yy - R(alpha)*R(x) - R(beta)
            if u == 0 or u.precision_absolute() - u.valuation() < 3:
                # близко к T: используем v = y + alpha x + beta, delta = v^{-1}
                v = yy + R(alpha)*R(x) + R(beta)
                if v == 0 or v.precision_absolute() - v.valuation() < 3: continue
                c = cube_class_Qp(v, p)
                c = tuple((-ci) % 3 for ci in c) if p != 3 else (( -c[0]) % 3, (-c[1]) % 3)
            else:
                c = cube_class_Qp(u, p)
            G.add(c)
        if G.dim() >= want_dim: break
    return G

def local_image_Ep(Ep, tang, p, want_dim, extra=0):
    """образ delta': E'(Q_p) -> (K_v*/K_v*^3)^-."""
    x0, s, m1, n1 = tang
    a1,a2,a3,a4,a6 = Ep.ainvs()
    R = Qp(p, PREC)
    KV = KvCubes(p)
    G = Subgroup3(p, mod9=(p == 3))
    for x in sample_xs(p, 400 + extra):
        fx = R(x)**3 + R(a2)*R(x)**2 + R(a4)*R(x) + R(a6)
        if fx == 0: continue
        if not fx.is_square(): continue
        y = fx.sqrt()
        for yy in (y, -y):
            # delta' = y - (m1 x + n1) w  : s-часть = y, t-часть = -(m1 x + n1)
            sv = yy; tv = -(R(m1)*R(x) + R(n1))
            try:
                c = KV.cls(sv, tv)
            except AssertionError as e:
                continue
            G.add(c)
        if G.dim() >= want_dim: break
    return G

def bad_primes(E):
    S = set(E.conductor().prime_factors()); S.add(3)
    return sorted(S)

def global_generators_phihat(S):
    return list(S)   # d = prod p^{e_p}

def pi_for_prime(p):
    """p = 1 mod 3: pi = u + v*w с нормой p (u^2+3v^2 = p) или (u+v w)/2 с u^2+3v^2 = 4p. Возвращаем (u,v) рациональные."""
    K = QuadraticField(-3, 'w')
    P = K.ideal(p).prime_factors()[0]
    g = P.gens_reduced()[0]
    u, v = g.list()   # g = u + v*w
    return QQ(u), QQ(v)

def selmer_3isog(a, b, verbose=False):
    E, Ep, phi, alpha, beta = curve_data(a, b)
    S = bad_primes(E)
    tang = tangent_Tprime(Ep)
    # --- локальные образы ---
    imE = {}; imEp = {}
    for p in S:
        n = h1_mu3_dim(p)
        extra = 0
        while True:
            G = local_image_E(E, alpha, beta, p, n, extra)
            Gp = local_image_Ep(Ep, tang, p, n, extra)
            if G.dim() + Gp.dim() == n: break
            if G.dim() + Gp.dim() > n:
                raise RuntimeError(f"({a},{b}) p={p}: сумма размерностей {G.dim()}+{Gp.dim()} > {n} — ОШИБКА ТЕОРИИ/КОДА")
            extra += 2000
            if extra > 20000:
                raise RuntimeError(f"({a},{b}) p={p}: не удалось сертифицировать образы ({G.dim()}+{Gp.dim()}<{n})")
        imE[p] = G; imEp[p] = Gp
    # --- Sel^phihat: d = prod_{p in S} p^{e_p} ---
    sel_hat = []
    for es in itertools.product(range(3), repeat=len(S)):
        d = prod(QQ(p)**e for p,e in zip(S, es))
        ok = all(imE[q].contains(cube_class_rational_Qp(d, q)) for q in S)
        if ok: sel_hat.append(d)
    # --- Sel^phi: theta = zeta3^{e0} prod (pi_p/pibar_p)^{e_p}, p in S, p = 1 mod 3 ---
    P1 = [p for p in S if p % 3 == 1]
    gens = [(QQ(-1)/2, QQ(1)/2)]  # zeta3 = (-1 + w)/2
    for p in P1:
        u, v = pi_for_prime(p)
        # pi/pibar = pi^2 / p = (u^2 - 3 v^2 + 2 u v w)/p
        gens.append(((u*u - 3*v*v)/p, 2*u*v/p))
    KVs = {q: KvCubes(q) for q in S}
    def kmul(x, y):
        return (x[0]*y[0] - 3*x[1]*y[1], x[0]*y[1] + x[1]*y[0])
    sel_phi = []
    for es in itertools.product(range(3), repeat=len(gens)):
        th = (QQ(1), QQ(0))
        for g, e in zip(gens, es):
            for _ in range(e): th = kmul(th, g)
        ok = True
        for q in S:
            c = KVs[q].cls(th[0], th[1])
            if not imEp[q].contains(c): ok = False; break
        if ok: sel_phi.append(es)
    dh = ZZ(len(sel_hat)).log(3); dp = ZZ(len(sel_phi)).log(3)
    assert 3**dh == len(sel_hat) and 3**dp == len(sel_phi)
    bound = dh + dp - 1
    # положительный контроль: delta(T) = 4ab(a+b) in Sel^phihat
    dT = QQ(4*a*b*(a+b))
    # приведём к представителю prod p^{e}, e in {0,1,2}
    dTr = QQ(1)
    for p, e in factor(dT):
        dTr *= QQ(p)**(e % 3)
    assert dTr in sel_hat, ("delta(T) не в Sel^phihat!", a, b, dTr, sel_hat)
    if verbose:
        print(f"({a},{b}) S={S} dims: hat={dh} phi={dp} bound={bound} sel_hat={sel_hat}")
    return dict(a=a, b=b, S=[int(p) for p in S], dim_hat=int(dh), dim_phi=int(dp), bound=int(bound),
                sel_hat=[str(d) for d in sel_hat], sel_phi=[list(map(int,es)) for es in sel_phi],
                imE={int(p): sorted([list(map(int,e)) for e in imE[p].elems if e is not None]) for p in S},
                imEp={int(p): sorted([list(map(int,e)) for e in imEp[p].elems if e is not None]) for p in S})

if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[1] == "test":
    for (a,b) in [(1,2),(2,3),(3,5),(7,2),(1,3),(19,60),(2,5)]:
        r = selmer_3isog(a, b, verbose=True)
