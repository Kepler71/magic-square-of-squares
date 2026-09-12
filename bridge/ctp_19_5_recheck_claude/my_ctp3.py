#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ЧАСТЬ 3.
 (A) ВТОРАЯ, структурно ДРУГАЯ реализация: L = Q[phi]/(phi^3-3I phi+J) как кольцо
     многочленов (а не как Q^3 с интерполяцией). Работает и для НЕрасщепимой L.
 (B) Контроль на СОБСТВЕННОМ примере Фишера 3.4 (571a1, L - кубическое ПОЛЕ):
     сверка m, gamma_1, места, значение спаривания с ОПУБЛИКОВАННЫМИ в статье.
 (C) Пересчёт (19,5) в этой второй реализации.
"""
from fractions import Fraction as F
import math, random
random.seed(1234)
_lib = open('/home/kep/magicKube/bridge/ctp_19_5_recheck_claude/my_ctp_lib.py').read()
exec(_lib.split("# ---- конкретная кривая")[0])

def content_of(q):
    dens = [c.denominator for c in q if c != 0]
    Lc = 1
    for d_ in dens: Lc = Lc*d_//math.gcd(Lc, d_)
    nums = [int(c*Lc) for c in q]
    gg = 0
    for v in nums: gg = math.gcd(gg, abs(v))
    return F(gg, Lc)

def find_local_point(gq, gamq, p, avoid=None):
    cands = [(1,0)] + [(x,1) for x in range(-250, 251)]
    for den in range(2, 14):
        for num in range(-60, 61):
            if math.gcd(abs(num), den) == 1: cands.append((num, den))
    random.shuffle(cands)
    if avoid is not None: cands = [c for c in cands if c != avoid]
    for (x, z) in cands:
        gv = evalq(gq, F(x), F(z))
        if gv == 0: continue
        gm = evalquad(gamq, F(x), F(z))
        if gm == 0: continue
        if is_square_local(gv, p): return (x, z, gv, gm)
    return None

# ---------- арифметика в L = Q[phi]/(phi^3 - 3I phi + J) ----------
class Lalg:
    def __init__(self, I, J):
        self.I = I; self.J = J
        # phi^3 = 3I phi - J
    def red(self, c):
        """c - список коэфф. при phi^0..phi^k -> привести к длине 3"""
        c = [F(t) for t in c]
        while len(c) > 3:
            k = len(c)-1
            top = c.pop()
            # phi^k = phi^(k-3) * phi^3 = phi^(k-3)*(3I phi - J)
            c[k-3] += -self.J*top
            c[k-2] += 3*self.I*top
        while len(c) < 3: c.append(F(0))
        return c
    def mul(self, a, b):
        r = [F(0)]*5
        for i in range(3):
            for j in range(3): r[i+j] += a[i]*b[j]
        return self.red(r)
    def add(self, a, b): return [a[i]+b[i] for i in range(3)]
    def scal(self, k, a): return [F(k)*t for t in a]
    def inv(self, a):
        """обратный в L через расширенный Евклид по многочленам"""
        f = [self.J, -3*self.I, F(0), F(1)]   # phi^3 - 3I phi + J  (коэфф. 0..3)
        def pdeg(p):
            d = len(p)-1
            while d > 0 and p[d] == 0: d -= 1
            return d
        def pdivmod(A, B):
            A = [F(t) for t in A]; db = pdeg(B)
            Q = [F(0)]*(max(len(A)-db, 1))
            while pdeg(A) >= db and any(t != 0 for t in A):
                da = pdeg(A)
                if da < db: break
                co = A[da]/B[db]
                Q[da-db] += co
                for i in range(db+1): A[da-db+i] -= co*B[i]
                if all(t == 0 for t in A): break
            return Q, A
        r0, r1 = f[:], [F(t) for t in a]+[F(0)]
        s0, s1 = [F(1)], [F(0)]
        t0, t1 = [F(0)], [F(1)]
        def psub(A, B):
            n = max(len(A), len(B)); R = [F(0)]*n
            for i,v in enumerate(A): R[i] += v
            for i,v in enumerate(B): R[i] -= v
            return R
        def pmulp(A, B):
            R = [F(0)]*(len(A)+len(B)-1)
            for i,u in enumerate(A):
                for j,v in enumerate(B): R[i+j] += u*v
            return R
        while any(t != 0 for t in r1):
            q, rem = pdivmod(r0[:], r1)
            r0, r1 = r1, rem
            s0, s1 = s1, psub(s0, pmulp(q, s1))
            t0, t1 = t1, psub(t0, pmulp(q, t1))
        c = r0[pdeg(r0)]
        assert pdeg(r0) == 0, "необратимый элемент L"
        return self.red([v/c for v in t0])

def z_L(g, LA):
    a,b,c,d,e = g
    return [F(3*b*b-8*a*c, 3), F(4*a, 3), F(0)]
def G_L(g, LA):
    """G(x,z): список 5 коэфф., каждый - элемент L"""
    h = hessian(g)
    return [[F(h[k],3), F(4*g[k],3), F(0)] for k in range(5)]
def H_L(g, LA):
    G = G_L(g, LA)
    I = LA.I
    # (2/9)(I - phi^2)
    extra = [F(2*I,9), F(0), F(-2,9)]
    c0 = G[0]
    c1 = [v/2 for v in G[1]]
    c2 = LA.add([v/6 for v in G[2]], extra)
    return [c0, c1, c2]

def gamma_L(g1, g2, g3, LA, m):
    """gamma_1 = коэфф. при phi^2 в (z(g2)z(g3)/m) H_1"""
    z2 = z_L(g2, LA); z3 = z_L(g3, LA)
    fac = LA.mul(LA.mul(z2, z3), LA.inv(m))
    H = H_L(g1, LA)
    out = []
    for k in range(3):
        out.append(LA.mul(fac, H[k])[2])
    return out   # [коэфф. при x^2, xz, z^2]

# ================= (B) КОНТРОЛЬ: пример Фишера 3.4, 571a1 =================
print("="*78)
print("КОНТРОЛЬ: СОБСТВЕННЫЙ пример Фишера, Example 3.4 (E = 571a1)")
g1 = [-11, 68, -52, -164, -64]
g2 = [-4, -60, -232, -52, -3]
g3 = [-31, -78, 32, 102, -53]
for nm, g in (("g1",g1),("g2",g2),("g3",g3)):
    print("  %s: I=%d J=%d" % (nm, inv_I(g), inv_J(g)))
I571, J571 = inv_I(g1), inv_J(g1)
print("  статья: I = 44608, J = 18842960 -> совпало:", (I571, J571) == (44608, 18842960))
assert (I571, J571) == (44608, 18842960)
D571 = F(16,27)*(4*I571**3 - J571**2)
print("  Delta =", D571, "   статья: -2^12*571 =", -2**12*571, " -> совпало:", D571 == -2**12*571)
assert D571 == -2**12*571
LA = Lalg(I571, J571)
# кубика X^3-3IX+J неприводима? проверим рациональные корни
def rat_roots(I, J):
    out = []
    N = abs(J)
    ds = [d for d in range(1, int(N**0.5)+2) if N % d == 0]
    cands = set()
    for d in ds:
        cands |= {d, N//d, -d, -(N//d)}
    cands.add(0)
    for c in cands:
        if c**3 - 3*I*c + J == 0: out.append(c)
    return sorted(out)
rr = rat_roots(I571, J571)
print("  рациональные корни X^3-3IX+J:", rr, " -> L", "поле (неприводима)" if not rr else "расщепима")
# m из статьи
m571 = [F(936032,9), F(-8656,9), F(20,9)]
# проверка: m^2 = z(g1)z(g2)z(g3)
lhs = LA.mul(m571, m571)
rhs = LA.mul(LA.mul(z_L(g1,LA), z_L(g2,LA)), z_L(g3,LA))
print("  m из статьи: m^2 == z1 z2 z3 ?", lhs == rhs)
assert lhs == rhs
# тождество G(1,0)G = H^2
def check_GH_L(g, LA):
    G = G_L(g, LA); H = H_L(g, LA)
    lhs = [LA.mul(G[0], G[k]) for k in range(5)]
    A_,B_,C_ = H
    rhs = [LA.mul(A_,A_), LA.scal(2, LA.mul(A_,B_)),
           LA.add(LA.mul(B_,B_), LA.scal(2, LA.mul(A_,C_))),
           LA.scal(2, LA.mul(B_,C_)), LA.mul(C_,C_)]
    return lhs == rhs
print("  тождество G(1,0)G = H^2 в L:", check_GH_L(g1, LA))
assert check_GH_L(g1, LA)
gam571 = gamma_L(g1, g2, g3, LA, m571)
print("  МОЁ gamma_1 =", [str(v) for v in gam571])
print("  СТАТЬЯ:      (4/9)(5x^2-16xz-12z^2) = [20/9, -64/9, -16/3]")
ok571 = gam571 == [F(20,9), F(-64,9), F(-16,3)]
print("  СОВПАЛО:", ok571)
assert ok571

# локальная часть
a571 = F(g2[0])
print("  a = g2(1,0) =", a571)
Sset = set([2,3,5,7]) | primes_of(D571) | primes_of(a571) | primes_of(content_of(gam571))
for c in gam571:
    if c != 0: Sset |= primes_of(c)
print("  мой набор мест:", sorted(Sset), "+ real")
tot = 1; mins = []
for p in sorted(Sset) + ['inf']:
    res = find_local_point(g1, gam571, p)
    assert res is not None, p
    x, z, gv, gm = res
    hh = hilbert(a571, gm, p)
    tot *= hh
    if hh == -1: mins.append(str(p))
    print("    v=%-5s (%s:%s) gamma=%-14s (a,gamma)_v=%+d" % (str(p), x, z, str(gm), hh))
print("  <g1,g2> =", tot, "   минусы:", mins)
print("  СТАТЬЯ: вклад ТОЛЬКО от вещественного места, спаривание НЕТРИВИАЛЬНО -> -1")
print("  локальная точка статьи (15,4): g1(15,4) =", evalq(g1,F(15),F(4)),
      " gamma(15,4) =", evalquad(gam571,F(15),F(4)))
print("  (a,gamma)_real при (15,4) =", hilbert(a571, evalquad(gam571,F(15),F(4)), 'inf'))
print("  => мой результат совпал со статьёй:", tot == -1 and mins == ['inf'])

# ================= (C) (19,5) во ВТОРОЙ реализации =================
print()
print("="*78)
print("(19,5): пересчёт gamma во ВТОРОЙ реализации (L как кольцо многочленов)")
m_, n_ = 19, 5
s_ = F(m_*m_+n_*n_,2); b_ = s_*m_*m_*n_*n_
er = [-b_, -s_*m_**4, -s_*n_**4]
sh = int(-sum(er)/3)
r = [int((e+sh)/16) for e in er]
A_ = r[0]*r[1]+r[0]*r[2]+r[1]*r[2]; B_ = -r[0]*r[1]*r[2]
I = -48*A_; J = -1728*B_
phi = [-12*t for t in r]
LA2 = Lalg(I, J)
G1 = [680904, 3128916, 1361808, -3128916, 680904]
G2 = [718704, -3036516, 1311408, 3288516, 592704]
G3 = [135072, 2334948, 8283996, 4669896, 540288]
# m: строим из компонент (L расщепима), переводим в базис (1,phi,phi^2) интерполяцией
def comp(el):  # элемент L (базис 1,phi,phi^2) -> три компоненты
    return [el[0]+el[1]*p+el[2]*p*p for p in phi]
def interp(v):  # три значения -> базис (1,phi,phi^2)
    Dd = [(phi[0]-phi[1])*(phi[0]-phi[2]), (phi[1]-phi[0])*(phi[1]-phi[2]), (phi[2]-phi[0])*(phi[2]-phi[1])]
    c = [F(0),F(0),F(0)]
    for i in range(3):
        j, k = [t for t in range(3) if t != i]
        # (X-phi_j)(X-phi_k) = X^2 - (pj+pk)X + pj pk
        pj, pk = phi[j], phi[k]
        coef = F(v[i], 1)/Dd[i] if not isinstance(v[i], F) else v[i]/Dd[i]
        c[2] += coef; c[1] += coef*(-(pj+pk)); c[0] += coef*(pj*pk)
    return c
zc = [ [F(4*g[0]*p + 3*g[1]**2 - 8*g[0]*g[2], 3) for p in phi] for g in (G1,G2,G3)]
prod = [zc[0][i]*zc[1][i]*zc[2][i] for i in range(3)]
def sqrtF2(x):
    a = sqrtint(x.numerator); b = sqrtint(x.denominator)
    assert a is not None and b is not None
    return F(a,b)
mcomp = [sqrtF2(v) for v in prod]
mL = interp(mcomp)
print("  m в базисе (1,phi,phi^2) =", [str(v) for v in mL])
assert LA2.mul(mL, mL) == LA2.mul(LA2.mul(z_L(G1,LA2), z_L(G2,LA2)), z_L(G3,LA2))
print("  проверка m^2 = z1z2z3 в L: ок")
print("  тождество G(1,0)G=H^2 в L:", check_GH_L(G1, LA2))
gamP = gamma_L(G1, G2, G3, LA2, mL)
print("  gamma (реализация 2, многочлены) =", [str(v) for v in gamP])
print("  gamma (реализация 1, интерполяция) = ['132076/3', '529592/3', '-775796/3']")
print("  СОВПАЛО:", [str(v) for v in gamP] == ['132076/3', '529592/3', '-775796/3'])

aV = F(G2[0])
Delta2 = F(16,27)*(4*I**3 - J*J)
Sset2 = set([2,3,5,7]) | primes_of(Delta2) | primes_of(aV) | primes_of(content_of(gamP))
for c in gamP:
    if c != 0: Sset2 |= primes_of(c)
tot2 = 1; mins2 = []
for p in sorted(Sset2) + ['inf']:
    res = find_local_point(G1, gamP, p)
    assert res is not None, p
    x, z, gv, gm = res
    hh = hilbert(aV, gm, p)
    tot2 *= hh
    if hh == -1: mins2.append(str(p))
print("  <g1,g2> (реализация 2) =", tot2, "  минусы:", mins2)
