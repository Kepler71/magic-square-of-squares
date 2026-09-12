#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""КОНТРОЛЬ на ложные пропуски/срабатывания: кривые с ПОЛНЫМ 2-кручением,
рангом 0 и НЕЗАВИСИМО известной Sha = 4 (Cremona/BSD, ранг 0 по Колывагину).
Тогда dim Sel^2 = 0 + 2 + 2 = 4, Sha[2]=(Z/2)^2, и спаривание Касселса-Тейта
ОБЯЗАНО иметь ранг 2: ядро = образ 2-кручения.
Моя реализация формулы Фишера (та же, что для (19,5)) применяется к квартикам
из ell2cover. Sage использован ТОЛЬКО чтобы получить квартики и эталонные rank/Sha.
"""
from fractions import Fraction as F
import math, random, json, itertools
random.seed(11)

L = open('/home/kep/magicKube/bridge/ctp_19_5_recheck_claude/my_ctp_lib.py').read()
# оставляем только generic-часть (до «конкретная кривая»)
generic = L.split("# ---- конкретная кривая")[0]
exec(generic)

class Fisher:
    def __init__(self, I, J, phi):
        self.I = I; self.J = J; self.phi = phi
        self.D = [(phi[0]-phi[1])*(phi[0]-phi[2]),
                  (phi[1]-phi[0])*(phi[1]-phi[2]),
                  (phi[2]-phi[0])*(phi[2]-phi[1])]
        self.Delta = F(16,27)*(4*I**3 - J*J)
        for p in phi: assert p**3 - 3*I*p + J == 0
    def z(self, g):
        a,b,c,d,e = g
        return [F(4*a*ph + 3*b*b - 8*a*c, 3) for ph in self.phi]
    def G(self, g, i):
        h = hessian(g)
        return [F(4*self.phi[i]*g[k] + h[k], 3) for k in range(5)]
    def H(self, g, i):
        G = self.G(g, i)
        return [G[0], F(G[1],2), F(G[2],6) + F(2,9)*(self.I - self.phi[i]**2)]
    def gamma(self, gA, gB, gC):
        zA = self.z(gA); zB = self.z(gB); zC = self.z(gC)
        mm = [sqrtF(zA[i]*zB[i]*zC[i]) for i in range(3)]
        out = [F(0),F(0),F(0)]
        for i in range(3):
            Hf = self.H(gA, i); coef = mm[i]/zA[i]
            for k in range(3): out[k] += coef*Hf[k]/self.D[i]
        return out
    def places(self, gA, gamq, aval):
        S = set([2,3,5,7]) | primes_of(self.Delta) | primes_of(aval)
        for c in gA:
            if c != 0: S |= set(factor(F(c).denominator))
        for c in gamq:
            if c != 0: S |= set(factor(c.denominator))
        S |= primes_of(content_of(gamq))
        for c in gamq:
            if c != 0: S |= primes_of(c)
        return sorted(S)
    def pair(self, gA, gB, gC):
        gamq = self.gamma(gA, gB, gC)
        aval = F(gB[0])
        if aval == 0: return None, None
        S = self.places(gA, gamq, aval)
        tot = 1; mn = []
        for p in S + ['inf']:
            res = find_local_point(gA, gamq, p)
            if res is None: return None, None
            x, z, gv, gm = res
            hh = hilbert(aval, gm, p)
            tot *= hh
            if hh == -1: mn.append(str(p))
        return tot, mn

def sqrtF(x):
    a = sqrtint(x.numerator); b = sqrtint(x.denominator)
    assert a is not None and b is not None, ("не квадрат в L", x)
    return F(a, b)

def content_of(q):
    dens = [c.denominator for c in q if c != 0]
    Lc = 1
    for d_ in dens: Lc = Lc*d_//math.gcd(Lc, d_)
    nums = [int(c*Lc) for c in q]
    gg = 0
    for v in nums: gg = math.gcd(gg, abs(v))
    return F(gg, Lc)

def rescale(g, I, J):
    """g -> mu^2 g с инвариантами ровно I,J; mu^2 обязан быть квадратом рационального"""
    Ig = inv_I(g); Jg = inv_J(g)
    if Ig == I and Jg == J: return g, F(1)
    t = F(I, Ig)                 # = mu^4
    r2 = None
    # mu^2 = sqrt(t)
    a = sqrtint(t.numerator); b = sqrtint(t.denominator)
    if a is None or b is None: return None, None
    mu2 = F(a, b)
    g2 = [mu2*c for c in g]
    if inv_I(g2) == I and inv_J(g2) == J:
        return [c for c in g2], mu2
    return None, None

data = json.load(open('/home/kep/magicKube/bridge/ctp_19_5_recheck_claude/ctrl_curves.json'))
for lbl, d in data.items():
    print("="*78)
    print("КРИВАЯ", lbl, " rank =", d['rank'], " Sha_an =", d['sha'])
    I = int(d['I']); J = int(d['J'])
    rts = [int(t) for t in d['roots']]
    phi = [-12*t for t in rts]
    Fi = Fisher(I, J, phi)
    print("  I,J =", I, J, "  корни", rts, " phi", phi)
    def act(g, M):
        a_,b_,c_,d_ = M
        def pmul(P,Q):
            R=[F(0)]*(len(P)+len(Q)-1)
            for i,p in enumerate(P):
                for j,q in enumerate(Q): R[i+j]+=p*q
            return R
        X=[F(a_),F(c_)]; Z=[F(b_),F(d_)]
        def powp(P,k):
            R=[F(1)]
            for _ in range(k): R=pmul(R,P)
            return R
        acc=[F(0)]*5
        for k in range(5):
            term=pmul(powp(X,4-k),powp(Z,k))
            term=[F(g[k])*t for t in term]
            for i,t in enumerate(term): acc[i]+=t
        return acc
    def make_unit(g):
        """Fisher, стр.3: заменой координат добиваемся, что z(g) - единица в L."""
        if all(v != 0 for v in Fi.z(g)): return g, (1,0,0,1)
        for _ in range(500):
            al = random.randint(-9,9); be = random.randint(-9,9)
            if math.gcd(abs(al),abs(be)) != 1: continue
            # дополняем (al,be) до матрицы с det 1
            gg,u,v = 1,0,0
            # расширенный Евклид: u*al+v*be=1
            def egcd(a,b):
                if b==0: return (a,1,0)
                g_,x_,y_=egcd(b,a%b); return (g_,y_,x_-(a//b)*y_)
            g_,u,v = egcd(al,be)
            if g_ != 1:
                if g_ == -1: u,v = -u,-v
                else: continue
            M = (al, -v, be, u)   # det = al*u + v*be = 1
            assert M[0]*M[3]-M[1]*M[2] == 1
            gn = act(g, M)
            if all(x != 0 for x in Fi.z(gn)) and inv_I(gn)==I and inv_J(gn)==J:
                return gn, M
        return None, None
    qs = []
    for qq in d['quartics']:
        g = [F(t) for t in qq]
        g2, mu2 = rescale(g, I, J)
        if g2 is None:
            print("  !! квартика не масштабируется квадратом:", qq); continue
        zz = Fi.z(g2)
        note = ""
        if any(v == 0 for v in zz):
            g2, M = make_unit(g2)
            note = "  [z был делителем нуля -> замена координат %s]" % (M,)
            assert g2 is not None
        qs.append(g2)
        print("   квартика", [str(c) for c in qq], " * ", mu2, " -> классы z =",
              [sqclass(v) for v in Fi.z(g2)], note)
    # образ 2-кручения
    def dmap(i):
        v = [None]*3
        for j in range(3):
            if j != i: v[j] = sqclass(rts[i]-rts[j])
        v[i] = sqclass((rts[i]-rts[(i+1)%3])*(rts[i]-rts[(i+2)%3]))
        return tuple(v)
    tors = {(1,1,1), dmap(0), dmap(1), dmap(2)}
    print("  образ E(Q)/2E(Q) (ранг 0 => только 2-кручение):", sorted(tors))
    cls = [tuple(sqclass(v) for v in Fi.z(g)) for g in qs]
    print("  классы квартик:", cls)
    def mulcl(A, B): return tuple(sqclass(A[i]*B[i]) for i in range(3))
    # достроим недостающие классы: для пары (i,j) нужна квартика класса cls_i*cls_j.
    # берём собственно эквивалентные варианты уже имеющихся (класс не меняется) -
    # если нужного класса нет, пара пропускается. Проверим заодно структуру группы.
    grp = set()
    allc = set(cls) | tors
    changed = True
    while changed:
        changed = False
        for A_ in list(allc):
            for B_ in list(allc):
                C_ = mulcl(A_, B_)
                if C_ not in allc: allc.add(C_); changed = True
    print("  порождённая классами группа: |G| =", len(allc), " (dim Sel^2 ожидается 4 => 16)")
    # спаривания по всем парам, для которых найдётся третья квартика нужного класса
    print("  --- спаривания ---")
    nz = 0
    for i in range(len(qs)):
        for j in range(len(qs)):
            need = mulcl(cls[i], cls[j])
            k = None
            for t in range(len(qs)):
                if cls[t] == need: k = t; break
            if k is None: continue
            try:
                v, mn = Fi.pair(qs[i], qs[j], qs[k])
            except AssertionError as ex:
                print("    <%d,%d>: z1z2z3 не квадрат -> %s" % (i, j, ex)); continue
            if v is None: continue
            tag = ""
            if cls[i] in tors or cls[j] in tors: tag = "  [один из аргументов в образе E(Q)/2E(Q) -> ОБЯЗАНО +1]"
            print("    <g%d,g%d> = %+d  минусы %s%s" % (i, j, v, mn, tag))
            if v == -1: nz += 1
            if tag and v != 1:
                print("      *** ПРОТИВОРЕЧИЕ: ненулевое спаривание с элементом образа ***")
    print("  ненулевых значений:", nz, " (ОБЯЗАНО быть > 0, т.к. Sha[2]=(Z/2)^2 и CTP невырождено)")
