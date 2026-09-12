# -*- coding: utf-8 -*-
"""общая библиотека: всё выведено заново, только fractions/math"""
from fractions import Fraction as F
import math, random

def sqrtint(n):
    if n < 0: return None
    rr = math.isqrt(n)
    return rr if rr*rr == n else None

def legendre(a, p):
    a %= p
    if a == 0: return 0
    rr = pow(a, (p-1)//2, p)
    return 1 if rr == 1 else -1

def val_unit(x, p):
    assert x != 0
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0: num //= p; v += 1
    while den % p == 0: den //= p; v -= 1
    return v, F(num, den)

def unit_mod(u, p, k):
    M = p**k
    return (u.numerator % M) * pow(u.denominator % M, -1, M) % M

def hilbert(a, b, p):
    a = F(a); b = F(b)
    assert a != 0 and b != 0
    if p == 'inf':
        return -1 if (a < 0 and b < 0) else 1
    if p == 2:
        al, u = val_unit(a, 2); be, v = val_unit(b, 2)
        um = unit_mod(u, 2, 3); vm = unit_mod(v, 2, 3)
        eu = ((um-1)//2) % 2; ev = ((vm-1)//2) % 2
        wu = ((um*um-1)//8) % 2; wv = ((vm*vm-1)//8) % 2
        return -1 if (eu*ev + al*wv + be*wu) % 2 else 1
    al, u = val_unit(a, p); be, v = val_unit(b, p)
    um = unit_mod(u, p, 1); vm = unit_mod(v, p, 1)
    s = 1
    if (al*be) % 2 == 1 and (p % 4) == 3: s = -s
    if be % 2 == 1: s *= legendre(um, p)
    if al % 2 == 1: s *= legendre(vm, p)
    return s

def is_square_local(x, p):
    if x == 0: return False
    if p == 'inf': return x > 0
    v, u = val_unit(x, p)
    if v % 2: return False
    if p == 2: return unit_mod(u, 2, 3) == 1
    return legendre(unit_mod(u, p, 1), p) == 1

def factor(n):
    n = abs(int(n)); f = {}
    d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n,0)+1
    return f

def sqclass(x):
    x = F(x)
    nn = x.numerator*x.denominator
    out = -1 if nn < 0 else 1
    for p, k in factor(nn).items():
        if k % 2: out *= p
    return out

def primes_of(x):
    x = F(x)
    return set(factor(x.numerator)) | set(factor(x.denominator))

def inv_I(g):
    a,b,c,d,e = g; return 12*a*e - 3*b*d + c*c
def inv_J(g):
    a,b,c,d,e = g; return 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c*c*c
def hessian(g):
    a,b,c,d,e = g
    return [3*b*b-8*a*c, 4*(b*c-6*a*d), 2*(2*c*c-24*a*e-3*b*d), 4*(c*d-6*b*e), 3*d*d-8*c*e]
def evalq(g, x, z):
    a,b,c,d,e = g
    return a*x**4 + b*x**3*z + c*x*x*z*z + d*x*z**3 + e*z**4
def evalquad(q, x, z):
    return q[0]*x*x + q[1]*x*z + q[2]*z*z

# ---- конкретная кривая (19,5), выводится заново ----
m_, n_ = 19, 5
s_ = F(m_*m_+n_*n_, 2)
b_ = s_*m_*m_*n_*n_
e_roots = [-b_, -s_*m_**4, -s_*n_**4]
shift = int(-sum(e_roots)/3)
r = [int((e+shift)/16) for e in e_roots]
A = r[0]*r[1]+r[0]*r[2]+r[1]*r[2]
B = -r[0]*r[1]*r[2]
I = -48*A
J = -1728*B
phi = [-12*x for x in r]
Delta = int(F(16,27)*(4*I**3 - J*J))
D = [(phi[0]-phi[1])*(phi[0]-phi[2]), (phi[1]-phi[0])*(phi[1]-phi[2]), (phi[2]-phi[0])*(phi[2]-phi[1])]

def z_inv(g):
    a,bb,c,d,e = g
    return [F(4*a*ph + 3*bb*bb - 8*a*c, 3) for ph in phi]
def G_form_g(g, i):
    h = hessian(g)
    return [F(4*phi[i]*g[k] + h[k], 3) for k in range(5)]
def H_form_g(g, i):
    G = G_form_g(g, i)
    return [G[0], F(G[1],2), F(G[2],6) + F(2,9)*(I - phi[i]**2)]
def sqrtF(x):
    a = sqrtint(x.numerator); bq = sqrtint(x.denominator)
    assert a is not None and bq is not None, ("не квадрат", x)
    return F(a, bq)
def gamma_form_g(gA, gB, gC, signs=(1,1,1)):
    zA = z_inv(gA); zB = z_inv(gB); zC = z_inv(gC)
    mm = [signs[i]*sqrtF(zA[i]*zB[i]*zC[i]) for i in range(3)]
    out = [F(0),F(0),F(0)]
    for i in range(3):
        H = H_form_g(gA, i); coef = mm[i]/zA[i]
        for k in range(3): out[k] += coef*H[k]/D[i]
    return out, mm
def content_of(q):
    dens = [c.denominator for c in q if c != 0]
    L_ = 1
    for d_ in dens: L_ = L_*d_//math.gcd(L_, d_)
    nums = [int(c*L_) for c in q]
    gg = 0
    for v in nums: gg = math.gcd(gg, abs(v))
    return F(gg, L_)
def place_set_g(gA, gamq, aval):
    S = set([2,3,5,7]) | primes_of(Delta) | primes_of(aval)
    for c in gA:
        if c != 0: S |= set(factor(F(c).denominator))
    for c in gamq:
        if c != 0: S |= set(factor(c.denominator))
    S |= primes_of(content_of(gamq))
    return sorted(S)
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
