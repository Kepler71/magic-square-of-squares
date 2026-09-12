#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
НЕЗАВИСИМЫЙ пересчёт спаривания Кассельса-Тейта для G1 (m,n)=(19,5).
Всё выведено заново: s,b, корни E, модель M, I,J, phi, Delta, z(g), H, gamma,
локальные точки (СВОИ, найденные собственным поиском), символы Гильберта (своя реализация),
набор мест (свой, по Remark 3.3).
Никакого Sage, PARI, ни кода проекта. Только fractions/math.
"""
from fractions import Fraction as F
import math, random, sys

random.seed(20260912)

def sqrtint(n):
    if n < 0: return None
    r = math.isqrt(n)
    return r if r*r == n else None

# ---------------------------------------------------------------- Гильберт
def legendre(a, p):
    a %= p
    if a == 0: return 0
    r = pow(a, (p-1)//2, p)
    return 1 if r == 1 else -1

def val_unit(x, p):
    """x - Fraction != 0. вернуть (v, u) с x = p^v * u, u - p-адическая единица (Fraction)."""
    assert x != 0
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0:
        num //= p; v += 1
    while den % p == 0:
        den //= p; v -= 1
    return v, F(num, den)

def unit_mod(u, p, k):
    """u - Fraction, p-адическая единица. вернуть u mod p^k как целое."""
    M = p**k
    return (u.numerator % M) * pow(u.denominator % M, -1, M) % M

def hilbert(a, b, p):
    """символ Гильберта (a,b)_p, p - простое или строка 'inf'. a,b - Fraction !=0."""
    a = F(a); b = F(b)
    assert a != 0 and b != 0
    if p == 'inf':
        return -1 if (a < 0 and b < 0) else 1
    if p == 2:
        al, u = val_unit(a, 2)
        be, v = val_unit(b, 2)
        um = unit_mod(u, 2, 3)   # mod 8
        vm = unit_mod(v, 2, 3)
        eu = ((um - 1)//2) % 2
        ev = ((vm - 1)//2) % 2
        wu = ((um*um - 1)//8) % 2
        wv = ((vm*vm - 1)//8) % 2
        e = (eu*ev + al*wv + be*wu) % 2
        return -1 if e else 1
    al, u = val_unit(a, p)
    be, v = val_unit(b, p)
    um = unit_mod(u, p, 1)
    vm = unit_mod(v, p, 1)
    s = 1
    if (al*be) % 2 == 1 and (p % 4) == 3:
        s = -s
    if be % 2 == 1:
        s *= legendre(um, p)
    if al % 2 == 1:
        s *= legendre(vm, p)
    return s

def is_square_local(x, p):
    """x - Fraction. является ли x квадратом в Q_p (x!=0)."""
    if x == 0: return False
    if p == 'inf': return x > 0
    v, u = val_unit(x, p)
    if v % 2: return False
    if p == 2:
        return unit_mod(u, 2, 3) == 1
    return legendre(unit_mod(u, p, 1), p) == 1

# --- самотест символа Гильберта ---
def factor(n):
    n = abs(n); f = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0)+1; n//=d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n, 0)+1
    return f

def selftest_hilbert(trials=400):
    bad = 0
    for _ in range(trials):
        a = F(random.randint(-200, 200) or 3, random.randint(1, 40))
        b = F(random.randint(-200, 200) or 5, random.randint(1, 40))
        if a == 0 or b == 0: continue
        ps = set([2]) | set(factor(a.numerator)) | set(factor(a.denominator)) \
             | set(factor(b.numerator)) | set(factor(b.denominator))
        prod = hilbert(a, b, 'inf')
        for p in ps: prod *= hilbert(a, b, p)
        if prod != 1: bad += 1; print("  RECIPROCITY FAIL", a, b)
        for p in list(ps)+['inf']:
            if hilbert(a, b, p) != hilbert(b, a, p): bad += 1; print("  SYM FAIL", a, b, p)
            if hilbert(a, -a, p) != 1: bad += 1; print("  (a,-a) FAIL", a, p)
            if a != 1 and hilbert(a, 1-a, p) != 1: bad += 1; print("  (a,1-a) FAIL", a, p)
    return bad

# ---------------------------------------------------------------- инварианты квартик
def inv_I(g):
    a,b,c,d,e = g
    return 12*a*e - 3*b*d + c*c

def inv_J(g):
    a,b,c,d,e = g
    return 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c*c*c

def hessian(g):
    a,b,c,d,e = g
    return [3*b*b - 8*a*c,
            4*(b*c - 6*a*d),
            2*(2*c*c - 24*a*e - 3*b*d),
            4*(c*d - 6*b*e),
            3*d*d - 8*c*e]

def evalq(g, x, z):
    a,b,c,d,e = g
    return a*x**4 + b*x**3*z + c*x*x*z*z + d*x*z**3 + e*z**4

def evalquad(q, x, z):
    """q = [A,B,C] -> A x^2 + B xz + C z^2"""
    return q[0]*x*x + q[1]*x*z + q[2]*z*z

# ---------------------------------------------------------------- G1-семейство, пара (19,5)
m, n = 19, 5
assert math.gcd(m, n) == 1
s = F(m*m + n*n, 2)
b = s * m*m * n*n
print("=== 1. Семейство G1, (m,n) = (19,5) ===")
print("  s =", s, "  b =", b)
assert s == 193 and b == 1741825
e_roots = [-b, -s*m**4, -s*n**4]
print("  корни E (e1,e2,e3) =", [str(x) for x in e_roots])

# необходимый квадратный класс: X-e1=(mn)^2 F4, X-e2 = s m^2 F0, X-e3 = s n^2 F8, X = b t^2
# ПРОВЕРКА КАК ТОЖДЕСТВ МНОГОЧЛЕНОВ ОТ t (коэффициенты при t^0, t^2)
F0 = (m*m, n*n)      # m^2 + n^2 t^2  -> (const, t^2)
F4 = (s, s)          # s(1+t^2)
F8 = (n*n, m*m)
X  = (F(0), b)       # b t^2
def polysub(P, c): return (P[0]-c, P[1])
def polymul(k, P):  return (k*P[0], k*P[1])
assert polysub(X, e_roots[0]) == polymul(F((m*n)**2), F4), "X-e1"
assert polysub(X, e_roots[1]) == polymul(s*m*m, F0), "X-e2"
assert polysub(X, e_roots[2]) == polymul(s*n*n, F8), "X-e3"
print("  тождества X-e_i проверены как многочлены от t  -> delta = (1, s, s) = (1,193,193)")
assert sqrtint(193) is None
delta = [1, 193, 193]

# ---------------------------------------------------------------- модель M
print("=== 2. Модель M (вывожу сам) ===")
S1 = sum(e_roots)
shift = -S1/3
assert shift.denominator == 1
shift = int(shift)
print("  сдвиг X -> X + %d  (центрирование)" % shift)
u2 = 16
r = [(e + shift)/u2 for e in e_roots]
assert all(x.denominator == 1 for x in r)
r = [int(x) for x in r]
print("  корни M:", r, " сумма =", sum(r))
assert sum(r) == 0
A = r[0]*r[1] + r[0]*r[2] + r[1]*r[2]
B = -r[0]*r[1]*r[2]
print("  M: y^2 = x^3 + (%d) x + (%d)" % (A, B))
assert A == -766425627513 and B == 254371088344251312
I = -48*A
J = -1728*B
print("  I = c4 =", I)
print("  J = 2c6 =", J)
assert I == 36788430120624 and J == -439553240658866267136
phi = [-12*x for x in r]
print("  phi =", phi)
for p_ in phi:
    assert p_**3 - 3*I*p_ + J == 0, "phi не корень X^3-3IX+J"
print("  все phi_i - корни X^3 - 3 I X + J: ок  => L = Q x Q x Q")
Delta = F(16, 27)*(4*I**3 - J*J)
assert Delta.denominator == 1
Delta = int(Delta)
print("  Delta =", Delta)
print("  факторизация Delta:", factor(Delta))

D = [ (phi[0]-phi[1])*(phi[0]-phi[2]),
      (phi[1]-phi[0])*(phi[1]-phi[2]),
      (phi[2]-phi[0])*(phi[2]-phi[1]) ]

# ---------------------------------------------------------------- квартики (ДАННЫЕ сертификата)
g1 = [680904, 3128916, 1361808, -3128916, 680904]
g2 = [718704, -3036516, 1311408, 3288516, 592704]
g3 = [135072, 2334948, 8283996, 4669896, 540288]
print("=== 3. Квартики сертификата: проверка инвариантов ===")
for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
    print("  %s: I=%s  J=%s   совпало: %s / %s" % (nm, inv_I(g), inv_J(g), inv_I(g) == I, inv_J(g) == J))
    assert inv_I(g) == I and inv_J(g) == J

# ---------------------------------------------------------------- z(g), H(g)
def z_inv(g):
    """z(g) в L = Q^3, компонента i: (4 a phi_i + 3b^2 - 8ac)/3"""
    a,bb,c,d,e = g
    return [F(4*a*ph + 3*bb*bb - 8*a*c, 3) for ph in phi]

def G_form(g, i):
    """G = (1/3)(4 phi g + h) в компоненте i, как список 5 коэффициентов"""
    h = hessian(g)
    return [F(4*phi[i]*g[k] + h[k], 3) for k in range(5)]

def H_form(g, i):
    G = G_form(g, i)
    return [G[0], F(G[1], 2), F(G[2], 6) + F(2, 9)*(I - phi[i]**2)]

def check_GH(g):
    """тождество G(1,0) G = H^2 как многочленов"""
    for i in range(3):
        G = G_form(g, i); H = H_form(g, i)
        lhs = [G[0]*c for c in G]
        # H^2
        A_,B_,C_ = H
        rhs = [A_*A_, 2*A_*B_, B_*B_+2*A_*C_, 2*B_*C_, C_*C_]
        if lhs != rhs:
            return False
    return True

print("=== 4. z(g) и тождество Гессиана ===")
for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
    ok = check_GH(g)
    print("  %s: G(1,0)*G == H^2 : %s" % (nm, ok))
    assert ok

def sqclass(x):
    """бесквадратное ядро рационального x != 0"""
    x = F(x)
    nn = x.numerator * x.denominator
    sgn = -1 if nn < 0 else 1
    f = factor(nn)
    out = sgn
    for p, k in f.items():
        if k % 2: out *= p
    return out

z1 = z_inv(g1); z2 = z_inv(g2); z3 = z_inv(g3)
print("  z(g1) =", [str(v) for v in z1], " классы:", [sqclass(v) for v in z1])
print("  z(g2) =", [str(v) for v in z2], " классы:", [sqclass(v) for v in z2])
print("  z(g3) =", [str(v) for v in z3], " классы:", [sqclass(v) for v in z3])
cls1 = [sqclass(v) for v in z1]
print("  требуемый класс delta =", delta, " -> z(g1) в нужном классе:", cls1 == delta)
prodz = [z1[i]*z2[i]*z3[i] for i in range(3)]
print("  z1*z2*z3 классы:", [sqclass(v) for v in prodz], " (должны быть все 1)")
assert all(sqclass(v) == 1 for v in prodz), "z1z2z3 не квадрат в L"

def sqrtF(x):
    """квадратный корень Fraction (должен быть точным)"""
    a = sqrtint(x.numerator); bq = sqrtint(x.denominator)
    assert a is not None and bq is not None, ("не квадрат", x)
    return F(a, bq)

# ---------------------------------------------------------------- gamma
def gamma_form(gA, gB, gC, signs=(1,1,1)):
    """gamma = коэффициент при phi^2 в (m/z(gA)) * H_{gA},
       где m^2 = z(gA)z(gB)z(gC). Возвращает [A,B,C] (квадратичная форма в x,z)."""
    zA = z_inv(gA); zB = z_inv(gB); zC = z_inv(gC)
    mm = [signs[i]*sqrtF(zA[i]*zB[i]*zC[i]) for i in range(3)]
    out = [F(0), F(0), F(0)]
    for i in range(3):
        H = H_form(gA, i)
        coef = mm[i]/zA[i]
        for k in range(3):
            out[k] += coef*H[k]/D[i]
    return out, mm

gam, mvals = gamma_form(g1, g2, g3)
print("=== 5. gamma (моя интерполяция) ===")
print("  m =", [str(v) for v in mvals])
print("  gamma =", [str(v) for v in gam])

# ---------------------------------------------------------------- набор мест по Remark 3.3
def primes_of(x):
    x = F(x)
    return set(factor(x.numerator)) | set(factor(x.denominator))

def content_primes(q):
    """простые, делящие содержание квадратичной формы q (с рац. коэфф.)"""
    from functools import reduce
    dens = [c.denominator for c in q if c != 0]
    L_ = 1
    for d_ in dens: L_ = L_*d_//math.gcd(L_, d_)
    nums = [int(c*L_) for c in q]
    gg = 0
    for v in nums: gg = math.gcd(gg, abs(v))
    cont = F(gg, L_)
    return primes_of(cont), cont

def place_set(gA, gamq, aval):
    S = set([2, 3, 5, 7])            # N(v) < 11
    S |= primes_of(Delta)            # v | Delta(g1)
    S |= primes_of(aval)             # g2(1,0) должно быть v-единицей
    for c in gA:
        if c != 0: S |= set(factor(F(c).denominator))
    for c in gamq:
        if c != 0: S |= set(factor(c.denominator))
    cp, cont = content_primes(gamq)
    S |= cp
    return sorted(S), cont

a_val = F(g2[0])
places, cont = place_set(g1, gam, a_val)
print("=== 6. Набор мест (мой, по Remark 3.3) ===")
print("  primes(Delta) =", sorted(primes_of(Delta)))
print("  primes(a=g2(1,0)=%s) =" % a_val, sorted(primes_of(a_val)))
print("  content(gamma) =", cont, " простые:", sorted(primes_of(cont)))
print("  НЕОБХОДИМЫЙ набор:", places, "+ real")

# консервативный сверх-набор: добавим числители gamma и g1
S_big = set(places)
for c in gam:
    if c != 0: S_big |= primes_of(c)
for c in g1:
    if c != 0: S_big |= primes_of(c)
S_big = sorted(S_big)
print("  сверх-набор:", S_big, "+ real")

# ---------------------------------------------------------------- локальные точки (мой поиск)
def find_local_point(gq, gamq, p, avoid=None):
    """ищу (x:z) из Q с g(x,z) квадрат !=0 в Q_p и gamma(x,z)!=0"""
    cands = []
    # z=0 лист
    cands.append((1, 0))
    R = 200 if p != 'inf' else 60
    for x in range(-R, R+1):
        cands.append((x, 1))
    for den in range(2, 12):
        for num in range(-40, 41):
            if math.gcd(abs(num), den) == 1:
                cands.append((num, den))
    random.shuffle(cands)
    if avoid is not None:
        cands = [c for c in cands if c != avoid]
    for (x, z) in cands:
        gv = evalq(gq, F(x), F(z))
        if gv == 0: continue
        gm = evalquad(gamq, F(x), F(z))
        if gm == 0: continue
        if is_square_local(gv, p):
            return (x, z, gv, gm)
    return None

def run_pairing(gq, gamq, aval, S, label, avoid_map=None, verbose=True):
    tot = 1
    rows = []
    for p in list(S) + ['inf']:
        av = None if avoid_map is None else avoid_map.get(p)
        res = find_local_point(gq, gamq, p, avoid=av)
        if res is None:
            print("  !!! НЕТ локальной точки в месте", p, " -> квартика НЕ ELS?")
            return None, None
        x, z, gv, gm = res
        h = hilbert(aval, gm, p)
        tot *= h
        rows.append((p, x, z, gv, gm, h))
        if verbose:
            print("    v=%-6s (x:z)=(%s:%s)  gamma=%-24s  (a,gamma)_v = %+d" % (str(p), x, z, str(gm), h))
    if verbose:
        print("  ПРОИЗВЕДЕНИЕ [%s] = %+d   минусы в: %s" % (label, tot, [str(rr[0]) for rr in rows if rr[5] == -1]))
    return tot, rows

print()
print("=== 7. Самотест символа Гильберта ===")
bad = selftest_hilbert(400)
print("  неудач:", bad)
assert bad == 0

print()
print("=== 8. ELS всех трёх квартик (мой поиск свидетелей) ===")
els_places = sorted(set([2, 3]) | primes_of(Delta))
for nm, g in (("g1", g1), ("g2", g2), ("g3", g3)):
    ok = True
    wit = []
    for p in els_places + ['inf']:
        found = None
        for (x, z) in [(1, 0)] + [(x, 1) for x in range(-300, 301)] + \
                      [(nu, de) for de in range(2, 10) for nu in range(-60, 61) if math.gcd(abs(nu), de) == 1]:
            gv = evalq(g, F(x), F(z))
            if gv != 0 and is_square_local(gv, p):
                found = (x, z, gv); break
        if found is None:
            ok = False; wit.append((p, None))
        else:
            wit.append((p, found[:2]))
    print("  %s ELS на %s : %s   свидетели: %s" % (nm, els_places+['inf'], ok, wit))
    assert ok

print()
print("=== 9. ОСНОВНОЙ ПЕРЕСЧЁТ <g1,g2> (сырая gamma, мои локальные точки) ===")
tot1, rows1 = run_pairing(g1, gam, a_val, S_big, "<g1,g2> сырая gamma, сверх-набор")

print()
print("=== 9b. Тот же расчёт по МИНИМАЛЬНОМУ набору мест (Remark 3.3) ===")
tot1b, _ = run_pairing(g1, gam, a_val, places, "<g1,g2> сырая gamma, минимальный набор")

print()
print("=== 9c. Другие локальные точки ===")
avoid = {rr[0]: (rr[1], rr[2]) for rr in rows1}
tot1c, _ = run_pairing(g1, gam, a_val, S_big, "<g1,g2> другие точки", avoid_map=avoid)

print()
print("=== 10. Нормированная gamma (деление на content) ===")
gam_norm = [c/cont for c in gam]
print("  gamma_norm =", [str(v) for v in gam_norm])
pl2, cont2 = place_set(g1, gam_norm, a_val)
S2 = set(pl2)
for c in gam_norm:
    if c != 0: S2 |= primes_of(c)
for c in g1:
    if c != 0: S2 |= primes_of(c)
S2 = sorted(S2)
print("  набор:", S2)
tot2, _ = run_pairing(g1, gam_norm, a_val, S2, "<g1,g2> нормированная gamma")

print()
print("=== 11. Все 8 выборов знаков m ===")
for sg in [(1,1,1),(1,1,-1),(1,-1,1),(-1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1),(-1,-1,-1)]:
    gg, _ = gamma_form(g1, g2, g3, signs=sg)
    plS, _c = place_set(g1, gg, a_val)
    SS = set(plS)
    for c in gg:
        if c != 0: SS |= primes_of(c)
    for c in g1:
        if c != 0: SS |= primes_of(c)
    t, _ = run_pairing(g1, gg, a_val, sorted(SS), "signs=%s" % (sg,), verbose=False)
    print("  знаки %s -> %+d" % (sg, t))

print()
print("=== 12. Обратный порядок <g2,g1> ===")
gamR, mR = gamma_form(g2, g1, g3)
aR = F(g1[0])
print("  gamma(по g2) =", [str(v) for v in gamR], "  a = g1(1,0) =", aR)
plR, contR = place_set(g2, gamR, aR)
SR = set(plR)
for c in gamR:
    if c != 0: SR |= primes_of(c)
for c in g2:
    if c != 0: SR |= primes_of(c)
totR, _ = run_pairing(g2, gamR, aR, sorted(SR), "<g2,g1>")

print()
print("=== 13. Fisher 3.2(v): <g1,g3> должно равняться <g1,g2> ===")
gam13, _ = gamma_form(g1, g3, g2)
a3 = F(g3[0])
pl3, _c = place_set(g1, gam13, a3)
S3 = set(pl3)
for c in gam13:
    if c != 0: S3 |= primes_of(c)
for c in g1:
    if c != 0: S3 |= primes_of(c)
tot13, _ = run_pairing(g1, gam13, a3, sorted(S3), "<g1,g3>")

print()
print("=== ИТОГ ===")
print("  <g1,g2> сырая/сверх-набор  =", tot1)
print("  <g1,g2> сырая/мин-набор    =", tot1b)
print("  <g1,g2> другие точки       =", tot1c)
print("  <g1,g2> нормированная gamma=", tot2)
print("  <g2,g1>                    =", totR)
print("  <g1,g3>                    =", tot13)
