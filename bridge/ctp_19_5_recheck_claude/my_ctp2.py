#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ЧАСТЬ 2: атака на слабые места сертификата.
 A. СВОЙ вывод связи «квартика -> класс спуска» через сизигию бинарной квартики
    (без опоры на T_i / X_num сертификата).
 B. Контроли: диагональ, тривиальный класс, образ 2-кручения.
 C. Независимость от представителя: SL2(Z)-замены g1 и g2 -> спаривание обязано не меняться.
"""
from fractions import Fraction as F
import math, random
random.seed(7)
exec(open('/home/kep/magicKube/bridge/ctp_19_5_recheck_claude/my_ctp_lib.py').read())

print("=== A. СВОЙ вывод: класс спуска квартики = класс z(g) ===")
# Работаем с многочленами от (x,z) как со списками коэффициентов (по убыванию степени x).
def pmul(P, Q):
    R = [0]*(len(P)+len(Q)-1)
    for i, p in enumerate(P):
        for j, q in enumerate(Q):
            R[i+j] += p*q
    return R
def padd(P, Q):
    n = max(len(P), len(Q)); R = [0]*n
    for i, p in enumerate(P): R[n-len(P)+i] += p
    for i, q in enumerate(Q): R[n-len(Q)+i] += q
    return R
def pscal(k, P): return [k*c for c in P]
def ptrim(P):
    i = 0
    while i < len(P)-1 and P[i] == 0: i += 1
    return P[i:]

g1 = [680904, 3128916, 1361808, -3128916, 680904]
g2 = [718704, -3036516, 1311408, 3288516, 592704]
g3 = [135072, 2334948, 8283996, 4669896, 540288]

h1 = hessian(g1)
# сизигия: h^3 - 48 I g^2 h - 64 J g^3 = ? * t^2
lhs = padd(padd(pmul(pmul(h1, h1), h1), pscal(-48*I, pmul(pmul(g1, g1), h1))),
           pscal(-64*J, pmul(pmul(g1, g1), g1)))
lhs = ptrim(lhs)
print("  deg(h^3 - 48I g^2 h - 64J g^3) =", len(lhs)-1, "(ожидается 12)")
# ищем c и многочлен t степени 6 с lhs = c * t^2
def poly_sqrt(P):
    """точный квадратный корень многочлена (если есть), P степени 2k"""
    n = len(P)-1
    assert n % 2 == 0
    k = n//2
    lead = P[0]
    r0 = sqrtint(lead) if isinstance(lead, int) else None
    if r0 is None: return None
    T = [F(0)]*(k+1); T[0] = F(r0)
    for i in range(1, k+1):
        ssum = sum(T[j]*T[i-j] for j in range(1, i))
        T[i] = (F(P[i]) - ssum)/(2*T[0])
    # проверка
    return T if pmul(T, T) == [F(c) for c in P] else None

found = None
for c in (27, -27, 1, -1, 729, -729, 108, -108):
    P = [F(v, c) for v in lhs]
    # привести к целому лидеру
    lead = P[0]
    if lead.denominator != 1: continue
    T = poly_sqrt([int(v) if v.denominator == 1 else v for v in P]) if all(v.denominator == 1 for v in P) else None
    if T is not None:
        found = (c, T); break
print("  сизигия: h^3 - 48I g^2 h - 64J g^3 = %s * t(x,z)^2 ; найдено c = %s" % ("c", found[0] if found else None))
assert found is not None
csyz, tcov = found
print("  t(x,z) степени", len(tcov)-1)

# Теперь ПРОВЕРЯЮ КАК ТОЖДЕСТВА МНОГОЧЛЕНОВ:
#   X(x,z) := 3 h(x,z) / (4 g(x,z))
#   X + 3 phi_i = 9 G_i /(4 g)  = 9 H_i^2 / (4 z_i g)
#   X^3 - 27 I X - 27 J = 729 * c/27 * t^2 / (64 g^3)
# 1) X + 3phi = 9G/(4g) <=> 3h + 12 phi g = 9G = 3(4 phi g + h): тождество очевидно, проверим:
for i in range(3):
    G = [F(4*phi[i]*g1[k] + hessian(g1)[k], 3) for k in range(5)]
    lhs2 = padd(pscal(F(3), h1), pscal(F(12*phi[i]), g1))     # 3h + 12 phi g
    rhs2 = pscal(F(9), G)
    assert [F(v) for v in lhs2] == [F(v) for v in rhs2], "X+3phi identity"
print("  тождество  3h + 12*phi_i*g = 9*G_i  проверено для i=1,2,3")
# 2) G_i = H_i^2 / z_i  -- это уже проверенное G(1,0)G=H^2
for i in range(3):
    G = [F(4*phi[i]*g1[k] + hessian(g1)[k], 3) for k in range(5)]
    H = H_form_g(g1, i)
    assert pmul([F(G[0])], [F(v) for v in G]) == pmul([F(v) for v in H], [F(v) for v in H])
print("  тождество  z_i * G_i = H_i^2  проверено")
# 3) сизигия -> точка лежит на E_{I,J}
print("  сизигия проверена как тождество многочленов степени 12")

print()
print("  СЛЕДСТВИЕ (мой вывод, не из сертификата):")
print("   если g1(x,z) = y^2 != 0, то X = 3h/(4g) даёт точку E_{I,J}: Y^2 = X^3-27IX-27J,")
print("   и X - (-3 phi_i) = 9 H_i(x,z)^2 / (4 z_i y^2)  =>  класс = класс z(g1)_i.")
z1 = z_inv(g1)
print("   класс z(g1) =", [sqclass(v) for v in z1])
print("   2-кручение E_{I,J} имеет x = -3 phi_i =", [-3*p for p in phi])
print("   E_{I,J} <-> M через X = 36 x (36 - квадрат): -3phi_i = 36 r_i, r =", r)
for i in range(3):
    assert -3*phi[i] == 36*r[i]
print("   => класс (x_M - r_i) точки M(Q), приходящей из y^2=g1, равен", [sqclass(v) for v in z1])
print("   ТРЕБУЕМЫЙ класс delta =", [1,193,193], " -> СОВПАДАЕТ:", [sqclass(v) for v in z1] == [1,193,193])

print()
print("=== A2. Обратная сторона: g1 действительно ИЗОМОРФНО C_delta (степень покрытия) ===")
# Проверяю, что отображение (x:z) -> X = 3h/(4g) имеет степень 4 (2-накрытие), т.е. h и g не пропорциональны
# и что ker: h*g' - h'*g не тождественный ноль.
ratio_ok = not all(h1[k]*g1[0] == g1[k]*h1[0] for k in range(5))
print("  h и g не пропорциональны:", ratio_ok, " => X непостоянно, deg = 4  => это 2-накрытие")
assert ratio_ok

print()
print("=== B1. Образ 2-кручения M(Q)[2] при спуске (не должен содержать delta) ===")
d12 = r[0]-r[1]; d13 = r[0]-r[2]; d23 = r[1]-r[2]
tors = {
    "O":  (1,1,1),
    "P1": (sqclass(d12*d13), sqclass(d12), sqclass(d13)),
    "P2": (sqclass(-d12), sqclass(d12*d23*(-1)*(-1)), sqclass(d23)),
    "P3": (sqclass(-d13), sqclass(-d23), sqclass(d13*d23)),
}
# аккуратно: delta(P_i) = (r_i-r_j для j!=i, а i-я компонента = произведение остальных)
def dmap(i):
    v = [None]*3
    for j in range(3):
        if j != i: v[j] = sqclass(r[i]-r[j])
    v[i] = sqclass((r[i]-r[(i+1)%3])*(r[i]-r[(i+2)%3]))
    return tuple(v)
tors = {"O": (1,1,1), "P1": dmap(0), "P2": dmap(1), "P3": dmap(2)}
for k, v in tors.items(): print("   delta(%s) = %s" % (k, v))
print("  delta=(1,193,193) среди образов 2-кручения:", (1,193,193) in tors.values())
assert (1,193,193) not in tors.values()

print()
print("=== B2. Тривиальный класс: контроль на ЛОЖНОЕ СРАБАТЫВАНИЕ ===")
gtr = [1361808, 1540140, -2723616, -1540140, 1361808]
print("  g_triv инварианты совпали:", inv_I(gtr) == I, inv_J(gtr) == J)
ztr = z_inv(gtr)
print("  классы z(g_triv) =", [sqclass(v) for v in ztr], "(должно быть (1,1,1))")
# ищу рациональную точку на y^2 = g_triv  -> доказывает [g_triv]=0
pt = None
for de in range(1, 60):
    for nu in range(-400, 401):
        if math.gcd(abs(nu), de) != 1: continue
        v = evalq(gtr, F(nu), F(de))
        if v > 0:
            s_ = None
            if v.denominator == 1:
                s_ = sqrtint(v.numerator)
            else:
                a_ = sqrtint(v.numerator); b_ = sqrtint(v.denominator)
                s_ = a_ if (a_ is not None and b_ is not None) else None
            if s_ is not None:
                pt = (nu, de, v); break
    if pt: break
print("  рациональная точка на y^2=g_triv:", pt)
# спаривание <g1, g_triv>: нужна тройка с z(g1) z(g_triv) z(g?) = квадрат -> g? в классе g1
# берём (g1, g_triv, g1') где g1' - собственно эквивалентная g1 (тот же класс)
def act(g, M):
    """g -> g(a x + b z, c x + d z), M = (a,b,c,d)"""
    a_, b_, c_, d_ = M
    X = [a_, c_]   # a x + c z  (коэфф. при x, z)
    Z = [b_, d_]
    # строим как многочлен от (x,z)
    def lin(P): return P
    res = [0]*5
    # g = sum g_k * X^(4-k) * Z^k
    def powp(P, k):
        R = [1]
        for _ in range(k): R = pmul(R, P)
        return R
    acc = [0]
    for k in range(5):
        term = pscal(g[k], pmul(powp(X, 4-k), powp(Z, k)))
        acc = padd(acc, term)
    while len(acc) < 5: acc = [0]+acc
    return acc

g1p = act(g1, (1, 1, 0, 1))   # det=1, собственная эквивалентность (lambda=1)
print("  g1' = g1(x+z, z) =", g1p, " инварианты ок:", inv_I(g1p) == I, inv_J(g1p) == J)
print("  классы z(g1') =", [sqclass(v) for v in z_inv(g1p)])

def pairing(gA, gB, gC, verbose=False, extra=()):
    gamq, mm = gamma_form_g(gA, gB, gC)
    aval = F(gB[0])
    assert aval != 0
    S = place_set_g(gA, gamq, aval)
    SS = set(S) | set(extra)
    for c in gamq:
        if c != 0: SS |= primes_of(c)
    for c in gA:
        if c != 0: SS |= primes_of(F(c))
    tot = 1; minus = []
    for p in sorted(SS) + ['inf']:
        res = find_local_point(gA, gamq, p)
        if res is None: return None, None
        x, z, gv, gm = res
        hh = hilbert(aval, gm, p)
        tot *= hh
        if hh == -1: minus.append(str(p))
        if verbose: print("     v=%-6s (%s:%s) -> %+d" % (p, x, z, hh))
    return tot, minus

v, mn = pairing(g1, gtr, g1p)
print("  <g1, g_triv> =", v, " минусы:", mn, "   (ОБЯЗАНО быть +1)")

print()
print("=== B3. Диагональ <g1,g1> и <g2,g2> (вычисляется, не постулируется) ===")
v11, m11 = pairing(g1, g1, gtr)
print("  <g1,g1> =", v11, " минусы:", m11, " (обязано +1)")
v22, m22 = pairing(g2, g2, gtr)
print("  <g2,g2> =", v22, " минусы:", m22, " (обязано +1)")
vtt, mtt = pairing(gtr, g1, g1p)
print("  <g_triv, g1> =", vtt, " минусы:", mtt, " (обязано +1)")

print()
print("=== C. Независимость от ПРЕДСТАВИТЕЛЯ класса (SL2(Z)-замены) ===")
mats = [(1,0,0,1), (1,1,0,1), (1,0,1,1), (2,1,1,1), (1,-3,0,1), (3,2,1,1), (0,-1,1,0), (1,0,5,1)]
for M1 in mats:
    gA = act(g1, M1)
    assert inv_I(gA) == I and inv_J(gA) == J
    assert [sqclass(v_) for v_ in z_inv(gA)] == [1,193,193], ("класс изменился!", M1, [sqclass(v_) for v_ in z_inv(gA)])
    v_, mn_ = pairing(gA, g2, g3)
    print("  g1 o %-12s : <g1,g2> = %+d   минусы %s" % (str(M1), v_, mn_))
print()
for M2 in mats:
    gB = act(g2, M2)
    assert inv_I(gB) == I and inv_J(gB) == J
    if gB[0] == 0:
        print("  g2 o %-12s : g2(1,0)=0, пропуск (Fisher 3.2(iii))" % str(M2)); continue
    v_, mn_ = pairing(g1, gB, g3)
    print("  g2 o %-12s : <g1,g2> = %+d   минусы %s" % (str(M2), v_, mn_))
print()
for M3 in mats[:5]:
    gC = act(g3, M3)
    v_, mn_ = pairing(g1, g2, gC)
    print("  g3 o %-12s : <g1,g2> = %+d   минусы %s" % (str(M3), v_, mn_))
