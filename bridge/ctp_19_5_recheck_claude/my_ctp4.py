#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""ЧАСТЬ 4: (a) сверка данных покрытия из сертификата (T_i, X_num, G_raw, lambda);
(b) тесты на подмену — проверяю, что МОЙ конвейер не вырожден."""
from fractions import Fraction as F
import math, random, json
random.seed(99)
_lib = open('/home/kep/magicKube/bridge/ctp_19_5_recheck_claude/my_ctp_lib.py').read()
exec(_lib)

cert = json.load(open('/home/kep/magicKube/bridge/ctp_cert_19_5.json'))
cv = cert['covering_g1']
T = [[F(t) for t in row] for row in cv['T']]
lam = F(cv['lambda'])
Dc = [F(t) for t in cv['D']]
Graw = [F(t) for t in cv['G_raw']]
Xnum = [F(t) for t in cv['X_num']]
d = [1, 193, 193]
g1 = [680904, 3128916, 1361808, -3128916, 680904]
g2 = [718704, -3036516, 1311408, 3288516, 592704]
g3 = [135072, 2334948, 8283996, 4669896, 540288]

print("=== 4a. Данные покрытия сертификата — проверяю КАК ТОЖДЕСТВА МНОГОЧЛЕНОВ ===")
def qsq(t):  # (A x^2+B xz+C z^2)^2 -> 5 коэфф.
    A,B,C = t
    return [A*A, 2*A*B, B*B+2*A*C, 2*B*C, C*C]
print("  D_i = prod_{j!=i}(r_i-r_j) (мои корни r):",
      [(r[i]-r[(i+1)%3])*(r[i]-r[(i+2)%3]) for i in range(3)])
print("  D_i из сертификата:                       ", [str(v) for v in Dc])
myD = [F((r[i]-r[(i+1)%3])*(r[i]-r[(i+2)%3])) for i in range(3)]
print("  совпало:", myD == Dc)
print("  lambda^2 * g1 == G_raw :", [lam*lam*F(c) for c in g1] == Graw)
sig = [d[i]*F(v)/Dc[i] for i, v in enumerate([1,1,1])]
acc = [F(0)]*5
for i in range(3):
    s = qsq(T[i])
    for k in range(5): acc[k] += d[i]*s[k]/Dc[i]
print("  sum_i d_i T_i^2 / D_i == 0 :", all(v == 0 for v in acc))
ok_all = True
for i in range(3):
    lhs = [Xnum[k] - F(r[i])*Graw[k] for k in range(5)]
    rhs = [d[i]*v for v in qsq(T[i])]
    print("  X_num - r_%d*G_raw == d_%d*T_%d^2 : %s" % (i+1, i+1, i+1, lhs == rhs))
    ok_all &= (lhs == rhs)
print("  ВСЕ тождества покрытия:", ok_all)
print("  => любая точка y^2=g1(x,z), y!=0, даёт x_M = X_num/G_raw на M с классами",
      d, "(независимо подтверждает мой вывод через сизигию)")

print()
print("=== 4b. Тесты на подмену в МОЁМ конвейере ===")
def gamma_and_pair(gA, gB, gC, aval=None, drop=None, gam_override=None):
    zA = z_inv(gA); zB = z_inv(gB); zC = z_inv(gC)
    mm = [sqrtF(zA[i]*zB[i]*zC[i]) for i in range(3)]
    gam = [F(0)]*3
    for i in range(3):
        H = H_form_g(gA, i); co = mm[i]/zA[i]
        for k in range(3): gam[k] += co*H[k]/D[i]
    if gam_override is not None: gam = gam_override
    a = F(gB[0]) if aval is None else F(aval)
    S = set([2,3,5,7]) | primes_of(Delta) | primes_of(a) | primes_of(content_of(gam))
    for c in gam:
        if c != 0: S |= primes_of(c)
    S = sorted(S)
    if drop is not None: S = [p for p in S if p != drop]
    tot = 1; mn = []
    for p in S + ['inf']:
        res = find_local_point(gA, gam, p)
        if res is None: return None, None
        x, z, gv, gm = res
        hh = hilbert(a, gm, p); tot *= hh
        if hh == -1: mn.append(str(p))
    return tot, mn

base, bm = gamma_and_pair(g1, g2, g3)
print("  эталон           : %+d  минусы %s" % (base, bm))
v, mn = gamma_and_pair(g1, g2, g3, drop=23)
print("  выброшено место 23: %+d  минусы %s   -> ЛОВИТСЯ (знак сменился): %s" % (v, mn, v != base))
v, mn = gamma_and_pair(g1, g2, g3, drop=3)
print("  выброшено место 3 : %+d  минусы %s   -> ЛОВИТСЯ: %s" % (v, mn, v != base))
bad = [F(132076,3)+1, F(529592,3), F(-775796,3)]
v, mn = gamma_and_pair(g1, g2, g3, gam_override=bad)
print("  испорчен gamma[0] : %+d  минусы %s   -> знак изменился: %s" % (v, mn, v != base))
gbad = list(g1); gbad[2] += 1
print("  испорчен g1[2]    : инварианты I,J совпадают?", inv_I(gbad) == I and inv_J(gbad) == J, " -> ЛОВИТСЯ проверкой инвариантов")
# подмена целевого класса: беру квартику класса (193,193,1) вместо (1,193,193)
print("  подмена g1 -> g2 (класс (193,193,1)): классы z =", [sqclass(x) for x in z_inv(g2)],
      " != требуемого", [1,193,193])
v, mn = gamma_and_pair(g2, g3, g1)
print("    <g2,g3> = %+d  минусы %s  (тоже -1: спаривание невырождено на клейновой четвёрке)" % (v, mn))

print()
print("=== 4c. Нет ли минуса в вещественном месте? ===")
print("  a = g2(1,0) =", g2[0], "> 0  =>  (a, *)_real = +1 тождественно:",
      all(hilbert(F(g2[0]), F(t), 'inf') == 1 for t in (-7, -1, 1, 13, -1000000)))
print("  C_(19,5) имеет вещественные точки: F0=361+25t^2>0, F4=193(1+t^2)>0, F8=25+361t^2>0 при всех t")
print("  => препятствие чисто глобальное, вещественное место ни при чём")

print()
print("=== 4d. Контроль содержания квартик (обоснование пропуска мест) ===")
for nm, g in (("g1",g1),("g2",g2),("g3",g3)):
    c = 0
    for t in g: c = math.gcd(c, abs(t))
    print("  content(%s) = %d, простые %s  (подмножество supp(Delta) = %s): %s"
          % (nm, c, sorted(factor(c)), sorted(primes_of(Delta)),
             set(factor(c)) <= primes_of(Delta)))
