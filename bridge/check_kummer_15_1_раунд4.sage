# -*- coding: utf-8 -*-
# Раунд 4: закрытие двух оставшихся дыр.
#  P1. p=17 и p=113 — те самые p, где мой строгий критерий (все три F_i —
#      НЕНУЛЕВЫЕ квадраты mod p) не срабатывает. Проверяем чётность
#      оценки + квадратность единицы, т.е. настоящий критерий Q_p.
#  P2. Замкнутость образа delta: перебор aG+cT1+dT2, |a| <= 12.

import sys
from sage.all import *

def hdr(s):
    print("\n" + "="*72); print(s); print("="*72); sys.stdout.flush()

Rt = PolynomialRing(QQ,'t'); t = Rt.gen()
f0 = 225 + t**2; f4 = 113*(1+t**2); f8 = 1 + 225*t**2

def sq_in_Qp(xv, p):
    xv = QQ(xv)
    if xv == 0: return None
    v = xv.valuation(p)
    if v % 2: return False
    u = xv/QQ(p)**v
    num = ZZ(u.numerator()); den = ZZ(u.denominator())
    if p == 2:
        return (num*inverse_mod(den,8)) % 8 == 1
    r = (num*inverse_mod(den,p)) % p
    return r != 0 and kronecker(r,p) == 1

hdr("P1. КРИТЕРИЙ Q_p ДЛЯ p = 17 И p = 113 (и контрольно p=2)")
for p, tv in [(17, 8), (113, 15), (2, 0), (3, 1), (7, 0)]:
    vals = [f0(tv), f4(tv), f8(tv)]
    res = [sq_in_Qp(v, p) for v in vals]
    vlt = [QQ(v).valuation(p) for v in vals]
    print("  p=%-4s t=%-3s : F=(%s, %s, %s)" % (p, tv, vals[0], vals[1], vals[2]))
    print("           оценки v_p = %s ; квадраты в Q_p = %s" % (vlt, res))
    print("           ВСЕ ТРИ квадраты в Q_p: %s" % all(res))

hdr("P1b. ЕСТЬ ЛИ ВООБЩЕ t mod 17 со всеми тремя ненулевыми квадратами?")
cnt = 0
for tv in range(17):
    a0 = (225+tv*tv) % 17; a4 = (113*(1+tv*tv)) % 17; a8 = (1+225*tv*tv) % 17
    if a0 and a4 and a8 and kronecker(a0,17)==1 and kronecker(a4,17)==1 and kronecker(a8,17)==1:
        cnt += 1
print("  таких t mod 17: %d  => мой строгий скан 'все три ненулевые квадраты'" % cnt)
print("  действительно не может найти точку при p=17; но точка Q_17 существует")
print("  через t=8, где F0 = 289 = 17^2 (чётная оценка, единичная часть 1).")
print("  Соответствующая F_17-точка кривой: (t,u0,u4,u8)=(8,0,u4,u8) — она ГЛАДКАЯ,")
print("  т.к. F0'(8) = 16 != 0 mod 17, поэтому поднимается по Гензелю.")
u4sq = (113*(1+64)) % 17; u8sq = (1+225*64) % 17
print("  u4^2 = %s (квадрат: %s), u8^2 = %s (квадрат: %s)"
      % (u4sq, kronecker(u4sq,17)==1, u8sq, kronecker(u8sq,17)==1))

hdr("P1c. p=113: F4 = 113(1+t^2) ВСЕГДА делится на 113")
print("  значит нужно 113 | 1+t^2, т.е. t = ±r mod 113, r^2 = -1")
r0 = ZZ(Mod(-1,113).sqrt())
print("  r = %s ; 15^2+1 = 226 = 2*113 => t=15 годится (15 = ±r mod 113? %s)"
      % (r0, (15 % 113) in [r0 % 113, (-r0) % 113]))
print("  2 — квадрат mod 113 (113 = 1 mod 8): %s" % (kronecker(2,113)==1))

hdr("P2. ЗАМКНУТОСТЬ ОБРАЗА delta (перебор |a| <= 12)")
S = QQ(113); M = ZZ(15); N = ZZ(1); B = S*M**2*N**2
e = [QQ(-B), QQ(-S*M**4), QQ(-S*N**4)]
Rx = PolynomialRing(QQ,'x'); x = Rx.gen()
f = (x-e[0])*(x-e[1])*(x-e[2]); c = f.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]])
def sqclass(v):
    v = QQ(v); z = ZZ(v.numerator()*v.denominator())
    return ZZ(z.sign())*ZZ(z.abs().squarefree_part())
def delta(P):
    if P.is_zero(): return (ZZ(1),ZZ(1),ZZ(1))
    d = [QQ(P[0]-e[i]) for i in range(3)]
    for i in range(3):
        if d[i]==0:
            j,k = [q for q in range(3) if q != i]
            d[i] = QQ((e[i]-e[j])*(e[i]-e[k]))
    return tuple(sqclass(v) for v in d)
G = E.gens()[0]; T1 = E(e[0],0); T2 = E(e[1],0)
REQ = (ZZ(1), ZZ(113), ZZ(113))
seen = set(); bad = []
for a in range(-12, 13):
    aG = a*G
    for cc in [0,1]:
        for dd in [0,1]:
            P = aG + cc*T1 + dd*T2
            if P.is_zero(): continue
            v = delta(P); seen.add(v)
            if v == REQ: bad.append((a,cc,dd))
print("  различных классов среди 100 точек: %d (ожидание 8)" % len(seen))
print("  точек с требуемым классом (1,113,113): %d %s" % (len(bad), bad))
print("  классы: %s" % sorted(seen))
print("\nГОТОВО")
