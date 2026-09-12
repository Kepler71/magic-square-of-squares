# -*- coding: utf-8 -*-
# Раунд 3 атаки на (15,1).
#  N1. Исправленная проверка констант торсора D + явная 2-адическая точка.
#  N2. Граница ранга БЕЗ второго спуска: первый спуск на изогенных кривых.
#  N3. Локальная разрешимость C: прямой счёт по t в F_p для всех p < 5000.
#  N4. Плохие простые покрытия C (где сливаются точки ветвления).
#  N5. Расширенный поиск точек на C другой параметризацией.

import sys
from sage.all import *

def hdr(s):
    print("\n" + "=" * 72); print(s); print("=" * 72); sys.stdout.flush()

S = QQ(113); M = ZZ(15); N = ZZ(1); B = S*M**2*N**2
e = [QQ(-B), QQ(-S*M**4), QQ(-S*N**4)]
Rx = PolynomialRing(QQ,'x'); x = Rx.gen()
f = (x-e[0])*(x-e[1])*(x-e[2]); c = f.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]]); Emin = E.minimal_model()

hdr("N1. ТОРСОР D: константы и 2-адическая точка")
print("  X-e1 = r^2, X-e2 = 113 z2^2, X-e3 = 113 z3^2")
print("  => r^2 - 113 z2^2 = e2-e1 = %s" % (e[1]-e[0]))
print("  => r^2 - 113 z3^2 = e3-e1 = %s" % (e[2]-e[0]))
ok = (e[1]-e[0] == -5695200) and (e[2]-e[0] == 25312)
print("  совпадает с торсором Codex (-5695200, 25312): %s" % ok)

# точка над Q_2 из t=0 на C: F0=225, F4=113, F8=1
print("\n  t=0 на C над Q_2:  F0=225=15^2, F4=113, F8=1")
print("  113 mod 8 = %s  => 113 квадрат в Q_2" % (113 % 8))
print("  X = 0;  X-e1 = 25425, 25425 mod 8 = %s => квадрат в Q_2" % (25425 % 8))
print("  (X-e2)/113 = %s = %s^2" % ((0-e[1])/113, sqrt((0-e[1])/113)))
print("  (X-e3)/113 = %s" % ((0-e[2])/113))
print("  => D(Q_2) != пусто; мой целочисленный перебор в раунде 2 просто не")
print("     видел r = sqrt(25425) (иррационален над Q, но квадрат в Q_2).")

hdr("N2. ГРАНИЦА РАНГА БЕЗ ВТОРОГО СПУСКА: первый спуск на изогенных кривых")
cl = Emin.isogeny_class()
print("  изогенный класс: %d кривых, матрица степеней:" % len(cl.curves))
print(cl.matrix())
best = None
for i, Ei in enumerate(cl.curves):
    try:
        mwi = Ei.mwrank_curve()
        mwi.two_descent(second_descent=False, verbose=False)
        rb = mwi.rank_bound(); cert = mwi.certain(); sr = mwi.selmer_rank()
        print("  кривая %d %s : первый спуск rank_bound=%s certain=%s selmer_rank=%s"
              % (i, Ei.ainvs(), rb, cert, sr))
        if best is None or rb < best:
            best = rb
    except Exception as ex:
        print("  кривая %d: ошибка %s" % (i, ex))
print("\n  ЛУЧШАЯ граница ранга от ОДНОГО ТОЛЬКО первого 2-спуска по классу: %s" % best)
if best is not None and best <= 1:
    print("  => граница rank<=1 получается БЕЗ второго спуска и без Касселса-Тейта.")
else:
    print("  => первого спуска нигде не хватает; граница 1 держится ТОЛЬКО на")
    print("     втором спуске (eclib) / спаривании Касселса-Тейта (PARI).")

hdr("N3. ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ C: прямой счёт по t в F_p, p < 5000")
badp = []
for p in prime_range(3, 5000):
    p = int(p)
    found = False
    for tv in range(p):
        t2 = (tv*tv) % p
        a0 = (225 + t2) % p
        if a0 == 0 or kronecker(a0, p) != 1: continue
        a4 = (113*(1 + t2)) % p
        if a4 == 0 or kronecker(a4, p) != 1: continue
        a8 = (1 + 225*t2) % p
        if a8 == 0 or kronecker(a8, p) != 1: continue
        found = True
        break
    if not found:
        badp.append(p)
print("  нечётные p<5000 без гладкой F_p-точки со всеми тремя квадратами: %s" % badp)
print("  (для таких p по Гензелю точка над F_p поднимается в Q_p при p нечётном)")
# p = 2 отдельно
print("  p=2: t=0 даёт F0=225 (225 mod 8 = %s), F4=113 (mod 8 = %s), F8=1 -> все квадраты в Q_2"
      % (225 % 8, 113 % 8))
print("  над R: t=0, все F_i > 0")

hdr("N4. ПЛОХИЕ ПРОСТЫЕ НАКРЫТИЯ C")
Rt = PolynomialRing(QQ,'t'); t = Rt.gen()
f0 = 225 + t**2; f4 = 113*(1+t**2); f8 = 1 + 225*t**2
prod = f0*f4*f8
print("  disc(F0) = %s" % factor(f0.discriminant()))
print("  disc(F4) = %s" % factor(f4.discriminant()))
print("  disc(F8) = %s" % factor(f8.discriminant()))
print("  res(F0,F4) = %s" % factor(f0.resultant(f4)))
print("  res(F0,F8) = %s" % factor(f0.resultant(f8)))
print("  res(F4,F8) = %s" % factor(f4.resultant(f8)))
print("  ведущие коэффициенты: F0 -> 1, F4 -> 113, F8 -> 225 = 3^2*5^2")
print("  => плохие простые накрытия содержатся в {2,3,5,7,113}")

hdr("N5. РАСШИРЕННЫЙ ПОИСК ТОЧЕК НА C_{15,1}")
# t = p/q; нужно 113(p^2+q^2)=квадрат, т.е. p^2+q^2 = 113*k^2
from math import isqrt
def issq(v):
    if v < 0: return False
    r = isqrt(v); return r*r == v
# p^2 + q^2 = 113 k^2  =>  p^2 = -q^2 mod 113  =>  p = ±r0*q mod 113
r0 = int(Mod(-1, 113).sqrt())
NB = 20000
cands = 0
hits = []
for q in range(1, NB+1):
    q2 = q*q
    for sgn in (1, -1):
        p0 = (sgn*r0*q) % 113
        p = p0
        while p <= NB:
            if p > 0 or q > 0:
                v = p*p + q2
                if v % 113 == 0 and issq(v // 113) and gcd(p, q) == 1:
                    cands += 1
                    if issq(225*q2 + p*p) and issq(q2 + 225*p*p):
                        hits.append((p, q))
            p += 113
    if q % 5000 == 0:
        print("    ... q = %d, кандидатов %d" % (q, cands)); sys.stdout.flush()
print("  t=p/q, 0<=p<=%d, 1<=q<=%d, gcd=1: кандидатов по F4 = %d, ПОЛНЫХ ТОЧЕК = %d %s"
      % (NB, NB, cands, len(hits), hits[:5]))

print("\nГОТОВО")
