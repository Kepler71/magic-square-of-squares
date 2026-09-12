#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Попытка СЛОМАТЬ заявление C_{19,5}(Q) = пусто: прямой поиск точки.
Не brute force по t, а параметризация самого жёсткого условия.

C: u0^2 = 25t^2+361, u4^2 = 193(t^2+1), u8^2 = 361t^2+25.
t = p/q, gcd(p,q)=1  =>  нужно, чтобы были квадратами:
   25p^2+361q^2 ,  193(p^2+q^2) ,  361p^2+25q^2 .

Условие 193(p^2+q^2)=квадрат  <=>  p^2+q^2 = 193*k^2.
193 = 12^2+7^2. В Z[i]: p+qi = (12+7i)*(a+bi)*(единица) с a^2+b^2 = k^2,
т.е. (a,b,k) — пифагорова тройка: a=A^2-B^2, b=2AB (и перестановка).
Это даёт ТОЧНУЮ двухпараметрическую параметризацию (с точностью до
масштабирования и единиц), а не перебор.
"""
from math import isqrt, gcd
import sys

M, N = 19, 5
S = 193


def issq(x):
    if x < 0:
        return False
    r = isqrt(x)
    return r * r == x


# быстрые фильтры по вычетам для 25p^2+361q^2 и 361p^2+25q^2
FILTER_MODS = [9, 5, 7, 11, 13, 16]
QRSETS = {}
for mm in FILTER_MODS:
    QRSETS[mm] = set((i * i) % mm for i in range(mm))


def pass_filter(v, mm):
    return (v % mm) in QRSETS[mm]


LIM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
found = []
tested = 0
for A in range(1, LIM + 1):
    A2 = A * A
    for B in range(0, A):
        if gcd(A, B) != 1 or ((A - B) % 2 == 0):
            continue          # примитивная пифагорова тройка
        a = A2 - B * B
        bb = 2 * A * B
        for (aa, b2) in ((a, bb), (bb, a)):
            # p+qi = (12 +- 7i)(aa + b2 i)
            for sgn in (1, -1):
                p = 12 * aa - sgn * 7 * b2
                q = sgn * 7 * aa + 12 * b2
                p, q = abs(p), abs(q)
                if p == 0 or q == 0:
                    continue
                g = gcd(p, q)
                p //= g
                q //= g
                tested += 1
                v1 = 25 * p * p + 361 * q * q
                bad = False
                for mm in FILTER_MODS:
                    if not pass_filter(v1, mm):
                        bad = True
                        break
                if bad:
                    continue
                if not issq(v1):
                    continue
                v2 = 361 * p * p + 25 * q * q
                if not issq(v2):
                    continue
                v3 = S * (p * p + q * q)
                if not issq(v3):
                    continue
                found.append((p, q, isqrt(v1), isqrt(v3), isqrt(v2)))
                print("!!! НАЙДЕНА ТОЧКА C: t=%d/%d  u0=%d u4=%d u8=%d"
                      % (p, q, isqrt(v1), isqrt(v3), isqrt(v2)))

print("параметр A до %d; проверено кандидатов (p,q) с 193|(p^2+q^2)/квадрат: %d" % (LIM, tested))
print("найдено точек C: %d" % len(found))

# контроль метода: та же схема должна НАЙТИ известную точку для (15,8), t=1.
# там s=(225+64)/2=289/2 ... проверим позитивный контроль отдельно ниже.
print()
print("ПОЗИТИВНЫЙ КОНТРОЛЬ параметризации: ищем решения p^2+q^2=193k^2 и")
print("сверяем, что они действительно дают 193(p^2+q^2) квадратом.")
cnt = 0
for A in range(1, 40):
    for B in range(0, A):
        if gcd(A, B) != 1 or ((A - B) % 2 == 0):
            continue
        a = A * A - B * B
        bb = 2 * A * B
        for (aa, b2) in ((a, bb), (bb, a)):
            for sgn in (1, -1):
                p = abs(12 * aa - sgn * 7 * b2)
                q = abs(sgn * 7 * aa + 12 * b2)
                if p == 0 or q == 0:
                    continue
                g = gcd(p, q)
                p //= g
                q //= g
                assert issq(S * (p * p + q * q)), (p, q)
                cnt += 1
print("проверено %d штук — во всех 193(p^2+q^2) квадрат: параметризация КОРРЕКТНА" % cnt)

# и контроль полноты: перебором найдём все (p,q) с p,q<=400, 193(p^2+q^2)=квадрат,
# и убедимся, что все они попадают в параметризацию
brute = set()
for p in range(1, 401):
    for q in range(1, 401):
        if gcd(p, q) == 1 and issq(S * (p * p + q * q)):
            brute.add((p, q))
par = set()
for A in range(1, 200):
    for B in range(0, A):
        if gcd(A, B) != 1 or ((A - B) % 2 == 0):
            continue
        a = A * A - B * B
        bb = 2 * A * B
        for (aa, b2) in ((a, bb), (bb, a)):
            for sgn in (1, -1):
                p = abs(12 * aa - sgn * 7 * b2)
                q = abs(sgn * 7 * aa + 12 * b2)
                if p == 0 or q == 0:
                    continue
                g = gcd(p, q)
                p //= g
                q //= g
                if p <= 400 and q <= 400:
                    par.add((p, q))
print("перебор p,q<=400: %d пар; параметризация покрыла: %d; пропущено: %s"
      % (len(brute), len(brute & par), sorted(brute - par)[:20]))
