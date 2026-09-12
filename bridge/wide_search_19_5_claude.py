#!/usr/bin/env python3
# Широкий независимый поиск рациональных точек на C_{19,5}.
# t = p/q, gcd(p,q)=1, требуется одновременно:
#   A = 361 q^2 +  25 p^2  — квадрат
#   B =  25 q^2 + 361 p^2  — квадрат
#   D = 193 (q^2 + p^2)    — квадрат
# Сильнейшее условие — третье: p^2+q^2 = 193 k^2. Используем его как основной фильтр:
# p^2 = -q^2 mod 193  =>  p = +-c*q mod 193, где c^2 = -1 mod 193.
from math import isqrt, gcd
import sys

P = 193
# c с c^2 = -1 mod 193
c = None
for x in range(1, P):
    if (x * x + 1) % P == 0:
        c = x
        break
assert c is not None
print("c =", c, " c^2 mod 193 =", (c * c) % P)

def issq(n):
    if n < 0:
        return False
    r = isqrt(n)
    return r * r == n

# быстрые фильтры по квадратичным вычетам
FMODS = [64, 63, 65, 11, 13, 17, 19, 23, 29, 31]
QRS = []
for M in FMODS:
    QRS.append(bytearray(M))
    for i in range(M):
        QRS[-1][(i * i) % M] = 1

def maybe_sq(n):
    for M, q in zip(FMODS, QRS):
        if not q[n % M]:
            return False
    return True

QMAX = 30000
hits = []
cands = 0
for q in range(1, QMAX + 1):
    q2 = q * q
    r1 = (c * q) % P
    r2 = (-c * q) % P
    starts = {r1, r2}
    for st in starts:
        p = st
        if p == 0:
            p = P
        while p <= QMAX:
            if gcd(p, q) == 1:
                p2 = p * p
                D = P * (q2 + p2)
                if maybe_sq(D) and issq(D):
                    cands += 1
                    A = 361 * q2 + 25 * p2
                    if maybe_sq(A) and issq(A):
                        B = 25 * q2 + 361 * p2
                        if maybe_sq(B) and issq(B):
                            hits.append((p, q))
                            print("HIT t = %d/%d" % (p, q))
            p += P
    if q % 5000 == 0:
        print("  q =", q, " кандидатов с квадратным F4:", cands, flush=True)

print("QMAX =", QMAX)
print("кандидатов, у которых F4 — квадрат:", cands)
print("полных попаданий (все три квадрата):", hits)
