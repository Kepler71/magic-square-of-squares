from fractions import Fraction as F
from math import isqrt, gcd

def v2(n):
    n = abs(n); e = 0
    while n % 2 == 0:
        n //= 2; e += 1
    return e

def ps(x):
    return x >= 0 and isqrt(x) ** 2 == x

def sqQ2(x):
    if x == 0:
        return True
    a, b = x.numerator, x.denominator
    neg = a < 0
    if neg:
        a = -a
    e = v2(a) - v2(b)
    if e % 2:
        return False
    ao = a >> v2(a); bo = b >> v2(b)
    r = ((-ao if neg else ao) * bo) % 8
    return r == 1

def cells(m, n, t):
    t = F(t); s = F(m * m + n * n, 2)
    return [m * m + n * n * t * t, (m * t + n) ** 2, s * (1 + t * t) - 2 * m * n * t,
            (m * t - n) ** 2, s * (1 + t * t), (m + n * t) ** 2,
            s * (1 + t * t) + 2 * m * n * t, (m - n * t) ** 2, n * n + m * m * t * t]

LINES = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]

print("=== 10. explicit RATIONAL points with BOTH red cells square, at even parameter ===")
for (m, n, t) in [(13, 8, F(27, 11)), (15, 8, F(289, 41)), (16, 5, F(-7077, 14773)),
                  (7, 2, F(95, 7)), (11, 2, F(31, 25))]:
    c = cells(m, n, t)
    q = t.denominator
    sc = [x * q * q for x in c]
    assert all(x.denominator == 1 for x in sc), (m, n, t)
    sc = [int(x) for x in sc]
    nsq = sum(1 for x in sc if ps(x))
    magic = all(sc[i] + sc[j] + sc[k] == 3 * sc[4] for i, j, k in LINES)
    print(f"({m},{n}) t={t}: c2={sc[2]} {'= '+str(isqrt(sc[2]))+'^2' if ps(sc[2]) else 'NOT a square'}; "
          f"c6={sc[6]} {'= '+str(isqrt(sc[6]))+'^2' if ps(sc[6]) else 'NOT a square'}; "
          f"squares among nine = {nsq}; distinct cells = {len(set(sc))}; magic = {magic}")

print()
print("=== 11. is the 2-adic witness forced to be the degenerate t=1? (NO) ===")
for (m, n) in [(13, 8), (15, 8), (19, 16), (16, 5)]:
    found = []
    for k in range(1, 2000):
        t = F(1 + 2 ** 7 * k, 1)
        c = cells(m, n, t)
        if all(sqQ2(c[i]) for i in (0, 2, 4, 6, 8)) and len(set(c)) == 9:
            found.append(t)
        if len(found) >= 3:
            break
    print(f"({m},{n}): non-degenerate t, all five cells squares in Q2: {found}")

print()
print("=== 12. Wesolowski chain mod 16: finite-state proof that s alternates 1,5 mod 8 ===")
seen = {}; st = (3, 1, 1); i = 0; seq = []
while st not in seen:
    seen[st] = i
    x, y, z = st
    seq.append((i, st, (x * x + y * y) % 16, (((x * x + y * y) % 32) // 2) % 8, z % 2))
    st = ((3 * x + 2 * z) % 16, x % 16, (x + z) % 16); i += 1
print(f"  state (x,y,z) mod 16: first repeat at i={i}, previous occurrence i={seen[st]}, period={i - seen[st]}")
print("  i | (x,y,z) mod 16 | x^2+y^2 mod 16 | s mod 8 | z mod 2")
for row in seq:
    print(f"  {row[0]:2d} | {row[1]} | {row[2]:3d} | {row[3]} | {row[4]}")
ok = all(r[3] == (5 if r[0] % 2 == 0 else 1) for r in seq) and all(r[4] == (1 if r[0] % 2 == 0 else 0) for r in seq)
print(f"  pattern 's=5 mod 8 and z odd exactly at even i' holds on the whole cycle: {ok}")
