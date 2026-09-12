from fractions import Fraction as F
from math import isqrt, gcd

def ps(x):
    return x >= 0 and isqrt(x) ** 2 == x

def red(m, n, u, v):
    # (4 v^2)-scaled red cells: H2 = 2[(mu-nv)^2+(mv-nu)^2], H6 = 2[(mu+nv)^2+(mv+nu)^2]
    return 2 * ((m * u - n * v) ** 2 + (m * v - n * u) ** 2), 2 * ((m * u + n * v) ** 2 + (m * v + n * u) ** 2)

print("search for t=u/v, |u|<=B, 1<=v<=B, with BOTH red cells rational squares")
B = 400
for (m, n) in [(19, 16), (13, 8), (15, 8), (16, 5), (7, 4), (3, 2)]:
    hits = []
    for v in range(1, B + 1):
        for u in range(-B, B + 1):
            if gcd(abs(u), v) != 1:
                continue
            a, b = red(m, n, u, v)
            if ps(a) and ps(b):
                hits.append(F(u, v))
    deg = [h for h in hits if h in (F(1), F(-1), F(0))]
    nondeg = [h for h in hits if h not in (F(1), F(-1), F(0))]
    print(f"({m},{n}): total {len(hits)} hits; degenerate t in {{0,1,-1}}: {sorted(set(deg))}; "
          f"non-degenerate: {sorted(set(nondeg))[:6]}{' ...' if len(set(nondeg))>6 else ''}")
