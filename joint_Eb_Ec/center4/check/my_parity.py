# Проверяющий: утверждение 9. (i) 1±u — квадраты в Q_2 <=> v_2(u) >= 3 (перебор рациональных u);
# (ii) наборы (w_b, w_c, w_{b+c}, w_{b-c}) при B, C = 0 mod 8, |B|,|C| <= 480; w(E_N)=+1 <=> N mod 8 in {1,2,3}.
import json
from fractions import Fraction
from math import gcd
def v2(n):
    k = 0
    while n % 2 == 0: n //= 2; k += 1
    return k
def q2square(x):
    x = Fraction(x)
    if x == 0: return False
    p, q = x.numerator, x.denominator
    v = v2(abs(p)) - v2(q)
    if v % 2: return False
    up = p // 2**v2(abs(p)); uq = q // 2**v2(q)
    return (up * uq) % 8 == 1
bad = []
for p in range(-300, 301):
    for q in range(1, 301):
        if p == 0 or gcd(p, q) != 1: continue
        u = Fraction(p, q)
        vu = v2(abs(p)) - v2(q)
        lhs = q2square(1 + u) and q2square(1 - u)
        if lhs != (vu >= 3): bad.append(str(u))
def sqf(n):
    n = abs(n); r = 1; d = 2
    while d * d <= n:
        while n % (d * d) == 0: n //= d * d
        if n % d == 0: r *= d; n //= d
        d += 1
    return r * n
def w(t):
    return 1 if sqf(t) % 8 in (1, 2, 3) else -1
pats = {}
for B in range(-480, 481, 8):
    for C in range(-480, 481, 8):
        if B == 0 or C == 0 or B == C or B == -C: continue
        pat = (w(B), w(C), w(B + C), w(B - C))
        pats.setdefault(pat, (B, C))
out = dict(q2_mismatch=bad[:20], n_mismatch=len(bad), patterns=len(pats), examples={str(k): v for k, v in pats.items()})
print(out); json.dump(out, open('my_parity.json', 'w'), indent=1)
