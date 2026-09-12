"""Проверка выживших пар: (1) CPU-повтор модульного фильтра (numpy, независимо от GPU); (2) точная арифметика:
z = x + y < 1 и 1 + z, 1 - z — рациональные квадраты (целочисленный isqrt). Для найденных точек — четвёртая пара."""
import sys, json
from math import gcd, isqrt
import numpy as np
from fractions import Fraction as F

B = int(sys.argv[1])
ab = np.load(f'gen_{B}_ab.npy'); A_s, B_s = ab[0].tolist(), ab[1].tolist()
n = len(A_s)
pairs = np.fromfile(f'pairs_{B}.bin', dtype=np.uint32).reshape(-1, 2)
res = np.fromfile(f'gen_{B}.bin', dtype=np.uint8).reshape(n, -1)
lines = open(f'gen_{B}_meta.txt').read().split('\n')
PR = list(map(int, lines[1].split()))
mw = [list(map(int, lines[2 + l].split())) for l in range(len(PR))]
mlo = np.array([row[0::2] for row in mw], dtype=np.uint64); mhi = np.array([row[1::2] for row in mw], dtype=np.uint64)
# (1) модульный повтор для всех выживших
I, J = pairs[:, 0].astype(np.int64), pairs[:, 1].astype(np.int64)
assert np.all(I < J)
ok = np.ones(len(pairs), dtype=bool)
for l in range(len(PR)):
    rj = res[J, l].astype(np.uint64)
    lo = mlo[l][res[I, l]]; hi = mhi[l][res[I, l]]
    w = np.where(rj < 64, lo, hi); sh = np.where(rj < 64, rj, rj - np.uint64(64))
    ok &= ((w >> sh) & np.uint64(1)) == 1
assert ok.all(), "GPU survivor fails CPU modular reference"
assert len({(int(i), int(j)) for i, j in pairs}) == len(pairs), "duplicate pairs"
def ND(i):
    a, b = A_s[i], B_s[i]; N = 4*a*b*(b*b - a*a); D = (a*a + b*b)**2; g = gcd(N, D); return N//g, D//g
def is_sq_frac(num, den):
    if num <= 0: return False
    t = num * den; r = isqrt(t); return r*r == t
valid = 0; one = 0; prod_sq = 0; points = []
for i, j in pairs.tolist():
    n1, d1 = ND(i); n2, d2 = ND(j)
    zn, zd = n1*d2 + n2*d1, d1*d2           # z = zn/zd
    if zn >= zd: continue
    valid += 1
    p_ = is_sq_frac(zd + zn, zd); m_ = is_sq_frac(zd - zn, zd)
    one += p_ or m_
    prod_sq += is_sq_frac(zd*zd - zn*zn, zd*zd)
    if p_ and m_:
        points.append({'t1': f'{A_s[i]}/{B_s[i]}', 't2': f'{A_s[j]}/{B_s[j]}', 'z': str(F(zn, zd))})
out = {'bound': B, 'values': n, 'survivors': len(pairs), 'cpu_modular_recheck': 'all pass',
       'positive_radical_survivors': valid, 'at_least_one_rational_radical': one,
       'product_square': prod_sq, 'six_line_points': points}
json.dump(out, open(f'verify_{B}.json', 'w'), indent=1)
print(json.dumps(out))
