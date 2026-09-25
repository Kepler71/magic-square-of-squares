# Claude, 25.09: проверка идеи — при rank E_n(Q) = 1 три точки 2n_i G с x в арифм. прогрессии невозможны
# (примитивный делитель знаменателя старшего члена). Здесь: (1) у D_{2m} есть примитивный делитель,
# m = 1..M; (2) прямой поиск прогрессий среди x(2mG), |m| <= M.
from sage.all import *
import itertools, json
M = 10
out = []
for N in [5, 6, 7, 13, 14, 15, 21, 22, 23, 29, 30, 31, 37, 38, 39, 41, 46, 47]:
    E = EllipticCurve([0, 0, 0, -N**2, 0])
    rk = E.rank(only_use_mwrank=False)
    if rk != 1: out.append((N, rk, 'пропуск')); continue
    G = E.gens()[0]
    xs = {}; dens = {}
    for m in range(1, M + 1):
        P = (2 * m) * G
        xs[m] = P[0]; dens[m] = ZZ(P[0].denominator())
    # примитивные делители D_{2m}
    noprim = []
    for m in range(1, M + 1):
        r = dens[m]
        for k in range(1, m):
            g = gcd(r, dens[k])
            while g > 1:
                r //= g; g = gcd(r, dens[k])
        if r == 1: noprim.append(m)
    # прогрессии: x(2aG) + x(2cG) = 2 x(2bG), a,b,c различны по модулю (x(-P) = x(P))
    aps = [(a, b, c) for a, b, c in itertools.permutations(range(1, M + 1), 3) if a < c and xs[a] + xs[c] == 2 * xs[b]]
    out.append((N, 1, noprim, aps))
    print(f"N={N}: ранг 1, без примитивного делителя: {noprim}, прогрессий: {aps}", flush=True)
