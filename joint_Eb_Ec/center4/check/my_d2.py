# Независимый пересчёт утверждения 6 (доля центров, отсекаемых D2). Метод: DFS по произведениям простых p = 1 (mod 4).
import sys, json
from math import isqrt
N = int(sys.argv[1])
sv = bytearray([1]) * (N + 1); sv[0] = sv[1] = 0
for i in range(2, isqrt(N) + 1):
    if sv[i]: sv[i*i::i] = bytearray(len(sv[i*i::i]))
P1 = [p for p in range(5, N + 1, 4) if sv[p]]
tot = exc = exc_pp = 0; exc_bigp = 0; examples = []
def dfs(idx, g, fac):
    global tot, exc, exc_pp, exc_bigp
    if g > 1:
        d = 1
        for p, e in fac: d *= 2 * e + 1
        if d >= 9:
            tot += 1
            if any(p ** (2 * e) > 2 * g for p, e in fac):
                exc += 1
                if len(fac) == 1: exc_pp += 1
                pmax = max(p for p, e in fac)
                if pmax * pmax > 2 * g: exc_bigp += 1
                if g in (65,) or (len(examples) < 5 and len(fac) >= 3): examples.append((g, fac[:]))
    for j in range(idx, len(P1)):
        p = P1[j]
        if g * p > N: break
        gg, e = g, 0
        while gg * p <= N:
            gg *= p; e += 1
            fac.append((p, e)); dfs(j + 1, gg, fac); fac.pop()
sys.setrecursionlimit(10000)
dfs(0, 1, [])
r = dict(N=N, total=tot, excluded=exc, frac=exc / tot, excluded_prime_powers=exc_pp, excluded_with_pmax2_gt_2g=exc_bigp,
         examples=examples)
print(r); json.dump(r, open(f'my_d2_N{N}.json', 'w'), indent=1)
