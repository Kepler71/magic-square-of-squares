# -*- coding: utf-8 -*-
"""Насколько теорема D (p^{2 v_p(γ)} <= 2γ) и теорема B (γ не степень простого) сужают перебор центров:
γ <= N, все простые делители ≡ 1 (mod 4) (Rabern), не менее 4 треугольников с гипотенузой γ: prod(2e_p+1) >= 9."""
import sys, json
from sympy import primerange
N = int(sys.argv[1]) if len(sys.argv) > 1 else 10**7
P = [p for p in primerange(5, N + 1) if p % 4 == 1]
cnt = dict(N=N, candidates=0, prime_powers=0, D2_excluded=0, D2_excluded_not_prime_power=0, examples=[])
def rec(i, g, fac):
    # fac: список (p, e)
    if fac:
        prod = 1
        for p, e in fac: prod *= 2 * e + 1
        if prod >= 9:
            cnt['candidates'] += 1
            if len(fac) == 1: cnt['prime_powers'] += 1
            mx = max(p ** (2 * e) for p, e in fac)
            if mx > 2 * g:
                cnt['D2_excluded'] += 1
                if len(fac) > 1:
                    cnt['D2_excluded_not_prime_power'] += 1
                    if len(cnt['examples']) < 10: cnt['examples'].append((g, fac))
    for j in range(i, len(P)):
        p = P[j]
        if g * p > N: break
        q, e = g * p, 1
        while q <= N:
            rec(j + 1, q, fac + [(p, e)])
            q *= p; e += 1
rec(0, 1, [])
cnt['share_excluded'] = round(cnt['D2_excluded'] / max(1, cnt['candidates']), 4)
print(cnt)
json.dump(cnt, open(f'd2_pruning_N{N}.json', 'w'), indent=1)
