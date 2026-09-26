# -*- coding: utf-8 -*-
"""t2: следствие 73. Наименьшее d>1, бесквадратное, d=1 mod 24, все простые делители = 1 mod 4.
Перебор d <= DMAX (ограниченный цикл) двумя способами: своё решето и Sage factor.
Также: числа = 1 mod 24 меньше 73 и их разложения; первые 20 допустимых d.
"""
import json, os
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
DMAX = 200000

spf = list(range(DMAX + 1))
for i in range(2, int(DMAX ** 0.5) + 1):
    if spf[i] == i:
        for j in range(i * i, DMAX + 1, i):
            if spf[j] == j:
                spf[j] = i

def ok_sieve(d):
    n = d; last = 0
    while n > 1:
        p = spf[n]; n //= p
        if p == last or p % 4 != 1:   # повтор = не бесквадратно
            return False
        last = p
    return d % 24 == 1

A = [d for d in range(2, DMAX + 1) if ok_sieve(d)]
from sage.all import factor, is_squarefree
B = [d for d in range(2, DMAX + 1) if d % 24 == 1 and is_squarefree(d)
     and all(p % 4 == 1 for p, _ in factor(d))]
below73 = [(d, str(factor(d))) for d in range(1, 73) if d % 24 == 1]
res = dict(DMAX=DMAX, min_A=A[0], min_B=B[0], lists_equal=(A == B), count=len(A),
           first20=A[:20], below73=below73,
           no_mod4_condition_first=[d for d in range(2, 400) if d % 24 == 1 and is_squarefree(d)][:12])
print(json.dumps(res, indent=1))
json.dump(res, open("t2_73.json", "w"), indent=1)
