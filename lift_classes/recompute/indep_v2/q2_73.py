# -*- coding: utf-8 -*-
"""q2: наименьший нетривиальный положительный бесквадратный d = 1 mod 24 со всеми простыми = 1 mod 4.
Собственная реализация пробным делением, граница перебора явная (d <= 10000).
Плюс: перечень первых таких d и проверка, что 73 -- первое после 1."""
import json, os

def fac(n):
    out = {}; p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out

good = []
for d in range(1, 10001):           # явная граница
    if d % 24 != 1:
        continue
    f = fac(d)
    if all(e == 1 and p % 4 == 1 for p, e in f.items()):
        good.append(d)
below73 = [d for d in range(1, 73) if d % 24 == 1]
res = {"d=1 mod 24 below 73": {d: fac(d) for d in below73},
       "admissible (first 20)": good[:20],
       "smallest nontrivial": good[1],
       "count<=10000": len(good)}
print(res)
assert good[0] == 1 and good[1] == 73
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "q2_73.json"), "w"),
          indent=1, default=str)
