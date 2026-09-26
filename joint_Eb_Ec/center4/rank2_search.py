# -*- coding: utf-8 -*-
"""Поиск конфигураций «три прямые с общим центром + линейная связь» РАНГА ≤ 2 в S^1(Q).
Точки окружности w = eps * (pi/pib)^r * (kap/kapb)^s, eps in {1,i,-1,-i}, |r|<=E, |s|<=F,
pi | p, kap | q — гауссовы простые (p, q = 1 mod 4, p<q<=PMAX).
Класс w по модулю квадратов в S^1(Q): (r mod 2, s mod 2, [eps = ±i]) (class(i)=class(1+i)=2).
Три прямые с общим центром a <=> три w одного класса (класс = класс a); u_j = a*Im(w_j).
Ищем связи типов 1–3 (как в search_twisted.py) между Im(w) внутри одного класса.
Нулевой класс (все w — квадраты) = квадратный центр (песочные часы/тройка линий магического квадрата).
Запуск: python3 rank2_search.py PMAX E"""
import sys, json, time
from itertools import combinations
from circle_lib import gauss_prime, split_primes
from sympy import primerange

PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
E = int(sys.argv[2]) if len(sys.argv) > 2 else 8

def gmul(z, w):
    return (z[0]*w[0] - z[1]*w[1], z[0]*w[1] + z[1]*w[0])

def gpow(z, k):
    r = (1, 0)
    for _ in range(k):
        r = gmul(r, z)
    return r

t0 = time.time()
primes = [p for p in primerange(5, PMAX + 1) if p % 4 == 1]
hits = []
checked = 0
for p, q in combinations(primes, 2):
    pi = gauss_prime(p); pib = (pi[0], -pi[1])
    ka = gauss_prime(q); kab = (ka[0], -ka[1])
    # w * p^E q^E = eps * pi^(E+r) pib^(E-r) kap^(E+s) kapb^(E-s)  -> гауссово целое G
    PP = [[gmul(gpow(pi, E + r), gpow(pib, E - r)) for r in range(-E, E + 1)]]
    Pp = [gmul(gpow(pi, E + r), gpow(pib, E - r)) for r in range(-E, E + 1)]
    Qq = [gmul(gpow(ka, E + s), gpow(kab, E - s)) for s in range(-E, E + 1)]
    byclass = {}
    for ir, r in enumerate(range(-E, E + 1)):
        for js, s in enumerate(range(-E, E + 1)):
            G = gmul(Pp[ir], Qq[js])
            for k, eps in enumerate(((1, 0), (0, 1), (-1, 0), (0, -1))):
                W = gmul(eps, G)
                cl = (r % 2, s % 2, k % 2)
                U = W[1]            # = p^E q^E * Im(w)
                if U == 0:
                    continue
                byclass.setdefault(cl, {})
                byclass[cl].setdefault(abs(U), (r, s, k))
    for cl, D in byclass.items():
        vals = sorted(D)
        S = set(vals)
        checked += len(vals)
        for i in range(len(vals)):
            x = vals[i]
            for j in range(i + 1, len(vals)):
                y = vals[j]
                rel = None
                if x + y in S: rel = (1, (x, y, x + y))
                elif (x + y) % 2 == 0 and (x + y) // 2 in S: rel = (2, (x, y, (x + y) // 2))
                elif (y - x) % 2 == 0 and (y - x) // 2 in S and (y - x) // 2 not in (x, y): rel = (3, ((y - x) // 2, x, y))
                elif x + 2 * y in S: rel = (3, (y, x, x + 2 * y))
                elif y + 2 * x in S: rel = (3, (x, y, y + 2 * x))
                if rel:
                    hits.append(dict(p=p, q=q, cls=cl, type=rel[0], U=rel[1], exps=[D[v] for v in rel[1]]))
print(json.dumps(dict(PMAX=PMAX, E=E, pairs=len(primes)*(len(primes)-1)//2, values=checked,
                      hits=len(hits), square_class_hits=[h for h in hits if h['cls'] == (0, 0, 0)],
                      first_hits=hits[:20], time_s=round(time.time()-t0, 1)), default=str, indent=1))
json.dump(dict(PMAX=PMAX, E=E, hits=hits), open(f'rank2_search_P{PMAX}_E{E}.json', 'w'), default=str)
