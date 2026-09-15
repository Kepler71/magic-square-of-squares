# Fable 15.09. Контроль структуры кратностей на случайных наклонах: (i) внутри генерического класса все кривые
# ℚ-изоморфны; (ii) разные классы не изогенны (подпись a_p по p<300 + is_isogenous для совпадающих подписей).
from sage.all import *
import sys, json, random
from classes import CLASSES, curve_of

random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
N = int(sys.argv[2]) if len(sys.argv) > 2 else 30
SMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 3000
res = []
while len(res) < N:
    s = random.randint(3, SMAX); r = random.randint(1, s - 1)
    if gcd(r, s) != 1 or 2 * r == s: continue
    reps = {}; ok_iso = True; bad = []
    for c, lst in CLASSES.items():
        Es = []
        for S in lst:
            E, L, c0, T = curve_of(S, r, s)
            if E.discriminant() == 0: bad.append(T); continue
            Es.append(E.minimal_model())
        for E in Es[1:]:
            if not Es[0].is_isomorphic(E): ok_iso = False
        reps[c] = Es[0]
    sig = {}
    for c, E in reps.items():
        N0 = E.conductor()
        sig.setdefault(tuple(E.ap(p) if N0 % p else None for p in primes(3, 300)), []).append(c)
    merges = []
    for k, cs in sig.items():
        if len(cs) > 1:
            for c2 in cs[1:]:
                if reps[cs[0]].is_isogenous(reps[c2]): merges.append([len(CLASSES[cs[0]]), len(CLASSES[c2])])
    res.append(dict(slope=f'{r}/{s}', iso_in_class=ok_iso, classes=len(reps), merges=merges, degenerate=bad))
    print(res[-1], flush=True)
json.dump(res, open('/home/kep/magicKube/isogeny_mult/fable/mult_control.json', 'w'))
print('наклонов:', len(res), 'изоморфизм внутри классов всюду:', all(d['iso_in_class'] for d in res),
      'слияний классов:', sum(len(d['merges']) for d in res))
