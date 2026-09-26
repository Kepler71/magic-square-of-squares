# -*- coding: utf-8 -*-
"""(a) Перебор в исходных переменных: r, s mod p^k (не оба делятся на p), z = p^e * w, w единица mod p^k,
       e в [-4, k]. Утверждение: если e < THR, шар никогда не PASS; если e >= THR, шар PASS.
   (b) Арифметика следствий: наименьший нетривиальный d (squarefree, d>0, d=1 mod 24, p|d => p=1 mod 4);
       списки D_-, D_+ для трёх наклонов в обеих ориентациях и со сменой знака."""
import json, time
import p23lib as L
from sympy import factorint

out = {}
t0 = time.time()
for p, k in ((2, 6), (3, 4)):
    pk = p ** k
    tal = {}
    for e in range(-4, k + 1):                   # явные границы
        for r in range(pk):
            for s in range(pk):
                if r % p == 0 and s % p == 0:
                    continue
                for w in range(1, pk):
                    if w % p == 0:
                        continue
                    if e >= 0:
                        m, K = 0, k + e
                        B = r * w * p ** e; C = s * w * p ** e
                    else:
                        m, K = -e, k
                        B = r * w; C = s * w
                    cat = L.classify_naive(L.ball_cells(B, C, m, p, K), p)
                    key = f'e={e}:{cat}'
                    tal[key] = tal.get(key, 0) + 1
        print(p, e, {kk: v for kk, v in tal.items() if kk.startswith(f'e={e}:')}, f'{time.time()-t0:.0f}s', flush=True)
    bad = [kk for kk in tal if (int(kk.split(':')[0][2:]) < L.THR[p] and kk.endswith('PASS')) or
           (int(kk.split(':')[0][2:]) >= L.THR[p] and not kk.endswith('PASS'))]
    out[f'rsz_p{p}_k{k}'] = {'tally': tal, 'violations': bad}

def ok_d(d):
    if d <= 0 or d % 24 != 1:
        return False
    f = factorint(d)
    return all(v == 1 and q % 4 == 1 for q, v in f.items())

small = [d for d in range(1, 2000) if ok_d(d)]
out['admissible_d_below_2000'] = small[:20]

def Dset(n):
    n = abs(n)
    primes = [q for q in factorint(n)]
    res = set()
    for mask in range(1 << len(primes)):     # конечное: 2^#простых
        d = 1
        for t, q in enumerate(primes):
            if mask >> t & 1:
                d *= q
        if ok_d(d):
            res.add(d)
    return sorted(res)

slopes = {}
for (a, b) in ((126, 451), (73, 362), (265, 298)):
    for (r, s) in ((a, b), (b, a), (-a, b), (a, -b)):
        slopes[f'r={r},s={s}'] = {'r-s': r - s, 'r+s': r + s,
                                  'factor(r-s)': str(factorint(abs(r - s))),
                                  'factor(r+s)': str(factorint(abs(r + s))),
                                  'D_minus(d_T)': Dset(r - s), 'D_plus(d_L)': Dset(r + s)}
out['slopes'] = slopes
print(json.dumps(out, indent=1))
json.dump(out, open('rsz_and_arith.json', 'w'), indent=1)
