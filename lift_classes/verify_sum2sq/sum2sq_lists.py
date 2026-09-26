# Независимая проверка шагов 4-5 и следствий теоремы об остаточных классах (Codex, 26.09).
# A: тождества пар (сумма 2), клетки T и L;  B: лемма о сумме двух квадратов при p=3 mod 4
#    (полный перебор по F_p) и локальная независимая альтернатива шагу 5 (нет полюсов при p=3 mod 4);
# C: списки D_-, D_+ для трёх наклонов своим кодом, с условием p=1 mod 4 и без него;
# D: перебор всех взаимно простых 1<=|r|,|s|<=RMAX: следствие <73, gcd(d_T,d_L)=1 (вкл. r,s нечётные).
# Все циклы с явными конечными границами.
import sys, time, json, random
from math import gcd
from itertools import combinations
import sympy as sp

RMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 300
out = {}
# ---- A
r, s, z = sp.symbols('r s z')
f = {(i, j): 1 + (i*r + j*s)*z for i in (-1, 0, 1) for j in (-1, 0, 1)}
assert all(sp.expand(f[(i, j)] + f[(-i, -j)] - 2) == 0 for (i, j) in f)
T = (1 - r*z)*(1 + (r+s)*z)*(1 - s*z); L = (1 - r*z)*(1 + (r-s)*z)*(1 + s*z)
assert sp.expand(T - f[(-1, 0)]*f[(1, 1)]*f[(0, -1)]) == 0
assert sp.expand(L - f[(-1, 0)]*f[(1, -1)]*f[(0, 1)]) == 0
# клетки (i,j) линии суммируются в (0,0) => это магические линии
assert (-1+1+0, 0+1-1) == (0, 0) and (-1+1+0, 0-1+1) == (0, 0)
print('A: opposite pairs sum to 2; T=f(-1,0)f(1,1)f(0,-1), L=f(-1,0)f(1,-1)f(0,1): OK', flush=True)
# ---- B
P = list(sp.primerange(3, 1000))
bad = [p for p in P if p % 4 == 3 and any((u*u + v*v) % p == 0 for u in range(1, p) for v in range(p))]
sharp = [p for p in P if p % 4 == 1 and not any((u*u + 1) % p == 0 for u in range(1, p))]
assert not bad and not sharp
print('B1: p=3 mod 4 (<1000): u^2+v^2=0 mod p only trivially; p=1 mod 4: -1 is QR (sharp): OK', flush=True)
# локальная альтернатива шагу 5: при p=3 mod 4 полюс b запрещён (b^-2 - 1 = -1 mod p невычет)
rng = random.Random(20260926); cnt = 0
for p in [q for q in P if q % 4 == 3][:30]:
    for _ in range(2000):                        # явная граница
        n = rng.randint(1, 6); u = rng.randint(1, 10**6)
        if u % p == 0: continue
        w = rng.randint(1, 10**6)
        if w % p == 0: continue
        # b = u/(w p^n): 1-b^2 = (w^2 p^{2n} - u^2)/(w^2 p^{2n}); числитель — единица с остатком -u^2
        num = w*w*p**(2*n) - u*u
        assert num % p != 0 and pow(num % p, (p-1)//2, p) == p-1
        cnt += 1
print(f'B2: {cnt} random poles at p=3 mod 4: 1-b^2 never a local square: OK', flush=True)

# ---- C
def pf(n):
    n = abs(n); o = []; q = 2
    while q*q <= n:
        if n % q == 0:
            o.append(q)
            while n % q == 0: n //= q
        q += 1
    if n > 1: o.append(n)
    return o
def D(n, need14=True):
    ps = [q for q in pf(n) if (q % 4 == 1 or not need14)]
    ds = []
    for k in range(len(ps) + 1):
        for c in combinations(ps, k):
            d = 1
            for q in c: d *= q
            if d % 24 == 1: ds.append(d)
    return sorted(ds)
tab = {}
for (a, b) in [(126, 451), (73, 362), (265, 298)]:
    tab[f'{a}/{b}'] = dict(r_minus_s=a-b, r_plus_s=a+b, Dminus=D(a-b), Dplus=D(a+b),
                           Dminus_no_mod4=D(a-b, False), Dplus_no_mod4=D(a+b, False))
    print('C:', a, b, tab[f'{a}/{b}'], flush=True)
assert tab['126/451']['Dminus'] == [1] and tab['126/451']['Dplus'] == [1, 577]
assert tab['73/362']['Dminus'] == [1] and tab['73/362']['Dplus'] == [1, 145]
assert tab['265/298']['Dminus'] == [1] and tab['265/298']['Dplus'] == [1]
same = all(v['Dminus'] == v['Dminus_no_mod4'] and v['Dplus'] == v['Dplus_no_mod4'] for v in tab.values())
print('C: lists without the p=1 mod 4 condition coincide (step 5 not load-bearing here):', same, flush=True)
out['table'] = tab; out['table_mod4_not_needed'] = same
# минимальный нетривиальный squarefree d=1 mod 24 (с условием и без условия p=1 mod 4)
sqf = lambda n: all(n % (q*q) for q in range(2, int(n**0.5) + 1))
m_all = min(d for d in range(2, 10**4) if d % 24 == 1 and sqf(d))
m_14 = min(d for d in range(2, 10**4) if d % 24 == 1 and sqf(d) and all(q % 4 == 1 for q in pf(d)))
print('C: min nontrivial squarefree d=1 mod24:', m_all, '; with p=1 mod4:', m_14, flush=True)
assert m_all == 73 and m_14 == 73

# ---- D
t0 = time.time(); npairs = 0; nboth_odd = 0; viol = []
for a in range(-RMAX, RMAX + 1):
    if a == 0: continue
    for b in range(1, RMAX + 1):                    # s>0 без потери общности (общий знак z)
        if gcd(a, b) != 1 or a == b or a == -b: continue
        npairs += 1
        Dm, Dp = D(a - b), D(a + b)
        if a % 2 and b % 2: nboth_odd += 1
        if abs(a - b) < 73 and Dm != [1]: viol.append(('cor-', a, b, Dm))
        if abs(a + b) < 73 and Dp != [1]: viol.append(('cor+', a, b, Dp))
        for dT in Dm:
            if (a - b) % dT: viol.append(('div-', a, b, dT))
            for dL in Dp:
                if gcd(dT, dL) != 1: viol.append(('gcd', a, b, dT, dL))
        for dL in Dp:
            if (a + b) % dL: viol.append(('div+', a, b, dL))
    if (a + RMAX) % (RMAX // 2) == 0:
        print(f'  D: r={a} pairs={npairs} viol={len(viol)} {time.time()-t0:.0f}s', flush=True)
print(f'D: {npairs} coprime pairs (|r|,s<={RMAX}), {nboth_odd} with r,s odd (2|r-s, 2|r+s); violations: {len(viol)}', flush=True)
out['D'] = dict(RMAX=RMAX, pairs=npairs, both_odd=nboth_odd, violations=viol[:20])
json.dump(out, open('/home/kep/magicKube/lift_classes/verify_sum2sq/sum2sq_lists.json', 'w'), indent=1)
assert not viol
print('ALL CHECKS PASSED')
