# Проверка звеньев шага «подгруппа ранга 1 => P_j = 2 m_j G» и применения примитивных делителей
# на кривых E_N ранга >= 2, для НЕбазисных примитивных направлений G = sum c_i G_i (+ T из E[2]).
# (1) у B_{2M}(G), M = 2..MMAX, есть примитивный простой (относительно всех B_k, k < 2M);
# (2) слабая форма, которая реально нужна: простой l | B_{2M}(G), l не делит B_{2m}(G) при m < M;
# (3) среди x(2mG), 1 <= m <= MMAX, нет трёх в арифметической прогрессии;
# (4) независимость от подъёма: x(2m(G+T)) = x(2mG);
# (5) решёточная лемма: 2R лежит на прямой ZG + E[2] => R = mG + T (перебор малых R).
# Это численный контроль, не доказательство.
from sage.all import *
import sys, itertools, json, time
MMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 8
CMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 3
NLIST = [int(t) for t in sys.argv[3].split(',')] if len(sys.argv) > 3 else [34, 41, 65, 137, 210, 1254]

def has_prim(Bn, prev):
    r = Bn
    for b in prev:
        g = gcd(r, b)
        while g > 1:
            r //= g; g = gcd(r, g)
            g = gcd(r, b)
    return r > 1

summary = {}
t0 = time.time()
for Nv in NLIST:
    E = EllipticCurve([0,0,0,-Nv**2,0])
    gens = E.gens()
    r = len(gens)
    T2 = [E(0), E(0,0), E(Nv,0), E(-Nv,0)]
    stats = dict(rank=r, directions=0, prim_fail=[], weak_fail=[], ap_found=[], lift_fail=0, lattice_fail=0)
    rng = range(-CMAX, CMAX+1)
    for c in itertools.product(rng, repeat=r):
        if all(ci == 0 for ci in c): continue
        if gcd(list(c)) != 1: continue
        # одно направление из пары ±
        first = next(ci for ci in c if ci != 0)
        if first < 0: continue
        G = sum((ci*g for ci, g in zip(c, gens)), E(0))
        stats['directions'] += 1
        B = {}
        X = {}
        for k in range(1, 2*MMAX+1):
            xk = (k*G)[0]
            B[k] = xk.denominator(); X[k] = xk
        for M in range(2, MMAX+1):
            n = 2*M
            if not has_prim(B[n], [B[k] for k in range(1, n)]):
                stats['prim_fail'].append((c, n))
            if not has_prim(B[n], [B[2*m] for m in range(1, M)]):
                stats['weak_fail'].append((c, n))
        xs = [X[2*m] for m in range(1, MMAX+1)]
        for i, j, k in itertools.permutations(range(MMAX), 3):
            if i < k and xs[i] + xs[k] == 2*xs[j]:
                stats['ap_found'].append((c, i+1, j+1, k+1))
        for T in T2[1:]:
            for m in (1, 2, 3):
                if (2*m*(G+T))[0] != X[2*m]: stats['lift_fail'] += 1
    # (5) решёточная лемма на малых R: если 2R = kG + T для направления G, то k чётно и R - (k/2)G в E[2]
    if r >= 2:
        G = gens[0] + 2*gens[1]  # примитивное направление (1,2,0,...)
        for c in itertools.product(range(-4, 5), repeat=r):
            R = sum((ci*g for ci, g in zip(c, gens)), E(0))
            # 2R на прямой <=> 2c пропорционален (1,2,0..)
            v = [2*ci for ci in c]
            base = [1, 2] + [0]*(r-2)
            if v[1] == 2*v[0] and all(vi == 0 for vi in v[2:]):
                k = v[0]
                if k % 2 != 0: stats['lattice_fail'] += 1; continue
                D = R - (k//2)*G
                if 2*D != E(0): stats['lattice_fail'] += 1
    stats['time'] = float(time.time() - t0)
    summary[Nv] = stats
    print("N=%d rank=%d направлений=%d prim_fail=%s weak_fail=%s AP=%s lift_fail=%d lattice_fail=%d t=%.1fs" % (
        Nv, r, stats['directions'], stats['prim_fail'], stats['weak_fail'], stats['ap_found'],
        stats['lift_fail'], stats['lattice_fail'], time.time() - t0), flush=True)
json.dump({str(k): {kk: (vv if not isinstance(vv, list) else [str(z) for z in vv]) for kk, vv in v.items()} for k, v in summary.items()},
          open('lattice_prim_M%d_C%d.json' % (MMAX, CMAX), 'w'), indent=1)
