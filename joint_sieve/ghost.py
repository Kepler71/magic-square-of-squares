from jsieve import *
from run_slope import build_factors
import json, sys
r, s = int(sys.argv[1]), int(sys.argv[2])
d = json.load(open(f'res_{r}_{s}.json'))
facs, info = build_factors(r, s)
byname = {f.name: f for f in facs}
F = [byname[','.join(map(str, x['T']))] for x in d['factors']]
N = d['N']
ex = d['extra_list']
print('невырожденных классов', len(ex), 'пример', ex[0])
# разбираем первый класс: для каждого множителя ранга 0 находим точку кручения и её z
def split(t):
    out = []; p = 0
    for f in F:
        out.append(t[p:p + f.dim]); p += f.dim
    return out
for t in ex[:3]:
    parts = split(t); desc = []
    for f, c in zip(F, parts):
        if f.rank == 0:
            P = f.E(0)
            for k, g in zip(c, f.tgens): P = P + k * g
            desc.append(f'{f.name}: z={f.z_of(P)}')
        else:
            desc.append(f'{f.name}: n={c[:f.rank]} t={c[f.rank:]}')
    print(' | '.join(desc[:8]))
print('--- множители ранга 1: класс призрака vs класс P_T(0); z малых представителей')
t = ex[0]; parts = split(t)
for f, c, Ni in zip(F, parts, N):
    if f.rank == 0: continue
    p0 = f.decompose(f.point(0, 1))
    zs = []
    for n in (c[0], c[0] - Ni):
        P = n * f.gens[0]
        for k, g in zip(c[1:], f.tgens): P = P + k * g
        zs.append((n, str(f.z_of(P)), round(float(P.height()), 3)))
    print(f'{f.name}: N={Ni}, призрак {c}, P0 {p0}, sat_index {f.sat_index}, представители {zs}')
