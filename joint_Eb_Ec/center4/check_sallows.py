# -*- coding: utf-8 -*-
"""Контроль на семилинейных квадратах (Саллоус и 2P/3P из FABLE_SYMMETRY) и на AB1 Бремнера–Саллоуса.
У семилинейных квадратов шесть клеток a±c, a±(b+c), a±(b−c) — квадраты при НЕквадратном a:
это три линии с общим (неквадратным) центром и связью (b+c) − (b−c) − 2c = 0 (N_L = 4).
Теоремы B и D доказаны для любого общего центра (квадратного или нет) — здесь их гипотезы выполнены."""
import json
from fractions import Fraction
from math import isqrt, gcd
from circle_lib import *

out = {}
F = json.load(open('/home/kep/magicKube/fable_symmetry/s1_sallows.json'))
a, b, c, lam = F['a'], F['b'], F['c'], F['lam']
print('Саллоус: a,b,c,lam =', a, b, c, lam)
# проверка семилинейного квадрата
sq = [[a-b+lam, a+b+c, a-c], [a+b-c, a+lam, a-b+c], [a+c, a-b-c, a+b+lam]]
print('клетки:', sq, 'все квадраты:', all(is_sq(x) for r in sq for x in r))
sq2P = [[Fraction(x) for x in r] for r in F['new_seven_line_squares'][0]]
sq3P = [[Fraction(x) for x in r] for r in F['new_seven_line_squares'][1]]
for name, S in (('2P', sq2P), ('3P', sq3P)):
    lamS = S[1][1] - a
    six = [S[0][1], S[0][2], S[1][0], S[1][2], S[2][0], S[2][1]]
    same = six == [a+b+c, a-c, a+b-c, a-b+c, a+c, a-b-c]
    print(name, ': шесть клеток вне главной диагонали те же:', same, '; lam =', lamS, '(≠0)')
    out[name] = dict(six_cells_same=same, lam=str(lamS))

def analyse(tag, a, us, ns):
    NL = sum(abs(n) for n in ns)
    assert sum(n*u for n, u in zip(us, ns)) == 0
    res = {}
    res['a'] = a; res['us'] = us; res['ns'] = ns
    res['a_is_square'] = is_sq(a)
    res['factor_a'] = {str(p): e for p, e in factorint(a).items()}
    res['lines'] = {str(u): {k: (list(v) if isinstance(v, tuple) else v) for k, v in line_data(a, u).items()} for u in us}
    prof = local_profile(a, us)
    res['profile'] = {str(p): {'e': d['e'], 'pi': d['pi'],
                     'lines': {str(u): d['lines'][u] for u in us}} for p, d in prof.items()}
    rk, M = circle_rank(a, us)
    res['circle_rank'] = rk; res['valuation_matrix'] = M
    ok1, bad1 = min_twice_ok(a, us); res['lemma_min_twice'] = ok1
    ok2, cls = balance_check(a, us, NL)
    res['thmD_ok'] = ok2
    res['thmD_classes'] = [(list(s), ps, Pi, str(r)) for s, ps, Pi, r in cls]
    Zs = []
    for s, ps, Pi, r in cls:
        Z, NZ, Pi2, div = relation_Z(a, us, ns, ps)
        Zs.append(dict(cls=ps, Z=Z, NZ=NZ, Pi2=Pi2, divisible=div, NZ_ge_Pi2=NZ >= Pi2, NZ_le_bound=NZ <= (NL*a//Pi)**2 if Pi else None))
    res['Z_checks'] = Zs
    print(f'--- {tag}: a={a} ({"квадрат" if is_sq(a) else "не квадрат"}) = {factorint(a)}, u={us}, связь n={ns}')
    for p, d in prof.items():
        print(f'   p={p} e={d["e"]}: ', {u: (d["lines"][u]["r"], d["lines"][u]["rb"], "чист" if d["lines"][u]["pure"] else f't={d["lines"][u]["t"]}', d['lines'][u]['sigma']) for u in us})
    print(f'   ранг <w_u> в S^1(Q): {rk}; матрица v_pi(w_u): {M}')
    print(f'   лемма (min t достигается ≥2 раз): {ok1} {bad1}')
    print(f'   теорема D: {ok2}; классы (ориентация, простые, Pi, Pi^2/(N_L a)):', [(s, ps, Pi, str(r)) for s, ps, Pi, r in cls])
    for z in Zs: print('   Z-контроль:', z)
    return res

# шесть клеток Саллоуса в исходном масштабе и после деления на 4
us = [c, b+c, b-c]; ns = [-2, 1, -1]
out['sallows_6cells_scale1'] = analyse('Саллоус, 6 клеток, a=4420', a, us, ns)
g = 4
out['sallows_6cells_primitive'] = analyse('Саллоус, 6 клеток /4, a=1105', a//g, [u//g for u in us], ns)

# Прогноз теоремы D для a'=1105: допустимы только классы с Pi^2 <= 4*1105=4420
print('Прогноз D для a=1105: допустимые классы Pi: 5,13,17,65 (65^2=4225<=4420); запрещены 85,221,1105')

# AB1 Бремнера–Саллоуса: центр 425^2, полные линии через центр: c и b+c (две линии; связи нет)
A, B, C = 425**2, 41496, -138600
lines = [u for u in (B, C, B+C, B-C) if is_sq(A-u) and is_sq(A+u)]
print('AB1: полные линии через квадратный центр:', lines, '(две — гипотезы теорем B/D не выполнены)')
rk, M = circle_rank(A, lines)
print('AB1: ранг двух точек окружности', rk, M)
out['AB1'] = dict(full_lines=lines, rank=rk, M=M)
json.dump(out, open('check_sallows.json', 'w'), indent=1, default=str)
