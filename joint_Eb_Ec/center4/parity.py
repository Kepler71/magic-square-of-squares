# -*- coding: utf-8 -*-
"""Корневые числа четырёх кривых E_b, E_c, E_{b+c}, E_{b−c} и 2-адика общего центра.
(1) 2-адический образ: при a=1 условия 1 ± u ∈ Q_2^{*2} <=> v_2(u) >= 3 (перебор по модулю 2^12).
(2) w(E_t) зависит только от знака t и 2-адического класса t; все 16 наборов знаков
    (w_b, w_c, w_{b+c}, w_{b−c}) реализуются парами целых B, C ≡ 0 (mod 8) => локальных (R, Q_2)
    ограничений на корневые числа нет.
(3) Данные Саллоуса: классы и корневые числа, ранги (Sage/mwrank) — иллюстрация того, что при
    НЕквадратном центре ранговые условия на E_{b±c} не обязаны выполняться.
Запуск: env DOT_SAGE=<scratch> python3 parity.py"""
import json
from sympy import factorint
out = {}

def sqfree(n):
    n = abs(n); s = 1
    for p, e in factorint(n).items():
        if e % 2: s *= p
    return s

def w_cong(t):
    N = sqfree(t)
    return 1 if N % 8 in (1, 2, 3) else -1

# (1) 2-адический образ при a = 1
def is_sq_Q2(x, prec=12):
    """x — целое != 0, проверка квадратности в Q_2 (x = 2^v * u)"""
    v = 0
    while x % 2 == 0:
        x //= 2; v += 1
    return v % 2 == 0 and x % 8 == 1

# u = 2^v * odd/1 и u = odd/2^k: проверяем 1 ± u
res = {}
for v in range(-4, 7):
    okv = False
    for odd in range(1, 256, 2):
        for sgn in (1, -1):
            if v >= 0:
                u_num, u_den = sgn * odd * 2 ** v, 1
            else:
                u_num, u_den = sgn * odd, 2 ** (-v)
            # 1 ± u = (u_den ± u_num)/u_den; u_den = 2^{-v} — квадрат по модулю квадратов, если -v чётно
            ok = True
            for s in (1, -1):
                num = u_den + s * u_num
                if num == 0: ok = False; break
                # квадратность num/u_den в Q_2: num * u_den — квадрат
                if not is_sq_Q2(num * u_den): ok = False; break
            okv = okv or ok
    res[v] = okv
out['Q2_image_1pm_u_square_by_v2'] = res
print('1±u — квадраты в Q_2 возможны при v_2(u) =', [v for v, o in res.items() if o])

# (2) все 16 наборов корневых чисел при B, C ≡ 0 mod 8
pats = {}
R = 60
for B in range(-8 * R, 8 * R + 1, 8):
    for C in range(-8 * R, 8 * R + 1, 8):
        if B == 0 or C == 0 or B == C or B == -C:
            continue
        pat = (w_cong(B), w_cong(C), w_cong(B + C), w_cong(B - C))
        if pat not in pats:
            pats[pat] = (B, C)
out['root_number_patterns_found'] = {str(k): v for k, v in pats.items()}
print('реализовано наборов (w_b,w_c,w_{b+c},w_{b-c}):', len(pats), 'из 16')

# (3) Саллоус
F = json.load(open('/home/kep/magicKube/fable_symmetry/s1_sallows.json'))
a, b, c = F['a'], F['b'], F['c']
data = {}
for name, t in (('b', b), ('c', c), ('b+c', b + c), ('b-c', b - c), ('a*c', a * c), ('a*(b+c)', a * (b + c)), ('a*(b-c)', a * (b - c))):
    data[name] = dict(t=t, N=sqfree(t), N_mod8=sqfree(t) % 8, w=w_cong(t))
try:
    from sage.all import EllipticCurve
    for name in data:
        N = data[name]['N']
        E = EllipticCurve([-N * N, 0])
        data[name]['rank'] = int(E.rank(only_use_mwrank=True))
except Exception as ex:
    out['sage_error'] = str(ex)
out['sallows_curves'] = data
for k, v in data.items(): print('Саллоус', k, v)
json.dump(out, open('parity.json', 'w'), indent=1, default=str)
