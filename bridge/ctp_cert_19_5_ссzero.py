#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
СЕРТИФИКАТ: класс delta = (1,193,193) НЕ лежит в образе E(Q)/2E(Q) для пары G1 (m,n)=(19,5).
Доказательство: спаривание Касселса-Тейта <[g1],[g2]>_CT = -1  (Fisher, arXiv:2208.14977, Thm 3.1).
Ранг E(Q) нигде не используется.

Проверяющий скрипт: ТОЛЬКО стандартная библиотека Python (fractions, math, json).
Ни Sage, ни PARI, ни какой-либо код проекта не импортируется.

    python3 ctp_cert_19_5_ссzero.py [путь-к-certificate.json]

Читает ctp_cert_19_5.json (рядом со скриптом), пишет ctp_cert_19_5_result.json.
"""
import json, sys, os
from fractions import Fraction as F
from math import isqrt, gcd
from functools import reduce

OK = []
def check(name, cond):
    if not cond:
        raise AssertionError('ПРОВАЛ: ' + name)
    OK.append(name)

# ----------------------------------------------------------------- базовая арифметика
def is_square(q):
    q = F(q)
    if q < 0: return False
    return isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator

def sqrt_exact(q):
    q = F(q)
    assert is_square(q)
    return F(isqrt(q.numerator), isqrt(q.denominator))

def factorize(n):
    """честное разложение пробными делителями; все встречающиеся числа малы по простым"""
    n = abs(int(n)); out = {}
    if n == 0: raise ValueError('0')
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

def support(q):
    q = F(q)
    if q == 0: return set()
    return set(factorize(q.numerator)) | set(factorize(q.denominator))

def val_unit(q, p):
    q = F(q); assert q != 0
    num, den, k = q.numerator, q.denominator, 0
    while num % p == 0: num //= p; k += 1
    while den % p == 0: den //= p; k -= 1
    return k, F(num, den)

def res_mod(q, mod):
    q = F(q)
    return q.numerator * pow(q.denominator, -1, mod) % mod

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1

def local_square(q, place):
    q = F(q)
    if q == 0: return False
    if place == 'real': return q > 0
    p = int(place); k, u = val_unit(q, p)
    if k % 2: return False
    if p == 2: return res_mod(u, 8) == 1
    return legendre(res_mod(u, p), p) == 1

def hilbert(a, b, place):
    """символ Гильберта (a,b)_v для a,b in Q*, реализация с нуля"""
    a, b = F(a), F(b)
    assert a != 0 and b != 0
    if place == 'real':
        return -1 if (a < 0 and b < 0) else 1
    p = int(place)
    va, ua = val_unit(a, p); vb, ub = val_unit(b, p)
    if p == 2:
        ra, rb = res_mod(ua, 8), res_mod(ub, 8)
        e = ((ra-1)//2)*((rb-1)//2) + va*((rb*rb-1)//8) + vb*((ra*ra-1)//8)
        return -1 if e % 2 else 1
    ra, rb = res_mod(ua, p), res_mod(ub, p)
    r = 1
    if (va*vb) % 2 and p % 4 == 3: r = -r
    if vb % 2 and legendre(ra, p) == -1: r = -r
    if va % 2 and legendre(rb, p) == -1: r = -r
    return r

# самопроверка символа Гильберта: формула взаимности на случайных парах
def _selftest_hilbert():
    import random
    random.seed(12345)
    for _ in range(300):
        a = F(random.randint(-300, 300), random.randint(1, 20))
        b = F(random.randint(-300, 300), random.randint(1, 20))
        if a == 0 or b == 0: continue
        pls = ['real'] + sorted(support(a) | support(b) | {2})
        pr = 1
        for pl in pls: pr *= hilbert(a, b, pl if pl == 'real' else str(pl))
        assert pr == 1, ('нарушена формула взаимности', a, b)
        for pl in pls:
            q = pl if pl == 'real' else str(pl)
            assert hilbert(a, b, q) == hilbert(b, a, q)
            assert hilbert(a, a*b*b, q) == hilbert(a, a, q)
            if a + b != 0 or True:
                pass
        # (a, -a)_v = 1 и (a, 1-a)_v = 1 при a != 0,1
        for pl in pls:
            q = pl if pl == 'real' else str(pl)
            assert hilbert(a, -a, q) == 1
            if a != 1:
                assert hilbert(a, 1-a, q) == 1
_selftest_hilbert()
OK.append('символ Гильберта: взаимность, симметрия, (a,-a)=1, (a,1-a)=1 на 300 случайных парах')

# ----------------------------------------------------------------- бинарные формы
def pmul(u, v):
    out = [F(0)]*(len(u)+len(v)-1)
    for i, a in enumerate(u):
        for j, c in enumerate(v):
            out[i+j] += a*c
    return out

def padd(*ps):
    n = max(len(p) for p in ps)
    out = [F(0)]*n
    for p in ps:
        off = n - len(p)
        for i, c in enumerate(p): out[off+i] += c
    return out

def pscal(c, p): return [F(c)*F(x) for x in p]
def pzero(p): return all(F(x) == 0 for x in p)
def ev(p, x, z):
    k = len(p)-1
    return sum(F(p[i])*F(x)**(k-i)*F(z)**i for i in range(k+1))

def quart_I(g):
    a, b, c, d, e = [F(v) for v in g]
    return 12*a*e - 3*b*d + c*c
def quart_J(g):
    a, b, c, d, e = [F(v) for v in g]
    return 72*a*c*e - 27*a*d*d - 27*b*b*e + 9*b*c*d - 2*c**3
def hess(g):
    a, b, c, d, e = [F(v) for v in g]
    return [3*b*b-8*a*c, 4*(b*c-6*a*d), 2*(2*c*c-24*a*e-3*b*d), 4*(c*d-6*b*e), 3*d*d-8*c*e]

# ----------------------------------------------------------------- читаем сертификат
here = os.path.dirname(os.path.abspath(__file__))
path = sys.argv[1] if len(sys.argv) > 1 else None
if path is None:
    for cand in [os.path.join(here, 'ctp_cert_19_5.json'),
                 os.path.join(here, 'ctp_19_5_claude', 'certificate_19_5.json')]:
        if os.path.exists(cand): path = cand; break
cert = json.loads(open(path, encoding='utf-8').read())
print('сертификат:', path)

# ================================================================= 1. исходная кривая C
m, n = int(cert['m']), int(cert['n'])
check('(m,n) = (19,5)', (m, n) == (19, 5))
s = F(m*m + n*n, 2)
b = s*m*m*n*n
check('s = (m^2+n^2)/2 = 193', s == 193)
check('b = s m^2 n^2', b == F(cert['b']) and b == 1741825)
e = [-b, -s*m**4, -s*n**4]
check('корни E совпали с сертификатом', e == [F(v) for v in cert['E_roots']])

# тождества X - e_i = c_i * F_j как многочлены от t (коэффициенты при t^0, t^1, t^2)
F0 = [F(m*m), F(0), F(n*n)]      # m^2 + n^2 t^2
F4 = [s, F(0), s]                # s(1+t^2)
F8 = [F(n*n), F(0), F(m*m)]      # n^2 + m^2 t^2
Xt = [F(0), F(0), b]             # X = b t^2
def sub(u, v): return [a-c for a, c in zip(u, v)]
check('X - e1 = (mn)^2 * F4', sub(Xt, [e[0], F(0), F(0)]) == [F(m*n)**2*c for c in F4])
check('X - e2 = s m^2 * F0', sub(Xt, [e[1], F(0), F(0)]) == [s*m*m*c for c in F0])
check('X - e3 = s n^2 * F8', sub(Xt, [e[2], F(0), F(0)]) == [s*n*n*c for c in F8])
delta = [F(v) for v in cert['target_class']]
check('требуемый класс delta = (1, s, s) = (1,193,193)', delta == [F(1), F(193), F(193)])
check('s = 193 НЕ квадрат -> точек C над t=oo нет', not is_square(s))
check('F0,F4,F8 = (pos) + 0*t + (pos)*t^2 -> при любом рациональном t все F_i > 0',
      all(p[0] > 0 and p[1] == 0 and p[2] > 0 for p in (F0, F4, F8)))

# ================================================================= 2. модель M и (I,J)
shift = F(cert['shift']); u2 = F(cert['u2'])
check('shift = -(e1+e2+e3)/3', shift == -(e[0]+e[1]+e[2])/3)
check('u2 = 16 — квадрат (квадратный класс сохраняется)', is_square(u2) and u2 == 16)
r = [(ei + shift)/u2 for ei in e]
check('корни модели целые и совпали', r == [F(v) for v in cert['model_roots']]
      and all(x.denominator == 1 for x in r))
check('сумма корней модели = 0', sum(r) == 0)
A = r[0]*r[1] + r[0]*r[2] + r[1]*r[2]
B = -r[0]*r[1]*r[2]
check('A,B модели y^2=x^3+Ax+B', A == F(cert['model_A']) and B == F(cert['model_B']))
I, J = F(cert['I']), F(cert['J'])
check('I = c4 = -48A', I == -48*A)
check('J = 2c6 = -1728B', J == -1728*B)
phi = [F(v) for v in cert['phi_roots']]
check('phi_i = -12 r_i', phi == [-12*x for x in r])
check('phi_i — корни X^3-3IX+J', all(p**3 - 3*I*p + J == 0 for p in phi))
check('phi попарно различны', len(set(phi)) == 3)
DISC = F(cert['discriminant'])
check('Disc = 16(4I^3-J^2)/27 != 0', DISC == F(16*(4*I**3 - J**2), 27) and DISC != 0)

# ================================================================= 3. g1 — это 2-накрытие C_delta
D = [(r[i]-r[(i+1) % 3])*(r[i]-r[(i+2) % 3]) for i in range(3)]
check('D_i совпали', D == [F(v) for v in cert['covering_g1']['D']])
Lc2 = [F(1)/D[i] for i in range(3)]
Lc1 = [-(sum(r) - r[i])/D[i] for i in range(3)]
Lc0 = [(r[(i+1) % 3]*r[(i+2) % 3])/D[i] for i in range(3)]
# контроль интерполяции: Theta^k = sum_i r_i^k L_i(Theta)
for k in range(3):
    c2 = sum(r[i]**k*Lc2[i] for i in range(3))
    c1 = sum(r[i]**k*Lc1[i] for i in range(3))
    c0 = sum(r[i]**k*Lc0[i] for i in range(3))
    want = {0: (F(0), F(0), F(1)), 1: (F(0), F(1), F(0)), 2: (F(1), F(0), F(0))}[k]
    check('интерполяция Theta^%d' % k, (c2, c1, c0) == want)

cov = cert['covering_g1']
T = [[F(v) for v in q] for q in cov['T']]      # T_i = A x^2 + B xz + C z^2
lam = F(cov['lambda'])
G_raw = [F(v) for v in cov['G_raw']]
X_num = [F(v) for v in cov['X_num']]
gs = [[F(c) for c in q] for q in cert['quartics']]
g1, g2, g3 = gs
Tsq = [pmul(t, t) for t in T]
check('коника: [delta t^2]_2 = 0 тождественно',
      pzero(padd(*[pscal(delta[i]*Lc2[i], Tsq[i]) for i in range(3)])))
check('G_raw = -[delta t^2]_1',
      G_raw == pscal(-1, padd(*[pscal(delta[i]*Lc1[i], Tsq[i]) for i in range(3)])))
check('X_num = [delta t^2]_0',
      X_num == padd(*[pscal(delta[i]*Lc0[i], Tsq[i]) for i in range(3)]))
for i in range(3):
    check('накрытие: X_num - r_%d*G_raw = d_%d * T_%d^2' % (i+1, i+1, i+1),
          padd(X_num, pscal(-r[i], G_raw)) == pscal(delta[i], Tsq[i]))
check('G_raw = lambda^2 * g1, lambda рационально', G_raw == pscal(lam*lam, g1) and lam != 0)
check('d1 d2 d3 — квадрат', is_square(delta[0]*delta[1]*delta[2]))

# ================================================================= 4. три квартики, z и m
for k, g in enumerate(gs):
    check('I(g%d) = I' % (k+1), quart_I(g) == I)
    check('J(g%d) = J' % (k+1), quart_J(g) == J)
def z_of(g):
    a, bb, c = F(g[0]), F(g[1]), F(g[2])
    return [(4*a*p + 3*bb*bb - 8*a*c)/3 for p in phi]
zs = [z_of(g) for g in gs]
check('z(g_i) совпали с сертификатом', zs == [[F(v) for v in row] for row in cert['z_values']])
check('z(g1) — единица в L (все компоненты != 0)', all(v != 0 for v in zs[0]))
check('КЛАСС z(g1) = (1,193,193) = требуемый delta',
      all(is_square(zs[0][i]/delta[i]) for i in range(3)))
mv = [F(v) for v in cert['m_values']]
check('m^2 = z(g1)z(g2)z(g3) покомпонентно',
      all(mv[i]**2 == zs[0][i]*zs[1][i]*zs[2][i] for i in range(3)))

# ================================================================= 5. gamma по формуле Фишера
def H_of(g):
    h = hess(g); out = []
    for p in phi:
        G = [(4*p*F(g[j]) + h[j])/3 for j in range(5)]
        H = [G[0], G[1]/2, G[2]/6 + F(2, 9)*(I - p*p)]
        check_sq = pmul(H, H)
        if pmul(H, H) != pscal(G[0], G):
            raise AssertionError('ПРОВАЛ: тождество G(1,0)G = H^2')
        out.append(H)
    return out

def gamma_of(gA, zA, mvals):
    Hs = H_of(gA); out = []
    for j in range(3):
        out.append(sum((mvals[i]/zA[i])*Hs[i][j] /
                       ((phi[i]-phi[(i+1) % 3])*(phi[i]-phi[(i+2) % 3])) for i in range(3)))
    return out
OK.append('тождество Гессиана G(1,0)G = H^2 для всех трёх компонент L')

gam = gamma_of(g1, zs[0], mv)
check('gamma пересчитана и совпала', gam == [F(v) for v in cert['gamma']])
gam_n = [F(v) for v in cert['gamma_normalized']]
ratios = set(gam_n[i]/gam[i] for i in range(3) if gam[i] != 0)
check('нормированная gamma = рациональное кратное gamma', len(ratios) == 1)

a = F(cert['a'])
check('a = g2(1,0) != 0 (Fisher 3.2(iii))', a == F(g2[0]) and a != 0)

# ================================================================= 6. полнота набора мест
def needed_places(gA, gamma, aa):
    """Fisher Remark 3.3: вклад v тривиален, если N(v)>=11, g1 и gamma v-целые,
    v не делит Disc*content(gamma), a — v-единица, v не делит 2."""
    S = {2, 3, 5, 7} | support(DISC) | support(aa)
    for c in list(gA) + list(gamma):
        if F(c) != 0:
            S |= support(F(F(c).denominator))
    nz = [F(c) for c in gamma if F(c) != 0]
    cnum = reduce(gcd, [F(c).numerator for c in nz])
    cden = reduce(lambda x, y: x*y//gcd(x, y), [F(c).denominator for c in nz])
    cont = F(cnum, cden)
    S |= support(F(cont.numerator)) | support(F(cont.denominator))
    return S

def verify_run(gA, gamma, aa, run, tag):
    need = needed_places(gA, gamma, aa)
    have = {v['place'] for v in run}
    check(tag + ': набор мест покрывает Remark 3.3', {str(p) for p in need} | {'real'} <= have)
    check(tag + ': места не повторяются', len(have) == len(run))
    check(tag + ': вещественное место присутствует', 'real' in have)
    pr = 1; neg = []
    for v in run:
        x, z = F(v['x']), F(v['z'])
        check(tag + ' ' + v['place'] + ': (x:z) != (0:0)', x != 0 or z != 0)
        q = ev(gA, x, z); gv = ev(gamma, x, z)
        check(tag + ' ' + v['place'] + ': g(x,z) совпало', q == F(v['g']))
        check(tag + ' ' + v['place'] + ': gamma(x,z) совпало и != 0', gv == F(v['gamma']) and gv != 0)
        check(tag + ' ' + v['place'] + ': g(x,z) — ненулевой квадрат в Q_v',
              local_square(q, v['place']))
        hb = hilbert(aa, gv, v['place'])
        check(tag + ' ' + v['place'] + ': символ Гильберта совпал', hb == v['hilbert'])
        pr *= hb
        if hb == -1: neg.append(v['place'])
    return pr, neg

results = {}
for k, run in enumerate(cert['pair_runs']):
    pr, neg = verify_run(g1, gam, a, run, 'прогон %d' % k)
    results['forward_run_%d' % k] = {'product': pr, 'negative_places': neg}
    check('<g1,g2>_CT = -1 (прогон %d)' % k, pr == -1)

for k, run in enumerate(cert['pair_run_gamma_normalized']):
    pr, neg = verify_run(g1, gam_n, a, run, 'нормир. gamma %d' % k)
    results['normalized_run_%d' % k] = {'product': pr, 'negative_places': neg}
    check('нормировка gamma не меняет произведение', pr == -1)

# обратный порядок аргументов: gamma строится по g2, a = g1(1,0)
a_rev = F(cert['reverse']['a'])
check('a_rev = g1(1,0) != 0', a_rev == F(g1[0]) and a_rev != 0)
mv_rev = [F(v) for v in cert['reverse']['m_values']]
check('m для обратного порядка: m^2 = z(g2)z(g1)z(g3)',
      all(mv_rev[i]**2 == zs[1][i]*zs[0][i]*zs[2][i] for i in range(3)))
gam_rev = gamma_of(g2, zs[1], mv_rev)
check('gamma обратного порядка пересчитана и совпала',
      gam_rev == [F(v) for v in cert['reverse']['gamma']])
for k, run in enumerate(cert['reverse']['runs']):
    pr, neg = verify_run(g2, gam_rev, a_rev, run, 'обратный прогон %d' % k)
    results['reverse_run_%d' % k] = {'product': pr, 'negative_places': neg}
    check('<g2,g1>_CT = -1 (обратный прогон %d)' % k, pr == -1)
check('СИММЕТРИЯ: <g1,g2> = <g2,g1>', True)

# ================================================================= 7. ELS всех трёх квартик
els_need = {2, 3} | support(DISC)
for k, (g, loc) in enumerate(zip(gs, cert['local_solubility'])):
    extra = set()
    for c in g:
        if F(c) != 0: extra |= support(F(F(c).denominator))
    have = {v['place'] for v in loc}
    check('ELS g%d: места покрывают {2,3} u supp(Disc) u знаменатели + real' % (k+1),
          {str(p) for p in (els_need | extra)} | {'real'} <= have)
    for v in loc:
        x, z = F(v['x']), F(v['z'])
        q = ev(g, x, z)
        check('ELS g%d @ %s' % (k+1, v['place']),
              (x != 0 or z != 0) and q == F(v['value']) and q != 0 and local_square(q, v['place']))
# Обоснование пропуска остальных мест: при p нечётном, p не делит Disc, квартика p-целая
# с p-единичным содержанием (иначе p | Disc), редукция — гладкая кривая рода 1 над F_p,
# #C(F_p) >= p+1-2sqrt(p) > 0, гладкая точка поднимается по Гензелю (p нечётно).
check('все p, делящие Disc, а также 2 и 3, перечислены явно',
      els_need == {2, 3, 5, 7, 19, 193})

# ================================================================= 8. вывод
pairing = results['forward_run_0']['product']
res = {
  'pair': '(m,n)=(19,5)',
  'verifier': 'чистый Python 3 (fractions, math, json); Sage/PARI/код проекта не импортируются',
  'target_class': '(1,193,193)',
  'pairing_value': pairing,
  'reverse_pairing_value': results['reverse_run_0']['product'],
  'runs': results,
  'checks_passed': len(OK),
  'rank_used': False,
  'conclusion': ('<[g1],[g2]>_CT = -1 != 0  =>  класс (1,193,193) не лежит в образе '
                 'E(Q)/2E(Q)  =>  на C_(19,5) нет конечных рациональных точек; '
                 'бесконечность исключена тем, что s=193 не квадрат  =>  C_(19,5)(Q) = 0'),
}
open(os.path.join(here, 'ctp_cert_19_5_result.json'), 'w', encoding='utf-8').write(
    json.dumps(res, indent=2, ensure_ascii=False) + '\n')

print()
print('пройдено проверок:', len(OK))
for k, v in results.items():
    print('  %-22s произведение = %+d   места с -1: %s' % (k, v['product'], v['negative_places']))
print()
print('<[g1],[g2]>_CT =', pairing, '   <[g2],[g1]>_CT =', results['reverse_run_0']['product'])
print(res['conclusion'])
