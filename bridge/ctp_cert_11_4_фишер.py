"""Независимая проверка сертификата Касселса-Тейта для пары G1 (m,n) = (11,4).

СТАНДАРТНЫЙ Python. Не импортирует ни Sage, ни PARI, ни код проекта.
Читает bridge/ctp_11_4/certificate_11_4.json и заново проверяет ВСЁ:

  1.  минимальная модель M, её корни, I = -48A, J = -1728B, phi = -12 e;
  2.  инварианты I, J каждой из трёх бинарных квартик;
  3.  z(g1) лежит в требуемом классе (1, 274, 274) в (Q*/Q*^2)^3;
  4.  z(g1) z(g2) z(g3) = m^2 покомпонентно;
  5.  тождество Гессиана G(1,0) G = H^2 как многочленов;
  6.  gamma_1 = коэффициент при phi^2 в (m/z1) H_1, через интерполяцию Лагранжа;
  7.  всюду локальная разрешимость всех трёх квартик (явные свидетели);
  8.  полнота множества мест по Замечанию 3.3 Фишера (через факторизацию);
  9.  каждый символ Гильберта и произведение по местам = -1;
  10. связь с исходной кривой C рода 5: X = b t^2, X - e_i = (квадрат) * s^eps,
      то есть необходимый класс подъёма действительно (1, s, s) = (1, 274, 274).

ГРАНИЦА РАНГА E(Q) НЕ ИСПОЛЬЗУЕТСЯ НИГДЕ.

Запуск:  python3 ctp_cert_11_4_фишер.py
Вывод:   bridge/ctp_11_4/independent_result_11_4.json
"""
from fractions import Fraction as Q
from math import isqrt, prod, gcd
from functools import reduce
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
CERT = HERE / 'ctp_11_4' / 'certificate_11_4.json'
r = json.loads(CERT.read_text())

# ---------------------------------------------------------------- элементарная арифметика


def square(a):
    a = Q(a)
    return a >= 0 and isqrt(a.numerator) ** 2 == a.numerator and isqrt(a.denominator) ** 2 == a.denominator


def fac(n):
    n = abs(int(n))
    ans = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            ans[p] = ans.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        ans[n] = ans.get(n, 0) + 1
    return ans


def support(a):
    a = Q(a)
    return set(fac(a.numerator)) | set(fac(a.denominator)) if a else set()


def unit(a, p):
    """a = p^k * u, u — p-единица; возвращает (k, u)."""
    a = Q(a)
    u, v = a.numerator, a.denominator
    k = 0
    assert u
    while u % p == 0:
        u //= p
        k += 1
    while v % p == 0:
        v //= p
        k -= 1
    return k, Q(u, v)


def residue(a, p):
    return a.numerator * pow(a.denominator, -1, p) % p


def local_square(a, pl):
    """квадрат ли a в Q_v (v = pl)."""
    a = Q(a)
    if not a:
        return True
    if pl == 'real':
        return a > 0
    p = int(pl)
    k, u = unit(a, p)
    if k % 2:
        return False
    return residue(u, 8) == 1 if p == 2 else pow(residue(u, p), (p - 1) // 2, p) == 1


def hilbert(a, b, pl):
    """символ Гильберта (a,b)_v, точно."""
    if pl == 'real':
        return -1 if a < 0 and b < 0 else 1
    p = int(pl)
    v, u = unit(a, p)
    w, t = unit(b, p)
    if p == 2:
        uu, tt = residue(u, 8), residue(t, 8)
        e = ((uu - 1) // 2) * ((tt - 1) // 2) + v * ((tt * tt - 1) // 8) + w * ((uu * uu - 1) // 8)
        return -1 if e % 2 else 1
    lu = pow(residue(u, p), (p - 1) // 2, p)
    lt = pow(residue(t, p), (p - 1) // 2, p)
    return ((-1 if (v * w * ((p - 1) // 2)) % 2 else 1)
            * (-1 if lu == p - 1 and w % 2 else 1)
            * (-1 if lt == p - 1 and v % 2 else 1))


# ---------------------------------------------------------------- инварианты бинарных квартик


def ij(g):
    a, b, c, d, e = g
    return (12 * a * e - 3 * b * d + c * c,
            72 * a * c * e - 27 * a * d * d - 27 * b * b * e + 9 * b * c * d - 2 * c ** 3)


def ev(g, x, z):
    n = len(g) - 1
    return sum(c * x ** (n - i) * z ** i for i, c in enumerate(g))


def hess(g):
    a, b, c, d, e = g
    return [3 * b * b - 8 * a * c, 4 * (b * c - 6 * a * d),
            2 * (2 * c * c - 24 * a * e - 3 * b * d), 4 * (c * d - 6 * b * e), 3 * d * d - 8 * c * e]


def zval(g, ph):
    a, b, c, d, e = g
    return (4 * a * ph + 3 * b * b - 8 * a * c) / 3


REPORT = []


def check(cond, msg):
    assert cond, 'ПРОВАЛ: ' + msg
    REPORT.append(msg)


# ---------------------------------------------------------------- 1. модель и (I,J)

I, J = Q(r['I']), Q(r['J'])
gs = [list(map(Q, g)) for g in r['quartics']]
phis = list(map(Q, r['phi_roots']))
roots = list(map(Q, r['minimal_roots']))

A, B = Q(r['M_ainvs'][3]), Q(r['M_ainvs'][4])
check([Q(c) for c in r['M_ainvs'][:3]] == [0, 0, 0], 'модель M имеет вид y^2 = x^3 + A x + B')
check(A == -4644297081900 and B == 3655579703416830000, 'M: y^2 = x^3 - 4644297081900 x + 3655579703416830000')
check(len(set(roots)) == 3 and all(e ** 3 + A * e + B == 0 for e in roots), 'три различных рациональных корня M')
check(sum(roots) == 0, 'сумма корней M равна нулю')
check(I == -48 * A and J == -1728 * B, 'I = c4(M) = -48A,  J = 2 c6(M) = -1728B')
check(phis == [-12 * e for e in roots], 'phi_i = -(12 e_i + b2),  b2 = 0')
check(all(ph ** 3 - 3 * I * ph + J == 0 for ph in phis), 'phi_i — корни X^3 - 3IX + J')
check(len(set(phis)) == 3, 'phi_i попарно различны (L = Q x Q x Q)')

# ---------------------------------------------------------------- 2. инварианты трёх квартик

for i, g in enumerate(gs):
    check(ij(g) == (I, J), f'квартика g{i+1} имеет инварианты ровно (I, J)')
disc = Q(16, 27) * (4 * I ** 3 - J ** 2)
check(disc != 0 and disc == Q(r['discriminant']), 'Delta = 16(4I^3 - J^2)/27 != 0, совпала с сертификатом')

# ---------------------------------------------------------------- 3. класс z(g1)

zs = [[zval(g, ph) for ph in phis] for g in gs]
check(zs == [list(map(Q, z)) for z in r['z_values']], 'значения z(g_i) в трёх компонентах L пересчитаны')
TRIP = [Q(t) for t in r['target_squareclass']]
check(TRIP == [1, 274, 274], 'целевой класс в сертификате — (1, 274, 274)')
check(all(z != 0 for z in zs[0]), 'z(g1) — единица в L (не делитель нуля)')
check(all(square(z / d) for z, d in zip(zs[0], TRIP)),
      'z(g1) лежит РОВНО в классе (1, 274, 274) в (Q*/Q*^2)^3')

# ---------------------------------------------------------------- 4. корень m

ms = list(map(Q, r['m_values']))
check(all(m * m == zs[0][i] * zs[1][i] * zs[2][i] for i, m in enumerate(ms)),
      'z(g1) z(g2) z(g3) = m^2 в L (сумма трёх классов равна нулю)')

# ---------------------------------------------------------------- 5-6. Гессиан и gamma_1

Hs = []
for ph in phis:
    h = hess(gs[0])
    G = [(4 * ph * c + d) / 3 for c, d in zip(gs[0], h)]
    H = [G[0], G[1] / 2, G[2] / 6 + Q(2, 9) * (I - ph * ph)]
    Hsq = [sum(H[i] * H[j] for i in range(3) for j in range(3) if i + j == k) for k in range(5)]
    check(Hsq == [G[0] * v for v in G], f'тождество G(1,0) G = H^2 в компоненте phi = {ph}')
    Hs.append(H)

gam = []
for j in range(3):
    gam.append(sum((ms[i] / zs[0][i]) * Hs[i][j] / prod(phis[i] - phis[k] for k in range(3) if i != k)
                   for i in range(3)))
check(gam == list(map(Q, r['gamma_raw'])),
      'gamma_1 = коэффициент при phi^2 в (z(g2)z(g3)/m) H_1, интерполяция Лагранжа')
scale = Q(r['gamma_scale'])
gam = [g * scale for g in gam]
check(gam == list(map(Q, r['gamma'])), 'глобальная нормировка gamma_1 (не меняет произведение по формуле взаимности)')
check(any(c != 0 for c in gam), 'gamma_1 не тождественный ноль')

a = Q(r['a'])
check(a == gs[1][0] and a != 0, 'a = g2(1,0) буквально из теоремы 3.1 (без замен и сокращений)')

# ---------------------------------------------------------------- 7-9. места и символы Гильберта
# Замечание 3.3 Фишера: вклад v тривиален, если N(v) >= 11, g1 и gamma_1 v-целые,
# v не делит Delta(g1)*content(gamma_1), a — v-единица и v не делит 2.
# Дополнение к этому набору обязано целиком входить в набор мест прогона.

need = {2, 3, 5, 7} | support(disc) | support(a)
for c in gs[0] + gam:
    need |= support(Q(c.denominator))
content = reduce(gcd, [abs(c.numerator) for c in gam if c])
need |= support(Q(content))
check(need == {2, 3, 5, 7, 11, 137, 2531},
      f'обязательный набор мест по Замечанию 3.3: {sorted(need)} + real')

products = []
for k, run in enumerate(r['pair_runs']):
    places = {v['place'] for v in run}
    check({str(p) for p in need} | {'real'} <= places, f'прогон {k}: все обязательные места присутствуют')
    check(len(places) == len(run), f'прогон {k}: места не повторяются')
    for v in run:
        x, z = Q(v['x']), Q(v['z'])
        check(x != 0 or z != 0, f'прогон {k}, место {v["place"]}: точка (x:z) корректна')
        q = ev(gs[0], x, z)
        gv = ev(gam, x, z)
        check(q == Q(v['g']) and gv == Q(v['gamma']), f'прогон {k}, место {v["place"]}: значения пересчитаны')
        check(gv != 0, f'прогон {k}, место {v["place"]}: gamma_1(x_v,z_v) != 0')
        check(local_square(q, v['place']), f'прогон {k}, место {v["place"]}: g1(x_v,z_v) — квадрат в Q_v')
        check(hilbert(a, gv, v['place']) == v['hilbert'], f'прогон {k}, место {v["place"]}: символ Гильберта пересчитан')
    pr = prod(v['hilbert'] for v in run)
    products.append(pr)
    check(pr == -1, f'прогон {k}: произведение символов Гильберта = -1')

minus_places = [[v['place'] for v in run if v['hilbert'] == -1] for run in r['pair_runs']]
check(all(mp == minus_places[0] for mp in minus_places), 'набор мест с -1 одинаков во всех прогонах')

# ---------------------------------------------------------------- 7'. ELS всех трёх квартик

els_places = []
for i, (g, local) in enumerate(zip(gs, r['local_solubility'])):
    ps = {2, 3} | support(disc)
    for c in g:
        ps |= support(Q(c.denominator))
    have = {v['place'] for v in local}
    check({str(p) for p in ps} | {'real'} <= have, f'g{i+1}: проверены все места плохой редукции {sorted(ps)} и real')
    for v in local:
        x, z = Q(v['x']), Q(v['z'])
        check(x != 0 or z != 0, f'g{i+1}, место {v["place"]}: свидетель корректен')
        q = ev(g, x, z)
        check(q == Q(v['value']) and q != 0, f'g{i+1}, место {v["place"]}: значение пересчитано')
        check(local_square(q, v['place']), f'g{i+1}, место {v["place"]}: y^2 = g имеет Q_v-точку')
    els_places.append(sorted(have))
    # В остальных местах p: g p-целая, content(g) и Delta(g) — p-единицы (иначе p попало бы в ps),
    # значит гладкая кривая рода 1 y^2 = g(x,z) имеет ХОРОШУЮ редукцию в p, p >= 5.
    # Оценка Хассе #C(F_p) >= p + 1 - 2 sqrt(p) > 0 даёт гладкую F_p-точку, Гензель её поднимает.

# ---------------------------------------------------------------- 10. связь с исходной кривой C

m, n = 11, 4
s = Q(137, 2)
b = s * m * m * n * n
check(b == 132616, 'b = s m^2 n^2 = 132616')
er = [-b, -s * m ** 4, -s * n ** 4]
check(er == [Q(-132616), Q(-2005817, 2), Q(-17536)], 'корни E: e1 = -b, e2 = -s m^4, e3 = -s n^4')
shift = roots[0] - 4 * er[0]
check(shift == 1537414, 'сдвиг минимальной модели = 1537414')
check(all(4 * e + shift == f for e, f in zip(er, roots)),
      'x_M = 4X + 1537414 переводит корни E в корни M В ТОМ ЖЕ ПОРЯДКЕ')
check(4 ** 3 == 8 ** 2, 'y_M = 8V согласовано с x_M = 4X + shift (изоморфизм, u = 1/2)')
check(square(Q(4)) and square(Q(36)), 'масштабы 4 (E->M) и 36 (M->E_{I,J}) — квадраты: класс delta сохраняется')
check(all(-3 * ph == 36 * e for ph, e in zip(phis, roots)),
      'Theta_i = -3 phi_i = 36 e_i: x-координаты 2-кручения E_{I,J}')

# F0 = m^2 + n^2 t^2,  F4 = s(1 + t^2),  F8 = n^2 + m^2 t^2;  X = b t^2
# [-e_i, b] — коэффициенты X - e_i как многочлена от t^2 (свободный, при t^2).
FF = [[Q(m * m), Q(n * n)], [s, s], [Q(n * n), Q(m * m)]]
for i, c, j, name in [(0, Q(m * m * n * n), 1, 'X - e1 = (mn)^2 u4^2'),
                      (1, s * m * m, 0, 'X - e2 = s m^2 u0^2'),
                      (2, s * n * n, 2, 'X - e3 = s n^2 u8^2')]:
    check([-er[i], b] == [c * v for v in FF[j]], name)
check(square(s / Q(274)) and not square(s), 's = 137/2 имеет класс 274 и НЕ является квадратом')
check(not square(Q(274)), '274 не квадрат: точка C над бесконечностью потребовала бы (u4/t)^2 = s — невозможно')

# ---------------------------------------------------------------- итог

fin_all_plus = all(v['hilbert'] == 1 for v in r['pair_runs'][0] if v['place'] != 'real')
res = {
    'pair': '(m,n) = (11,4)',
    'verifier': 'стандартный Python (fractions/math), без Sage, PARI и кода проекта',
    's': '137/2', 'b': 132616,
    'M': 'y^2 = x^3 - 4644297081900 x + 3655579703416830000',
    'I': str(I), 'J': str(J),
    'target_squareclass': [1, 274, 274],
    'z_g1_squareclass_verified': True,
    'invariants_verified': True,
    'Hessian_identities_verified': True,
    'Fisher_gamma_verified': True,
    'three_quartics_ELS_verified': True,
    'ELS_places_checked': els_places,
    'required_places_Remark_3_3': sorted(str(p) for p in need) + ['real'],
    'places_used': [v['place'] for v in r['pair_runs'][0]],
    'hilbert_symbols': {v['place']: v['hilbert'] for v in r['pair_runs'][0]},
    'places_with_minus_one': minus_places[0],
    'finite_hilbert_symbols_all_plus_one': fin_all_plus,
    'hilbert_products_per_run': products,
    'pairing_value': -1,
    'pairing_in_Q_mod_Z': '1/2',
    'reverse_pair_g2_g1': r.get('reverse_pair'),
    'pair_g1_g3': r.get('pair_g1_g3'),
    'self_pair_g1_g1': r.get('self_pair_g1_g1'),
    'control_image_of_E_Q_pairs_trivially': r.get('control_image_pairs'),
    'rank_bound_required': False,
    'checks_passed': len(REPORT),
    'conclusion_using_Fisher_Theorem_3_1':
        '<[g1],[g2]>_CT = -1 != 0, поэтому [g1] не лежит в образе E(Q)/2E(Q); '
        'значит класс (1,274,274) не достигается, конечных рациональных точек на C_(11,4) нет; '
        'бесконечность исключена тем, что s = 137/2 не квадрат. C_(11,4)(Q) = пусто.',
}
(HERE / 'ctp_11_4' / 'independent_result_11_4.json').write_text(json.dumps(res, indent=2, ensure_ascii=False) + '\n')

print(f'ПРОВЕРОК ПРОЙДЕНО: {len(REPORT)}')
for line in REPORT:
    print('  OK  ' + line)
print()
print(json.dumps(res, indent=2, ensure_ascii=False))
