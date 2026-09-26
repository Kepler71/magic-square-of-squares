#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lift_mass_direct.py — проверка массового утверждения §5 REVIEW_LIFT_CLASSES_2026-09-26:
  «у всех 2760 наклонов 0<r<s<=300 с нетривиальными списками Codex после фильтра по простым
   p = 3 (mod 4), 7<=p<=400, остаётся только (d_T, d_L) = (1, 1)».

Claude (подагент), 26.09.2026. Скрипт lift_classes/review_20260926/p7_filter_check.py НЕ открывался.
Черновик lm_common.py в этом каталоге (чужой, 19:54) не импортируется и не используется.

Способ (отличается от формулы I_p сводки): для каждого вычета z0 = z mod p напрямую вычисляются
классы ВСЕХ ДЕВЯТИ клеток f_ij = 1 + (i r + j s) z в Q_p*/Q_p*^2, затем проверяются восемь линий
АП-сетки, затем берутся произведения [T], [L] по клеткам. Обоснование (см. NOTE.md, §1):
  (a) при p = 3 mod 4 полюса нет: z в Z_p;
  (b) клетка с 1 + mu z0 != 0 mod p — единица, её класс (0, chi_p(1 + mu z0));
      клетка с 1 + mu z0 = 0 mod p: противоположная клетка 1 - mu z = 2 - (1 + mu z) — единица с вычетом 2,
      а центральная строка/столбец/диагонали (центр = 1) требуют равенства классов противоположных клеток,
      поэтому в любой конфигурации, где восемь произведений — квадраты, класс этой клетки равен (0, chi_p(2)).
  Итог: множество IMG_p(r, s) пар ([T], [L]), полученных так по всем z0, — ВЕРХНЯЯ оценка локального образа.
  Глобальный класс d_T (все простые = 1 mod 4, т.е. p не делит d_T) имеет в Q_p класс (0, chi_p(d_T)).
  Пара (d_T, d_L) исключается простым p, если (chi_p(d_T), chi_p(d_L)) не лежит в IMG_p(r, s).
Второй независимый код (метод B): клетка с нулевым вычетом получает ЛЮБОЙ из 4 классов Q_p*/Q_p*^2,
отбор только по восьми линиям. Образы A и B обязаны совпасть (проверяется для каждого вычисленного образа).

Все циклы имеют явные верхние границы; прогресс пишется в lift_mass_direct.log.
"""
import os, sys, time, json, random
from math import gcd
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
LOGF = open(os.path.join(HERE, 'lift_mass_direct.log'), 'w', encoding='utf-8')
T0 = time.time()


def log(*a):
    line = f'[{time.time() - T0:7.1f}s] ' + ' '.join(str(x) for x in a)
    print(line, flush=True)
    LOGF.write(line + '\n'); LOGF.flush()


S_MAX = 300
P_LO, P_HI = 7, 400
FAIL = []          # сюда пишутся все нарушенные внутренние проверки


def check(cond, msg):
    if not cond:
        FAIL.append(msg); log('!!! НАРУШЕНИЕ:', msg)


# ---------------------------------------------------------------- арифметика
def is_prime(n):
    if n < 2: return False
    d = 2
    while d * d <= n:                       # граница: d <= sqrt(n)
        if n % d == 0: return False
        d += 1
    return True


def factor(n):
    n = abs(n); f = {}; d = 2
    while d * d <= n:                       # граница: d <= sqrt(n) <= 25 при n <= 600
        while n % d == 0:
            f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f


def D_list(n):
    """Список Codex: d>0 бесквадратное, p|d => p|n, p = 1 mod 4, d = 1 mod 24 (lift/RESULT.md, теорема)."""
    ps = sorted(q for q in factor(n) if q % 4 == 1)
    out = []
    for k in range(len(ps) + 1):
        for comb in combinations(ps, k):
            d = 1
            for q in comb: d *= q
            if d % 24 == 1: out.append(d)
    return sorted(out)


def legendre_bit(a, p):
    """0, если a — ненулевой квадрат mod p, 1 — невычет. a не делится на p."""
    a %= p
    assert a != 0
    return 0 if pow(a, (p - 1) // 2, p) == 1 else 1


def chi_table(p):
    sq = {x * x % p for x in range(1, p)}
    t = [None] + [0 if x in sq else 1 for x in range(1, p)]
    # сверка с критерием Эйлера
    for x in range(1, p):
        assert t[x] == legendre_bit(x, p)
    return t


def vsplit(n, p):
    v = 0
    while n % p == 0:                        # граница: v <= log_p |n|
        n //= p; v += 1
    return v, n


def cls_exact(num, den, p):
    """Точный класс ненулевого рационального num/den в Q_p*/Q_p*^2 (p нечётно): код 2*(v mod 2) + chi(ед. часть)."""
    assert num != 0 and den != 0
    a, nu = vsplit(num, p); b, de = vsplit(den, p)
    u = nu * pow(de, -1, p) % p
    return 2 * ((a - b) & 1) + legendre_bit(u, p)


PRIMES = [p for p in range(P_LO, P_HI + 1) if is_prime(p) and p % 4 == 3]

# ---------------------------------------------------------------- сетка
CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
K = {c: n for n, c in enumerate(CELLS)}
OPP = [K[(-i, -j)] for (i, j) in CELLS]
ROWS = [tuple(K[(i, j)] for j in (-1, 0, 1)) for i in (-1, 0, 1)]
COLS = [tuple(K[(i, j)] for i in (-1, 0, 1)) for j in (-1, 0, 1)]
DIAG = [(K[(-1, -1)], K[(0, 0)], K[(1, 1)]), (K[(-1, 1)], K[(0, 0)], K[(1, -1)])]
LINES8 = ROWS + COLS + DIAG
LINES6 = ROWS + COLS                       # только для отрицательного контроля
T_IDX = (K[(-1, 0)], K[(1, 1)], K[(0, -1)])   # T = (1 - r z)(1 + (r+s) z)(1 - s z)
L_IDX = (K[(-1, 0)], K[(1, -1)], K[(0, 1)])   # L = (1 - r z)(1 + (r-s) z)(1 + s z)


def mus(r, s):
    return [i * r + j * s for (i, j) in CELLS]


def selftest_grid():
    from fractions import Fraction as F
    rng = random.Random(1)
    for _ in range(200):
        r = rng.randint(-50, 50); s = rng.randint(-50, 50); z = F(rng.randint(-99, 99), rng.randint(1, 99))
        cell = [1 + m * z for m in mus(r, s)]
        T = (1 - r * z) * (1 + (r + s) * z) * (1 - s * z)
        L = (1 - r * z) * (1 + (r - s) * z) * (1 + s * z)
        check(cell[T_IDX[0]] * cell[T_IDX[1]] * cell[T_IDX[2]] == T, 'T по клеткам')
        check(cell[L_IDX[0]] * cell[L_IDX[1]] * cell[L_IDX[2]] == L, 'L по клеткам')
        check(cell[K[(0, 0)]] == 1, 'центр = 1')
        for k in range(9):
            check(cell[k] + cell[OPP[k]] == 2, 'сумма противоположных = 2')
    # в каждой линии — три разные клетки; всего 8 линий; через центр проходят ровно 4
    check(len(LINES8) == 8 and all(len(set(l)) == 3 for l in LINES8), '8 линий')
    check(sum(K[(0, 0)] in l for l in LINES8) == 4, '4 центральные линии')


# ---------------------------------------------------------------- метод A: прямые классы девяти клеток
def residue_configs(r0, s0, p, chi, lines=LINES8):
    """Для каждого z0 в F_p: (z0, классы 9 клеток, есть ли нулевая клетка, линии выполнены).
    Код класса 2e + c; здесь всегда e = 0 (см. обоснование в шапке)."""
    mu = [m % p for m in mus(r0, s0)]
    out = []
    for z0 in range(p):                      # граница: p
        vals = [(1 + m * z0) % p for m in mu]
        cl = [0] * 9; hz = False
        for k in range(9):
            v = vals[k]
            if v:
                cl[k] = chi[v]
            else:
                hz = True
                vo = vals[OPP[k]]
                if vo != 2 % p:
                    FAIL.append(f'противоположная к нулевой клетке не = 2 mod p: {(r0, s0, p, z0, k)}')
                cl[k] = chi[vo]
        ok = all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in lines)
        out.append((z0, cl, hz, ok))
    return out


def image_A(r0, s0, p, chi, lines=LINES8, both=False):
    """Образ ([T],[L]) по всем z0; при both=True ещё и образ только по z0 без нулевых клеток."""
    img = set(); img_nz = set()
    for z0, cl, hz, ok in residue_configs(r0, s0, p, chi, lines):
        if ok:
            t = cl[T_IDX[0]] ^ cl[T_IDX[1]] ^ cl[T_IDX[2]]
            l = cl[L_IDX[0]] ^ cl[L_IDX[1]] ^ cl[L_IDX[2]]
            img.add((t, l))
            if not hz: img_nz.add((t, l))
    return (img, img_nz) if both else img


# ---------------------------------------------------------------- метод B: нулевые клетки — любой из 4 классов
def image_B(r0, s0, p, chi):
    mu = [m % p for m in mus(r0, s0)]
    img = set()
    for z0 in range(p):                      # граница: p
        known = [0] * 9; unk = []
        for k in range(9):
            v = (1 + mu[k] * z0) % p
            if v: known[k] = chi[v]
            else: unk.append(k)
        for combo in product(range(4), repeat=len(unk)):   # граница: 4^(<=3)
            cl = known[:]
            for k, c in zip(unk, combo): cl[k] = c
            if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8):
                t = cl[T_IDX[0]] ^ cl[T_IDX[1]] ^ cl[T_IDX[2]]
                l = cl[L_IDX[0]] ^ cl[L_IDX[1]] ^ cl[L_IDX[2]]
                img.add((t, l))
    return img


CHI = {}
IMG_CACHE = {}
STATS = dict(images=0, zero_cell_adds=0)


def get_img(r, s, p):
    key = (r % p, s % p, p)
    if key not in IMG_CACHE:
        chi = CHI.setdefault(p, chi_table(p))
        a, a_nz = image_A(r, s, p, chi, both=True)
        b = image_B(r, s, p, chi)
        # в методе B классы T, L — полные коды; оценка должна быть чётной
        check(all(t < 2 and l < 2 for (t, l) in b), f'нечётная оценка T/L в методе B {key}')
        check(a == b, f'методы A и B разошлись {key}: {a} vs {b}')
        check((0, 0) in a, f'(0,0) не в образе {key}')
        if a_nz != a: STATS['zero_cell_adds'] += 1
        STATS['images'] += 1
        IMG_CACHE[key] = a
    return IMG_CACHE[key]


# ================================================================= основная часть
def main():
    log('Проверка массового утверждения §5; простые p = 3 mod 4 в [7,400]:', len(PRIMES), 'шт.,', PRIMES[0], '...', PRIMES[-1])
    selftest_grid()
    log('самопроверка сетки: нарушений', len(FAIL))

    slopes = [(r, s) for s in range(2, S_MAX + 1) for r in range(1, s) if gcd(r, s) == 1]
    nontriv = []
    for (r, s) in slopes:
        Dm, Dp = D_list(r - s), D_list(r + s)
        if Dm != [1] or Dp != [1]:
            nontriv.append((r, s, Dm, Dp))
    log(f'наклонов 0<r<s<={S_MAX}, gcd=1: {len(slopes)}; с нетривиальными списками: {len(nontriv)}')
    npairs = sum(len(Dm) * len(Dp) - 1 for (r, s, Dm, Dp) in nontriv)
    log('всего нетривиальных пар (d_T, d_L):', npairs,
        '; max |D_-|, |D_+|:', max(len(x[2]) for x in nontriv), max(len(x[3]) for x in nontriv))

    records = []
    survivors_total = 0
    p7_insufficient = []
    needed_hist = {}
    for idx, (r, s, Dm, Dp) in enumerate(nontriv):          # граница: 2760 (или сколько получится)
        pairs = [(a, b) for a in Dm for b in Dp if (a, b) != (1, 1)]
        excl = {pr: [] for pr in pairs}
        for p in PRIMES:                                     # граница: len(PRIMES)
            img = get_img(r, s, p)
            for (a, b) in pairs:
                if (legendre_bit(a, p), legendre_bit(b, p)) not in img:
                    excl[(a, b)].append(p)
        surv = [pr for pr in pairs if not excl[pr]]
        survivors_total += len(surv)
        first = {pr: (excl[pr][0] if excl[pr] else None) for pr in pairs}
        needed = None if surv else max(first.values())
        needed_hist[needed] = needed_hist.get(needed, 0) + 1
        if any(7 not in excl[pr] for pr in pairs):
            p7_insufficient.append((r, s))
        records.append(dict(r=r, s=s, D_minus=Dm, D_plus=Dp,
                            pairs=[dict(dT=a, dL=b, first_excluding_p=first[(a, b)],
                                        n_excluding_p=len(excl[(a, b)])) for (a, b) in pairs],
                            survivors=surv, needed_p=needed))
        if (idx + 1) % 250 == 0 or idx + 1 == len(nontriv):
            log(f'  наклонов обработано {idx + 1}/{len(nontriv)}; выживших пар пока {survivors_total}; '
                f'образов вычислено {STATS["images"]}; нарушений {len(FAIL)}')

    log('ИТОГ ФИЛЬТРА: выживших нетривиальных пар:', survivors_total,
        '; наклонов, где остаётся что-то кроме (1,1):', sum(1 for x in records if x['survivors']))
    log('наклонов, где одного p = 7 не хватает:', len(p7_insufficient))
    log('распределение «нужного» простого (max по парам от наименьшего исключающего p):',
        sorted(needed_hist.items(), key=lambda kv: (kv[0] is None, kv[0] or 0)))
    worst = [x for x in records if x['needed_p'] == max(v for v in needed_hist if v is not None)]
    log('наклоны с наибольшим нужным p:', [(x['r'], x['s'], x['D_minus'], x['D_plus']) for x in worst[:10]])
    minexcl = min(pp['n_excluding_p'] for x in records for pp in x['pairs'])
    log('наименьшее число исключающих простых у одной пары:', minexcl)
    log('образов, в которые клетки с нулевым вычетом добавили новые пары:', STATS['zero_cell_adds'], 'из', STATS['images'])

    # ---- именованные наклоны (126/451 и 73/362 вне диапазона s<=300)
    for (r, s) in [(126, 451), (73, 362), (265, 298)]:
        Dm, Dp = D_list(r - s), D_list(r + s)
        pairs = [(a, b) for a in Dm for b in Dp if (a, b) != (1, 1)]
        res = {}
        for (a, b) in pairs:
            res[(a, b)] = [p for p in PRIMES if (legendre_bit(a, p), legendre_bit(b, p)) not in get_img(r, s, p)][:6]
        log(f'  наклон {r}/{s}: D_-={Dm}, D_+={Dp}, исключающие простые (первые 6):', res)

    # ---- отрицательные контроли фильтра
    # (1) только p = 7 (должно не хватать — в сводке 504)
    # (2) роли T и L переставлены (d_T из списка r+s, d_L из списка r-s)
    swapped_surv = 0; swapped_slopes = 0
    for (r, s, Dm, Dp) in nontriv:
        pairs = [(a, b) for a in Dp for b in Dm if (a, b) != (1, 1)]
        ss = 0
        for (a, b) in pairs:
            if all((legendre_bit(a, p), legendre_bit(b, p)) in get_img(r, s, p) for p in PRIMES):
                ss += 1
        swapped_surv += ss; swapped_slopes += (ss > 0)
    log('КОНТРОЛЬ «T и L переставлены»: выживших пар', swapped_surv, 'на', swapped_slopes, 'наклонах')
    # (3) без диагоналей (6 линий): образ шире — сколько пар выживает
    chi_cache = CHI
    six_surv = 0; six_slopes = 0; img6 = {}
    for (r, s, Dm, Dp) in nontriv:
        pairs = [(a, b) for a in Dm for b in Dp if (a, b) != (1, 1)]
        ss = 0
        for (a, b) in pairs:
            alive = True
            for p in PRIMES:
                key = (r % p, s % p, p)
                if key not in img6:
                    img6[key] = image_A(r, s, p, chi_cache.setdefault(p, chi_table(p)), lines=LINES6)
                if (legendre_bit(a, p), legendre_bit(b, p)) not in img6[key]:
                    alive = False; break
            ss += alive
        six_surv += ss; six_slopes += (ss > 0)
    log('КОНТРОЛЬ «без двух диагоналей (6 линий)»: выживших пар', six_surv, 'на', six_slopes, 'наклонах')

    out = dict(S_MAX=S_MAX, primes=PRIMES, n_slopes=len(slopes), n_nontrivial=len(nontriv),
               n_nontrivial_pairs=npairs, survivors_total=survivors_total,
               p7_insufficient=p7_insufficient, needed_hist={str(k): v for k, v in needed_hist.items()},
               swapped_survivors=swapped_surv, six_line_survivors=six_surv,
               zero_cell_adds=STATS['zero_cell_adds'], images=STATS['images'],
               records=records)
    with open(os.path.join(HERE, 'lift_mass_direct.json'), 'w', encoding='utf-8') as f:
        json.dump(out, f, ensure_ascii=False)
    log('записано lift_mass_direct.json')
    return nontriv


# ================================================================= контроли
def union_residue_stats(p):
    """По всем (A0, B0) в F_p^2 (A = r z, B = s z): число вычетов без нулевых клеток с выполненными 8 линиями,
    из них — с нетривиальным классом хотя бы одной клетки; объединённый образ ([T],[L])."""
    chi = CHI.setdefault(p, chi_table(p))
    n_ok_nz = 0; n_nontriv_nz = 0; union_img = set(); any_cell_nontriv = 0
    for A0 in range(p):
        for B0 in range(p):
            # z0 = 1, (r0, s0) = (A0, B0): одна конфигурация
            vals = [(1 + i * A0 + j * B0) % p for (i, j) in CELLS]
            cl = [0] * 9; hz = False
            for k in range(9):
                if vals[k]: cl[k] = chi[vals[k]]
                else: hz = True; cl[k] = chi[vals[OPP[k]]]
            if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8):
                union_img.add((cl[T_IDX[0]] ^ cl[T_IDX[1]] ^ cl[T_IDX[2]], cl[L_IDX[0]] ^ cl[L_IDX[1]] ^ cl[L_IDX[2]]))
                if any(cl): any_cell_nontriv += 1
                if not hz:
                    n_ok_nz += 1
                    n_nontriv_nz += any(cl)
    return n_ok_nz, n_nontriv_nz, union_img, any_cell_nontriv


def exact_grid_count(p, M, poles=()):
    """Точный перебор A, B в [0, p^M) (целый z) и полюсов A = a/p^n, B = b/p^n (a, b в [0, p^M), не оба = 0 mod p).
    Класс клетки — точно, по целому числителю; нулевая клетка или числитель = 0 mod p^M — «не определено» (пропуск).
    Возвращает (eight_ok, из них не все 9 квадраты) для целого случая и для каждого n из poles."""
    Q = p ** M
    res = {}
    for n in (0,) + tuple(poles):
        sh = p ** n
        eight = 0; bad = 0; undet = 0
        tab = {}
        for A in range(Q):                     # граница: p^M
            for B in range(Q):                 # граница: p^M
                if n and A % p == 0 and B % p == 0: continue
                cl = []
                skip = False
                for (i, j) in CELLS:
                    num = sh + i * A + j * B
                    if num % Q == 0: skip = True; break
                    c = tab.get(num)
                    if c is None:
                        c = cls_exact(num, sh, p); tab[num] = c
                    cl.append(c)
                if skip: undet += 1; continue
                if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8):
                    eight += 1
                    if any(cl): bad += 1
        res[n] = (eight, bad, undet)
    return res


def controls(nontriv):
    log('=== КОНТРОЛЬ C1: p = 7, автоматический подъём (метод A по всем (A0,B0) в F_7^2) ===')
    n_ok, n_bad, uimg, anyc = union_residue_stats(7)
    log(f'  p=7: вычетов без нулевых клеток с 8 линиями {n_ok}, из них нетривиальных {n_bad}; '
        f'объединённый образ {sorted(uimg)}; конфигураций (вкл. нулевые клетки) с нетривиальной клеткой: {anyc}')
    check(uimg == {(0, 0)} and anyc == 0, 'p=7: подъём не автоматический')
    # объединение образов по всем (r0, s0) != (0, 0)
    chi7 = CHI.setdefault(7, chi_table(7))
    u = set()
    for r0 in range(7):
        for s0 in range(7):
            if (r0, s0) != (0, 0): u |= image_A(r0, s0, 7, chi7)
    log('  p=7: объединение IMG_7(r0,s0) по всем (r0,s0) != (0,0):', sorted(u))
    check(u == {(0, 0)}, 'p=7: объединение образов не {(0,0)}')

    log('=== КОНТРОЛЬ C2: p = 29, 37 — контрпримеры к подъёму обязаны найтись ===')
    for p, exp in [(29, (17, 4)), (37, (29, 12))]:
        n_ok, n_bad, uimg, anyc = union_residue_stats(p)
        log(f'  p={p}: вычетов (A0,B0) без нулевых клеток с 8 линиями {n_ok}, из них нетривиальных {n_bad}; '
            f'ожидалось из claude_p7 (счёт по mod p^2, делённый на p^2): {exp}; объединённый образ {sorted(uimg)}')
        check((n_ok, n_bad) == exp, f'p={p}: не совпало с claude_p7')
        check(n_bad > 0, f'p={p}: контрпримеров нет — контроль пуст')
    log('  точный перебор по A, B mod p^M своим кодом cls_exact (сравнение с таблицей claude_p7):')
    for p, M, poles, exp in [(7, 3, (1, 2), 'всего 2413, контрпримеров 0'),
                             (29, 2, (), 'целый: 14297 / 3364'),
                             (37, 2, (), 'целый: 39701 / 16428')]:
        t = time.time()
        res = exact_grid_count(p, M, poles)
        log(f'  p={p}, M={M}: (8 линий, из них не все 9 квадраты, не определено) по n=0 (целый)/полюсам: {res}; '
            f'claude_p7: {exp}; {time.time() - t:.1f}s')
        if p == 7:
            check(sum(v[0] for v in res.values()) == 2413 and sum(v[1] for v in res.values()) == 0, 'p=7 точный счёт')
            check(res[1][0] == 0 and res[2][0] == 0, 'p=7: полюс прошёл 8 линий')
        if p == 29: check(res[0][:2] == (14297, 3364), 'p=29 точный счёт')
        if p == 37: check(res[0][:2] == (39701, 16428), 'p=37 точный счёт')
    # явный рациональный контрпример при p = 29 и 37: целые r, s, z
    for p in (29, 37):
        chi = CHI[p]; found = None
        for A0 in range(1, p):
            for B0 in range(1, p):
                vals = [(1 + i * A0 + j * B0) % p for (i, j) in CELLS]
                if 0 in vals: continue
                cl = [chi[v] for v in vals]
                if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8) and any(cl):
                    found = (A0, B0); break
            if found: break
        A0, B0 = found
        r, s, z = A0, B0, 1
        k = 0
        while gcd(r, s) != 1 or r == s:        # граница: k < p
            k += 1; s = B0 + k * p
            if k > p: break
        cl = [cls_exact(1 + m * z, 1, p) for m in mus(r, s)]
        lines_ok = all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8)
        tT = cl[T_IDX[0]] ^ cl[T_IDX[1]] ^ cl[T_IDX[2]]; tL = cl[L_IDX[0]] ^ cl[L_IDX[1]] ^ cl[L_IDX[2]]
        log(f'  явный пример p={p}: r={r}, s={s}, z=1: классы клеток {cl}; 8 линий в Q_p: {lines_ok}; '
            f'[T]={tT}, [L]={tL} (коды 2e+c)')
        check(lines_ok and any(cl), f'p={p}: явный пример не сработал')

    log('=== КОНТРОЛЬ C3: полюс при p = 3 mod 4 ===')
    # (i) лемма: при v_p(b) = -n < 0 число 1 - b^2 не квадрат в Q_p; b = a / p^n, a — единица, a mod p^2, n = 1..3
    viol = 0; tests = 0
    for p in PRIMES:
        for n in (1, 2, 3):
            for a in range(1, p * p):          # граница: p^2
                if a % p == 0: continue
                tests += 1
                c = cls_exact(p ** (2 * n) - a * a, p ** (2 * n), p)
                if c == 0: viol += 1
    log(f'  (i) 1 - (a/p^n)^2 — квадрат в Q_p: {viol} из {tests} (p = 3 mod 4 в [7,400], n = 1,2,3, a mod p^2)')
    check(viol == 0, 'полюс: 1 - b^2 оказался квадратом')
    # (ii) случайные рациональные z с v_p(z) < 0 и случайные наклоны из списка: 8 линий в Q_p не выполняются
    rng = random.Random(20260926)
    passed = 0; tests2 = 0
    for _ in range(20000):                     # граница: 20000
        r, s, Dm, Dp = nontriv[rng.randrange(len(nontriv))]
        p = PRIMES[rng.randrange(len(PRIMES))]
        n = rng.randint(1, 3)
        num = rng.randint(-10 ** 6, 10 ** 6)
        if num % p == 0: continue
        den = p ** n * rng.randint(1, 50)
        cells = [(den + m * num) for m in mus(r, s)]
        if any(c == 0 for c in cells): continue
        tests2 += 1
        cl = [cls_exact(c, den, p) for c in cells]
        if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8): passed += 1
    log(f'  (ii) случайные полюсы z (v_p(z) = -1..-3) на наклонах из списка: 8 линий выполнены в {passed} из {tests2}')
    check(passed == 0, 'полюс прошёл 8 линий')
    # (iii) отрицательный контроль: при p = 1 mod 4 полюсы с 8 линиями существуют (тест не пустой)
    for p in (5, 13):
        res = exact_grid_count(p, 2, poles=(1,))
        log(f'  (iii) p={p} (= 1 mod 4), M=2, полюс n=1: (8 линий, не все 9 квадраты, не опр.) = {res[1]}')
        check(res[1][0] > 0, f'p={p}: полюсов с 8 линиями нет — отрицательный контроль пуст')
    res = exact_grid_count(11, 2, poles=(1, 2))
    log(f'  (iii) p=11 (= 3 mod 4), M=2: целый {res[0]}, полюсы n=1: {res[1]}, n=2: {res[2]}')
    check(res[1][0] == 0 and res[2][0] == 0, 'p=11: полюс прошёл 8 линий')

    log('=== КОНТРОЛЬ C4: точные классы по целым z против метода A (корректность верхней оценки) ===')
    rng = random.Random(4242)
    sample = rng.sample(nontriv, 25)
    n_z = 0; n_ok = 0; mism = 0; realized = {}; zero_ok = 0
    for (r, s, Dm, Dp) in sample:              # граница: 25 наклонов
        for p in PRIMES:                       # граница: len(PRIMES)
            chi = CHI.setdefault(p, chi_table(p))
            conf = {z0: (cl, hz, ok) for (z0, cl, hz, ok) in residue_configs(r, s, p, chi)}
            img = get_img(r, s, p)
            zs = []
            for z0 in range(p):
                zs += [z0 + p * k for k in range(4)]         # граница: 4 подъёма на вычет
                if conf[z0][1]:
                    # подъёмы, где нулевая клетка делится на p^2: z = z* + p^2 k
                    for m in mus(r, s):
                        if (1 + m * z0) % p == 0 and m % p:
                            zst = (-pow(m, -1, p * p)) % (p * p)
                            zs += [zst + p * p * k for k in range(8)]   # граница: 8
            for z in zs:
                cells = [1 + m * z for m in mus(r, s)]
                if any(c == 0 for c in cells): continue
                n_z += 1
                cl = [cls_exact(c, 1, p) for c in cells]
                if all((cl[a] ^ cl[b] ^ cl[c]) == 0 for (a, b, c) in LINES8):
                    n_ok += 1
                    rcl, hz, ok = conf[z % p]
                    if not ok or cl != rcl: mism += 1
                    if hz: zero_ok += 1
                    tl = (cl[T_IDX[0]] ^ cl[T_IDX[1]] ^ cl[T_IDX[2]], cl[L_IDX[0]] ^ cl[L_IDX[1]] ^ cl[L_IDX[2]])
                    if tl not in img: mism += 1
                    realized.setdefault((r, s, p), set()).add(tl)
            if len(img) > 1 and realized.get((r, s, p), set()) != img:
                realized.setdefault(('не реализовано', r, s, p), set()).update(img - realized.get((r, s, p), set()))
    notreal = [k for k in realized if k[0] == 'не реализовано']
    log(f'  целых z проверено {n_z}; с 8 линиями в Q_p {n_ok} (из них с клеткой, делящейся на p: {zero_ok}); '
        f'расхождений с методом A (классы 9 клеток или ([T],[L]) вне образа): {mism}')
    log(f'  образов с элементами, не реализованными этими z (оценка сверху могла быть не точной): {len(notreal)}')
    check(mism == 0, 'точные классы разошлись с методом A')


if __name__ == '__main__':
    nt = main()
    controls(nt)
    log('ВСЕГО НАРУШЕНИЙ ВНУТРЕННИХ ПРОВЕРОК:', len(FAIL))
    for m in FAIL[:20]: log('   ', m)
    log('готово')
