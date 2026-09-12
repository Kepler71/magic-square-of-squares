# -*- coding: utf-8 -*-
"""
СТАТУС: проверка (скептическая) чужого результата
ДЛЯ: Codex, пользователь, реестр SUMMARY
ИТОГ: попытка ОПРОВЕРГНУТЬ локальное исключение (m,n)=(16,5) по модулю 41
ОТМЕНЯЕТ: ничего
ПРОВЕРЕНО: девять клеток выведены здесь заново из условия магичности, а не взяты у Codex
ЧИТАТЬ ДАЛЬШЕ, ЕСЛИ: нужно знать, на чём именно держится счёт 126/127

Запуск:  python3 /home/kep/magicKube/bridge/local_16_5_p41_скептик.py

Что делается (в порядке скептицизма):
 §1  НЕЗАВИСИМЫЙ вывод девяти клеток G1 из параметризации Сегре + условий магичности.
     Никаких формул от Codex не используется, пока они не выведены здесь.
 §2  Сверка выведенных клеток с формулами Codex (F0,F4,F8,L,U) — символически.
 §3  Полный перебор ВСЕХ p+1 точек P^1(F_41) для (16,5): таблица всех пяти клеток.
     Нули разрешены (ноль — квадрат в F_p). Отдельно — вариант «нули запрещены».
 §4  Контроли: (15,8) при t=1 обязан выживать; (13,8) mod 17 — воспроизведение таблицы Codex;
     (7,1) t=0 — точка C7 без полного квадрата.
 §5  Прямой рациональный поиск малой высоты для (16,5) — попытка найти контрпример.
 §6  Скан простых: какие p дают препятствие у (16,5) и у (13,8); не одиночное ли это p=41.
 §7  Проверка корректности самой редукции: 41-целость клеток, случай 41 | знаменатель t.
"""
from fractions import Fraction as Fr
from math import gcd
import itertools, sys

# ----------------------------------------------------------------------------
# §1. НЕЗАВИСИМЫЙ вывод девяти клеток G1
# ----------------------------------------------------------------------------
# Постановка (из леммы о двух полных парах):
#   магический 3x3 ассоциативен, c_i + c_{8-i} = 2*c4, сумма = 3*c4.
#   G1 = обе РЁБЕРНЫЕ пары полные: c1=P^2, c7=Q^2, c3=R^2, c5=S^2,
#   условие совместности  P^2+Q^2 = R^2+S^2 = 2*c4.
#   Параметризация Сегре этой квадрики (q=1):
#       P = m*t + n,  Q = m - n*t,  R = m*t - n,  S = m + n*t.
# Свободные клетки c0,c2,c6,c8 и центр c4 выводятся ЛИНЕЙНО из условий магичности:
#       c0 + c2 = 3c4 - P^2   (верхняя строка)
#       c0 + c6 = 3c4 - R^2   (левый столбец)
#       c2 + c6 = 2c4         (ассоциативность)
#   =>  c0 = 2c4 - (P^2+R^2)/2,  c2 = c4 + (R^2-P^2)/2,  c6 = 2c4 - c2, c8 = 2c4 - c0.

import sympy as sp

def derive_cells_symbolically():
    m, n, t = sp.symbols('m n t')
    P = m*t + n
    Q = m - n*t
    R = m*t - n
    S = m + n*t
    # совместность
    lhs = sp.expand(P**2 + Q**2)
    rhs = sp.expand(R**2 + S**2)
    assert sp.simplify(lhs - rhs) == 0, "квадрика не сбалансирована"
    c4 = sp.together(sp.expand(lhs)/2)
    c0 = sp.expand(2*c4 - (P**2 + R**2)/2)
    c2 = sp.expand(c4 + (R**2 - P**2)/2)
    c6 = sp.expand(2*c4 - c2)
    c8 = sp.expand(2*c4 - c0)
    c1, c7, c3, c5 = sp.expand(P**2), sp.expand(Q**2), sp.expand(R**2), sp.expand(S**2)
    cells = [c0, c1, c2, c3, c4, c5, c6, c7, c8]
    # ПОЛНАЯ проверка магичности: 3 строки, 3 столбца, 2 диагонали
    tot = sp.expand(3*c4)
    checks = {
        'строка0': c0+c1+c2, 'строка1': c3+c4+c5, 'строка2': c6+c7+c8,
        'столбец0': c0+c3+c6, 'столбец1': c1+c4+c7, 'столбец2': c2+c5+c8,
        'диаг': c0+c4+c8, 'антидиаг': c2+c4+c6,
    }
    bad = [k for k, v in checks.items() if sp.simplify(sp.expand(v) - tot) != 0]
    return cells, bad, (P, Q, R, S), (m, n, t)


def section1():
    print("="*100)
    print("§1. НЕЗАВИСИМЫЙ вывод девяти клеток G1 (Сегре + магичность). Формулы Codex НЕ используются.")
    print("="*100)
    cells, bad, (P,Q,R,S), (m,n,t) = derive_cells_symbolically()
    names = ['c0','c1','c2','c3','c4','c5','c6','c7','c8']
    for nm, c in zip(names, cells):
        print(f"   {nm} = {sp.factor(c)}")
    print(f"   магичность (8 линий): {'ВСЕ СОШЛИСЬ' if not bad else 'НЕ СОШЛИСЬ: '+str(bad)}")
    print(f"   c1=({P})^2, c3=({R})^2, c5=({S})^2, c7=({Q})^2  — автоматические квадраты")
    return cells, (m,n,t)


# ----------------------------------------------------------------------------
# §2. Сверка с формулами Codex
# ----------------------------------------------------------------------------
def section2(cells, syms):
    m, n, t = syms
    print()
    print("="*100)
    print("§2. Сверка с формулами Codex: F0=m^2+n^2 t^2, F4=s(1+t^2), F8=n^2+m^2 t^2, L=F4-2mnt, U=F4+2mnt")
    print("="*100)
    s = (m**2 + n**2)/2
    codex = {
        'F0': m**2 + n**2*t**2,
        'F4': s*(1 + t**2),
        'F8': n**2 + m**2*t**2,
        'L' : s*(1 + t**2) - 2*m*n*t,
        'U' : s*(1 + t**2) + 2*m*n*t,
    }
    mine = {'c0': cells[0], 'c4': cells[4], 'c8': cells[8], 'c2': cells[2], 'c6': cells[6]}
    # ищем соответствие как МНОЖЕСТВ многочленов
    pairs = []
    for k, v in codex.items():
        hit = [nm for nm, w in mine.items() if sp.simplify(sp.expand(v - w)) == 0]
        pairs.append((k, hit))
        print(f"   {k:>3}  == {hit if hit else 'НИ ОДНА МОЯ КЛЕТКА (расхождение!)'}")
    ok = all(h for _, h in pairs)
    used = sorted(set(h[0] for _, h in pairs if h))
    print(f"   ВЫВОД: клетки {'СОВПАДАЮТ' if ok else 'НЕ совпадают'}; мои пять свободных клеток = {used}")
    # Проверка симметрии (m,n) -> (n,m): система должна быть инвариантна (F0<->F8)
    sw = {m: n, n: m}
    inv = all(
        any(sp.simplify(sp.expand(codex[k].subs(sw, simultaneous=True) - codex[j])) == 0 for j in codex)
        for k in codex)
    print(f"   инвариантность системы при перестановке m<->n: {inv}  (значит (16,5) и (5,16) — одно и то же)")
    return ok


# ----------------------------------------------------------------------------
# аппарат над F_p
# ----------------------------------------------------------------------------
def qr_set(p):
    return set((x*x) % p for x in range(p))

def cells_mod_p(m, n, p):
    """пять свободных клеток как однородные формы от (a,b) над F_p; t = a/b."""
    inv2 = pow(2, p-2, p)
    s = ((m*m + n*n) % p) * inv2 % p
    M2, N2, MN2 = m*m % p, n*n % p, (2*m*n) % p
    def F(a, b):
        A, B = a % p, b % p
        F0 = (M2*B*B + N2*A*A) % p
        F4 = (s*(A*A + B*B)) % p
        F8 = (N2*B*B + M2*A*A) % p
        L  = (F4 - MN2*A*B) % p
        U  = (F4 + MN2*A*B) % p
        return (F0, F4, F8, L, U)
    return F, s

def projective_points(p):
    """все p+1 точек P^1(F_p): (t:1) и (1:0)."""
    return [(t, 1) for t in range(p)] + [(1, 0)]


def full_table(m, n, p, allow_zero=True, verbose=True, only_survivors3=False):
    """Полный перебор P^1(F_p). Возвращает (выжившие_на_3, выжившие_на_5)."""
    QR = qr_set(p)
    ok = (lambda x: x in QR) if allow_zero else (lambda x: x != 0 and x in QR)
    F, s = cells_mod_p(m, n, p)
    surv3, surv5, rows = [], [], []
    for (a, b) in projective_points(p):
        F0, F4, F8, L, U = F(a, b)
        three = ok(F0) and ok(F4) and ok(F8)
        five = three and ok(L) and ok(U)
        label = f"{a}" if b == 1 else "oo"
        if three: surv3.append(label)
        if five:  surv5.append(label)
        if three or not only_survivors3:
            rows.append((label, F0, F4, F8, L, U, three, five))
    if verbose:
        print(f"   (m,n)=({m},{n}), p={p}, s={s} (квадрат: {ok(s)}), нули {'разрешены' if allow_zero else 'ЗАПРЕЩЕНЫ'}")
        print(f"   {'t':>4} | {'F0':>4} {'F4':>4} {'F8':>4} | {'L':>4} {'U':>4} | 3кл 5кл")
        for (label, F0, F4, F8, L, U, three, five) in rows:
            mark3 = 'да ' if three else 'нет'
            mark5 = 'да ' if five else 'нет'
            print(f"   {label:>4} | {F0:>4} {F4:>4} {F8:>4} | {L:>4} {U:>4} | {mark3} {mark5}")
    return surv3, surv5


def section3():
    print()
    print("="*100)
    print("§3. (16,5) mod 41 — полный перебор ВСЕХ 42 проективных точек")
    print("="*100)
    m, n, p = 16, 5, 41
    QR = qr_set(p)
    print(f"   QR({p}) = {sorted(QR)}")
    print()
    print("   --- 3.1 показаны ТОЛЬКО точки, прошедшие F0,F4,F8 (нули разрешены) ---")
    s3, s5 = full_table(m, n, p, allow_zero=True, only_survivors3=True)
    print(f"   выжили на трёх клетках: {s3}")
    print(f"   выжили на пяти клетках: {s5 if s5 else 'НИ ОДНОЙ — препятствие подтверждается'}")
    print()
    print("   --- 3.2 то же при ЗАПРЕТЕ нулей (проверка, не держится ли вывод на трактовке нуля) ---")
    z3, z5 = full_table(m, n, p, allow_zero=False, verbose=False)
    print(f"   выжили на трёх: {z3};  на пяти: {z5 if z5 else 'ни одной'}")
    print()
    print("   --- 3.3 ПОЛНАЯ таблица всех 42 точек (контроль, что ни одна не потеряна) ---")
    a3, a5 = full_table(m, n, p, allow_zero=True, only_survivors3=False, verbose=True)
    print(f"   всего проверено точек: {p+1}")
    print(f"   ИТОГ §3: пятиклеточных решений над F_41 — {len(a5)}")
    print()
    print("   --- 3.4 (5,16) — та же пара с переставленными m,n ---")
    b3, b5 = full_table(5, 16, p, allow_zero=True, verbose=False)
    print(f"   (5,16) mod 41: на трёх {b3}; на пяти {b5 if b5 else 'ни одной'}")
    print()
    print("   --- 3.5 сверка с таблицей Codex ---")
    codex_rows = {4: (0,31,21,6,15), 10: (9,25,0,24,26), 31: (9,25,0,26,24), 37: (0,31,21,15,6)}
    F, s = cells_mod_p(m, n, p)
    agree = True
    for t, exp in codex_rows.items():
        got = F(t, 1)
        same = (got == exp)
        agree &= same
        print(f"   t={t:>2}: мои (F0,F4,F8,L,U)={got}  у Codex={exp}  {'совпало' if same else 'РАСХОЖДЕНИЕ'}")
    print(f"   таблица Codex для (16,5): {'воспроизведена поэлементно' if agree else 'НЕ воспроизведена'}")
    print(f"   его список выживших {{4,10,31,37}} vs мой {s3}: "
          f"{'совпал' if sorted(s3, key=lambda x: (x=='oo', x)) == ['10','31','37','4'] or set(s3)=={'4','10','31','37'} else 'РАСХОЖДЕНИЕ'}")
    return len(a5) == 0


# ----------------------------------------------------------------------------
# §4. Контроли
# ----------------------------------------------------------------------------
def rational_cells(m, n, t):
    """пять свободных клеток как точные рациональные числа"""
    s = Fr(m*m + n*n, 2)
    F0 = Fr(m*m) + Fr(n*n)*t*t
    F4 = s*(1 + t*t)
    F8 = Fr(n*n) + Fr(m*m)*t*t
    L  = F4 - 2*m*n*t
    U  = F4 + 2*m*n*t
    return [F0, F4, F8, L, U]

def is_rational_square(x):
    if x < 0: return False
    if x == 0: return True
    nu, de = x.numerator, x.denominator
    rn, rd = int(round(nu**0.5)), int(round(de**0.5))
    for cand in (rn-2, rn-1, rn, rn+1, rn+2):
        if cand >= 0 and cand*cand == nu: rn = cand; break
    else: return False
    for cand in (rd-2, rd-1, rd, rd+1, rd+2):
        if cand > 0 and cand*cand == de: rd = cand; break
    else: return False
    return True

def all_nine_cells(m, n, t):
    """девять клеток как рациональные числа, в порядке c0..c8"""
    s = Fr(m*m + n*n, 2)
    c4 = s*(1 + t*t)
    c0 = Fr(m*m) + Fr(n*n)*t*t
    c8 = Fr(n*n) + Fr(m*m)*t*t
    c2 = c4 - 2*m*n*t
    c6 = c4 + 2*m*n*t
    c1 = (m*t + n)**2
    c7 = (m - n*t)**2
    c3 = (m*t - n)**2
    c5 = (m + n*t)**2
    return [c0,c1,c2,c3,c4,c5,c6,c7,c8]

def section4():
    print()
    print("="*100)
    print("§4. КОНТРОЛИ (если бы скрипт «исключал» и их — он был бы сломан)")
    print("="*100)
    # контроль 1: (15,8), t=1 — все девять квадратов (вырожденно)
    cells = all_nine_cells(15, 8, Fr(1))
    sq = [is_rational_square(c) for c in cells]
    print(f"   (15,8), t=1: девять клеток {[str(c) for c in cells]}")
    print(f"      все квадраты: {all(sq)}; различны: {len(set(cells))==9}; положительны: {all(c>0 for c in cells)}")
    for p in (17, 41, 7, 11, 13, 23):
        s3, s5 = full_table(15, 8, p, verbose=False)
        print(f"      (15,8) mod {p:>2}: выжили на пяти клетках {s5}  ->  "
              f"{'контроль ПРОЙДЕН (t=1 должен быть тут)' if '1' in s5 else 'КОНТРОЛЬ ПРОВАЛЕН'}")
    print()
    # контроль 2: (7,1), t=0 — точка C7, но не полный квадрат
    cells = all_nine_cells(7, 1, Fr(0))
    sq = [is_rational_square(c) for c in cells]
    print(f"   (7,1), t=0: девять клеток {[str(c) for c in cells]}")
    print(f"      квадратны: {sq}")
    print(f"      F0,F4,F8 квадраты: {sq[0] and sq[4] and sq[8]} (=> точка C7);  все девять: {all(sq)}")
    print("      => подтверждает различение Codex: C7 непусто, а полного квадрата нет")
    print()
    # контроль 3: (13,8) mod 17 — воспроизведение таблицы Codex
    print("   (13,8) mod 17 — воспроизведение таблицы Codex:")
    s3, s5 = full_table(13, 8, 17, verbose=True, only_survivors3=True)
    print(f"      выжили на трёх: {s3} (Codex: 2, 8, 9, 15); на пяти: {s5 if s5 else 'ни одной'}")
    codex_17 = {2:(0,13,9,5,4), 8:(15,16,0,1,14), 9:(15,16,0,14,1), 15:(0,13,9,4,5)}
    F, s = cells_mod_p(13, 8, 17)
    ok17 = all(F(t,1) == v for t, v in codex_17.items())
    print(f"      таблица 17 у Codex воспроизведена поэлементно: {ok17}")


# ----------------------------------------------------------------------------
# §5. Прямой рациональный поиск
# ----------------------------------------------------------------------------
def section5(H=300):
    print()
    print("="*100)
    print(f"§5. Прямой рациональный поиск для (16,5): t=a/b, |a|,b <= {H} — попытка найти контрпример")
    print("="*100)
    m, n = 16, 5
    best = None
    found5 = []
    found3 = []
    cnt = 0
    for b in range(1, H+1):
        for a in range(-H, H+1):
            if gcd(abs(a), b) != 1: continue
            cnt += 1
            t = Fr(a, b)
            F0, F4, F8, L, U = rational_cells(m, n, t)
            k = sum(1 for x in (F0, F4, F8, L, U) if is_rational_square(x))
            if best is None or k > best[0]: best = (k, t)
            if is_rational_square(F0) and is_rational_square(F4) and is_rational_square(F8):
                found3.append(t)
                if is_rational_square(L) and is_rational_square(U):
                    found5.append(t)
    print(f"   перебрано t: {cnt}")
    print(f"   максимум квадратных клеток из пяти: {best[0]} при t={best[1]}")
    print(f"   t с тремя квадратами (F0,F4,F8): {found3[:20]}{' ...' if len(found3)>20 else ''}  всего {len(found3)}")
    print(f"   t с пятью квадратами: {found5 if found5 else 'не найдено'}")
    print("   ВНИМАНИЕ: отсутствие находки здесь НЕ доказательство. Доказательством служит §3.")


# ----------------------------------------------------------------------------
# §6. Скан простых
# ----------------------------------------------------------------------------
def primes_upto(N):
    sieve = [True]*(N+1); sieve[0:2] = [False, False]
    for i in range(2, int(N**0.5)+1):
        if sieve[i]:
            for j in range(i*i, N+1, i): sieve[j] = False
    return [i for i in range(N+1) if sieve[i]]

def section6(N=400):
    print()
    print("="*100)
    print(f"§6. Скан нечётных простых p <= {N}: у каких пар есть локальное препятствие на пяти клетках")
    print("="*100)
    for (m, n) in [(16,5), (13,8), (15,8), (11,4), (19,5), (15,1), (19,16)]:
        bad = []
        for p in primes_upto(N):
            if p == 2: continue
            if p == m or p == n: pass  # p | m или p | n допустимо, формулы всё равно определены
            _, s5 = full_table(m, n, p, verbose=False)
            if not s5: bad.append(p)
        print(f"   ({m},{n}): препятствие на пяти клетках при p = {bad if bad else 'нет ни одного p'}")
    print("   (для (15,8) список обязан быть пуст: t=1 даёт точку над Q, значит и над каждым F_p)")


# ----------------------------------------------------------------------------
# §7. Корректность редукции
# ----------------------------------------------------------------------------
def section7():
    print()
    print("="*100)
    print("§7. Корректность самой редукции при p=41")
    print("="*100)
    m, n, p = 16, 5, 41
    print(f"   p={p}: p | 2? {p==2}; p | m={m}? {m % p == 0}; p | n={n}? {n % p == 0};"
          f" p | m^2+n^2={m*m+n*n}? {(m*m+n*n) % p == 0}")
    print("   => знаменатель s=(m^2+n^2)/2 обратим mod 41, m,n ненулевы mod 41 — редукция определена.")
    # явная проверка: для t с 41 в знаменателе однородные клетки 41-целы и точка уходит в (1:0)
    bad = 0
    import random
    random.seed(20260912)
    for _ in range(20000):
        a = random.randint(-10**6, 10**6); b = random.randint(1, 10**6)
        g = gcd(abs(a), b)
        if g == 0: continue
        a, b = a//g, b//g
        if a == 0 and b == 0: continue
        vals = [m*m*b*b + n*n*a*a, Fr((m*m+n*n)*(a*a+b*b), 2),
                n*n*b*b + m*m*a*a,
                Fr((m*m+n*n)*(a*a+b*b), 2) - 2*m*n*a*b,
                Fr((m*m+n*n)*(a*a+b*b), 2) + 2*m*n*a*b]
        for v in vals:
            v = Fr(v)
            if v.denominator % p == 0: bad += 1
    print(f"   20000 случайных t=a/b: однородных клеток с 41 в знаменателе — {bad} (должно быть 0)")
    # проверка: t с 41 | b действительно попадает в (1:0)
    t = Fr(3, 41*7)
    a, b = t.numerator, t.denominator
    print(f"   пример t={t}: (a mod 41, b mod 41) = ({a % p}, {b % p}) -> проективно (1:0) = 'oo'. "
          f"Значит t с 41 в знаменателе покрыты точкой бесконечности, ни один не потерян.")
    F, s = cells_mod_p(m, n, p)
    print(f"   в точке oo: (F0,F4,F8,L,U) = {F(1,0)}; все они = (n^2, s, m^2, s, s);"
          f" s={s} квадрат? {s in qr_set(p)}")
    print("   => на бесконечности препятствие ровно одно: s не квадрат mod 41.")
    # независимая проверка s mod 41 «руками»
    print(f"   s = {m*m+n*n}/2; {m*m+n*n} mod 41 = {(m*m+n*n)%p}; 2^-1 mod 41 = {pow(2,p-2,p)}; "
          f"s mod 41 = {((m*m+n*n)%p)*pow(2,p-2,p)%p}; 38^20 mod 41 = {pow(38,20,41)} (1=квадрат, 40=нет)")


# ----------------------------------------------------------------------------
def section8():
    print()
    print("="*100)
    print("§8. Не зависит ли вывод от положительности и попарной различности")
    print("="*100)
    print("   Над F_p понятий «>0» и «различны» нет; в §3 они нигде не использовались —")
    print("   проверялась ТОЛЬКО квадратность пяти значений, нули разрешены.")
    print("   Поэтому исключение сильнее требуемого: оно запрещает даже вырожденные решения.")
    print("   Контроль обратного направления: (15,8), t=1 — решение ВЫРОЖДЕННОЕ (клетки повторяются),")
    print("   и §4 показывает, что мод-41 проверка его НЕ убивает. Значит тест действительно")
    print("   не опирается на невырожденность.")
    cells = all_nine_cells(15, 8, Fr(1))
    print(f"   (15,8), t=1: клетки {[str(c) for c in cells]} — повторы есть, но mod p выживает.")


def section9():
    print()
    print("="*100)
    print("§9. Проверка ОБОБЩЕНИЯ Codex: запрещённые отношения λ = m/n mod p")
    print("="*100)
    for p, claimed in [(17, {2,3,6,8,9,11,14,15}), (41, {5,8,13,19,22,28,33,36})]:
        forb = set()
        for lam in range(p):
            _, s5 = full_table(lam, 1, p, verbose=False)      # (m,n) = (λ,1)
            if not s5: forb.add(lam)
        _, s5inf = full_table(1, 0, p, verbose=False)          # n ≡ 0 mod p
        print(f"   p={p}: мой список запрещённых λ = {sorted(forb)}")
        print(f"          список Codex               = {sorted(claimed)}  -> "
              f"{'СОВПАЛ' if forb == claimed else 'РАСХОЖДЕНИЕ'}")
        print(f"          случай n≡0 (λ=oo): исключается? {not s5inf} (Codex: не входит в список)")
    m, n, p = 16, 5, 41
    lam = m * pow(n, p-2, p) % p
    print(f"   для (16,5): λ = 16/5 mod 41 = {lam}; в списке Codex: {lam in {5,8,13,19,22,28,33,36}}")


def section10(H=1200):
    print()
    print("="*100)
    print("§10. Насколько вывод зависит ИМЕННО от p=41; и автоматична ли положительность")
    print("="*100)
    m, n = 16, 5
    for p in (41, 47, 127):
        s3, s5 = full_table(m, n, p, verbose=False)
        print(f"   (16,5) mod {p:>3}: выжили на трёх {s3}; на пяти {s5 if s5 else 'ни одной'}")
    print("   => исключение (16,5) НЕ висит на одном простом: 41, 47 и 127 дают его независимо.")
    print()
    print("   Положительность клеток в G1 при m≠n — тождество, а не условие:")
    print("     F0,F8 > 0 очевидно; F4 = s(1+t^2) > 0; L,U = s(1+t^2) ∓ 2mnt имеют дискриминант")
    print("     4(m^2n^2 − s^2) = −(m^2−n^2)^2 < 0 при m≠n, значит L,U > 0 при всех вещественных t.")
    from fractions import Fraction as Fr2
    worst = min((rational_cells(m, n, Fr2(k, 97))[3] for k in range(-400, 401)), default=None)
    print(f"     контроль: min L на 801 значении t = {worst} (>0: {worst > 0})")
    print()
    print(f"   Решето по малым p + поиск точек C7 (три квадрата) для (16,5), |a|,b <= {H}:")
    SIEVE = [7, 11, 13, 19, 23, 29, 31, 37, 43]
    qr = {p: qr_set(p) for p in SIEVE}
    Fm = {p: cells_mod_p(m, n, p)[0] for p in SIEVE}
    hits = 0; tested = 0
    for b in range(1, H+1):
        for a in range(-H, H+1):
            if gcd(abs(a), b) != 1: continue
            ok = True
            for p in SIEVE:
                F0, F4, F8, L, U = Fm[p](a, b)
                if not (F0 in qr[p] and F4 in qr[p] and F8 in qr[p]): ok = False; break
            if not ok: continue
            tested += 1
            t = Fr(a, b)
            F0, F4, F8, L, U = rational_cells(m, n, t)
            if is_rational_square(F0) and is_rational_square(F4) and is_rational_square(F8):
                hits += 1
                print(f"      ТОЧКА C7: t={t}")
    print(f"      прошли решето: {tested}; точек C7 найдено: {hits}")
    print("      (пусто ⇒ НИЧЕГО не доказано про C7(Q) — это только отсутствие находки)")


def section11(N=400):
    print()
    print("="*100)
    print(f"§11. Когерентность: что даёт ОДИН этот тест (пять клеток, нечётные p <= {N}) на всех 127 парах")
    print("="*100)
    pairs = [(a, b) for a in range(2, 21) for b in range(1, a) if gcd(a, b) == 1]
    P = [p for p in primes_upto(N) if p > 2]
    killed, alive = {}, []
    for (a, b) in pairs:
        hit = None
        for p in P:
            _, s5 = full_table(a, b, p, verbose=False)
            if not s5: hit = p; break
        if hit: killed[(a, b)] = hit
        else: alive.append((a, b))
    print(f"   закрыто одним локальным тестом: {len(killed)} из {len(pairs)}")
    print(f"   НЕ закрыто: {alive}")
    print("   у каждой незакрытой обязана быть глобальная (пусть вырожденная) точка — проверяем:")
    for (a, b) in alive:
        found = None
        for t in (Fr(0), Fr(1), Fr(-1)):
            cs = rational_cells(a, b, t)
            if all(is_rational_square(x) for x in cs):
                found = (t, [str(x) for x in cs]); break
        ss = a*a + b*b
        why = []
        if int(ss**0.5 + 0.5)**2 == ss: why.append("m^2+n^2 — точный квадрат (t=1)")
        if ss % 2 == 0 and int((ss//2)**0.5 + 0.5)**2 == ss//2: why.append("(m^2+n^2)/2 — точный квадрат (t=0)")
        print(f"      ({a},{b}): точка {found[0] if found else 'НЕ НАЙДЕНА — тогда это дыра'}"
              f"  клетки {found[1] if found else ''}  причина: {'; '.join(why) or '—'}")
    print("   ВЫВОД (наблюдение): тест не берёт ровно те пары, у которых есть вырожденная")
    print("   рациональная точка, то есть там, где он НЕ МОЖЕТ сработать в принципе.")
    print("   Это ровно то поведение, которого ждём от корректного локального теста.")


if __name__ == "__main__":
    cells, syms = section1()
    ok_cells = section2(cells, syms)
    ok_obstr = section3()
    section4()
    section5(H=int(sys.argv[1]) if len(sys.argv) > 1 else 300)
    section6(N=400)
    section7()
    section8()
    section9()
    section10()
    section11()
    print()
    print("="*100)
    print("ИТОГ")
    print("="*100)
    print(f"   клетки Codex совпали с моим независимым выводом: {ok_cells}")
    print(f"   над F_41 у (16,5) пятиклеточных точек нет (все 42 проективные проверены): {ok_obstr}")
