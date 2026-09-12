# -*- coding: utf-8 -*-
"""
СКЕПТИЧЕСКАЯ ПРОВЕРКА заявленного локального исключения (m,n) = (13,8) по модулю 17.

Источник проверяемого утверждения: PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md, §«(13,8), p=17».
Утверждение Codex: полная система девяти клеток G1 для (13,8) не имеет решений уже над F_17,
причём вывод не опирается ни на ранги, ни на положительность/различность клеток.

Задача скрипта — ПОПЫТКА ОПРОВЕРЖЕНИЯ. Он:
  A. выводит девять клеток G1 НЕЗАВИСИМО (из аксиом магического квадрата, а не из чужой формулы)
     и проверяет их символически (sympy);
  B. сверяет полученные клетки с клетками, которые использует Codex;
  C. перебирает ВСЕ p+1 точек P^1(F_p), включая бесконечность и точки с нулями,
     при четырёх разных соглашениях (ноль — квадрат/не квадрат; пять клеток/все девять;
     плюс нетривиальный квадратичный твист c);
  D. прогоняет положительные контроли: (15,8) (есть рациональная точка t=±1),
     (7,1) (есть рациональная точка t=0) — машина обязана их НЕ исключать;
  E. проверяет, что над Q_17 (а не только над F_17) точек тоже нет — явным подъёмом Гензеля;
  F. независимо перебирает рациональные t = a/b и считает, сколько клеток вышло квадратами;
  G. проверяет опубликованный список «запрещённых λ = m/n mod 17».

Запуск: python3 /home/kep/magicKube/bridge/local_13_8_p17_скептик.py
"""

from fractions import Fraction as Fr
from math import gcd, isqrt
import itertools
import sys

import sympy as sp

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


# ----------------------------------------------------------------------------
# A. НЕЗАВИСИМЫЙ ВЫВОД ДЕВЯТИ КЛЕТОК
# ----------------------------------------------------------------------------
# Общий магический квадрат 3x3 с центром c и двумя свободными параметрами A, B:
#
#     [ c+A     c-A-B   c+B   ]
#     [ c-A+B   c       c+A-B ]
#     [ c-B     c+A+B   c-A   ]
#
# (все восемь сумм равны 3c тождественно — проверяется ниже).
#
# Семейство G1 фиксируется требованиями:
#   центр            c  = s(1+t^2),  s = (m^2+n^2)/2,
#   левый верх       c+A = F0 = m^2 + n^2 t^2   =>  A = F0 - c,
#   правый низ       c-A = F8 = n^2 + m^2 t^2   (выполняется автоматически, т.к. F0+F8 = 2c),
#   четыре «автоматические» клетки c-A-B, c-A+B, c+A-B, c+A+B — полные квадраты.
# Последнее требование решается подстановкой B = -2mn t (проверяется ниже символически).
# Тогда две оставшиеся клетки: L = c+B = c - 2mn t и U = c-B = c + 2mn t — это «красные» c2, c6.

m, n, t, c, A, B = sp.symbols('m n t c A B')


def derive_cells():
    s = (m**2 + n**2) / 2
    c_val = s * (1 + t**2)
    F0 = m**2 + n**2 * t**2
    F8 = n**2 + m**2 * t**2
    A_val = sp.expand(F0 - c_val)
    B_val = -2 * m * n * t

    M_gen = sp.Matrix([[c + A, c - A - B, c + B],
                       [c - A + B, c, c + A - B],
                       [c - B, c + A + B, c - A]])

    # все восемь сумм тождественно равны 3c
    sums = [sum(M_gen.row(i)) for i in range(3)] + [sum(M_gen.col(j)) for j in range(3)]
    sums.append(M_gen[0, 0] + M_gen[1, 1] + M_gen[2, 2])
    sums.append(M_gen[0, 2] + M_gen[1, 1] + M_gen[2, 0])
    magic_ok = all(sp.simplify(x - 3 * c) == 0 for x in sums)

    M = M_gen.subs({c: c_val, A: A_val, B: B_val}).applyfunc(sp.expand)
    # контроль: F0 и F8 на диагонали
    diag_ok = (sp.simplify(M[0, 0] - F0) == 0) and (sp.simplify(M[2, 2] - F8) == 0)

    auto = {
        (0, 1): (m * t + n)**2,
        (1, 0): (m * t - n)**2,
        (1, 2): (m + n * t)**2,
        (2, 1): (m - n * t)**2,
    }
    auto_ok = all(sp.simplify(M[i, j] - sp.expand(v)) == 0 for (i, j), v in auto.items())
    return M, magic_ok, diag_ok, auto_ok, sp.expand(c_val)


def section_A():
    say("=" * 78)
    say("A. НЕЗАВИСИМЫЙ ВЫВОД ДЕВЯТИ КЛЕТОК G1 (из аксиом магического квадрата)")
    say("=" * 78)
    M, magic_ok, diag_ok, auto_ok, c_val = derive_cells()
    say("Матрица (мой вывод):")
    for i in range(3):
        say("   [ " + " | ".join(str(sp.factor(M[i, j])) for j in range(3)) + " ]")
    say(f"  все 8 магических сумм = 3*центр тождественно : {magic_ok}")
    say(f"  на диагонали стоят F0 = m^2+n^2 t^2 и F8 = n^2+m^2 t^2 : {diag_ok}")
    say(f"  четыре клетки — полные квадраты (mt+n)^2,(mt-n)^2,(m+nt)^2,(m-nt)^2 : {auto_ok}")
    say(f"  центр = {sp.factor(c_val)}")
    say("  ВЫВОД: неавтоматических клеток ровно пять: F0, F4=центр, F8, L=F4-2mnt, U=F4+2mnt.")
    say("         L стоит в позиции (0,2), U — в позиции (2,0); это и есть «красные» c2, c6.")
    return M


# ----------------------------------------------------------------------------
# B. СВЕРКА С КЛЕТКАМИ CODEX
# ----------------------------------------------------------------------------
def section_B(M):
    say()
    say("=" * 78)
    say("B. СВЕРКА С ФОРМУЛАМИ CODEX (PARALLEL_RESULT_FULL_NINE_126_2026-09-12.md)")
    say("=" * 78)
    a, b = sp.symbols('a b')
    s = (m**2 + n**2) / 2
    # однородные формы Codex, t = a/b
    cod = {
        'F0': m**2 * b**2 + n**2 * a**2,
        'F4': s * (a**2 + b**2),
        'F8': n**2 * b**2 + m**2 * a**2,
        'L': s * (a**2 + b**2) - 2 * m * n * a * b,
        'U': s * (a**2 + b**2) + 2 * m * n * a * b,
    }
    mine = {'F0': M[0, 0], 'F4': M[1, 1], 'F8': M[2, 2], 'L': M[0, 2], 'U': M[2, 0]}
    allok = True
    for k in ['F0', 'F4', 'F8', 'L', 'U']:
        diff = sp.simplify(sp.expand(cod[k] / b**2 - mine[k].subs(t, a / b)))
        ok = (diff == 0)
        allok &= ok
        say(f"  {k:2s}: Codex(a,b)/b^2 == моё(t=a/b) : {ok}")
    say(f"  ИТОГ СВЕРКИ КЛЕТОК: {'совпадают полностью' if allok else 'РАСХОЖДЕНИЕ!'}")
    # также проверяем, что обезразмеривание b^2 — квадрат (не меняет класс квадратов)
    say("  Переход неоднородное -> однородное есть умножение на b^2 (квадрат): класс квадратов сохранён.")
    return allok


# ----------------------------------------------------------------------------
# C. ПЕРЕБОР P^1(F_p)
# ----------------------------------------------------------------------------
def qr_set(p):
    return set((x * x) % p for x in range(p))


def cells_mod(mm, nn, a, b, p):
    """пять неавтоматических клеток как однородные формы, приведённые в F_p.
       s = (m^2+n^2)/2 — при нечётном p делитель 2 обратим."""
    inv2 = pow(2, p - 2, p)
    s = ((mm * mm + nn * nn) % p) * inv2 % p
    F0 = (mm * mm * b * b + nn * nn * a * a) % p
    F4 = s * (a * a + b * b) % p
    F8 = (nn * nn * b * b + mm * mm * a * a) % p
    L = (F4 - 2 * mm * nn * a * b) % p
    U = (F4 + 2 * mm * nn * a * b) % p
    return {'F0': F0, 'F4': F4, 'F8': F8, 'L': L, 'U': U}


def auto_cells_mod(mm, nn, a, b, p):
    """четыре автоматические клетки, однородные: (ma+nb)^2 и т.д."""
    return {
        '(mt+n)^2': pow((mm * a + nn * b) % p, 2, p),
        '(mt-n)^2': pow((mm * a - nn * b) % p, 2, p),
        '(m+nt)^2': pow((mm * b + nn * a) % p, 2, p),
        '(m-nt)^2': pow((mm * b - nn * a) % p, 2, p),
    }


def projective_points(p):
    """все p+1 точек P^1(F_p) как пары (a,b): (t:1) и (1:0)."""
    pts = [(x, 1) for x in range(p)]
    pts.append((1, 0))
    return pts


def scan(mm, nn, p, zero_ok=True, nine=False, twist=1, verbose=False):
    """возвращает список выживших точек (a,b)."""
    QR = qr_set(p)
    QRnz = QR - {0}
    survivors = []
    rows = []
    for (a, b) in projective_points(p):
        vals = cells_mod(mm, nn, a, b, p)
        if nine:
            vals = dict(vals, **auto_cells_mod(mm, nn, a, b, p))
        bad = []
        for k, v in vals.items():
            w = v * twist % p
            good = (w in QR) if zero_ok else (w in QRnz)
            if not good:
                bad.append(k)
        rows.append(((a, b), vals, bad))
        if not bad:
            survivors.append((a, b))
    if verbose:
        for (ab, vals, bad) in rows:
            tt = f"t={ab[0]}" if ab[1] == 1 else "t=oo"
            vs = " ".join(f"{k}={v}" for k, v in vals.items())
            mark = "ВЫЖИЛА" if not bad else ("не кв.: " + ",".join(bad))
            say(f"    {tt:8s} {vs:60s} {mark}")
    return survivors, rows


def section_C():
    say()
    say("=" * 78)
    say("C. ПЕРЕБОР ВСЕХ ТОЧЕК P^1(F_17) ДЛЯ (m,n) = (13,8)")
    say("=" * 78)
    p = 17
    mm, nn = 13, 8
    inv2 = pow(2, p - 2, p)
    say(f"  константы mod {p}: m^2={mm*mm%p}, n^2={nn*nn%p}, "
        f"s={((mm*mm+nn*nn)%p)*inv2%p}, 2mn={2*mm*nn%p}")
    say(f"  QR_{p} (с нулём) = {sorted(qr_set(p))}")
    say(f"  проверка совпадения с таблицей Codex: m^2=16 n^2=13 s=6 2mn=4 -> "
        f"{(mm*mm%p, nn*nn%p, ((mm*mm+nn*nn)%p)*inv2%p, 2*mm*nn%p) == (16,13,6,4)}")
    say()
    say("  Полная таблица всех 18 проективных точек (пять клеток, ноль считается квадратом):")
    surv, rows = scan(mm, nn, p, zero_ok=True, nine=False, verbose=True)
    say(f"  ВЫЖИВШИХ при соглашении «ноль — квадрат, пять клеток»: {len(surv)} -> {surv}")

    # сколько точек проходит первые три квадратности (проверка фразы Codex про «ровно четыре»)
    QR = qr_set(p)
    first3 = [(a, b) for (a, b) in projective_points(p)
              if all(cells_mod(mm, nn, a, b, p)[k] in QR for k in ('F0', 'F4', 'F8'))]
    say(f"  проходят только F0,F4,F8 (кривая C7): {len(first3)} точек -> "
        f"{[(a if b==1 else 'oo') for (a,b) in first3]}")
    say("  (у Codex заявлено ровно четыре конечных t: 2, 8, 9, 15, и отказ на бесконечности)")

    say()
    say("  Другие соглашения (ищем дыру):")
    for zero_ok in (True, False):
        for nine in (False, True):
            s1, _ = scan(mm, nn, p, zero_ok=zero_ok, nine=nine)
            say(f"    ноль-квадрат={zero_ok!s:5s} клеток={9 if nine else 5}: выживших {len(s1)} {s1}")
    # квадратичный твист: cells = c * F_i, c — неквадрат
    nonres = min(x for x in range(2, p) if x not in qr_set(p))
    for nine in (False, True):
        s1, _ = scan(mm, nn, p, zero_ok=True, nine=nine, twist=nonres)
        say(f"    твист c={nonres} (неквадрат), клеток={9 if nine else 5}: выживших {len(s1)} {s1}")
    say("    (твист c нужен только для контроля: над Q множитель c обязан быть квадратом,")
    say("     иначе четыре автоматические клетки c*(mt±n)^2 не квадраты — а все четыре линейные")
    say("     формы одновременно в нуль обратиться не могут.)")
    return surv


# ----------------------------------------------------------------------------
# D. ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ
# ----------------------------------------------------------------------------
def rational_cells(mm, nn, t):
    s = Fr(mm * mm + nn * nn, 2)
    F4 = s * (1 + t * t)
    return {
        'F0': Fr(mm * mm) + Fr(nn * nn) * t * t,
        'F4': F4,
        'F8': Fr(nn * nn) + Fr(mm * mm) * t * t,
        'L': F4 - 2 * mm * nn * t,
        'U': F4 + 2 * mm * nn * t,
        '(mt+n)^2': (mm * t + nn) ** 2,
        '(mt-n)^2': (mm * t - nn) ** 2,
        '(m+nt)^2': (mm + nn * t) ** 2,
        '(m-nt)^2': (mm - nn * t) ** 2,
    }


def is_rat_square(x: Fr) -> bool:
    if x < 0:
        return False
    a, b = x.numerator, x.denominator
    return isqrt(a) ** 2 == a and isqrt(b) ** 2 == b


def section_D():
    say()
    say("=" * 78)
    say("D. ПОЛОЖИТЕЛЬНЫЕ КОНТРОЛИ: машина обязана НЕ исключать пары с реальной точкой")
    say("=" * 78)
    for (mm, nn, t0) in [(15, 8, Fr(1)), (15, 8, Fr(-1)), (7, 1, Fr(0))]:
        vals = rational_cells(mm, nn, t0)
        sq = {k: is_rat_square(v) for k, v in vals.items()}
        say(f"  (m,n)=({mm},{nn}), t={t0}: все девять клеток рациональные квадраты: {all(sq.values())}")
        say(f"     клетки: {[str(v) for v in vals.values()]}")
        if not all(sq.values()):
            say(f"     не квадраты: {[k for k,v in sq.items() if not v]}")
    for (mm, nn) in [(15, 8), (7, 1)]:
        for p in (17, 41, 5, 7, 11, 13, 19, 23):
            s1, _ = scan(mm, nn, p, zero_ok=True, nine=True)
            s5, _ = scan(mm, nn, p, zero_ok=True, nine=False)
            if len(s5) == 0:
                say(f"  !!! ТРЕВОГА: ({mm},{nn}) mod {p} пусто, хотя рациональная точка есть")
        say(f"  ({mm},{nn}): ни при одном из простых 5,7,11,13,17,19,23,41 пустоты нет — "
            f"механизм не исключает пары с точкой. ОК")
    # контроль на самой (13,8): при других простых выжившие должны быть
    say()
    say("  Контроль на (13,8) при других простых (если код систематически врёт, пусто будет везде):")
    for p in (5, 7, 11, 13, 19, 23, 29, 31, 37, 41, 43):
        s1, _ = scan(13, 8, p, zero_ok=True, nine=False)
        say(f"    p={p:3d}: выживших {len(s1):3d}")


# ----------------------------------------------------------------------------
# E. УРОВЕНЬ Q_17 (явный Гензель) — строго сильнее, чем F_17
# ----------------------------------------------------------------------------
def vp(x: Fr, p):
    if x == 0:
        return None
    num, den = x.numerator, x.denominator
    v = 0
    while num % p == 0:
        num //= p
        v += 1
    while den % p == 0:
        den //= p
        v -= 1
    return v


def is_sq_qp(x: Fr, p):
    v = vp(x, p)
    if v is None:
        return True
    if v % 2:
        return False
    y = x / Fr(p) ** v
    u = (y.numerator * pow(y.denominator, p - 2, p)) % p
    return pow(u, (p - 1) // 2, p) == 1


def q17_check(mm, nn, p, KMAX=6):
    """есть ли t в Z_p (и в окрестности бесконечности) со всеми пятью клетками — квадратами в Q_p.
       Возвращает True / False / None(не решено)."""
    def fs(a, b):
        s = Fr(mm * mm + nn * nn, 2)
        F4 = s * (Fr(a) ** 2 + Fr(b) ** 2)
        return [Fr(mm * mm * b * b + nn * nn * a * a), F4,
                Fr(nn * nn * b * b + mm * mm * a * a),
                F4 - 2 * mm * nn * a * b, F4 + 2 * mm * nn * a * b]

    def disk(r, k, inf):
        vals = fs(1, r) if inf else fs(r, 1)
        undec = False
        for x in vals:
            v = vp(x, p)
            if v is not None and v < k:
                if not is_sq_qp(x, p):
                    return False
            else:
                undec = True
        if not undec:
            return True
        if k >= KMAX:
            return None
        res = False
        for z in range(p):
            out = disk(r + z * p ** k, k + 1, inf)
            if out is True:
                return True
            if out is None:
                res = None
        return res

    r1 = disk(0, 0, False)
    if r1 is True:
        return True
    r2 = disk(0, 1, True)   # t = 1/u, u in p Z_p
    if r2 is True:
        return True
    if r1 is False and r2 is False:
        return False
    return None


def section_E():
    say()
    say("=" * 78)
    say("E. УРОВЕНЬ Q_17: явный подъём Гензеля (независимо от перебора по F_17)")
    say("=" * 78)
    for (mm, nn, p) in [(13, 8, 17), (15, 8, 17), (7, 1, 17), (16, 5, 41), (15, 8, 41)]:
        r = q17_check(mm, nn, p)
        word = {True: "ЕСТЬ точка", False: "ТОЧЕК НЕТ", None: "НЕ РЕШЕНО (не доказательство!)"}[r]
        say(f"  ({mm},{nn}) над Q_{p}: {word}")


# ----------------------------------------------------------------------------
# F. НЕЗАВИСИМЫЙ ПЕРЕБОР РАЦИОНАЛЬНЫХ t
# ----------------------------------------------------------------------------
def section_F(BOUND=260):
    say()
    say("=" * 78)
    say(f"F. ПРЯМОЙ ПЕРЕБОР РАЦИОНАЛЬНЫХ t = a/b, |a|,|b| <= {BOUND}, для (13,8)")
    say("=" * 78)
    best = (-1, None, None)
    cnt = 0
    for b in range(1, BOUND + 1):
        for a in range(-BOUND, BOUND + 1):
            if gcd(abs(a), b) != 1:
                continue
            cnt += 1
            t0 = Fr(a, b)
            vals = rational_cells(13, 8, t0)
            k = sum(1 for v in vals.values() if is_rat_square(v))
            if k > best[0]:
                best = (k, t0, vals)
            if k == 9:
                say(f"  !!! ОПРОВЕРЖЕНИЕ: t={t0} даёт девять квадратов")
    say(f"  перебрано несократимых t: {cnt}")
    say(f"  максимум квадратных клеток из девяти: {best[0]} при t={best[1]}")
    say(f"  (четыре клетки квадраты всегда, поэтому 4 — это «ничего не сошлось»)")
    say(f"  значения в лучшей точке: {[str(v) for v in best[2].values()]}")
    say("  ЭТО НЕ ДОКАЗАТЕЛЬСТВО ОТСУТСТВИЯ. Это лишь согласие с локальным выводом.")


# ----------------------------------------------------------------------------
# G. ПРОВЕРКА ЗАЯВЛЕННОГО СПИСКА ЗАПРЕЩЁННЫХ λ
# ----------------------------------------------------------------------------
def section_G():
    say()
    say("=" * 78)
    say("G. ПРОВЕРКА ОБОБЩЕНИЯ: запрещённые классы λ = m/n mod p")
    say("=" * 78)
    for p, claimed in [(17, {2, 3, 6, 8, 9, 11, 14, 15}), (41, {5, 8, 13, 19, 22, 28, 33, 36})]:
        found = set()
        for lam in range(p):
            s1, _ = scan(lam, 1, p, zero_ok=True, nine=False)
            if len(s1) == 0:
                found.add(lam)
        # n = 0 mod p  <=>  (m:n) = (1:0)
        s_inf, _ = scan(1, 0, p, zero_ok=True, nine=False)
        say(f"  p={p}: мой список запрещённых λ = {sorted(found)}")
        say(f"         заявлено Codex           = {sorted(claimed)}   совпадает: {found == claimed}")
        say(f"         случай n=0 mod p: выживших {len(s_inf)} (в список запретов не входит: "
            f"{len(s_inf) > 0})")
    lam_13_8 = 13 * pow(8, 15, 17) % 17
    say(f"  для (13,8): λ = 13/8 mod 17 = {lam_13_8}")


def main():
    M = section_A()
    okB = section_B(M)
    surv = section_C()
    section_D()
    section_E()
    section_F()
    section_G()

    say()
    say("=" * 78)
    say("ИТОГ СКЕПТИЧЕСКОЙ ПРОВЕРКИ")
    say("=" * 78)
    say(f"  клетки совпали с Codex: {okB}")
    say(f"  выживших точек P^1(F_17) для (13,8) при самом ПЕРМИССИВНОМ соглашении: {len(surv)}")
    say("  опровержения не найдено" if len(surv) == 0 else "  НАЙДЕНА ВЫЖИВШАЯ ТОЧКА — ОПРОВЕРЖЕНИЕ")

    with open('/home/kep/magicKube/bridge/local_13_8_p17_скептик.log', 'w', encoding='utf-8') as f:
        f.write("\n".join(OUT) + "\n")


if __name__ == '__main__':
    main()
