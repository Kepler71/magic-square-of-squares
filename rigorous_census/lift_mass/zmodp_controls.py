#!/usr/bin/env python3
"""zmodp_controls.py -- контроли к zmodp_mass.py.

C1  нет полюса при p = 3 mod 4: выборка рациональных z с v_p(z) < 0, точная арифметика Fraction;
    отрицательный контроль при p = 1 mod 4 (полюса с квадратными центральными линиями находятся).
C2  движок против РЕАЛИЗОВАННЫХ классов: для всех 2760 наклонов при p = 7 и p = 11 перебираются целые
    z in [0, p^3) плюс целые z в глубоких окрестностях нулей клеток (z = -1/mu mod p^k, k <= 6);
    классы считаются по настоящему целому 1+mu*z, восемь условий -- по классу целого произведения
    (без XOR-кодов движка). Требуется: реализованное множество векторов ⊆ образ движка (это то
    направление, на котором держится исключение), и отчёт, совпадают ли они.
    То же для 40 случайных нетривиальных наклонов при всех p in {19,23,31,43,47,59}.
C3  p = 3 и p = 7 по всем 27397 наклонам: образ движка = {нулевой вектор} (автоматический подъём;
    §4 Codex для p=3, §5 сводки и claude_p7 для p=7).
C4  воспроизведение чисел claude_p7 (перебор (A,B) mod p^2 без наклона, целый z) при p = 13,17,29,37
    и p = 7 по mod 7^3 -- другой, бесхитростный код.
C5  образ при p = 29, 37 (целая часть) нетривиален у части наклонов -- движок не «всегда тривиален».
C6  роли: при p = 1 mod 4 нечётная оценка T бывает только при p | (r-s), L -- только при p | (r+s)
    (выборка полюсов, точная арифметика) -- проверка, что списки D_- относятся к T, D_+ к L.
C7  отрицательные контроли фильтра: роли T/L переставлены; выброшена одна диагональ; только p=7.
Все циклы конечны (явные границы), прогресс -- в stdout.
"""
import random
import sys
import time
from fractions import Fraction
from math import gcd

from zmodp_common import (CELLS, IDX, LINES8, LINES8_IDX, T_IDX, L_IDX, codex_list, local_image_Zp,
                          primes_3mod4, qr_table, qp_class_rational, tl_of, unit_code, cell_mus)
import zmodp_common as ZC

T0 = time.time()
P = primes_3mod4(7, 400)
random.seed(20260926)


def log(*a):
    print("[%7.1f c]" % (time.time() - T0), *a, flush=True)


def int_class(n, p, qr):
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return (v & 1) | (0 if qr[n % p] else 2)


SLOPES = [(r, s) for s in range(2, 301) for r in range(1, s) if gcd(r, s) == 1]
NONTRIV = [(r, s) for (r, s) in SLOPES if len(codex_list(s - r)) > 1 or len(codex_list(r + s)) > 1]


# ---------------------------------------------------------------- C1
def c1():
    log("C1: полюс при p = 3 mod 4")
    bad = 0
    checked = 0
    for p in [3] + P:
        for trial in range(200):
            r, s = random.choice(SLOPES)
            n = random.randint(1, 3)
            u = random.randint(1, 10 ** 6)
            while u % p == 0:
                u = random.randint(1, 10 ** 6)
            w = random.randint(1, 10 ** 4)
            while w % p == 0:
                w = random.randint(1, 10 ** 4)
            z = Fraction(u * random.choice((1, -1)), p ** n * w)
            c_row = qp_class_rational(1 - s * s * z * z, p)   # средняя строка i=0: (1-sz)*1*(1+sz)
            c_col = qp_class_rational(1 - r * r * z * z, p)   # средний столбец j=0
            checked += 1
            # утверждение: для p-единицы lambda класс 1-lambda^2 z^2 равен ровно классу -1 (код 2)
            for lam, c in ((s, c_row), (r, c_col)):
                if lam % p != 0 and c != 2:
                    bad += 1
            if c_row == 0 and c_col == 0:
                bad += 1
    log("   p=3 и 39 простых из P, %d случайных полюсов: нарушений %d (ожидается 0)" % (checked, bad))
    # отрицательный контроль: p = 1 mod 4
    for p in (5, 13, 29, 37):
        found = 0
        for trial in range(2000):
            r, s = random.choice(SLOPES)
            u = random.randint(1, 10 ** 6)
            if u % p == 0:
                continue
            z = Fraction(u, p)
            if qp_class_rational(1 - s * s * z * z, p) == 0 and qp_class_rational(1 - r * r * z * z, p) == 0:
                found += 1
        log("   отр. контроль p=%d (1 mod 4): полюсов с квадратными средней строкой и столбцом: %d из 2000"
            % (p, found))
    return bad


# ---------------------------------------------------------------- C2
def realized(r, s, p, zs):
    qr = qr_table(p)
    mus = cell_mus(r, s)
    out = set()
    for z in zs:
        vals = [1 + m * z for m in mus]
        if 0 in vals:
            continue
        ok = True
        for ln in LINES8:
            a, b, c = (vals[IDX[x]] for x in ln)
            if int_class(a * b * c, p, qr) != 0:
                ok = False
                break
        if ok:
            out.add(tuple(int_class(v, p, qr) for v in vals))
    return out


def z_samples(r, s, p, M):
    zs = list(range(p ** M))
    for m in set(cell_mus(r, s)):
        if m % p == 0:
            continue
        for k in range(M, 7):
            pk = p ** k
            z0 = (-pow(m, -1, pk)) % pk
            zs.extend(z0 + t * pk for t in range(p))
    return zs


def c2():
    log("C2: движок против реализованных целыми z")
    bad = 0
    stats = {}
    for p in (7, 11):
        eq = sub = 0
        for n, (r, s) in enumerate(NONTRIV, 1):
            E = local_image_Zp(r, s, p)
            R = realized(r, s, p, z_samples(r, s, p, 3))
            if not R <= E:
                bad += 1
                log("   НАРУШЕНИЕ: реализовано вне образа движка", r, s, p, sorted(R - E)[:5])
            if R == E:
                eq += 1
            else:
                sub += 1
            if n % 690 == 0:
                log("   p=%d: %d/%d" % (p, n, len(NONTRIV)))
        stats[p] = (eq, sub)
        log("   p=%d: все 2760 наклонов: R ⊆ E везде? %s; R = E у %d, R ⊊ E у %d" % (p, bad == 0, eq, sub))
    sample = random.sample(NONTRIV, 40)
    for p in (19, 23, 31, 43, 47, 59):
        eq = sub = 0
        nontrivial_R = 0
        M = 2 if p > 23 else 3
        for (r, s) in sample:
            E = local_image_Zp(r, s, p)
            R = realized(r, s, p, z_samples(r, s, p, M))
            if not R <= E:
                bad += 1
                log("   НАРУШЕНИЕ", r, s, p, sorted(R - E)[:5])
            eq += R == E
            sub += R != E
            nontrivial_R += any(v != (0,) * 9 for v in R)
        log("   p=%d (40 наклонов, z < p^%d + окрестности нулей): R = E у %d, R ⊊ E у %d; наклонов с нетрив. "
            "реализованным вектором: %d" % (p, M, eq, sub, nontrivial_R))
    # явная проверка ключевого факта: наклоны, у которых p=7 недостаточно, убиты при p=11
    return bad


# ---------------------------------------------------------------- C3
def c3():
    log("C3: p=3 и p=7, все 27397 наклонов")
    bad = {3: 0, 7: 0}
    zero = (0,) * 9
    for p in (3, 7):
        for (r, s) in SLOPES:
            E = local_image_Zp(r, s, p)
            if E != {zero}:
                bad[p] += 1
        log("   p=%d: наклонов с образом != {все клетки квадраты}: %d (ожидается 0)" % (p, bad[p]))
    return bad[3] + bad[7]


# ---------------------------------------------------------------- C4
def c4():
    log("C4: слепое к наклону воспроизведение чисел claude_p7 (A=rz, B=sz произвольны mod p^M)")
    expect = {13: (845, None), 17: (1445, None), 29: (14297, 10933), 37: (39701, 23273), 7: (None, None)}
    res = {}
    for p, M in ((13, 2), (17, 2), (29, 2), (37, 2), (7, 3)):
        qr = qr_table(p)
        pm = p ** M
        n8 = n9 = 0
        for A in range(pm):
            for B in range(pm):
                vals = [(1 + i * A + j * B) % pm for (i, j) in CELLS]
                if 0 in vals:
                    continue
                cl = [int_class(v, p, qr) for v in vals]   # v < p^M, так что v_p < M -- корректно
                if all(cl[a] ^ cl[b] ^ cl[c] == 0 for (a, b, c) in LINES8_IDX):
                    n8 += 1
                    if not any(cl):
                        n9 += 1
        res[p] = (n8, n9)
        log("   p=%d mod p^%d: 8 квадратов: %d, из них все 9: %d, контрпримеров %d; claude_p7: %s"
            % (p, M, n8, n9, n8 - n9, expect[p]))
    return res


# ---------------------------------------------------------------- C5
def c5():
    log("C5: целая часть образа при p = 29, 37 (p = 1 mod 4) -- нетривиальна?")
    for p in (29, 37):
        cnt = 0
        for (r, s) in NONTRIV[:600]:
            E = local_image_Zp(r, s, p)
            if any(tl_of(v) != (0, 0) for v in E):
                cnt += 1
        log("   p=%d: из первых 600 нетривиальных наклонов с нетривиальным ([T],[L]) в Z_p-части: %d" % (p, cnt))


# ---------------------------------------------------------------- C6
def c6():
    log("C6: роли T <-> r-s, L <-> r+s (полюса при p = 1 mod 4, точная арифметика)")
    viol = 0
    seen = {"T_odd": 0, "L_odd": 0}
    cases = []
    for p in (5, 13, 17, 29):
        cases += [(r, s, p) for (r, s) in SLOPES if (s - r) % p == 0 or (r + s) % p == 0 or random.random() < 0.02]
    random.shuffle(cases)
    cases = cases[:3000]
    for (r, s, p) in cases:
        mus = cell_mus(r, s)
        for trial in range(60):
            n = random.randint(1, 3)
            u = random.randint(1, p ** 5)
            if u % p == 0:
                continue
            z = Fraction(u, p ** n)
            vals = [1 + m * z for m in mus]
            if any(v == 0 for v in vals):
                continue
            if any(qp_class_rational(vals[a] * vals[b] * vals[c], p) != 0 for (a, b, c) in LINES8_IDX):
                continue
            vT = qp_class_rational(vals[T_IDX[0]] * vals[T_IDX[1]] * vals[T_IDX[2]], p) & 1
            vL = qp_class_rational(vals[L_IDX[0]] * vals[L_IDX[1]] * vals[L_IDX[2]], p) & 1
            if vT:
                seen["T_odd"] += 1
                if (s - r) % p != 0:
                    viol += 1
            if vL:
                seen["L_odd"] += 1
                if (r + s) % p != 0:
                    viol += 1
    log("   нечётная v_p(T) встречена %d раз, v_p(L) -- %d раз; нарушений (не делит r-s / r+s): %d"
        % (seen["T_odd"], seen["L_odd"], viol))
    return viol


# ---------------------------------------------------------------- C7
def filter_survivors(TLfun, primes):
    """Сколько наклонов из 2760 сохраняют пару кроме (1,1)."""
    bad = 0
    for (r, s) in NONTRIV:
        Dm, Dp = codex_list(s - r), codex_list(r + s)
        TL = {p: TLfun(r, s, p) for p in primes}
        surv = [(a, b) for a in Dm for b in Dp
                if all((unit_code(a, p), unit_code(b, p)) in TL[p] for p in primes)]
        if surv != [(1, 1)]:
            bad += 1
    return bad


def c7():
    log("C7: отрицательные контроли фильтра (должны дать выживших)")
    std = lambda r, s, p: set(tl_of(v) for v in local_image_Zp(r, s, p))
    swapped = lambda r, s, p: set((l, t) for (t, l) in std(r, s, p))
    log("   только p=7: наклонов с выжившими: %d (ожидается 504)" % filter_survivors(std, [7]))
    log("   p in {7,11}: %d (по основному прогону ожидается 0)" % filter_survivors(std, [7, 11]))
    log("   роли T/L переставлены, все 39 p: %d" % filter_survivors(swapped, P))
    saved = ZC.LINES8_IDX
    ZC.LINES8_IDX = saved[:7]            # выброшена побочная диагональ
    try:
        log("   без побочной диагонали, все 39 p: %d" % filter_survivors(std, P))
    finally:
        ZC.LINES8_IDX = saved
    log("   контроль восстановления: все 39 p: %d (ожидается 0)" % filter_survivors(std, P))


def main():
    which = sys.argv[1:] or ["c1", "c2", "c3", "c4", "c5", "c6", "c7"]
    log("наклонов %d, нетривиальных %d" % (len(SLOPES), len(NONTRIV)))
    for w in which:
        globals()[w]()
    log("готово")


if __name__ == "__main__":
    main()
