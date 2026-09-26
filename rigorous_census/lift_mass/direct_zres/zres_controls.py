#!/usr/bin/env python3
"""zres_controls.py -- контроли к zres_mass.py (Claude, 26.09.2026).
Блоки (запуск: python3 zres_controls.py A|B|C|D|E; лог zres_controls_<X>.log):
 A  движок image_Zp против прямого перебора: целые z in [0,p^K) и случайные дроби a/b (p∤b);
    классы клеток -- точные (class_of_int от точного целого). Проверка: перебор ⊆ движок (корректность)
    и движок ⊆ перебор (полнота, т.е. каждая пара образа реализуется конкретным z).
 B  полюс при p = 3 mod 4: z = u/p^n точными дробями -- восемь условий не выполняются никогда;
    исчерпывающе для (A,B) = (a/p^n, b/p^n) при p=7,11 (n=1,2 / n=1); отрицательный контроль p = 1 mod 4.
 C  контроли claude_p7: p=7 -- автоматический подъём (все наклоны); p=29,37 -- контрпримеры существуют;
    воспроизведение таблицы NOTE_P7 перебором (A,B) mod p^2.
 D  отрицательный контроль фильтра: без двух диагоналей (6 условий) -- выжившие должны появиться.
 E  наклоны 126/451, 73/362 (s>300, вне массовой проверки): фильтр по тем же простым.
Все циклы с явными границами; время каждого блока печатается."""
import random
import sys
import time
from fractions import Fraction
from math import gcd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import zres_core as Z
from zres_core import (GRID, LINES, T_CELLS, L_CELLS, chi_table, class_of_int, codex_lists,
                       image_Zp, lines_ok, primes_upto, tl)

OUT = __file__.rsplit("/", 1)[0]
P3 = [p for p in primes_upto(400) if p % 4 == 3 and p >= 7]
SLOPES = [(r, s) for s in range(2, 301) for r in range(1, s) if gcd(r, s) == 1]
NONTRIV = [sl for sl in SLOPES if any(len(x) > 1 for x in codex_lists(*sl))]
logf = None


def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    logf.write(msg + "\n")
    logf.flush()


def classes_exact_int(z, r, s, p, chi):
    cls = {}
    for c in GRID:
        x = 1 + (c[0] * r + c[1] * s) * z
        if x == 0:
            return None
        cls[c] = class_of_int(x, p, chi)
    return cls


def classes_exact_frac(num, den, r, s, p, chi):
    """z = num/den: клетка = (den + mu*num)/den."""
    cls = {}
    cden = class_of_int(den, p, chi)
    for c in GRID:
        x = den + (c[0] * r + c[1] * s) * num
        if x == 0:
            return None
        cls[c] = class_of_int(x, p, chi) ^ cden
    return cls


def special(p, r, s):
    """p делит одно из r, s, r±s, 2r±s, r±2s -- тогда в дереве бывают деления узлов."""
    return any(v % p == 0 for v in (r, s, r + s, r - s, 2 * r + s, 2 * r - s, r + 2 * s, r - 2 * s))


def block_A():
    rng = random.Random(20260926)
    ps = [7, 11, 19, 23, 31, 43, 47, 59, 67, 71]
    tot_z = tot_pass = 0
    n_cases = n_sound_bad = n_complete_bad = n_flag_bad = 0
    for p in ps:
        chi = chi_table(p)
        spec = [sl for sl in NONTRIV if special(p, *sl)]
        rest = [sl for sl in NONTRIV if not special(p, *sl)]
        sample = spec[:60] + rng.sample(rest, 60)
        t = time.time()
        for (r, s) in sample:
            img, info = image_Zp(p, r, s, chi)
            K = info["maxdepth"] + 2
            while p ** K > 30000 and K > 2:
                K -= 1
            brute = set()
            bflag = False
            for z in range(p ** K):
                cls = classes_exact_int(z, r, s, p, chi)
                tot_z += 1
                if cls is None or not lines_ok(cls):
                    continue
                tot_pass += 1
                brute.add(tl(cls))
                bflag = bflag or any(cls.values())
            for _ in range(300):            # случайные дроби a/b, p ∤ b
                b = rng.randrange(1, 10 ** 6)
                if b % p == 0:
                    continue
                a = rng.randrange(-10 ** 6, 10 ** 6)
                cls = classes_exact_frac(a, b, r, s, p, chi)
                tot_z += 1
                if cls is None or not lines_ok(cls):
                    continue
                tot_pass += 1
                brute.add(tl(cls))
                bflag = bflag or any(cls.values())
            n_cases += 1
            if not brute <= img:
                n_sound_bad += 1
                log("  НАРУШЕНИЕ КОРРЕКТНОСТИ p=%d r=%d s=%d brute=%s img=%s" % (p, r, s, brute, img))
            if not img <= brute:
                n_complete_bad += 1
                log("  неполнота перебора (K=%d) p=%d r=%d s=%d img-brute=%s" % (K, p, r, s, img - brute))
            if bflag != info["nontriv_cell"]:
                n_flag_bad += 1
                log("  флаг нетривиальной клетки расходится p=%d r=%d s=%d" % (p, r, s))
        log("  p=%d: %d наклонов (особых %d), %.1f c" % (p, len(sample), len(spec[:60]), time.time() - t))
    log("A ИТОГ: случаев (p,наклон) %d; z проверено %d, из них прошли 8 условий %d" % (n_cases, tot_z, tot_pass))
    log("A ИТОГ: перебор ⊄ движок: %d; движок ⊄ перебор: %d; флаг клетки расходится: %d"
        % (n_sound_bad, n_complete_bad, n_flag_bad))


def block_B():
    rng = random.Random(7)
    t = time.time()
    # (1) случайные наклоны и полюса, все p = 3 mod 4 из фильтра
    tested = passed = 0
    for p in P3:
        chi = chi_table(p)
        sl = [x for x in SLOPES if (x[0] % p == 0 or x[1] % p == 0 or (x[0] + x[1]) % p == 0
                                     or (x[1] - x[0]) % p == 0)][:40] + rng.sample(SLOPES, 60)
        for (r, s) in sl:
            for n in (1, 2, 3):
                for _ in range(8):
                    u = rng.randrange(1, p ** 4)
                    if u % p == 0:
                        continue
                    cls = classes_exact_frac(u, p ** n, r, s, p, chi)
                    tested += 1
                    if cls is not None and lines_ok(cls):
                        passed += 1
                        log("  ПОЛЮС ПРОШЁЛ p=%d r=%d s=%d z=%d/%d^%d" % (p, r, s, u, p, n))
    log("B1 p=3 mod 4 (все 39): полюсов z=u/p^n (n=1..3) проверено %d, прошли 8 условий %d  [%.1f c]"
        % (tested, passed, time.time() - t))
    # (2) исчерпывающе по (A,B) = (a/p^n, b/p^n), хотя бы одно из a,b -- p-единица
    for p, ns, M in ((7, (1, 2), 7 ** 3), (11, (1,), 11 ** 2 * 3), (19, (1,), 19 ** 2)):
        chi = chi_table(p)
        t = time.time()
        tested = passed = 0
        for n in ns:
            den = p ** n
            cden = class_of_int(den, p, chi)
            for a in range(M):
                for b in range(M):
                    if a % p == 0 and b % p == 0:
                        continue
                    ok = True
                    cls = {}
                    for c in GRID:
                        x = den + c[0] * a + c[1] * b
                        if x == 0:
                            ok = False
                            break
                        cls[c] = class_of_int(x, p, chi) ^ cden
                    tested += 1
                    if ok and lines_ok(cls):
                        passed += 1
        log("B2 p=%d исчерпывающе, n=%s, a,b in [0,%d): пар %d, прошли 8 условий %d  [%.1f c]"
            % (p, ns, M, tested, passed, time.time() - t))
    # (3) отрицательный контроль: p = 1 mod 4 -- полюсные решения ДОЛЖНЫ находиться
    for p in (5, 13, 17, 29, 37):
        chi = chi_table(p)
        t = time.time()
        tested = passed = 0
        ex = None
        M = p ** 2
        for n in (1, 2):
            den = p ** n
            cden = class_of_int(den, p, chi)
            for a in range(M):
                for b in range(0, M, 1 if p < 20 else 3):
                    if a % p == 0 and b % p == 0:
                        continue
                    cls = {}
                    ok = True
                    for c in GRID:
                        x = den + c[0] * a + c[1] * b
                        if x == 0:
                            ok = False
                            break
                        cls[c] = class_of_int(x, p, chi) ^ cden
                    tested += 1
                    if ok and lines_ok(cls):
                        passed += 1
                        if ex is None:
                            ex = (a, b, n)
        log("B3 отриц. контроль p=%d (=1 mod 4): полюсных (A,B) %d, прошли 8 условий %d, пример (a,b,n)=%s  [%.1f c]"
            % (p, tested, passed, ex, time.time() - t))


def block_C():
    t = time.time()
    # (a) p=7: для всех 27397 наклонов образ {(0,0)}, нетривиальных клеток нет (повторено и в zres_mass)
    chi7 = chi_table(7)
    bad = sum(1 for (r, s) in SLOPES if image_Zp(7, r, s, chi7) != ({(0, 0)}, image_Zp(7, r, s, chi7)[1])
              or image_Zp(7, r, s, chi7)[1]["nontriv_cell"])
    log("C1 p=7: наклонов с нетривиальным образом или нетривиальной клеткой: %d из %d  [%.1f c]" % (bad, len(SLOPES), time.time() - t))
    # (b) p=29,37 (и 5,13,17 для сравнения): у скольких наклонов 8 условий НЕ дают 9 квадратов в Z_p
    for p in (5, 13, 17, 29, 37):
        chi = chi_table(p)
        t = time.time()
        nt = sum(1 for (r, s) in SLOPES if image_Zp(p, r, s, chi)[1]["nontriv_cell"])
        log("C2 p=%d: наклонов, где в Z_p есть z с 8 условиями, но не все 9 клеток квадраты: %d из %d  [%.1f c]"
            % (p, nt, len(SLOPES), time.time() - t))
    # (c) воспроизведение таблицы NOTE_P7 (целый z): все (A,B) mod p^2, клетки 1+iA+jB mod p^2,
    #     конфигурация пропускается, если клетка = 0 mod p^2; класс = (v, chi(ед. часть)), v in {0,1}
    for p in (13, 17, 29, 37):
        chi = chi_table(p)
        M = p * p
        t = time.time()
        n8 = n9 = skipped = 0
        for A in range(M):
            for B in range(M):
                cls = {}
                ok = True
                for c in GRID:
                    x = (1 + c[0] * A + c[1] * B) % M
                    if x == 0:
                        ok = False
                        break
                    cls[c] = class_of_int(x, p, chi)
                if not ok:
                    skipped += 1
                    continue
                if lines_ok(cls):
                    n8 += 1
                    if not any(cls.values()):
                        n9 += 1
        log("C3 p=%d, (A,B) mod p^2: с 8 квадратами %d, из них все 9 %d, контрпримеров %d; пропущено %d  [%.1f c]"
            % (p, n8, n9, n8 - n9, skipped, time.time() - t))


def filter_with_lines(lines, label):
    t = time.time()
    surv_slopes = 0
    surv_pairs = 0
    chis = {p: chi_table(p) for p in P3}
    for (r, s) in NONTRIV:
        DT, DL = codex_lists(r, s)
        pairs = [(a, b) for a in DT for b in DL if (a, b) != (1, 1)]
        alive = set(pairs)
        for p in P3:
            chi = chis[p]
            saved = Z.LINES
            Z.LINES = lines          # тот же движок с подменённым списком линий
            try:
                img, _ = image_Zp(p, r, s, chi)
            finally:
                Z.LINES = saved
            alive = {pr for pr in alive if (class_of_int(pr[0], p, chi), class_of_int(pr[1], p, chi)) in img}
            if not alive:
                break
        if alive:
            surv_slopes += 1
            surv_pairs += len(alive)
    log("D %s: наклонов с выжившей нетривиальной парой %d из %d, пар %d  [%.1f c]"
        % (label, surv_slopes, len(NONTRIV), surv_pairs, time.time() - t))


def block_D():
    # (i) без двух диагоналей (6 условий); (ii) только 4 линии через центр (средние строка/столбец, диагонали)
    filter_with_lines(LINES[:6], "без диагоналей (6 условий)")
    central = [LINES[1], LINES[4], LINES[6], LINES[7]]
    assert all((0, 0) in l for l in central)
    filter_with_lines(central, "только 4 центральные линии")


def block_E():
    for (r, s) in ((126, 451), (73, 362), (265, 298)):
        DT, DL = codex_lists(r, s)
        pairs = [(a, b) for a in DT for b in DL if (a, b) != (1, 1)]
        ex = {pr: [] for pr in pairs}
        for p in P3:
            chi = chi_table(p)
            img, _ = image_Zp(p, r, s, chi)
            for pr in pairs:
                if (class_of_int(pr[0], p, chi), class_of_int(pr[1], p, chi)) not in img:
                    ex[pr].append(p)
        log("E %d/%d: D_-=%s D_+=%s; исключающие простые по парам: %s" % (r, s, DT, DL, ex))


if __name__ == "__main__":
    which = sys.argv[1]
    logf = open(OUT + "/zres_controls_%s.log" % which, "w")
    log("zres_controls.py блок", which, time.strftime("%Y-%m-%d %H:%M:%S"))
    t0 = time.time()
    {"A": block_A, "B": block_B, "C": block_C, "D": block_D, "E": block_E}[which]()
    log("время блока %.1f c" % (time.time() - t0))
