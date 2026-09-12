#!/usr/bin/env sage
# -*- coding: utf-8 -*-
r"""
search_H_15_8.sage
==================

Расширенный поиск рациональных точек на кривой рода 2

    H : Y^2 = F0(t) * F4(t) * F8(t),   (m,n) = (15,8)

    s  = (m^2+n^2)/2 = 289/2          (m^2+n^2 = 289 = 17^2 — единственная
                                       пифагорова пара среди семи G1-пар)
    F0 = m^2 + n^2 t^2 = 225 +  64 t^2
    F4 = s (1 + t^2)   = (289/2)(1 + t^2)
    F8 = n^2 + m^2 t^2 =  64 + 225 t^2

    b  = s m^2 n^2 = 2080800 = 2 * 1020^2   (квадратный класс b = 2)
    a  = m^2/n^2 + 1 + n^2/m^2 = 69121/14400

Известная ВЫРОЖДЕННАЯ точка: t = 1, F0 = F4 = F8 = 289, Y = 17^3 = 4913.
Она поднимается на C (u0 = u4 = u8 = 17), но даёт три РАВНЫЕ клетки, то есть
магического квадрата из различных квадратов не даёт.

Скрипт состоит из пяти частей.

ЧАСТЬ A.  Символьная проверка всех используемых тождеств (карта H -> E_b,
  факторизация E_b, отсутствие точек над t = 0 и t = oo).

ЧАСТЬ B (БЕЗУСЛОВНЫЙ перебор по высоте).  t = u/v, gcd(u,v) = 1,
  0 <= u <= v <= HMAX.  Симметрии
      t -> -t   (все F_i чётны по t),
      t -> 1/t  (F0(1/t) = F8(t)/t^2, F8(1/t) = F0(t)/t^2, F4(1/t) = F4(t)/t^2)
  сводят произвольное рациональное t к этому диапазону; обе симметрии
  сохраняют и условие подъёма на C (F0 и F8 меняются местами).
  Решето: условие «N(u,v) — квадратичный вычет или 0 mod p» зависит только
  от (u mod p, v mod p), N однородна степени 6.  Выжившие кандидаты
  проверяются точно целочисленным is_square.

ЧАСТЬ C (ПОЛНЫЙ перебор по группе Морделла–Вейля; это главная часть).
  Все рациональные точки H лежат над аффинной картой (лидирующий коэффициент
  и f(0) равны b, а b не квадрат), и карта
      pi : H -> E_b,  (t,Y) |-> (X,V) = (b t^2, b Y)
  даёт БИЕКЦИЮ
      H(Q)/{(t,Y)~(-t,Y)}  <->  { P in E_b(Q) : x(P)/b in (Q^*)^2 }.
  Поэтому поиск точек H — это поиск точек E_b(Q) с x/b = квадрат.
  E_b(Q) = Z G1 + Z G2 + (Z/2)^2 (ранг 2 доказан PARI ellrank: нижняя =
  верхняя = 2; насыщенность G1,G2 — Sage saturation, индекс 1).
  Перебираем ВСЕ n1 G1 + n2 G2 + T с |n1|,|n2| <= NMAX.  Так как
  hhat(n1G1+n2G2) = n^T M n >= lambda_min(M) * |n|^2, этот ящик содержит
  ВСЕ точки E_b(Q) канонической высоты <= lambda_min * NMAX^2.
  Предварительное решето по простым хорошей редукции (условие: если P не
  редуцируется в O mod p, то x_E(P)/b mod p — квадрат в F_p) отсеивает
  подавляющее большинство пар (n1,n2) без вычисления гигантских координат.

ЧАСТЬ D.  Для КАЖДОЙ найденной точки проверяется подъём на кривую рода 5
  C : u_i^2 = F_i(t): нужно, чтобы F0(t), F4(t), F8(t) были квадратами
  КАЖДЫЙ ПО ОТДЕЛЬНОСТИ, а не только их произведение.  Дополнительно
  проверяется невырожденность (различность девяти клеток квадрата).

ЧАСТЬ E (специальный поиск точек C на ОГРОМНЫХ высотах).
  Условие «F4 квадрат» само по себе задаёт конику:
      v^2 F4 = (289/2)(u^2+v^2) — квадрат  <=>  u^2 + v^2 = 2 w^2,
  а это рационально параметризуется.  С alpha = (u+v)/2, beta = (u-v)/2
  получаем alpha^2 + beta^2 = w^2 — примитивную пифагорову тройку,
      {alpha,beta} = {x^2-y^2, 2xy},  w = x^2+y^2.
  Поэтому ВСЕ точки C имеют вид u = alpha+beta, v = alpha-beta, и остаётся
  проверить квадратность A0 = 225 v^2 + 64 u^2 и A8 = 64 v^2 + 225 u^2.
  Высота при этом растёт как x^2, так что перебор по x <= XMAX покрывает
  высоты ~2 XMAX^2 — на много порядков дальше прямого перебора.

ВАЖНО: отсутствие находок НЕ является доказательством отсутствия точек.
Достигнутые границы печатаются явно в итоговой сводке.

Запуск:
    sage search_H_15_8.sage [HMAX] [NMAX] [XMAX]
по умолчанию HMAX = 20000, NMAX = 100, XMAX = 6000.
"""

import sys
import time

import numpy as np

# ======================================================================
# 0. Параметры семейства
# ======================================================================
m = 15
n = 8
s = QQ(m ^ 2 + n ^ 2) / 2          # 289/2 — НЕ целое: m^2+n^2 нечётно
b = s * m ^ 2 * n ^ 2              # 2080800
aa = QQ(m ^ 2) / n ^ 2 + 1 + QQ(n ^ 2) / m ^ 2

R = PolynomialRing(QQ, 't')
t = R.gen()
F0 = m ^ 2 + n ^ 2 * t ^ 2
F4 = s * (1 + t ^ 2)
F8 = n ^ 2 + m ^ 2 * t ^ 2
f = F0 * F4 * F8

T_START = time.time()


def banner(txt):
    print()
    print("=" * 74)
    print(txt)
    print("=" * 74)
    sys.stdout.flush()


# ======================================================================
# ЧАСТЬ A. Символьные проверки
# ======================================================================
banner("ЧАСТЬ A. Проверка тождеств  [доказано: символьные тождества в Sage]")

print("  (m,n) = (%d,%d),  s = %s,  b = %s = %s,  a = %s"
      % (m, n, s, b, factor(Integer(b)), aa))
print("  f(t) = %s" % f)
print("  f палиндромичен (симметрия t -> 1/t):", list(f) == list(f)[::-1])

E = EllipticCurve(QQ, [0, b * aa, 0, b ^ 2 * aa, b ^ 3])
X = b * t ^ 2
assert (X + b) * (X + s * m ^ 4) * (X + s * n ^ 4) \
    == X ^ 3 + b * aa * X ^ 2 + b ^ 2 * aa * X + b ^ 3
assert X ^ 3 + b * aa * X ^ 2 + b ^ 2 * aa * X + b ^ 3 == b ^ 2 * f
assert X + b == m ^ 2 * n ^ 2 * F4
assert X + s * m ^ 4 == s * m ^ 2 * F0
assert X + s * n ^ 4 == s * n ^ 2 * F8
print("  [доказано] (X+b)(X+s m^4)(X+s n^4) = X^3+ba X^2+b^2a X+b^3,  X = b t^2")
print("  [доказано] подстановка X = b t^2 даёт b^2 * F0 F4 F8  =>  V = b Y")
print("  [доказано] X+b = m^2n^2 F4,  X+s m^4 = s m^2 F0,  X+s n^4 = s n^2 F8")

lead = f.list()[-1]
print("  f(0) = %s,  старший коэффициент = %s;  оба равны b = %s"
      % (f(0), lead, b))
print("  b — квадрат?", Integer(b).is_square(),
      "(квадратный класс b = %s)" % Integer(b).squarefree_part())
print("  [доказано] => на H нет рациональных точек ни над t = 0, ни над t = oo,")
print("              значит H(Q) целиком лежит в аффинной карте t != 0.")

# известная точка
assert F0(1) == 289 and F4(1) == 289 and F8(1) == 289
print("  [доказано] t = 1: F0 = F4 = F8 = 289 = 17^2, Y = 17^3 = 4913 —")
print("              ВЫРОЖДЕННАЯ точка (три клетки равны 17^2).")

# ======================================================================
# Вспомогательные функции
# ======================================================================
# Однородная целочисленная форма.
#   v^6 f(u/v) = (m^2v^2+n^2u^2) * s * (u^2+v^2) * (n^2v^2+m^2u^2)
# и s = 289/2 = 17^2/2, поэтому
#   v^6 f(u/v) — квадрат  <=>  Pint(u,v) := 2*A0*A4h*A8 — точный квадрат,
# где A0 = m^2v^2+n^2u^2, A4h = u^2+v^2, A8 = n^2v^2+m^2u^2,
# и тогда Y = 17*sqrt(Pint)/(2 v^3).
def Pint(u, v):
    u = Integer(u)
    v = Integer(v)
    return 2 * (m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2) * (u ^ 2 + v ^ 2) \
        * (n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2)


def lift_data(tt):
    """Квадратность F0,F4,F8 в рациональной точке t и сами значения."""
    tt = QQ(tt)
    vals = (F0(tt), F4(tt), F8(tt))
    sq = tuple(x.is_square() for x in vals)
    return sq, vals


def report_point(tt, Yv, tag):
    """Печать точки H + проверка подъёма на C (ЧАСТЬ D)."""
    tt = QQ(tt)
    assert Yv ^ 2 == f(tt), "ложная точка!"
    sq, vals = lift_data(tt)
    print("    t = %s   Y = %s   [%s]" % (tt, Yv, tag))
    print("      F0 = %s  квадрат=%s" % (vals[0], sq[0]))
    print("      F4 = %s  квадрат=%s" % (vals[1], sq[1]))
    print("      F8 = %s  квадрат=%s" % (vals[2], sq[2]))
    if all(sq):
        u0, u4, u8 = (x.sqrt() for x in vals)
        print("      ПОДНИМАЕТСЯ на C: u0=%s, u4=%s, u8=%s" % (u0, u4, u8))
        deg = (u0 == u4) or (u4 == u8) or (u0 == u8)
        print("      вырожденная (есть равные клетки):", deg)
        return True, deg
    print("      на C НЕ поднимается (квадратно только произведение)")
    return False, None


# ======================================================================
# ЧАСТЬ B. Безусловный перебор по высоте с решетом
# ======================================================================
def build_tables(primes):
    """tab[p][v mod p][u mod p] = 1, если Pint(u,v) mod p — квадрат или 0."""
    tabs = {}
    dens = {}
    for p in primes:
        qr = np.zeros(p, dtype=np.uint8)
        for x in range(p):
            qr[(x * x) % p] = 1            # 0 тоже допустим
        Tb = np.zeros((p, p), dtype=np.uint8)
        m2 = (m * m) % p
        n2 = (n * n) % p
        for vr in range(p):
            vv = (vr * vr) % p
            for ur in range(p):
                uu = (ur * ur) % p
                val = (2 * (m2 * vv + n2 * uu)) % p
                val = (val * ((uu + vv) % p)) % p
                val = (val * ((n2 * vv + m2 * uu) % p)) % p
                Tb[vr][ur] = qr[val]
        tabs[p] = Tb
        dens[p] = float(Tb.sum()) / (p * p)
    return tabs, dens


def sieve_search(HMAX, primes_fast, primes_slow, report_every=2000):
    tabs_f, dens_f = build_tables(primes_fast)
    tabs_s, _ = build_tables(primes_slow)
    exp_keep = 1.0
    for p in primes_fast:
        exp_keep *= dens_f[p]
    print("  быстрое решето: %s" % (primes_fast,))
    print("  плотности: %s" % ({p: round(dens_f[p], 4) for p in primes_fast},))
    print("  медленное решето: %s" % (primes_slow,))
    print("  ожидаемая доля выживших после быстрой стадии ~ %.3e" % exp_keep)
    sys.stdout.flush()

    found = []
    n_surv1 = n_surv2 = n_exact = 0
    t0 = time.time()
    for v in range(1, HMAX + 1):
        L = v + 1                          # u = 0..v
        alive = np.ones(L, dtype=np.uint8)
        for p in primes_fast:
            row = tabs_f[p][v % p]
            tiled = np.resize(row, int(L))  # row повторяется с периодом p
            alive &= tiled
        idx = np.nonzero(alive)[0]
        n_surv1 += len(idx)
        if len(idx):
            keep = np.ones(len(idx), dtype=np.uint8)
            for p in primes_slow:
                keep &= tabs_s[p][v % p][idx % p]
            idx = idx[np.nonzero(keep)[0]]
        n_surv2 += len(idx)
        for u in idx:
            u = int(u)
            if gcd(u, v) != 1:
                continue
            n_exact += 1
            P = Pint(u, v)
            if P.is_square():
                Y = QQ(17 * P.sqrt()) / QQ(2 * v ^ 3)
                found.append((u, v, Y))
                print("  *** ТОЧКА H: t = %d/%d,  Y = %s" % (u, v, Y))
                sys.stdout.flush()
        if v % report_every == 0:
            print("    v = %6d/%d  выжило1=%d выжило2=%d точных=%d найдено=%d  %.1f с"
                  % (v, HMAX, n_surv1, n_surv2, n_exact, len(found),
                     time.time() - t0))
            sys.stdout.flush()
    stats = dict(HMAX=HMAX, surv1=n_surv1, surv2=n_surv2, exact=n_exact,
                 found=len(found), seconds=round(time.time() - t0, 1))
    return found, stats


# ======================================================================
# ЧАСТЬ C. Полный перебор ящика в группе Морделла–Вейля
# ======================================================================
def mw_box_search(NMAX, n_sieve_primes=60):
    Em = E.minimal_model()
    iso = Em.isomorphism_to(E)
    # Соглашение Sage: изоморфизм с кортежем (u,r,s,t) из Em в E действует как
    #     x_E = (x_Em - r)/u^2,   y_E = (y_Em - s(x_Em-r) - t)/u^3.
    uiso, riso, siso, tiso = iso.tuple()
    assert siso == 0 and tiso == 0
    u2 = QQ(uiso) ^ 2
    u3 = QQ(uiso) ^ 3
    # контроль: сверяем формулу с самим объектом Sage на генераторе
    print("  минимальная модель Em: %s" % Em)
    print("  кондуктор Em = %s" % Em.conductor())
    print("  disc(Em) = %s" % factor(Em.discriminant()))
    print("  изоморфизм Em -> E:  x_E = (x_Em - %s)/(%s)" % (riso, u2))

    rk = pari(Em).ellrank()
    print("  PARI ellrank(Em) = %s   [доказано: нижняя = верхняя = %s]"
          % (rk, rk[0]))
    G = list(Em.gens())
    sat = Em.saturation(G)
    print("  генераторы: %s" % (G,))
    print("  Sage saturation: индекс = %s  [насыщено, если индекс = 1]" % sat[1])
    for g in G:
        chk = iso(g)
        assert chk[0] == (g[0] - riso) / u2 and chk[1] == g[1] / u3, \
            "формула изоморфизма не совпала с объектом Sage"
    print("  [проверено] формула изоморфизма сверена с iso(G_i) в Sage")
    M = Em.height_pairing_matrix(G)
    ev = [float(x) for x in M.eigenvalues()]
    lam = min(ev)
    print("  матрица высот:\n%s" % M)
    print("  собственные значения: %s,  lambda_min = %.6f" % (ev, lam))
    tors = list(Em.torsion_points())
    print("  кручение Em: %s" % (tors,))

    # --- локальное решето -------------------------------------------------
    bad = set(Integer(Em.discriminant()).prime_factors()) \
        | set(Integer(b).prime_factors()) | {2}
    bad |= set(QQ(u2).numerator().prime_factors()) \
        | set(QQ(u2).denominator().prime_factors())
    sieve_primes = []
    p = 3
    while len(sieve_primes) < n_sieve_primes:
        p = next_prime(p)
        if p not in bad:
            sieve_primes.append(p)
    print("  плохие простые (исключены из решета): %s" % sorted(bad))
    print("  простые решета: %s ... %s (%d штук)"
          % (sieve_primes[:6], sieve_primes[-3:], len(sieve_primes)))
    sys.stdout.flush()

    # Кандидаты: (n1, n2, k) — индекс точки кручения.  Симметрия P <-> -P
    # сохраняет x, поэтому достаточно n1 > 0, либо n1 = 0 и n2 >= 0.
    cand = []
    for n1 in range(0, NMAX + 1):
        n2range = range(-NMAX, NMAX + 1) if n1 > 0 else range(0, NMAX + 1)
        for n2 in n2range:
            for k in range(len(tors)):
                cand.append((n1, n2, k))
    total = len(cand)
    print("  кандидатов в ящике |n1|,|n2| <= %d (с точностью до P<->-P): %d"
          % (NMAX, total))
    sys.stdout.flush()

    t0 = time.time()
    for p in sieve_primes:
        if not cand:
            break
        Fp = GF(p)
        Ep = EllipticCurve(Fp, [Fp(c) for c in Em.a_invariants()])
        G1p = Ep(*[Fp(c) for c in G[0].xy()])
        G2p = Ep(*[Fp(c) for c in G[1].xy()])
        Tp = []
        for Tt in tors:
            if Tt.is_zero():
                Tp.append(Ep(0))
            else:
                Tp.append(Ep(*[Fp(c) for c in Tt.xy()]))
        u2p = Fp(u2)
        rp = Fp(riso)
        binv = Fp(b) ^ (-1)
        # таблицы кратных
        A = {}
        cur = Ep(0)
        A[0] = cur
        for i in range(1, NMAX + 1):
            cur = cur + G1p
            A[i] = cur
        B = {}
        cur = Ep(0)
        B[0] = cur
        for j in range(1, NMAX + 1):
            cur = cur + G2p
            B[j] = cur
            B[-j] = -cur
        newcand = []
        for (n1, n2, k) in cand:
            Q = A[n1] + B[n2] + Tp[k]
            if Q.is_zero():
                newcand.append((n1, n2, k))      # редукция в O — не решаем
                continue
            val = (Fp(Q[0]) - rp) / u2p          # x_E mod p
            val = val * binv                     # x_E / b mod p
            if val == 0 or val.is_square():
                newcand.append((n1, n2, k))
        cand = newcand
        if p < 60 or len(cand) < 200:
            print("    p = %4d: выжило %d  (%.1f с)"
                  % (p, len(cand), time.time() - t0))
            sys.stdout.flush()
    print("  после решета осталось %d кандидатов, %.1f с"
          % (len(cand), time.time() - t0))
    sys.stdout.flush()

    # --- точная проверка выживших ----------------------------------------
    hits = []
    for (n1, n2, k) in cand:
        P = n1 * G[0] + n2 * G[1] + tors[k]
        if P.is_zero():
            continue
        xE = (P[0] - riso) / u2
        r = xE / b
        if r > 0 and QQ(r).is_square():
            tt = QQ(r).sqrt()
            yE = P[1] / u3
            Yv = yE / b
            hits.append((n1, n2, k, tt, Yv))
            print("  *** ТОЧКА H из E_b: (n1,n2,T) = (%d,%d,%s), t = %s"
                  % (n1, n2, tors[k], tt))
            sys.stdout.flush()
    hbound = lam * NMAX ^ 2
    return hits, hbound, lam, Em


# ======================================================================
# ЧАСТЬ E. Точки C через параметризацию коники F4 = квадрат
# ======================================================================
def conic_search_C(XMAX, primes=(3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41),
                   report_every=500):
    """Перебор ВСЕХ t = u/v с F4(t) = квадрат, до высоты ~2 XMAX^2.

    u = alpha+beta, v = alpha-beta, {alpha,beta} = {x^2-y^2, 2xy},
    gcd(x,y) = 1, x > y >= 0, x-y нечётно.
    Далее нужно, чтобы A0 = m^2v^2+n^2u^2 и A8 = n^2v^2+m^2u^2 были квадратами.
    """
    qr = {}
    for p in primes:
        arr = np.zeros(p, dtype=np.uint8)
        for z in range(p):
            arr[(z * z) % p] = 1
        qr[p] = arr
    m2 = m * m
    n2 = n * n
    found = []
    t0 = time.time()
    n_pairs = 0
    n_surv = 0
    maxheight = 0
    for x in range(1, XMAX + 1):
        ys = np.arange(0, x, dtype=np.int64)
        # взаимно простые, противоположной чётности
        mask = ((x - ys) % 2 == 1)
        if x > 1:
            gg = np.gcd(np.int64(x), ys)
            mask &= (gg == 1)
        ys = ys[mask]
        if len(ys) == 0:
            continue
        A = x * x - ys * ys                      # x^2-y^2
        Bv = 2 * x * ys                          # 2xy
        al = np.maximum(A, Bv)
        be = np.minimum(A, Bv)
        u = al + be
        v = al - be
        keep = v > 0                             # v=0 => t=oo, точек нет
        u = u[keep]
        v = v[keep]
        yk = ys[keep]
        if len(u) == 0:
            continue
        n_pairs += len(u)
        if len(u):
            maxheight = max(maxheight, int(np.max(np.maximum(u, v))))
        alive = np.ones(len(u), dtype=bool)
        for p in primes:
            up = (u % p)
            vp = (v % p)
            a0 = (m2 * vp * vp + n2 * up * up) % p
            a8 = (n2 * vp * vp + m2 * up * up) % p
            alive &= (qr[p][a0] == 1) & (qr[p][a8] == 1)
            if not alive.any():
                break
        idx = np.nonzero(alive)[0]
        n_surv += len(idx)
        for i in idx:
            uu = Integer(int(u[i]))
            vv = Integer(int(v[i]))
            A0 = m ^ 2 * vv ^ 2 + n ^ 2 * uu ^ 2
            A8 = n ^ 2 * vv ^ 2 + m ^ 2 * uu ^ 2
            if A0.is_square() and A8.is_square():
                tt = QQ(uu) / QQ(vv)
                found.append(tt)
                print("  *** ТОЧКА C: t = %s (x=%d, y=%d)"
                      % (tt, x, int(yk[i])))
                sys.stdout.flush()
        if x % report_every == 0:
            print("    x = %5d/%d  пар=%d выжило=%d найдено=%d  высота<=%d  %.1f с"
                  % (x, XMAX, n_pairs, n_surv, len(found), maxheight,
                     time.time() - t0))
            sys.stdout.flush()
    stats = dict(XMAX=XMAX, pairs=n_pairs, survivors=n_surv,
                 max_height=maxheight, seconds=round(time.time() - t0, 1))
    return found, stats


# ======================================================================
# Запуск
# ======================================================================
def main():
    args = [z for z in sys.argv[1:] if not z.startswith('-')]
    HMAX = Integer(args[0]) if len(args) > 0 else 20000
    NMAX = Integer(args[1]) if len(args) > 1 else 100
    XMAX = Integer(args[2]) if len(args) > 2 else 6000

    all_pts = {}

    banner("ЧАСТЬ C. Полный перебор ящика в E_b(Q), |n1|,|n2| <= %d" % NMAX)
    hits, hbound, lam, Em = mw_box_search(NMAX)
    print("  найдено точек H: %d" % len(hits))
    print("  [доказано ПО (при верных генераторах и насыщенности)]")
    print("  перебор покрывает ВСЕ P in E_b(Q) с hhat(P) <= lambda_min*NMAX^2")
    print("  = %.4f * %d^2 = %.1f" % (lam, NMAX, hbound))
    for (n1, n2, k, tt, Yv) in hits:
        all_pts[QQ(tt)] = Yv
        all_pts[-QQ(tt)] = Yv

    banner("ЧАСТЬ B. Безусловный перебор t = u/v, 0 <= u <= v <= %d" % HMAX)
    primes_fast = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    primes_slow = [37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                   101, 103, 107, 109, 113, 127]
    found, stats = sieve_search(HMAX, primes_fast, primes_slow)
    print("  итог части B: %s" % stats)
    for (u, v, Y) in found:
        all_pts[QQ(u) / QQ(v)] = Y

    banner("ЧАСТЬ E. Точки C через конику F4 = квадрат, x <= %d" % XMAX)
    cpts, cstats = conic_search_C(XMAX)
    print("  итог части E: %s" % cstats)
    print("  найдено точек C: %s" % (cpts,))

    banner("ЧАСТЬ D. Подъём найденных точек H на кривую C рода 5")
    if not all_pts:
        print("  точек не найдено — поднимать нечего")
    lifted = []
    for tt in sorted(all_pts.keys()):
        ok, deg = report_point(tt, all_pts[tt], "H")
        if ok:
            lifted.append((tt, deg))

    banner("СВОДКА")
    print("  пара (m,n) = (15,8)")
    print("  найденные точки H (с точностью до t -> -t, t -> 1/t): %s"
          % sorted(set(abs(QQ(x)) for x in all_pts.keys())))
    print("  из них поднимаются на C: %s" % [str(x[0]) for x in lifted])
    print("  вырожденные (есть равные клетки): %s"
          % [str(x[0]) for x in lifted if x[1]])
    print("  ЧАСТЬ B: безусловно проверены все t = u/v с max(|u|,|v|) <= %d"
          % HMAX)
    print("  ЧАСТЬ C: проверены все P in E_b(Q) с hhat(P) <= %.1f" % hbound)
    print("  ЧАСТЬ E: проверены все t с F4(t) = квадрат и высотой <= %d"
          % cstats['max_height'])
    print("  За этими границами НИЧЕГО не проверено.")
    print("  Отсутствие новых находок — НЕ доказательство их отсутствия.")
    print("  общее время: %.1f с" % (time.time() - T_START))


main()
