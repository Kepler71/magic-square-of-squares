#!/usr/bin/env sage
# -*- coding: utf-8 -*-
r"""
search_H_19_5.sage      (Claude, 12.09.2026)
============================================

Расширенный поиск рациональных точек на кривой рода 2

    H : Y^2 = F0(t) * F4(t) * F8(t),   (m,n) = (19,5)

    s  = (m^2+n^2)/2 = 193
    F0 = m^2 + n^2 t^2 = 361 + 25 t^2
    F4 = s (1 + t^2)   = 193 (1 + t^2)
    F8 = n^2 + m^2 t^2 = 25 + 361 t^2

Запуск:   sage search_H_19_5.sage [HMAX] [NMAX]
          HMAX — граница высоты в части B (по умолчанию 2000)
          NMAX — граница кратности генератора в части C (по умолчанию 100000)


ЧАСТЬ B — безусловный перебор по высоте.
  t = u/v, gcd(u,v) = 1, 0 <= u <= v <= HMAX.  Симметрии
      t -> -t   (f чётна)                      сводят t к t >= 0,
      t -> 1/t  (f палиндромична, F0 <-> F8)   сводят t к 0 <= t <= 1,
  поэтому этот диапазон исчерпывает ВСЕ рациональные t с точностью до
  симметрии, и обе симметрии сохраняют также условие подъёма на C.
  Решето: N(u,v) = v^6 F0F4F8 однородна степени 6, значит условие
  "N — квадратичный вычет mod p" зависит только от (u:v) in P^1(F_p).
  Стадия 1: одна свёрнутая маска по модулю M = 11*17*37 = 6919
            (плотности 0.174, 0.336, 0.316 — три самых сильных простых).
  Стадия 2: ~30 дополнительных простых по выжившим индексам.
  Стадия 3: точный тест Integer.is_square.

ЧАСТЬ C — поиск по эллиптической кривой (условен: опирается на rank E_b = 1).
  Jac(H) ~ E_b^2; карта степени 2
      pi : H -> E_b,  (t,Y) |-> (X,V) = (b t^2, b Y),   b = s m^2 n^2.
  Обратно: P = (X,V) in E_b(Q) даёт точку H тогда и только тогда, когда
  X/b — квадрат рационального числа; тогда t = ±sqrt(X/b), Y = V/b.
  Значит H(Q) <-> { P in E_b(Q) : x(P)/b in (Q^*)^2 } — БИЕКЦИЯ с точностью
  до t -> -t.  E_b(Q) = Z*G + (Z/2)^2 (PARI ellrank = [1,1]; G насыщен).
  Так как кручение 2-е, x(nG+T) = x(-nG+T), достаточно n >= 1.
  Стадия 1: для набора простых p хорошей редукции строится орбита
            {x(nG+T) mod p} (период = порядок G в E(F_p)); n отбрасывается,
            если kronecker(x*b, p) = -1.  Маска по n тиражируется периодом.
  Стадия 2: точный пересчёт nG+T в Q для выживших n.
  Покрывает все P in E_b(Q) с hhat(P) <= hhat(G)*NMAX^2.

ЧАСТЬ D — подъём на кривую рода 5  C : u_i^2 = F_i(t).
  Нужно, чтобы F0(t), F4(t), F8(t) были квадратами КАЖДЫЙ ПО ОТДЕЛЬНОСТИ.

ОТСУТСТВИЕ НАХОДКИ НЕ ЕСТЬ ДОКАЗАТЕЛЬСТВО H(Q) = пусто.
"""

import time, sys
import numpy as np

# ----------------------------------------------------------------------
# 0. Параметры семейства
# ----------------------------------------------------------------------
m = Integer(19)
n = Integer(5)
s = (m ^ 2 + n ^ 2) // 2          # 193, целое именно для (19,5)
b = s * m ^ 2 * n ^ 2             # 1741825
aa = QQ(m ^ 2) / n ^ 2 + 1 + QQ(n ^ 2) / m ^ 2

R.<t> = QQ[]
F0 = m ^ 2 + n ^ 2 * t ^ 2
F4 = s * (1 + t ^ 2)
F8 = n ^ 2 + m ^ 2 * t ^ 2
f = F0 * F4 * F8

# целые коэффициенты E_b
A2 = Integer(b * aa)
A4 = Integer(b ^ 2 * aa)
A6 = Integer(b ^ 3)
E = EllipticCurve([0, A2, 0, A4, A6])

# генератор (см. ниже: получен PARI ellrank, насыщен Sage E.saturation)
GX = QQ(-28561963657) / 1369
GY = QQ(2089086742828800) / 50653


def banner():
    print("=" * 74)
    print("H : Y^2 = F0*F4*F8,   (m,n) = (%d,%d)" % (m, n))
    print("  s = %d,  b = s*m^2*n^2 = %d,  a = %s" % (s, b, aa))
    print("  f(t) =", f)
    print("  f палиндромичен:", f.list() == f.list()[::-1])
    print("  b квадрат? %s  =>  точек H над t=0 и t=oo нет"
          " (f(0) = старший коэф. = b = %d)" % (Integer(b).is_square(), b))
    print("=" * 74)


def identities():
    """[доказано] символьные тождества, на которых стоит часть C."""
    X = b * t ^ 2
    assert X ^ 3 + A2 * X ^ 2 + A4 * X + A6 == b ^ 2 * f
    assert X + b == m ^ 2 * n ^ 2 * F4
    assert X + s * m ^ 4 == s * m ^ 2 * F0
    assert X + s * n ^ 4 == s * n ^ 2 * F8
    assert E.two_torsion_rank() == 2
    assert f(-t) == f(t)
    assert (F0.subs(t=1 / t) * t ^ 2).numerator() == F8
    assert (F4.subs(t=1 / t) * t ^ 2).numerator() == F4
    print("[доказано] символьно проверено:")
    print("    X^3 + A2 X^2 + A4 X + A6 = b^2 f(t)   при X = b t^2")
    print("    X + b     = m^2 n^2 * F4 ;  X + s m^4 = s m^2 * F0 ;"
          "  X + s n^4 = s n^2 * F8")
    print("    f(-t) = f(t) ;  t^2 F0(1/t) = F8(t) ;  t^2 F4(1/t) = F4(t)")
    print("    E_b: y^2 = x^3 + %d x^2 + %d x + %d,  E_b[2](Q) = (Z/2)^2"
          % (A2, A4, A6))
    print()


# ----------------------------------------------------------------------
# 1. Точные проверки отдельных точек
# ----------------------------------------------------------------------
def Nint(u, v):
    """N(u,v) = v^6 * F0F4F8(u/v) — целое, однородное степени 6."""
    u = Integer(u); v = Integer(v)
    return ((m ^ 2 * v ^ 2 + n ^ 2 * u ^ 2) * s * (u ^ 2 + v ^ 2)
            * (n ^ 2 * v ^ 2 + m ^ 2 * u ^ 2))


def lift_to_C(tt):
    """Квадратичность F0,F4,F8 по отдельности в рациональной точке t."""
    tt = QQ(tt)
    vals = tuple(Fi(tt) for Fi in (F0, F4, F8))
    return tuple(x.is_square() for x in vals), vals


# ----------------------------------------------------------------------
# 2. ЧАСТЬ B: решето по высоте
# ----------------------------------------------------------------------
def qr_table(p):
    """T[vr][ur] = 1, если N(ur,vr) mod p — квадрат (0 считается квадратом)."""
    qr = np.zeros(p, dtype=np.uint8)
    for x in range(p):
        qr[(x * x) % p] = 1
    T = np.zeros((p, p), dtype=np.uint8)
    mm = int(m) ** 2 % p; nn = int(n) ** 2 % p; ss = int(s) % p
    for vr in range(p):
        for ur in range(p):
            val = (mm * vr * vr + nn * ur * ur) % p * ss % p
            val = val * ((ur * ur + vr * vr) % p) % p
            val = val * ((nn * vr * vr + mm * ur * ur) % p) % p
            T[vr][ur] = qr[val]
    return T


TILE_PRIMES = [int(11), int(17), int(37)]
FANCY_PRIMES = [int(q) for q in (79, 137, 113, 149, 89, 131, 167, 103, 97,
                109, 157, 127, 173, 107, 53, 3, 7, 19, 211, 233, 199, 41, 43,
                29, 59, 61, 67, 71, 73, 83, 101, 139, 151, 163, 179, 181,
                191, 197)]


def sieve_search(HMAX, report_every=None, logf=None):
    t0 = time.time()
    M = 1
    for p in TILE_PRIMES:
        M *= p
    tabs_t = {p: qr_table(p) for p in TILE_PRIMES}
    tabs_f = {p: qr_table(p) for p in FANCY_PRIMES}
    dens_t = {p: float(tabs_t[p].sum()) / p ** 2 for p in TILE_PRIMES}
    dens_f = {p: float(tabs_f[p].sum()) / p ** 2 for p in FANCY_PRIMES}
    # предтиражированные строки длины M
    tiled = {p: np.zeros((p, M), dtype=np.uint8) for p in TILE_PRIMES}
    for p in TILE_PRIMES:
        for r in range(p):
            tiled[p][r] = np.resize(tabs_t[p][r], int(M))

    keep1 = 1.0
    for p in TILE_PRIMES:
        keep1 *= dens_t[p]
    keep2 = keep1
    for p in FANCY_PRIMES:
        keep2 *= dens_f[p]
    print("  стадия 1: свёрнутая маска mod M = %d (простые %s), плотности %s"
          % (M, TILE_PRIMES, {p: round(dens_t[p], 4) for p in TILE_PRIMES}))
    print("  стадия 2: %d простых, %s"
          % (len(FANCY_PRIMES), {p: round(dens_f[p], 3) for p in FANCY_PRIMES}))
    print("  ожидаемая доля выживших: после ст.1 ~ %.3e, после ст.2 ~ %.3e"
          % (keep1, keep2))
    print("  построение таблиц: %.1f s" % (time.time() - t0))
    sys.stdout.flush()

    if report_every is None:
        report_every = max(1, HMAX // 20)

    found = []
    n_surv1 = 0
    n_surv2 = 0
    n_exact = 0
    t0 = time.time()
    p11, p17, p37 = TILE_PRIMES
    for v in range(int(1), int(HMAX) + int(1)):
        comb = tiled[p11][v % p11] & tiled[p17][v % p17]
        comb &= tiled[p37][v % p37]
        if not comb.any():
            continue                      # весь столбец v убит (напр. 11|v)
        L = int(v) + int(1)
        idx = np.nonzero(np.resize(comb, L) if L > int(M) else comb[:L])[0]
        n_surv1 += len(idx)
        for p in FANCY_PRIMES:
            if idx.size == 0:
                break
            idx = idx[np.nonzero(tabs_f[p][v % p][idx % p])[0]]
        n_surv2 += idx.size
        for uu in idx:
            uu = int(uu)
            if gcd(uu, v) != 1:
                continue
            n_exact += 1
            N = Nint(uu, v)
            if N.is_square():
                Y = QQ(N.sqrt()) / QQ(v) ** 3
                found.append((uu, v, Y))
                msg = "  *** ТОЧКА H: t = %d/%d, Y = %s" % (uu, v, Y)
                print(msg); sys.stdout.flush()
                if logf:
                    logf.write(msg + "\n"); logf.flush()
        if v % report_every == 0:
            print("    v = %7d / %d   ст1=%d ст2=%d точных=%d найдено=%d   %.1f s"
                  % (v, HMAX, n_surv1, n_surv2, n_exact, len(found),
                     time.time() - t0))
            sys.stdout.flush()
    stats = dict(HMAX=int(HMAX), surv1=n_surv1, surv2=n_surv2,
                 exact=n_exact, found=len(found),
                 seconds=round(time.time() - t0, 1))
    return found, stats


# ----------------------------------------------------------------------
# 3. ЧАСТЬ C: перебор по E_b(Q) с модулярным предфильтром
# ----------------------------------------------------------------------
def ec_add(P, Q, p, a2, a4):
    """Сложение на y^2 = x^3 + a2 x^2 + a4 x + a6 над F_p, P,Q — (x,y) или None."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2:
        if (y1 + y2) % p == 0:
            return None
        lam = (3 * x1 * x1 + 2 * a2 * x1 + a4) * pow(2 * y1, p - 2, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, p - 2, p) % p
    x3 = (lam * lam - a2 - x1 - x2) % p
    y3 = (lam * (x1 - x3) - y1) % p
    return (x3, y3)


def build_n_masks(p, Gp, Tp, a2, a4, bmod, qr):
    """Для простого p и КАЖДОГО класса кручения T отдельно: маска по
    n mod ord(G) — допустимо ли n (т.е. kronecker(x(nG+T)*b, p) != -1).
    Возвращает (список из len(Tp) масок, period) или None."""
    if Gp is None:
        return None
    orbit = []
    cur = Gp
    while cur is not None:
        orbit.append(cur)
        cur = ec_add(cur, Gp, p, a2, a4)
        if len(orbit) > 4 * p + 10:
            return None
    period = len(orbit) + 1               # ord(G); orbit[i] = (i+1)G
    masks = []
    for T in Tp:
        allowed = np.zeros(period, dtype=np.uint8)
        for i, P in enumerate(orbit):
            Q = ec_add(P, T, p, a2, a4)
            if Q is None:
                allowed[(i + 1) % period] = 1
                continue
            val = Q[0] * bmod % p
            if val == 0 or qr[val]:
                allowed[(i + 1) % period] = 1
        # n ≡ 0 (mod period): точка nG+T = T (или O)
        if T is None:
            allowed[0] = 1
        else:
            val = T[0] * bmod % p
            allowed[0] = 1 if (val == 0 or qr[val]) else 0
        masks.append(allowed)
    return masks, period


def elliptic_search(NMAX, n_primes=30, pstart=2000, NDIRECT=200, logf=None):
    NMAX = int(NMAX)
    NDIRECT = int(min(int(NDIRECT), int(NMAX)))
    G = E(GX, GY)
    tors = E.torsion_points()
    print("  E_b =", E)
    print("  генератор G = %s" % (G,))
    print("  hhat(G) = %.6f" % float(G.height()))
    print("  кручение E_b(Q)_tors = %s" % (tors,))
    print("  x(nG+T) = x(-nG+T) при 2-кручении T  =>  достаточно n >= 0")
    sys.stdout.flush()

    a2 = int(A2); a4 = int(A4); a6 = int(A6)
    disc = Integer(E.discriminant())
    # выбираем простые хорошей редукции
    ps = []
    p = next_prime(pstart)
    while len(ps) < n_primes:
        if disc % p != 0 and Integer(b) % p != 0:
            ps.append(int(p))
        p = next_prime(p)

    t0 = time.time()
    Tlist = [None] + [T for T in tors if not T.is_zero()]
    per_prime = []
    for p in ps:
        qr = np.zeros(p, dtype=np.uint8)
        for x in range(p):
            qr[(x * x) % p] = 1
        Gp = (int(GX.numerator() * pow(int(GX.denominator()), p - 2, p) % p),
              int(GY.numerator() * pow(int(GY.denominator()), p - 2, p) % p))
        x1, y1 = Gp                       # контроль: точка лежит на кривой
        assert (y1 * y1 - (x1 ** 3 + a2 * x1 * x1 + a4 * x1 + a6)) % p == 0
        Tp = [None] + [(int(T[0]) % p, int(T[1]) % p)
                       for T in tors if not T.is_zero()]
        res = build_n_masks(p, Gp, Tp, a2, a4, int(b) % p, qr)
        if res is None:
            continue
        masks, period = res
        per_prime.append((p, period, masks))
    dens = [1.0] * len(Tlist)
    for (p, period, masks) in per_prime:
        for j in range(len(Tlist)):
            dens[j] *= float(masks[j].sum()) / period
    print("  модулярный предфильтр: %d простых, периоды %s..."
          % (len(per_prime), [z[1] for z in per_prime[:6]]))
    print("  ожидаемая доля выживших n по 4 классам кручения: %s"
          % ["%.2e" % d for d in dens])
    print("  подготовка предфильтра: %.1f s" % (time.time() - t0))
    sys.stdout.flush()

    t0 = time.time()
    NP1 = int(NMAX) + int(1)
    cands = []
    for j, T in enumerate(Tlist):
        alive = np.ones(NP1, dtype=np.uint8)
        for (p, period, masks) in per_prime:
            alive &= np.resize(masks[j], NP1)
        cj = np.nonzero(alive)[0]
        cands.append(cj)
    tot = sum(len(c) for c in cands)
    print("  выживших пар (n,T): %d из %d  (%.1f s)"
          % (tot, 4 * NP1, time.time() - t0))
    sys.stdout.flush()

    # прямая точная проверка малых n — контроль, что предфильтр ничего
    # не прячет (считаем nG в Q без всякого решета)
    hits = []
    t0 = time.time()
    P = E(0)
    for k in range(int(1), int(NDIRECT) + int(1)):
        P = P + G
        for T in Tlist:
            Q = P if T is None else P + T
            if Q.is_zero():
                continue
            r = Q[0] / b
            if r > 0 and r.is_square():
                tt = QQ(r).sqrt()
                hits.append((k, T, tt, Q[1] / b))
                msg = "  *** ТОЧКА E->H (прямая): n=%d, T=%s, t=%s" % (k, T, tt)
                print(msg); sys.stdout.flush()
                if logf:
                    logf.write(msg + "\n"); logf.flush()
    print("  прямая точная проверка n <= %d (без решета): %d точек, %.1f s"
          % (NDIRECT, len(hits), time.time() - t0))
    print("    (для контроля: x(%dG) имеет %d цифр в числителе)"
          % (NDIRECT, len(str((NDIRECT * G)[0].numerator()))))
    sys.stdout.flush()

    nexact = 0
    for j, T in enumerate(Tlist):
        for nn in cands[j]:
            nn = int(nn)
            if nn == 0 and T is None:
                continue
            if nn <= int(NDIRECT):
                continue          # уже проверено точно выше
            nexact += 1
            P = nn * G if T is None else nn * G + T
            if P.is_zero():
                continue
            r = P[0] / b
            if r > 0 and r.is_square():
                tt = QQ(r).sqrt()
                hits.append((nn, T, tt, P[1] / b))
                msg = "  *** ТОЧКА E->H: n=%d, T=%s, t=%s" % (nn, T, tt)
                print(msg); sys.stdout.flush()
                if logf:
                    logf.write(msg + "\n"); logf.flush()
    print("  точных проверок в Q: %d" % nexact)
    hb = float(G.height()) * NMAX ** 2
    return hits, hb, tot


# ----------------------------------------------------------------------
# 4. main
# ----------------------------------------------------------------------
def main():
    args = [x for x in sys.argv[1:] if not x.startswith('-')]
    HMAX = Integer(args[0]) if len(args) > 0 else Integer(2000)
    NMAX = Integer(args[1]) if len(args) > 1 else Integer(100000)

    logf = open("/home/kep/magicKube/bridge/search_H_19_5_hits.log", "a")
    logf.write("\n=== run HMAX=%s NMAX=%s %s ===\n"
               % (HMAX, NMAX, time.strftime("%Y-%m-%d %H:%M:%S")))

    banner()
    identities()

    print("### ЧАСТЬ C: перебор E_b(Q) = Z*G + (Z/2)^2, n <= %d ###" % NMAX)
    hits, hbound, ncand = elliptic_search(NMAX, logf=logf)
    print("  ИТОГ C: найдено %d точек H." % len(hits))
    print("  Покрыты ВСЕ P in E_b(Q) с hhat(P) <= %.4g" % hbound)
    print("  (условно: rank E_b = 1 по PARI ellrank, G насыщен по"
          " E.saturation => E_b(Q) = Z*G + (Z/2)^2)")
    print()

    print("### ЧАСТЬ B: решето по t = u/v, 0 <= u <= v <= %d ###" % HMAX)
    found, stats = sieve_search(HMAX, logf=logf)
    print("  ИТОГ B:", stats)
    print()

    print("### ЧАСТЬ D: подъём на C рода 5 (F0,F4,F8 — квадраты отдельно) ###")
    pts = [(QQ(u) / QQ(v), Y) for (u, v, Y) in found]
    pts += [(tt, Yv) for (_, _, tt, Yv) in hits]
    if not pts:
        print("  рациональных точек H не найдено — поднимать нечего")
    for (tt, Yv) in pts:
        sq, vals = lift_to_C(tt)
        print("  t = %s: F0 кв=%s F4 кв=%s F8 кв=%s -> точка C: %s"
              % (tt, sq[0], sq[1], sq[2], all(sq)))
        logf.write("lift t=%s -> %s %s\n" % (tt, sq, all(sq)))
    print()
    print("СТАТУС: отсутствие находок — НЕ доказательство H(Q) = пусто.")
    logf.close()


main()
