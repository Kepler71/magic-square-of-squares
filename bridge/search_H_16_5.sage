#!/usr/bin/env sage
# -*- coding: utf-8 -*-
r"""
search_H_16_5.sage
==================

Расширенный поиск рациональных точек на кривой рода 2

    H : Y^2 = F0(t) * F4(t) * F8(t),   (m,n) = (16,5)

    s  = (m^2+n^2)/2 = 281/2                 (m,n разной чётности => s НЕ целое)
    F0 = m^2 + n^2 t^2 = 256 + 25 t^2
    F4 = s (1 + t^2)   = (281/2)(1 + t^2)
    F8 = n^2 + m^2 t^2 = 25 + 256 t^2

Целая модель:  W = 2Y,  W^2 = f(t) := 4 F0 F4 F8 = 562(256+25t^2)(1+t^2)(25+256t^2).

Четыре независимые части.

ЧАСТЬ B (безусловный перебор по высоте).  t = u/v, gcd(u,v)=1, 0 <= u <= v <= HMAX.
  Симметрии t -> -t и t -> 1/t сводят все рациональные t к этому диапазону:
  все F_i чётны по t (значит t -> -t ничего не меняет), а t -> 1/t даёт
  F0(1/t) = F8(t)/t^2, F8(1/t) = F0(t)/t^2, F4(1/t) = F4(t)/t^2, то есть
  и произведение, и каждый сомножитель по отдельности сохраняют квадратность.
  Решето: для каждого малого простого p условие "N(u,v) — квадратичный вычет
  или 0 mod p" зависит только от (u mod p, v mod p); N однородна степени 6.
  Кандидаты, прошедшие решето, проверяются точно (is_square целого числа).

ЧАСТЬ C (решето Морделла-Вейля на E_b).  Jac(H) ~ E_b^2, карта степени 2
      pi : H -> E,  (t,W) |-> (X,V) = (B t^2, B W),   B = 4 s m^2 n^2 = 3596800.
  E : V^2 = X^3 + B a X^2 + B^2 a X + B^3 = (X+4b)(X+4s m^4)(X+4s n^4),
      a = m^2/n^2 + 1 + n^2/m^2.
  Так как точек H над t=0 и t=oo нет (класс квадрата b=562 нетривиален),
      H(Q)/{t -> -t}  <->  { P in E(Q), P != O : x(P)/B in (Q^*)^2 }.
  E(Q) = Z*P1 + Z*P2 + (Z/2)^2 (PARI ellrank; насыщенность проверена Sage).
  Для КАЖДОЙ точки P = n1 P1 + n2 P2 + T условие x(P)/B in (Q^*)^2 влечёт
  локальные условия: x(P)/B — квадрат в R и в каждом Q_p.  Для простого p
  хорошей редукции, не делящего B, это условие на (n1 mod ord(P1 mod p),
  n2 mod ord(P2 mod p), T) — обычное решето Морделла-Вейля.  Прогоняем весь
  прямоугольник |n1|,|n2| <= NMAX и проверяем выживших точно.
  Это ПОЛНАЯ проверка всех точек H, чей образ лежит в этом прямоугольнике,
  то есть всех точек с ограниченной канонической высотой.

ЧАСТЬ E (полное перечисление C-кандидатов по коническому сечению).
  Для точки на кривой рода 5  C : u_i^2 = F_i(t)  нужно, чтобы КАЖДОЕ из
  F0, F4, F8 было квадратом.  Условие "F4(u/v) — квадрат" при gcd(u,v)=1
  равносильно  562(u^2+v^2) = квадрат, то есть
        u^2 + v^2 = 562 w^2,   w in Z.
  Это коника с рациональной точкой (21,11,1) (21^2+11^2=562).  Перечисляем
  ВСЕ её целые точки с w <= WMAX через гауссовы целые (полное разложение
  562 w^2 в сумму двух взаимно простых квадратов).  Так как
        u^2+v^2 = 562 w^2  =>  w = sqrt((u^2+v^2)/562) <= max(|u|,|v|)*sqrt(2/562),
  перебор w <= WMAX ПОЛНОСТЬЮ покрывает все t = u/v наивной высоты
        max(|u|,|v|) <= WMAX/sqrt(2/562) = 16.76... * WMAX.
  Для каждой такой t остаётся проверить, квадраты ли 25u^2+256v^2 (=v^2 F0)
  и 256u^2+25v^2 (=v^2 F8).  Заодно получаем все точки H с квадратным F4.

ЧАСТЬ D. Для каждой найденной точки H проверяется подъём на C.

Отсутствие находки НЕ является доказательством H(Q) = пусто и НЕ является
доказательством C(Q) = пусто.

Запуск:
    sage search_H_16_5.sage [HMAX] [NMAX] [WMAX]
по умолчанию HMAX=2000, NMAX=60, WMAX=60000.
"""

import time, sys
from math import isqrt
import numpy as np

# ----------------------------------------------------------------------
# 0. Параметры семейства и целая модель
# ----------------------------------------------------------------------
m = 16
n = 5
s = QQ(m ^ 2 + n ^ 2) / 2          # 281/2
b = s * m ^ 2 * n ^ 2              # 899200
B = ZZ(4 * b)                      # 3596800
aa = QQ(m ^ 2) / n ^ 2 + 1 + QQ(n ^ 2) / m ^ 2      # 72561/6400

R.<t> = QQ[]
F0 = m ^ 2 + n ^ 2 * t ^ 2
F4 = s * (1 + t ^ 2)
F8 = n ^ 2 + m ^ 2 * t ^ 2
Pprod = F0 * F4 * F8
Rz.<T> = ZZ[]
f = Rz(4 * Pprod)                  # W^2 = f(t),  W = 2Y

print("=" * 74)
print("H : Y^2 = F0*F4*F8,  (m,n) = (%d,%d)" % (m, n))
print("  s = %s,  b = s*m^2*n^2 = %s,  B = 4b = %s" % (s, b, B))
print("  a = %s,  класс квадрата b = %s" % (aa, QQ(b).squarefree_part()))
print("  F0 = %s;  F4 = %s;  F8 = %s" % (F0, F4, F8))
print("  целая модель: W = 2Y,  W^2 = f(t) =", f)
print("  f =", factor(f))
print("  f палиндромичен:", all(f[i] == f[6 - i] for i in range(7)))
print("=" * 74)

# ----------------------------------------------------------------------
# 1. Символьная проверка модели и карты H -> E
# ----------------------------------------------------------------------
E = EllipticCurve([0, B * aa, 0, B ^ 2 * aa, B ^ 3])
X = B * t ^ 2
assert X ^ 3 + B * aa * X ^ 2 + B ^ 2 * aa * X + B ^ 3 == B ^ 2 * f(t), \
    "карта (t,W) -> (B t^2, B W) неверна"
assert X + 4 * b == 4 * m ^ 2 * n ^ 2 * F4
assert X + 4 * s * m ^ 4 == 4 * s * m ^ 2 * F0
assert X + 4 * s * n ^ 4 == 4 * s * n ^ 2 * F8
assert E.a_invariants()[1] in ZZ and E.a_invariants()[3] in ZZ and E.a_invariants()[4] in ZZ
assert E.two_torsion_rank() == 2
print("[доказано] тождества проверены символьно:")
print("    (t,W) -> (X,V) = (B t^2, B W)  переводит W^2=f(t) в уравнение E")
print("    X + 4b      = 4 m^2 n^2 * F4")
print("    X + 4s m^4  = 4 s m^2  * F0")
print("    X + 4s n^4  = 4 s n^2  * F8")
print("    E =", E)
print("    E = (X+%d)(X+%d)(X+%d)" % (4 * b, 4 * s * m ^ 4, 4 * s * n ^ 4))
print("[доказано] точек H над t=0 и t=oo нет: f(0) = %s, старший коэф. = %s,"
      % (f(0), f.leading_coefficient()))
print("           оба имеют класс квадрата %s != 1."
      % QQ(f(0)).squarefree_part())
print()


# ----------------------------------------------------------------------
# 2. Вспомогательные функции
# ----------------------------------------------------------------------
def Nint(u, v):
    """N(u,v) = v^6 * f(u/v) = 562 (m^2v^2+n^2u^2)(u^2+v^2)(n^2v^2+m^2u^2).
    Точка H существует <=> N(u,v) — точный квадрат; тогда Y = sqrt(N)/(2 v^3)."""
    u = int(u); v = int(v)
    return (2 * (m ** 2 + n ** 2)
            * (m ** 2 * v * v + n ** 2 * u * u)
            * (u * u + v * v)
            * (n ** 2 * v * v + m ** 2 * u * u))


def is_sq_int(x):
    if x < 0:
        return False
    r = isqrt(x)
    return r * r == x


def lift_data(u, v):
    """v^2*F0, v^2*F4 (умноженное на 4 -> целое), v^2*F8 для t=u/v."""
    u = int(u); v = int(v)
    A0 = m ** 2 * v * v + n ** 2 * u * u          # v^2 F0
    A4 = 2 * (m ** 2 + n ** 2) * (u * u + v * v)  # 4 v^2 F4
    A8 = n ** 2 * v * v + m ** 2 * u * u          # v^2 F8
    return A0, A4, A8


def check_lift_to_C(u, v):
    """Квадратны ли F0, F4, F8 ПО ОТДЕЛЬНОСТИ при t = u/v.
    F0 квадрат <=> v^2F0 квадрат; F4 квадрат <=> 4v^2F4 квадрат."""
    A0, A4, A8 = lift_data(u, v)
    return (is_sq_int(A0), is_sq_int(A4), is_sq_int(A8)), (A0, A4, A8)


def check_lift_to_C_rational(tt):
    tt = QQ(tt)
    vals = tuple(Fi(tt) for Fi in (F0, F4, F8))
    return tuple(x.is_square() for x in vals), vals


# ----------------------------------------------------------------------
# 3. ЧАСТЬ B: решето + точный перебор t = u/v
# ----------------------------------------------------------------------
def build_tables(primes_list):
    """tab[p][vr][ur] = 1, если N(ur,vr) mod p — квадрат (0 тоже допустим)."""
    tabs = {}
    dens = {}
    for p in primes_list:
        qr = np.zeros(p, dtype=np.uint8)
        for x in range(p):
            qr[(x * x) % p] = 1
        ur = np.arange(p, dtype=np.int64)
        U2 = (ur * ur) % p
        Tb = np.zeros((p, p), dtype=np.uint8)
        c = (2 * (m ** 2 + n ** 2)) % p
        for vr in range(p):
            v2 = (vr * vr) % p
            val = (c * ((m ** 2 * v2 + n ** 2 * U2) % p)) % p
            val = (val * ((U2 + v2) % p)) % p
            val = (val * ((n ** 2 * v2 + m ** 2 * U2) % p)) % p
            Tb[vr] = qr[val]
        tabs[p] = Tb
        dens[p] = float(Tb.sum()) / (p * p)
    return tabs, dens


def sieve_search(HMAX, primes_fast, primes_slow, report_every=2000):
    tabs_f, dens_f = build_tables(primes_fast)
    tabs_s, dens_s = build_tables(primes_slow)
    exp_keep = 1.0
    for p in primes_fast:
        exp_keep *= dens_f[p]
    for p in primes_slow:
        exp_keep *= dens_s[p]
    print("  быстрые простые: %s" % (primes_fast,))
    print("  плотности: %s" % ({p: round(dens_f[p], 4) for p in primes_fast},))
    print("  медленные простые: %s" % (primes_slow,))
    print("  ожидаемая доля выживших после всего решета ~ %.3e" % exp_keep)
    sys.stdout.flush()

    found = []
    n_pairs = 0
    n_surv1 = 0
    n_surv2 = 0
    n_exact = 0
    t0 = time.time()

    for v in range(1, HMAX + 1):
        L = v + 1
        alive = np.ones(L, dtype=np.uint8)
        for p in primes_fast:
            row = tabs_f[p][v % p]
            tiled = np.tile(row, L // p + 1)[:L]
            np.logical_and(alive, tiled, out=alive)
        idx = np.nonzero(alive)[0]
        n_surv1 += len(idx)
        if len(idx):
            keep = np.ones(len(idx), dtype=np.uint8)
            for p in primes_slow:
                row = tabs_s[p][v % p]
                np.logical_and(keep, row[idx % p], out=keep)
            idx = idx[np.nonzero(keep)[0]]
        n_surv2 += len(idx)
        for u in idx:
            u = int(u)
            if gcd(u, v) != 1:
                continue
            n_exact += 1
            N = Nint(u, v)
            if is_sq_int(N):
                W = isqrt(N)
                Y = QQ(W) / (2 * QQ(v) ** 3)
                found.append((u, v, Y))
                print("  *** ТОЧКА H: t = %d/%d, Y = %s" % (u, v, Y))
                sys.stdout.flush()
        n_pairs += (2 if v == 1 else euler_phi(v))
        if v % report_every == 0:
            print("    v = %6d/%d  выжило1=%d выжило2=%d точных=%d найдено=%d  %.1f s"
                  % (v, HMAX, n_surv1, n_surv2, n_exact, len(found), time.time() - t0))
            sys.stdout.flush()

    stats = dict(HMAX=HMAX, pairs=n_pairs, surv1=n_surv1, surv2=n_surv2,
                 exact=n_exact, found=len(found), seconds=round(time.time() - t0, 1))
    return found, stats


# ----------------------------------------------------------------------
# 4. ЧАСТЬ C: решето Морделла-Вейля на E(Q)
# ----------------------------------------------------------------------
def mw_sieve(NMAX, P1, P2, maxprimes=40, ordlimit=12000, pmax=3000):
    """Перебор P = n1 P1 + n2 P2 + T, |n1|,|n2| <= NMAX, с локальным решетом.
    Возвращает список точно проверенных точек H."""
    tors = list(E.torsion_points())
    print("  E(Q) = Z*P1 + Z*P2 + (Z/2)^2")
    print("  P1 =", P1, " hhat =", float(P1.height()))
    print("  P2 =", P2, " hhat =", float(P2.height()))
    print("  кручение:", tors)
    print("  прямоугольник |n1|,|n2| <= %d, всего %d точек (с кручением)"
          % (NMAX, (2 * NMAX + 1) ** 2 * 4))
    sys.stdout.flush()

    L = 2 * NMAX + 1
    ns = np.arange(-NMAX, NMAX + 1, dtype=np.int64)
    # маска[ti][n1][n2]
    mask = [np.ones((L, L), dtype=bool) for _ in tors]

    # --- условие в бесконечном месте: x(P) должен быть > 0 ---
    # E(R) имеет две компоненты (disc>0). Компонентный характер eps: E(R)->Z/2
    # — гомоморфизм. Точки овала имеют x < 0, поэтому eps(P)=0 обязательно.
    e2 = sorted([QQ(r) for r, _ in E.division_polynomial(2).roots()])
    x_id = e2[-1]      # наибольший корень = начало компоненты единицы
    def eps(Pt):
        if Pt.is_zero():
            return 0
        return 0 if Pt[0] >= x_id else 1
    e1, e2p = eps(P1), eps(P2)
    print("  корни 2-кручения: %s; компонента единицы: x >= %s" % (e2, x_id))
    print("  eps(P1)=%d, eps(P2)=%d, eps(T)=%s"
          % (e1, e2p, [eps(Tt) for Tt in tors]))
    n1g, n2g = np.meshgrid(ns, ns, indexing='ij')
    for ti, Tt in enumerate(tors):
        par = (n1g * e1 + n2g * e2p + eps(Tt)) % 2
        mask[ti] &= (par == 0)
    alive = sum(int(mk.sum()) for mk in mask)
    print("  после вещественного места выжило %d из %d"
          % (alive, L * L * len(tors)))
    sys.stdout.flush()

    # корни 2-кручения (целые): x(P+T) = r + (r-r')(r-r'')/(x(P)-r)
    roots2 = [ZZ(4 * b), ZZ(4 * s * m ^ 4), ZZ(4 * s * n ^ 4)]
    roots2 = [-r for r in roots2]
    tors_idx = {}
    for ti, Tt in enumerate(tors):
        tors_idx[ti] = None if Tt.is_zero() else ZZ(Tt[0])

    badp = set(ZZ(p) for p, _ in factor(E.discriminant())) | {ZZ(2)}
    used = []
    t0 = time.time()
    NONE = -1                     # маркер точки O
    for p in primes(5, pmax):
        if len(used) >= maxprimes:
            break
        if p in badp or B % p == 0:
            continue
        Ep = E.change_ring(GF(p))
        try:
            r1 = Ep(P1); r2 = Ep(P2)
        except Exception:
            continue              # точка редуцируется плохо — пропускаем простое
        o1 = int(r1.order()); o2 = int(r2.order())
        if o1 * o2 > ordlimit:
            continue
        Bp = int(GF(p)(B))
        # сетка x-координат для P = i*r1 + j*r2
        A2list = []
        cur = Ep(0)
        for _ in range(o2):
            A2list.append(cur); cur = cur + r2
        xs = np.empty((o1, o2), dtype=np.int64)
        base = Ep(0)
        for i in range(o1):
            for j in range(o2):
                Q = base + A2list[j]
                xs[i, j] = NONE if Q.is_zero() else int(Q[0])
            base = base + r1
        # таблица квадратичных вычетов mod p
        qr = np.zeros(p, dtype=np.int8)
        for x in range(1, p):
            qr[(x * x) % p] = 1
        Bpinv = pow(Bp, p - 2, p)
        for ti in range(len(tors)):
            r = tors_idx[ti]
            if r is None:
                xt = xs
            else:
                rp = int(GF(p)(r))
                others = [int(GF(p)(rr)) for rr in roots2 if rr != r]
                c = (rp - others[0]) * (rp - others[1]) % p
                xt = np.empty_like(xs)
                flat = xs.ravel(); out = xt.ravel()
                for k in range(flat.size):
                    xv = int(flat[k])
                    if xv == NONE:
                        out[k] = rp                     # O + T = T
                    elif xv == rp:
                        out[k] = NONE                   # T + T = O
                    else:
                        out[k] = (rp + c * pow((xv - rp) % p, p - 2, p)) % p
            tab = np.zeros(xs.shape, dtype=bool)
            fl = xt.ravel(); tb = tab.ravel()
            for k in range(fl.size):
                xv = int(fl[k])
                if xv == NONE or xv == 0:
                    tb[k] = True          # редукция в O или v_p(x)>0: не ограничиваем
                else:
                    tb[k] = bool(qr[(xv * Bpinv) % p])
            i1 = (ns % o1).astype(np.int64)
            i2 = (ns % o2).astype(np.int64)
            mask[ti] &= tab[np.ix_(i1, i2)]
        used.append((p, o1, o2))
        alive = sum(int(mk.sum()) for mk in mask)
        print("    p=%4d (o1=%3d,o2=%3d): выжило %d   %.1f s"
              % (p, o1, o2, alive, time.time() - t0))
        sys.stdout.flush()
        if alive == 0:
            break

    survivors = []
    for ti in range(len(tors)):
        ii, jj = np.nonzero(mask[ti])
        for k in range(len(ii)):
            survivors.append((int(ns[ii[k]]), int(ns[jj[k]]), ti))
    print("  использовано простых: %d; выживших пар (n1,n2,T): %d"
          % (len(used), len(survivors)))
    sys.stdout.flush()

    hits = []
    if survivors:
        print("  точная проверка выживших...")
        for (k1, k2, ti) in survivors:
            Pt = k1 * P1 + k2 * P2 + tors[ti]
            if Pt.is_zero():
                continue
            r = QQ(Pt[0]) / QQ(B)
            if r > 0 and r.is_square():
                tt = r.sqrt()
                Yv = QQ(Pt[1]) / (2 * QQ(B))
                hits.append((k1, k2, ti, tt, Yv))
                print("  *** ТОЧКА H из E: n1=%d n2=%d T=%s, t=%s"
                      % (k1, k2, tors[ti], tt))
    # достигнутая граница канонической высоты
    h11 = float(P1.height()); h22 = float(P2.height())
    h12 = float(((P1 + P2).height() - h11 - h22) / 2)
    lam_min = ((h11 + h22) - sqrt((h11 - h22) ** 2 + 4 * h12 ** 2)) / 2
    hbound = float(lam_min) * (NMAX + 1) ** 2
    return hits, survivors, hbound, used


# ----------------------------------------------------------------------
# 5. ЧАСТЬ E: полное перечисление точек коники u^2+v^2 = 562 w^2
# ----------------------------------------------------------------------
def gauss_reps(M, plist):
    """Все примитивные представления M = u^2+v^2 (u,v>=0, gcd(u,v)=1),
    M = 2 * prod p_i^{e_i}, все p_i = 1 mod 4, 4 не делит M.
    plist = [(p, e)] для нечётных p."""
    zs = [(1, 1)]                      # (1+i)
    for (p, e) in plist:
        # гауссово простое над p
        aa_, bb_ = two_squares(p)
        pi = (int(aa_), int(bb_))
        # pi^e и conj(pi)^e
        def powc(z, k):
            r = (1, 0)
            for _ in range(k):
                r = (r[0] * z[0] - r[1] * z[1], r[0] * z[1] + r[1] * z[0])
            return r
        z1 = powc(pi, e)
        z2 = powc((pi[0], -pi[1]), e)
        new = []
        for z in zs:
            new.append((z[0] * z1[0] - z[1] * z1[1], z[0] * z1[1] + z[1] * z1[0]))
            new.append((z[0] * z2[0] - z[1] * z2[1], z[0] * z2[1] + z[1] * z2[0]))
        zs = new
    out = set()
    for (x, y) in zs:
        x = abs(x); y = abs(y)
        if x > y:
            x, y = y, x
        out.add((x, y))
    return out


def conic_search(WMAX, report_every=200000):
    """Полное перечисление t = u/v с квадратным F4 и w <= WMAX."""
    t0 = time.time()
    # 1) решето: w нечётно, все простые делители = 1 mod 4
    ok = np.ones(WMAX + 1, dtype=bool)
    ok[0] = False
    if WMAX >= 2:
        ok[2::2] = False
    spf = np.zeros(WMAX + 1, dtype=np.int32)
    sieve_lim = isqrt(WMAX) + 1
    for p in range(3, WMAX + 1, 2):
        if spf[p] == 0:
            spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
            if p % 4 == 3:
                ok[p::p] = False
    adm = np.nonzero(ok)[0]
    print("  допустимых w (нечётные, все простые делители = 1 mod 4) до %d: %d"
          % (WMAX, len(adm)))
    print("  (решето построено за %.1f s)" % (time.time() - t0))
    sys.stdout.flush()

    Hpts = []
    Cpts = []
    near = []
    cnt_points = 0
    for ii, w in enumerate(adm):
        w = int(w)
        # разложение w
        pl = {}
        x = w
        while x > 1:
            p = int(spf[x])
            e = 0
            while x % p == 0:
                x //= p; e += 1
            pl[p] = e
        # M = 2 * 281 * w^2
        fac = {281: 1}
        for p, e in pl.items():
            fac[p] = fac.get(p, 0) + 2 * e
        M = 2 * 281 * w * w
        for (u, v) in gauss_reps(M, sorted(fac.items())):
            if v == 0:
                continue
            if u * u + v * v != M:
                raise RuntimeError("плохое представление")
            if gcd(u, v) != 1:
                continue
            cnt_points += 1
            # t = u/v и t = v/u — обе получаются (симметрия t->1/t);
            # условия ниже симметричны относительно перестановки u<->v
            Q0 = 25 * u * u + 256 * v * v
            Q8 = 256 * u * u + 25 * v * v
            s0 = is_sq_int(Q0); s8 = is_sq_int(Q8)
            if s0 and s8:
                Cpts.append((u, v, w))
                print("  *** ТОЧКА C: t = %d/%d (w=%d)" % (u, v, w))
                sys.stdout.flush()
            elif s0 or s8:
                near.append((u, v, w, s0, s8))
                print("  [почти] t = %d/%d: F0 кв=%s, F8 кв=%s (F4 кв всегда)"
                      % (u, v, s0, s8))
                sys.stdout.flush()
            if is_sq_int(Q0 * Q8):
                Hpts.append((u, v, w))
                print("  *** ТОЧКА H (с квадратным F4): t = %d/%d" % (u, v))
                sys.stdout.flush()
        if report_every and (ii + 1) % report_every == 0:
            print("    обработано %d/%d допустимых w, точек коники %d, %.1f s"
                  % (ii + 1, len(adm), cnt_points, time.time() - t0))
            sys.stdout.flush()
    print("  всего примитивных точек коники: %d,  время %.1f s"
          % (cnt_points, time.time() - t0))
    return Hpts, Cpts, near, cnt_points


# ----------------------------------------------------------------------
# 6. Запуск
# ----------------------------------------------------------------------
def main():
    args = [x for x in sys.argv[1:] if not x.startswith('-')]
    HMAX = Integer(args[0]) if len(args) > 0 else 2000
    NMAX = Integer(args[1]) if len(args) > 1 else 60
    WMAX = Integer(args[2]) if len(args) > 2 else 60000

    allpts = []

    print("### ЧАСТЬ C: решето Морделла-Вейля на E(Q), |n1|,|n2| <= %d ###" % NMAX)
    P1 = E(QQ(-3853673901961) / 323761, QQ(9009520099488504981) / 184220009)
    P2 = E(QQ(-2933770124800) / 776161, QQ(3114315170322048000) / 683797841)
    hits, surv, hbound, used = mw_sieve(NMAX, P1, P2)
    print("  итог части C: найдено %d точек H." % len(hits))
    print("  ПОКРЫТО: все P in E(Q) с |n1|,|n2| <= %d, то есть все точки H," % NMAX)
    print("  чей образ имеет каноническую высоту <= %.4g (нижняя оценка)." % hbound)
    for (k1, k2, ti, tt, Yv) in hits:
        allpts.append((tt, Yv, "E(Q): %d*P1+%d*P2+T%d" % (k1, k2, ti)))
    print()

    print("### ЧАСТЬ B: решето по t = u/v, 0 <= u <= v <= %d ###" % HMAX)
    primes_fast = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    primes_slow = [37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                   101, 103, 107, 109, 113, 127]
    found, stats = sieve_search(HMAX, primes_fast, primes_slow)
    print("  итог части B:", stats)
    for (u, v, Y) in found:
        allpts.append((QQ(u) / QQ(v), Y, "перебор t=%d/%d" % (u, v)))
    print()

    print("### ЧАСТЬ E: коника u^2+v^2 = 562 w^2, w <= %d ###" % WMAX)
    print("  покрывает ВСЕ t с квадратным F4 и наивной высотой <= %d"
          % int(WMAX / sqrt(2.0 / 562)))
    Hp, Cp, near, ncon = conic_search(WMAX)
    print("  итог части E: точек коники %d; точек H %d; точек C %d; почти-C %d"
          % (ncon, len(Hp), len(Cp), len(near)))
    for (u, v, w) in Hp:
        allpts.append((QQ(u) / QQ(v), None, "коника w=%d" % w))
    print()

    print("### ЧАСТЬ D: подъём найденных точек H на C ###")
    if not allpts:
        print("  точек H не найдено — поднимать нечего")
    seen = set()
    for (tt, Yv, src) in allpts:
        if tt in seen:
            continue
        seen.add(tt)
        sq, vals = check_lift_to_C_rational(tt)
        print("  t = %s  [%s]" % (tt, src))
        print("    F0 = %s квадрат=%s" % (vals[0], sq[0]))
        print("    F4 = %s квадрат=%s" % (vals[1], sq[1]))
        print("    F8 = %s квадрат=%s" % (vals[2], sq[2]))
        print("    => поднимается на C:", all(sq))
    print()
    print("ВАЖНО: отсутствие находок — НЕ доказательство H(Q) = пусто")
    print("и НЕ доказательство C(Q) = пусто.")


main()
