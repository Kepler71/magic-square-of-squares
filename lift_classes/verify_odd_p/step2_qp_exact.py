#!/usr/bin/env python3
"""Точная Q_p-проверка Шага 2 (нечётное p) на рациональных z = u * p^e.
Только целочисленная арифметика: класс клетки в Q_p*/Q_p*^2 = (v_p mod 2, символ Лежандра единичной части).
Если все 8 произведений по линиям арифметической сетки — квадраты в Q_p, проверяем:
  A1 противоположные клетки одного полного Q_p-класса;
  A3 e>=0 (z p-целое) => все оценки чётны;
  A4 v_p(T) нечётна => p|(r-s);  v_p(L) нечётна => p|(r+s)   (T, L считаются напрямую);
  A5 полюс и (p|r или p|s) => n чётно и S=U=0;
  A6 полюс, r,s единицы, n нечётно => U=1 и ровно одна из r+-s делится на p, S соответствует;
  A7 полюс, r,s единицы, n чётно => S=U=0;
  A8 формула q_lambda (v<n: n-v mod 2; v>=n: 0) для обеих клеток каждой пары;
  A9 p = 3 mod 4 => нетривиальных S, S+U нет (Шаг 5 локально);
  A10 p = 3 => полюса нет, z in 3Z_3 (Шаг 4).
Аргументы: список простых. Все циклы конечны.
"""
import sys, time, math, random

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {c: k for k, c in enumerate(CELLS)}
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)] +
         [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)] +
         [[IDX[(-1, -1)], IDX[(0, 0)], IDX[(1, 1)]], [IDX[(-1, 1)], IDX[(0, 0)], IDX[(1, -1)]]])
OPP = [(IDX[(i, j)], IDX[(-i, -j)]) for (i, j) in CELLS]

def vp_int(N, p):
    v = 0
    while N % p == 0:          # N != 0 guaranteed by caller; loop bounded by log_p |N|
        N //= p
        v += 1
    return v, N

def vp_lam(lam, p):
    return vp_int(lam, p)[0]

def pairs_for(p, R):
    out = set()
    for r in range(-R, R + 1):
        for s in range(1, R + 1):
            if r == 0 or r == s or r == -s or math.gcd(r, s) != 1:
                continue
            out.add((r, s))
    # special pairs with high valuations of r, s, r+s, r-s
    for k in (1, 2, 3, 4):
        q = p ** k
        for m in (1, 2, 3):
            for t in (1, 2, 3, 4, 5, 7):
                for sg in (1, -1):
                    cand = [(sg * q * m, t), (t, sg * q * m), (t, sg * q * m - t), (t, t - sg * q * m),
                            (sg * q * m + t, t), (t - sg * q * m, t)]
                    for (r, s) in cand:
                        if r == 0 or s == 0 or r == s or r == -s or math.gcd(r, s) != 1:
                            continue
                        if s < 0:
                            r, s = -r, -s
                        out.add((r, s))
    return sorted(out)

def units_for(p, cap):
    M = p ** 3
    us = [u for u in range(1, M) if u % p]
    if len(us) > cap:
        rng = random.Random(1000 + p)
        us = sorted(rng.sample(us, cap))
    return us

def run_prime(p, R, ucap, EMIN=-4, EMAX=2):
    t0 = time.time()
    prs = pairs_for(p, R)
    us = units_for(p, ucap)
    stats = dict(tested=0, passed=0, pole_passed=0, S1=0, SU1=0, S1_n3plus=0, pr_or_ps_pole_passed=0,
                 eqval_pole_passed=0, highval_rs=0)
    for idx_pair, (r, s) in enumerate(prs):
        lam = [i * r + j * s for (i, j) in CELLS]
        vr, vs, vrp, vrm = vp_lam(r, p), vp_lam(s, p), vp_lam(r + s, p), vp_lam(r - s, p)
        for e in range(EMIN, EMAX + 1):
            n = -e
            for u in us:
                cls = []
                bad = False
                for L in lam:
                    if e >= 0:
                        N = 1 + L * u * p ** e
                        shift = 0
                    else:
                        N = p ** n + L * u
                        shift = n
                    if N == 0:
                        bad = True
                        break
                    v, unit = vp_int(N, p)
                    chi = 0 if pow(unit % p, (p - 1) // 2, p) == 1 else 1
                    cls.append(((v - shift), chi))
                if bad:
                    continue
                stats['tested'] += 1
                ok = True
                for Lz in LINES:
                    if sum(cls[k][0] for k in Lz) % 2 or sum(cls[k][1] for k in Lz) % 2:
                        ok = False
                        break
                if not ok:
                    continue
                stats['passed'] += 1
                par = [c[0] % 2 for c in cls]
                # A1
                for a, b in OPP:
                    assert cls[a][0] % 2 == cls[b][0] % 2 and cls[a][1] == cls[b][1], ("A1", p, r, s, e, u)
                S = par[IDX[(1, 1)]]
                U = par[IDX[(1, 0)]]
                SU = par[IDX[(1, -1)]]
                assert SU == (S + U) % 2 and par[IDX[(0, 1)]] == U
                # direct T, L parities
                vT = cls[IDX[(-1, 0)]][0] + cls[IDX[(1, 1)]][0] + cls[IDX[(0, -1)]][0]
                vL = cls[IDX[(-1, 0)]][0] + cls[IDX[(1, -1)]][0] + cls[IDX[(0, 1)]][0]
                assert vT % 2 == S and vL % 2 == SU
                # A3
                if e >= 0:
                    assert S == 0 and U == 0, ("A3", p, r, s, e, u)
                    if p == 3 and e == 0:
                        raise AssertionError(("A10 unit z at 3 passed", r, s, e, u))
                    continue
                stats['pole_passed'] += 1
                # A4
                if S:
                    stats['S1'] += 1
                    assert (r - s) % p == 0, ("A4T", p, r, s, e, u)
                    if n >= 3:
                        stats['S1_n3plus'] += 1
                if SU:
                    stats['SU1'] += 1
                    assert (r + s) % p == 0, ("A4L", p, r, s, e, u)
                # A5
                if vr > 0 or vs > 0:
                    stats['pr_or_ps_pole_passed'] += 1
                    assert n % 2 == 0 and S == 0 and U == 0, ("A5", p, r, s, e, u)
                else:
                    if n % 2:
                        # A6
                        assert U == 1, ("A6U", p, r, s, e, u)
                        assert (vrp > 0) != (vrm > 0), ("A6x", p, r, s, e, u)
                        assert S == (1 if vrm > 0 else 0), ("A6S", p, r, s, e, u)
                    else:
                        assert S == 0 and U == 0, ("A7", p, r, s, e, u)
                # A8
                for lamv, (c1, c2) in ((r, ((1, 0), (-1, 0))), (s, ((0, 1), (0, -1))),
                                       (r + s, ((1, 1), (-1, -1))), (r - s, ((1, -1), (-1, 1)))):
                    vl = vp_lam(lamv, p)
                    pred = (n - vl) % 2 if vl < n else 0
                    if vl == n:
                        stats['eqval_pole_passed'] += 1
                    if vl >= 2:
                        stats['highval_rs'] += 1
                    assert par[IDX[c1]] == pred and par[IDX[c2]] == pred, ("A8", p, r, s, e, u, lamv)
                # A9
                if p % 4 == 3:
                    assert S == 0 and SU == 0, ("A9", p, r, s, e, u)
                # A10
                if p == 3:
                    raise AssertionError(("A10 pole at 3 passed", r, s, e, u))
        if idx_pair % 200 == 0:
            print("  p=%d pair %d/%d tested=%d passed=%d pole_passed=%d S1=%d SU1=%d  %.0fs"
                  % (p, idx_pair, len(prs), stats['tested'], stats['passed'], stats['pole_passed'],
                     stats['S1'], stats['SU1'], time.time() - t0), flush=True)
    print("p=%d DONE pairs=%d units=%d e in [%d,%d]: %s  %.0fs" % (p, len(prs), len(us), EMIN, EMAX, stats,
                                                                  time.time() - t0), flush=True)
    return stats

if __name__ == "__main__":
    ps = [int(x) for x in sys.argv[1].split(',')]
    R = int(sys.argv[2]) if len(sys.argv) > 2 else 18
    ucap = int(sys.argv[3]) if len(sys.argv) > 3 else 120
    for p in ps:
        run_prime(p, R, ucap)
    print("ALL ASSERTIONS PASSED for", ps, flush=True)
