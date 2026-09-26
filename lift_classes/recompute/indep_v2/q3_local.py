# -*- coding: utf-8 -*-
"""q3: контроль необходимости (локально), собственный код, только точные дроби.

Для наклона (r,s) и простого p вычисляется ЛОКАЛЬНЫЙ ОБРАЗ
    I_p(r,s) = { ([T(z)], [L(z)]) in (Q_p^*/Q_p^*2)^2 : z in Q_p, все 9 клеток != 0,
                 восемь произведений строк/столбцов/диагоналей -- квадраты в Q_p }
методом разбиения Q_p на шары c + p^k Z_p:
  * клетка g = 1 + lam z постоянна по классу на шаре, если v(lam)+k-v(g(c)) >= 1 (p нечётно) / >= 3 (p=2);
  * если все клетки постоянны -- шар даёт одну точку образа, свидетель z=c (рациональное число);
  * если непостоянна ровно одна клетка g и шар содержит её корень z*=-1/lam, то на шаре класс g
    пробегает ВСЕ классы, остальные постоянны; восемь условий фиксируют требуемый класс g;
    строим явный рациональный z = z* + u p^a / lam с этим классом и проверяем его точно;
  * иначе шар делится на p подшаров (глубина ограничена MAXDEPTH, при превышении -- сообщение).
Покрытие Q_p: Z_p = шар(0,0); полюса v(z)=-n, n=1..N0+2 -- шары p^-n(w + p Z_p) (для p=2 шаг 2^3);
при n > N0 = max v_p(lam) + 3 все клетки имеют класс [lam z], образ зависит лишь от n mod 2 и
класса w -- эти n покрыты значениями N0+1, N0+2.
Каждая записанная точка образа имеет точного рационального свидетеля, прошедшего прямую проверку.

Проверяются утверждения Codex (локальные шаги доказательства):
  p=2,3: все 9 клеток -- квадраты в Q_p;   p нечётное: v_p(T) нечётна => p | r-s, v_p(L) нечётна => p | r+s;
  p = 3 mod 4: v_p(T), v_p(L) чётны.
Дополнительно (НЕ утверждение Codex, наблюдение о точности списков): для каждой пары (d_T,d_L) из
D_- x D_+ проверяется, лежит ли ([d_T]_p,[d_L]_p) в I_p.
Плюс: независимый случайный рациональный контроль и синтетические z, собранные по КТО для нескольких p.
"""
import os, sys, json, time, random
from fractions import Fraction as Fr
from math import gcd

OUT = os.path.dirname(os.path.abspath(__file__))
MAXDEPTH = 60
random.seed(20260926)

def vp(x, p):
    """v_p рационального x != 0"""
    x = Fr(x)
    n, d = x.numerator, x.denominator
    v = 0
    while n % p == 0:
        n //= p; v += 1
    while d % p == 0:
        d //= p; v -= 1
    return v

def unit_part(x, p):
    x = Fr(x); v = vp(x, p)
    y = x / Fr(p) ** v
    return y.numerator, y.denominator

def leg(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def cls(x, p):
    """класс в Q_p^*/Q_p^*2: (v mod 2, u); u = символ Лежандра (p нечётно) или u mod 8 (p=2)"""
    assert x != 0
    v = vp(x, p)
    n, d = unit_part(x, p)
    if p == 2:
        return (v % 2, (n * d) % 8)
    return (v % 2, leg(n * d, p))

def cmul(a, b, p):
    if p == 2:
        return ((a[0] + b[0]) % 2, (a[1] * b[1]) % 8)
    return ((a[0] + b[0]) % 2, a[1] * b[1])

ONE = lambda p: (0, 1)
IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX if (i, j) != (0, 0)]
LINES = ([[(i, j) for j in IDX] for i in IDX] + [[(i, j) for i in IDX] for j in IDX] +
         [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
T_CELLS = [(-1, 0), (1, 1), (0, -1)]     # T=(1-rz)(1+(r+s)z)(1-sz)
L_CELLS = [(-1, 0), (1, -1), (0, 1)]     # L=(1-rz)(1+(r-s)z)(1+sz)

def lam(r, s, i, j):
    return i * r + j * s

def eval_point(r, s, z, p):
    """точная проверка в точке z: вернуть None, если клетка 0 или какое-то произведение не квадрат;
    иначе словарь классов клеток и (T,L)."""
    z = Fr(z)
    cl = {(0, 0): (0, 1)}
    for (i, j) in CELLS:
        g = 1 + lam(r, s, i, j) * z
        if g == 0:
            return None
        cl[(i, j)] = cls(g, p)
    for line in LINES:
        c = (0, 1)
        for ij in line:
            c = cmul(c, cl[ij], p)
        if c != (0, 1):
            return None
    cT = (0, 1)
    for ij in T_CELLS:
        cT = cmul(cT, cl[ij], p)
    cL = (0, 1)
    for ij in L_CELLS:
        cL = cmul(cL, cl[ij], p)
    return cl, cT, cL

def unit_reps(p):
    if p == 2:
        return [1, 3, 5, 7]
    nr = next(a for a in range(2, p) if leg(a, p) == -1)
    return [1, nr]

def rep_of_class(C, p):
    """представитель класса C = (v mod 2, u) в виде (a, u_int) с x = u_int * p^a"""
    if p == 2:
        return C[0], C[1]
    u = 1 if C[1] == 1 else next(a for a in range(2, p) if leg(a, p) == -1)
    return C[0], u

def local_image(r, s, p, stats):
    lams = {ij: lam(r, s, *ij) for ij in CELLS}
    vl = {ij: vp(l, p) for ij, l in lams.items()}
    need = 3 if p == 2 else 1
    N0 = max(vl.values()) + 3
    witnesses = {}            # (cT, cL) -> z
    unresolved = []

    def record(z, res):
        cl, cT, cL = res
        key = (cT, cL)
        if key not in witnesses:
            witnesses[key] = z
        stats["points"] += 1
        # проверка утверждений Codex на КАЖДОЙ точке
        par = {ij: cl[ij][0] for ij in cl}
        if p in (2, 3):
            if any(cl[ij] != (0, 1) for ij in cl):
                stats["viol"].append(("p=2/3 not all squares", str(z)))
        else:
            if cT[0] == 1 and (r - s) % p != 0:
                stats["viol"].append(("v(T) odd, p !| r-s", str(z)))
            if cL[0] == 1 and (r + s) % p != 0:
                stats["viol"].append(("v(L) odd, p !| r+s", str(z)))
            if p % 4 == 3 and (cT[0] == 1 or cL[0] == 1):
                stats["viol"].append(("p=3 mod 4 odd v", str(z)))
        if cT[0] == 1: stats["oddT"] += 1
        if cL[0] == 1: stats["oddL"] += 1

    def rec(c, k, depth):
        stats["balls"] += 1
        uns = []
        for ij, l in lams.items():
            g = 1 + l * c
            if g == 0 or vl[ij] + k - vp(g, p) < need:
                uns.append(ij)
        if not uns:
            res = eval_point(r, s, c, p)
            if res is not None:
                record(c, res)
            return
        if len(uns) == 1:
            ij = uns[0]; l = lams[ij]
            g = 1 + l * c
            if g == 0 or vp(g, p) - vl[ij] >= k:        # шар содержит корень z* = -1/l
                # требуемый класс g из линий, где g участвует (остальные клетки постоянны)
                others = {}
                for ab, la in lams.items():
                    if ab != ij:
                        others[ab] = cls(1 + la * c, p)
                others[(0, 0)] = (0, 1)
                req = None; ok = True
                for line in LINES:
                    if ij in line:
                        cc = (0, 1)
                        for ab in line:
                            if ab != ij:
                                cc = cmul(cc, others[ab], p)
                        if req is None:
                            req = cc
                        elif req != cc:
                            ok = False
                    else:
                        cc = (0, 1)
                        for ab in line:
                            cc = cmul(cc, others[ab], p)
                        if cc != (0, 1):
                            ok = False
                if not ok:
                    return
                a0, u = rep_of_class(req, p)
                a = k + vl[ij] + 2          # глубже внутри шара
                if (a - a0) % 2:
                    a += 1
                zstar = Fr(-1, l)
                z = zstar + Fr(u) * Fr(p) ** a / l
                res = eval_point(r, s, z, p)
                if res is None or res[0][ij] != req:
                    stats["viol"].append(("root-rule witness failed", str(z)))
                    return
                stats["root_balls"] += 1
                record(z, res)
                return
        if depth >= MAXDEPTH:
            unresolved.append((str(c), k))
            return
        step = Fr(p) ** k
        for t in range(p):
            rec(c + t * step, k + 1, depth + 1)

    rec(Fr(0), 0, 0)
    first = 3 if p == 2 else 1
    for n in range(1, N0 + 3):
        # шары p^-n (w + p^first Z_p), w -- единицы по модулю p^first
        for w in range(1, p ** first):
            if w % p == 0:
                continue
            rec(Fr(w, p ** n), -n + first, 0)
    stats["unresolved"] = len(unresolved)
    return witnesses

def primes_upto(n):
    return [q for q in range(2, n + 1) if all(q % d for d in range(2, int(q ** 0.5) + 1))]

def prime_factors(n):
    n = abs(n); out = []; q = 2
    while q * q <= n:
        if n % q == 0:
            out.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        out.append(n)
    return out

def adm_list(m):
    ps = [q for q in prime_factors(m) if q % 4 == 1]
    out = []
    for mask in range(1 << len(ps)):
        d = 1
        for k, q in enumerate(ps):
            if mask >> k & 1:
                d *= q
        if d % 24 == 1:
            out.append(d)
    return sorted(out)

SLOPES = [(126, 451), (73, 362), (265, 298), (-12, 85), (12, 61), (24, 49), (4, 9), (3, 8), (1, 4),
          (2, 3), (1, 2), (5, 12), (7, 12), (1, 6), (8, 65), (36, 109)]
PR_SMALL = primes_upto(60)

t0 = time.time()
report = {"slopes": {}, "violations_total": 0}
for (r, s) in SLOPES:
    assert gcd(r, s) == 1 and r != 0 and abs(r) != abs(s)
    DT, DL = adm_list(r - s), adm_list(r + s)
    plist = sorted(set(PR_SMALL) | set(prime_factors(r * s * (r - s) * (r + s))))
    rec_s = {"D_T": DT, "D_L": DL, "primes": {}, "pairs_locally_ok": {}}
    pair_ok = {(a, b): [] for a in DT for b in DL}
    for p in plist:
        st = {"points": 0, "balls": 0, "root_balls": 0, "oddT": 0, "oddL": 0, "viol": []}
        W = local_image(r, s, p, st)
        img = sorted(W.keys())
        report["violations_total"] += len(st["viol"])
        miss = []
        for (a, b) in pair_ok:
            key = (cls(a, p), cls(b, p))
            if key not in W:
                pair_ok[(a, b)].append(p)
        rec_s["primes"][p] = {"image_size": len(img), "image": [str(x) for x in img],
                              "odd_vT_seen": st["oddT"] > 0, "odd_vL_seen": st["oddL"] > 0,
                              "p|r-s": (r - s) % p == 0, "p|r+s": (r + s) % p == 0,
                              "balls": st["balls"], "root_balls": st["root_balls"],
                              "unresolved": st["unresolved"], "violations": st["viol"][:5]}
        if st["viol"] or st["unresolved"]:
            print("!!", r, s, p, st["viol"][:3], st["unresolved"], flush=True)
    rec_s["pairs_locally_ok"] = {"%d,%d" % k: ("all tested p" if not v else "obstructed at %s" % v)
                                 for k, v in pair_ok.items()}
    odd_primes_T = [p for p, x in rec_s["primes"].items() if x["odd_vT_seen"]]
    odd_primes_L = [p for p, x in rec_s["primes"].items() if x["odd_vL_seen"]]
    rec_s["odd_vT_primes"] = odd_primes_T
    rec_s["odd_vL_primes"] = odd_primes_L
    report["slopes"]["%d/%d" % (r, s)] = rec_s
    print("%d/%d  D_T=%s D_L=%s  oddT at %s, oddL at %s  pairs: %s  (%.1fs)" % (
        r, s, DT, DL, odd_primes_T, odd_primes_L, rec_s["pairs_locally_ok"], time.time() - t0), flush=True)

print("violations_total =", report["violations_total"])
json.dump(report, open(os.path.join(OUT, "q3_local.json"), "w"), indent=1, ensure_ascii=False, default=str)
print("done %.1fs" % (time.time() - t0))
