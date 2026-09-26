# -*- coding: utf-8 -*-
"""t5: (A) СИНТЕТИЧЕСКИЕ РАЦИОНАЛЬНЫЕ z (КТО): для наклона (r,s) и набора простых P строим ОДНО рациональное z,
    при котором восемь произведений арифметической сетки -- квадраты в Q_p ОДНОВРЕМЕННО для всех p из P;
    локальные решения подбираются случайно (включая полюса и нечётные оценки), склеиваются по КТО.
    Квадратность и классы проверяются ТОЛЬКО через Sage (QQ.is_padic_square, valuation) -- другой путь, чем в t3.
    Проверяем утверждения Codex (P1 матрица классов, P2 полный подъём при 2 и 3, P3 носитель нечётных оценок,
    P4 p=3 mod 4 => оценки чётны) и свои уточнения (L7: T,L -- квадраты в Q_7; H: нечётность при p|r-s с чётной
    v_p(r-s) лишь при p=1 mod 8; Ip: для p=3 mod 4 пара символов ([T],[L]) лежит в образе по z mod p).
(B) Лемма p=7 на случайных рациональных z (Sage): 8 произведений -- квадраты в Q_7 => T, L -- квадраты в Q_7.
(C) Контроль метода ранга: столбец-тавтология B*C не повышает ранг, посторонняя форма a+b+2c повышает.
Все циклы ограничены.
"""
import os, sys, json, random, time
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
from math import gcd
from sage.all import QQ, ZZ, GF, matrix, CRT_list, factor, kronecker, prime_range
t0 = time.time()
IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]
LINES = ([[(i, j) for j in IDX] for i in IDX] + [[(i, j) for i in IDX] for j in IDX]
         + [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
def cells(r, s, z): return {(i, j): 1 + (i * r + j * s) * z for (i, j) in CELLS}
def T_of(f): return f[(-1, 0)] * f[(1, 1)] * f[(0, -1)]      # (1-rz)(1+(r+s)z)(1-sz)
def L_of(f): return f[(-1, 0)] * f[(1, -1)] * f[(0, 1)]      # (1-rz)(1+(r-s)z)(1+sz)
def eight_ok(f, p):
    return all(QQ(f[l[0]] * f[l[1]] * f[l[2]]).is_padic_square(p) for l in LINES)
def sq(x, p): return QQ(x).is_padic_square(p)
def same_class(x, y, p): return sq(QQ(x) * QQ(y), p)
def vpar(x, p): return int(QQ(x).valuation(p)) % 2
def chi_unit(x, p):
    x = QQ(x); u = x / QQ(p) ** x.valuation(p)
    return int(kronecker(u.numerator() * u.denominator(), p))

def c_of(mu, z0, p):
    a, b = kronecker(1 + mu * z0, p), kronecker(1 - mu * z0, p)
    if a and b and a != b: return None
    return int(a if a else b)
def image_Ip(r, s, p):
    img = set()
    for z0 in range(p):
        c = {}
        for mu in (r, s, r + s, r - s):
            v = c_of(mu, z0, p)
            if v is None: break
            c[mu] = v
        else:
            if c[r] == c[s] and c[r + s] * c[r - s] == c[r]:
                img.add((c[r + s], c[r - s]))
    return img

def check_point(r, s, z, p, Ip_cache):
    """возвращает список нарушений для z, прошедшего 8 условий в Q_p"""
    f = cells(r, s, z); bad = []
    S_rep, U_rep = f[(-1, -1)], f[(-1, 0)]
    # P1: матрица классов
    exp = {(-1, -1): [S_rep], (1, 1): [S_rep], (-1, 0): [U_rep], (1, 0): [U_rep], (0, -1): [U_rep], (0, 1): [U_rep],
           (-1, 1): [S_rep, U_rep], (1, -1): [S_rep, U_rep], (0, 0): [1]}
    for c, reps in exp.items():
        prod = 1
        for x in reps: prod *= x
        if not same_class(f[c], prod, p): bad.append("P1")
    T, L = T_of(f), L_of(f)
    if not same_class(T, S_rep, p) or not same_class(L, S_rep * U_rep, p): bad.append("P1-TL")
    vT, vL = vpar(T, p), vpar(L, p)
    if p in (2, 3):
        if not all(sq(x, p) for x in f.values()): bad.append("P2")
        if QQ(z).valuation(p) < (3 if p == 2 else 1): bad.append("P2v")
    else:
        if vT and not ((r - s) % p == 0 and p % 4 == 1): bad.append("P3T")
        if vL and not ((r + s) % p == 0 and p % 4 == 1): bad.append("P3L")
        if p % 4 == 3 and any(vpar(x, p) for x in f.values()): bad.append("P4")
        if p == 7 and not (sq(T, 7) and sq(L, 7)): bad.append("L7")
        if vT and not (ZZ(r - s).valuation(p) % 2 == 1 or p % 8 == 1): bad.append("H-T")
        if vL and not (ZZ(r + s).valuation(p) % 2 == 1 or p % 8 == 1): bad.append("H-L")
        if p % 4 == 3:
            if p not in Ip_cache: Ip_cache[p] = image_Ip(r, s, p)
            if (chi_unit(T, p), chi_unit(L, p)) not in Ip_cache[p]: bad.append("Ip")
    return bad, vT, vL

def _cls_int(n, p):
    """быстрый класс целого n != 0 в Q_p^*/Q_p^{*2}: маска (v mod 2, признак невычета / (a,b) при p=2)"""
    v = 0
    while n % p == 0:
        n //= p; v += 1
    if p == 2:
        u8 = n % 8
        return (v & 1) | ((1 if u8 in (3, 7) else 0) << 1) | ((1 if u8 in (3, 5) else 0) << 2)
    return (v & 1) | ((0 if pow(n % p, (p - 1) // 2, p) == 1 else 1) << 1)

def _fast_cells(r, s, p, e, a):
    out = {}
    for (i, j) in CELLS:
        lam = i * r + j * s
        if e >= 0:
            n = 1 + lam * a * p ** e
            if n == 0: return None
            out[(i, j)] = _cls_int(n, p)
        else:
            n = p ** (-e) + lam * a
            if n == 0: return None
            out[(i, j)] = _cls_int(n, p) ^ ((-e) & 1)
    return out

LOCAL_CACHE = {}
def local_solutions(r, s, p, rng, want, K=6):
    """случайные локальные решения z_p = p^e * a (подбор -- своей быстрой арифметикой; итог проверяет Sage)"""
    key = (r, s, p)
    if key in LOCAL_CACHE: return LOCAL_CACHE[key]
    sols = {}
    for _ in range(60000):                                   # ограничено
        e = rng.choice([-3, -2, -1, 0, 0, 1, 1, 2, 3, 4])
        a = rng.randrange(1, p ** K) * rng.choice([1, -1])
        if a % p == 0: continue
        C = _fast_cells(r, s, p, e, a)
        if C is None or any(C[l[0]] ^ C[l[1]] ^ C[l[2]] for l in LINES): continue
        cT = C[(-1, 0)] ^ C[(1, 1)] ^ C[(0, -1)]; cL = C[(-1, 0)] ^ C[(1, -1)] ^ C[(0, 1)]
        typ = (cT & 1, cL & 1, e < 0)
        sols.setdefault(typ, [])
        if len(sols[typ]) < 6: sols[typ].append((e, a))
        if sum(len(v) for v in sols.values()) >= want and len(sols) >= 1 and _ > 20000: break
    LOCAL_CACHE[key] = sols
    return sols

def build_global(r, s, P, rng, K=6):
    choice = {}
    for p in P:
        sols = local_solutions(r, s, p, rng, 24, K)
        if not sols: return None, None
        typs = sorted(sols)
        typ = rng.choice(typs)                              # равновероятно по типам -> чаще редкие
        choice[p] = (typ,) + rng.choice(sols[typ])
    D = 1
    for p, (typ, e, a) in choice.items():
        if e < 0: D *= p ** (-e)
    res, mods = [], []
    for p, (typ, e, a) in choice.items():
        zp = QQ(a) * QQ(p) ** e
        target = zp * D                                     # p-целое
        m = p ** (int(target.valuation(p)) + K + 4)
        num, den = target.numerator(), target.denominator()
        res.append(int(num * pow(int(den), -1, m) % m)); mods.append(m)
    x = CRT_list(res, mods)
    z = QQ(x) / D
    return z, choice

SLOPES = [(126, 451), (-126, 451), (73, 362), (-73, 362), (265, 298), (12, 85), (2, 3), (5, 8), (3, 10), (7, 18),
          (1, 24), (19, 30), (11, 14), (4, 9), (362, 363), (61, 64), (1, 97), (24, 49)]
BASE = list(prime_range(2, 44))
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 12
rng = random.Random(260926)
out = {"A": {}, "B": None, "C": None}
totals = dict(z_built=0, z_ok=0, checks=0, viol=0, odd_T=0, odd_L=0)
for (r, s) in SLOPES:
    assert gcd(r, s) == 1 and abs(r) != abs(s)
    P = sorted(set(BASE) | set(int(q) for q, _ in factor(abs(r - s))) | set(int(q) for q, _ in factor(abs(r + s))))
    rec = dict(P=P, built=0, ok=0, viol=[], oddT_primes=set(), oddL_primes=set(), types_used={})
    Ip_cache = {}
    for k in range(NZ):
        z, choice = build_global(r, s, P, rng)
        if z is None: continue
        rec["built"] += 1; totals["z_built"] += 1
        f = cells(r, s, z)
        if any(x == 0 for x in f.values()): continue
        okP = [p for p in P if eight_ok(f, p)]
        if len(okP) != len(P):
            rec["viol"].append(("CRT-precision", str(z)[:40], sorted(set(P) - set(okP)))); continue
        rec["ok"] += 1; totals["z_ok"] += 1
        for p in P:
            bad, vT, vL = check_point(r, s, z, p, Ip_cache)
            totals["checks"] += 1
            if vT: rec["oddT_primes"].add(p); totals["odd_T"] += 1
            if vL: rec["oddL_primes"].add(p); totals["odd_L"] += 1
            if bad:
                totals["viol"] += 1
                if len(rec["viol"]) < 5: rec["viol"].append((p, bad, str(z)[:60]))
        if k == 0:
            rec["example_z_height_digits"] = len(str(z))
    rec["oddT_primes"] = sorted(rec["oddT_primes"]); rec["oddL_primes"] = sorted(rec["oddL_primes"])
    out["A"]["%d/%d" % (r, s)] = rec
    print("(A) %d/%d: z построено %d, годны при всех %d простых: %d, нечётная v(T) при %s, v(L) при %s, нарушений %d  t=%.1fs"
          % (r, s, rec["built"], len(P), rec["ok"], rec["oddT_primes"], rec["oddL_primes"],
             len(rec["viol"]), time.time() - t0), flush=True)
print("(A) ИТОГ:", totals, flush=True)
out["A_totals"] = totals

# (B) лемма p=7 на случайных рациональных z
rngB = random.Random(77)
tested = passed = bad = 0
for _ in range(60000):
    s = rngB.randrange(2, 3000); r = rngB.randrange(-s + 1, s)
    if r == 0 or gcd(r, s) != 1: continue
    e = rngB.choice([-2, -1, 0, 1, 1, 2, 3])
    num = rngB.randrange(1, 10 ** 6) * rngB.choice([1, -1]); den = rngB.randrange(1, 10 ** 4)
    z = QQ(num) / QQ(den) * QQ(7) ** e
    f = cells(r, s, z)
    if any(x == 0 for x in f.values()): continue
    tested += 1
    if not eight_ok(f, 7): continue
    passed += 1
    if not (sq(T_of(f), 7) and sq(L_of(f), 7) and all(sq(x, 7) for x in f.values())): bad += 1
out["B"] = dict(tested=tested, passed=passed, nontrivial=bad)
print("(B) p=7, случайные рациональные z: проверено %d, 8 произведений -- квадраты в Q_7: %d, из них T или L или клетка не квадрат: %d  t=%.1fs"
      % (tested, passed, bad, time.time() - t0), flush=True)

# (C) контроль метода ранга (свежие точки, p=10007)
p = 10007; F = GF(p); rngC = random.Random(5)
def leg(x): return 0 if int(kronecker(x % p, p)) == 1 else 1
rowsT, rowsX, n = [], [], 0
for _ in range(100000):
    if len(rowsT) >= 200: break
    a, b, c = (rngC.randrange(1, p) for _ in range(3))
    cv = {(i, j): (a + i * b + j * c) % p for (i, j) in CELLS}
    if any(kronecker(v, p) != 1 for v in cv.values()) or (a + b + 2 * c) % p == 0: continue
    R = {k: (int(F(v).sqrt()) if rngC.random() < .5 else p - int(F(v).sqrt())) for k, v in cv.items()}
    Ls = [(R[(-1, j)], R[(0, j)], R[(1, j)]) for j in IDX] + [(R[(i, -1)], R[(i, 0)], R[(i, 1)]) for i in IDX] \
         + [(R[(-1, -1)], R[(0, 0)], R[(1, 1)]), (R[(-1, 1)], R[(0, 0)], R[(1, -1)])]
    vals = []
    for (al, ga, be) in Ls: vals += [(al + ga) * (ga + be), (al + ga) * (al + be)]
    vals += [b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c] + [R[k] for k in CELLS]
    if any(v % p == 0 for v in vals): continue
    al, ga, be = Ls[0]
    base = [leg(v) for v in vals] + [1]
    rowsT.append(base + [leg((ga + be) * (al + be))]); rowsX.append(base + [leg(a + b + 2 * c)])
rT = int(matrix(GF(2), rowsT).rank()); rX = int(matrix(GF(2), rowsX).rank())
out["C"] = dict(points=len(rowsT), rank_with_tautology_col=rT, rank_with_foreign_form=rX)
print("(C) 35 столбцов: с тавтологией BC ранг %d (ожидается 34), с посторонней формой a+b+2c ранг %d (ожидается 35)" % (rT, rX), flush=True)
out["time_s"] = round(time.time() - t0, 1)
json.dump(out, open("/home/kep/magicKube/lift_classes/recompute/v3_a51edd4e/t5_crt_local.json", "w"), indent=1, default=str)
print("время %.1f с" % (time.time() - t0))
