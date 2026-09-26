# -*- coding: utf-8 -*-
"""t3: локальный контроль необходимости теоремы Codex (частичный контроль, НЕ доказательство).

Для наклона (r,s) и простого p перебираем z = p^e * a (a -- целое, взаимно простое с p; e в [EMIN,EMAX]).
Клетки f_ij = 1 + (i r + j s) z. Квадратные классы в Q_p^*/Q_p^{*2} считаются ТОЧНО для рациональных чисел:
  класс кодируется битовой маской (XOR = умножение классов):
    нечётное p: бит0 = v mod 2, бит1 = [единичная часть -- невычет mod p];
    p = 2:      бит0 = v mod 2, (бит1,бит2) = (a,b) из u = (-1)^a 5^b mod 8.
Если восемь произведений (3 строки, 3 столбца, 2 диагонали арифметической сетки) -- квадраты в Q_p,
проверяем утверждения теоремы:
  P1 классы клеток образуют матрицу S U S+U / U 0 U / S+U U S, [T]=S, [L]=S+U;
  P2 p in {2,3}: все 9 клеток -- квадраты, v_2(z)>=3 / v_3(z)>=1;
  P3 p нечётно: v_p(T) нечётна => p|(r-s) и p=1 mod 4; v_p(L) нечётна => p|(r+s) и p=1 mod 4;
  P4 p=3 mod 4: все v_p(клеток) чётны.
Гипотеза H (своё уточнение, проверяется здесь же): v_p(T) нечётна => v_p(r-s) нечётна или p=1 mod 8
  (аналогично L с r+s); и обратно -- если это разрешено, пример находится.
Вещественное место: 8 произведений > 0  =>  все клетки > 0.
Второй реализацией (Sage QQ.is_padic_square) перепроверяется случайная подвыборка.
Все циклы ограничены явно.
"""
import os, sys, json, random, time
from fractions import Fraction
from math import gcd
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")

IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]
LINES = ([[(i, j) for j in IDX] for i in IDX] + [[(i, j) for i in IDX] for j in IDX]
         + [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
T_CELLS = [(-1, 0), (1, 1), (0, -1)]   # T=(1-rz)(1+(r+s)z)(1-sz)
L_CELLS = [(-1, 0), (1, -1), (0, 1)]   # L=(1-rz)(1+(r-s)z)(1+sz)

def vp(n, p):
    assert n != 0
    v = 0
    while n % p == 0:
        n //= p; v += 1
    return v

def cls_int(n, p):
    """класс ненулевого целого n в Q_p^*/Q_p^{*2} как битовая маска"""
    v = vp(n, p)
    u = n // p ** v
    if p == 2:
        u8 = u % 8
        a = 1 if u8 in (3, 7) else 0
        b = 1 if u8 in (3, 5) else 0
        return (v & 1) | (a << 1) | (b << 2)
    leg = pow(u % p, (p - 1) // 2, p)
    return (v & 1) | ((0 if leg == 1 else 1) << 1)

def cls_frac(x, p):
    x = Fraction(x)
    return cls_int(x.numerator, p) ^ cls_int(x.denominator, p)

def cell_classes(r, s, p, e, a):
    """классы 9 клеток при z = p^e a; None если какая-то клетка = 0"""
    out = {}
    for (i, j) in CELLS:
        lam = i * r + j * s
        if e >= 0:
            n = 1 + lam * a * p ** e
            if n == 0: return None
            out[(i, j)] = cls_int(n, p)
        else:
            k = -e
            n = p ** k + lam * a
            if n == 0: return None
            out[(i, j)] = cls_int(n, p) ^ (k & 1)
    return out

def eight_ok(C):
    for ln in LINES:
        if C[ln[0]] ^ C[ln[1]] ^ C[ln[2]]:
            return False
    return True

def allowed_odd(p, m):
    """по теореме Codex: может ли p делить класс (m = r-s для T, r+s для L)"""
    return p % 2 == 1 and p != 3 and m % p == 0 and p % 4 == 1

def allowed_odd_H(p, m):
    if not allowed_odd(p, m): return False
    return (vp(m, p) % 2 == 1) or (p % 8 == 1)

def prime_factors(n):
    n = abs(n); f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f

def samples_a(p, r, s, rng, NA):
    """набор единиц a: все единицы mod p^k (p^k<=NA) либо случайная выборка, + прицельные почти-корни"""
    k = 1
    while p ** (k + 1) <= NA: k += 1
    M = p ** k
    if M <= NA and M >= 2 * p:
        base = [a for a in range(1, M) if a % p]
    else:
        base = [x for x in (rng.randrange(1, p * p) for _ in range(NA)) if x % p]
    extra = []
    for lam in {r, s, r + s, r - s}:
        if lam % p == 0: continue
        for t in range(1, 8):
            mod = p ** t
            a0 = (-pow(lam, -1, mod)) % mod
            for w in range(0, 3):
                for sg in (1, -1):
                    extra.append(sg * (a0 + w * mod))
    # отрицательные a -- тоже
    base = base + [-a for a in base[: len(base) // 2]]
    return [a for a in base + extra if a % p]

def run_local(r, s, primes, NA, EMIN, EMAX, rng, sage_sub, log):
    st = {}
    for p in primes:
        rec = dict(tests=0, passed=0, viol=[], T_odd=0, L_odd=0, nontriv_unit=0,
                   e_passed=set(), witness_T=None, witness_L=None)
        A = samples_a(p, r, s, rng, NA)
        for e in range(EMIN, EMAX + 1):
            for a in A:
                C = cell_classes(r, s, p, e, a)
                if C is None: continue
                rec["tests"] += 1
                if not eight_ok(C): continue
                rec["passed"] += 1; rec["e_passed"].add(e)
                S = C[(-1, -1)]; U = C[(-1, 0)]
                pat = {(-1, -1): S, (-1, 0): U, (-1, 1): S ^ U, (0, -1): U, (0, 0): 0, (0, 1): U,
                       (1, -1): S ^ U, (1, 0): U, (1, 1): S}
                cT = C[T_CELLS[0]] ^ C[T_CELLS[1]] ^ C[T_CELLS[2]]
                cL = C[L_CELLS[0]] ^ C[L_CELLS[1]] ^ C[L_CELLS[2]]
                bad = []
                if C != pat: bad.append("P1-pattern")
                if cT != S or cL != (S ^ U): bad.append("P1-TL")
                if p in (2, 3):
                    if any(C.values()): bad.append("P2-nontrivial")
                    if e < (3 if p == 2 else 1): bad.append("P2-valuation")
                else:
                    if cT & 1 and not allowed_odd(p, r - s): bad.append("P3-T")
                    if cL & 1 and not allowed_odd(p, r + s): bad.append("P3-L")
                    if p % 4 == 3 and any(c & 1 for c in C.values()): bad.append("P4")
                    if cT & 1 and not allowed_odd_H(p, r - s): bad.append("H-T")
                    if cL & 1 and not allowed_odd_H(p, r + s): bad.append("H-L")
                if cT & 1:
                    rec["T_odd"] += 1
                    if rec["witness_T"] is None: rec["witness_T"] = (e, a)
                if cL & 1:
                    rec["L_odd"] += 1
                    if rec["witness_L"] is None: rec["witness_L"] = (e, a)
                if (S | U) & ~1: rec["nontriv_unit"] += 1
                if bad and len(rec["viol"]) < 5:
                    rec["viol"].append((e, a, bad))
                if bad: rec.setdefault("nviol", 0); rec["nviol"] = rec.get("nviol", 0) + 1
                if len(sage_sub) < 4000 and rng.random() < 0.02:
                    sage_sub.append((r, s, p, e, a, True))
                elif len(sage_sub) < 8000 and rng.random() < 0.0005:
                    sage_sub.append((r, s, p, e, a, None))
        # реализуемость, предсказанная гипотезой H
        predT = p % 2 == 1 and allowed_odd_H(p, r - s)
        predL = p % 2 == 1 and allowed_odd_H(p, r + s)
        rec["pred_T_odd_possible_H"] = predT; rec["pred_L_odd_possible_H"] = predL
        rec["codex_T_odd_allowed"] = allowed_odd(p, r - s); rec["codex_L_odd_allowed"] = allowed_odd(p, r + s)
        rec["H_realization_mismatch"] = ((rec["T_odd"] > 0) != predT) or ((rec["L_odd"] > 0) != predL)
        rec["e_passed"] = sorted(rec["e_passed"])
        st[p] = rec
        log("  %s p=%d tests=%d passed=%d T_odd=%d L_odd=%d viol=%d e=%s H_mismatch=%s" % (
            "%d/%d" % (r, s), p, rec["tests"], rec["passed"], rec["T_odd"], rec["L_odd"],
            rec.get("nviol", 0), rec["e_passed"][:6], rec["H_realization_mismatch"]))
    return st

def real_check(r, s, rng, N=40000):
    tested = pos8 = bad = 0
    for _ in range(N):
        den = rng.choice([1, 2, 3, 7, 10, 97, 1000, 10 ** 6])
        num = rng.randrange(-10 ** 7, 10 ** 7)
        z = Fraction(num, den * rng.choice([1, 1, 1000, abs(r) + abs(s)]))
        f = {c: 1 + (c[0] * r + c[1] * s) * z for c in CELLS}
        if any(v == 0 for v in f.values()): continue
        tested += 1
        if all(f[l[0]] * f[l[1]] * f[l[2]] > 0 for l in LINES):
            pos8 += 1
            if any(v < 0 for v in f.values()): bad += 1
    return dict(tested=tested, eight_positive=pos8, negative_cell=bad)

def sage_crosscheck(sub):
    from sage.all import QQ
    mism = 0
    for (r, s, p, e, a, passed) in sub:
        z = Fraction(a) * Fraction(p) ** e
        f = {c: 1 + (c[0] * r + c[1] * s) * z for c in CELLS}
        ok8 = all(QQ(f[l[0]] * f[l[1]] * f[l[2]]).is_padic_square(p) for l in LINES)
        # сверка поклеточных классов: клетка -- квадрат <=> маска 0
        C = cell_classes(r, s, p, e, a)
        for c in CELLS:
            if QQ(f[c]).is_padic_square(p) != (C[c] == 0): mism += 1
            if cls_frac(f[c], p) != C[c]: mism += 1
        if ok8 != eight_ok(C): mism += 1
        if passed and not ok8: mism += 1
    return mism

SLOPES_DEFAULT = [(126, 451), (-126, 451), (73, 362), (-73, 362), (265, 298), (12, 85), (362, 363),
                  (7, 18), (11, 14), (2, 3), (4, 9), (8, 9), (5, 8), (3, 10), (1, 24), (19, 30)]
BASE_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43]

if __name__ == "__main__":
    part = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    nparts = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    NA = int(sys.argv[3]) if len(sys.argv) > 3 else 2500
    slopes = SLOPES_DEFAULT[part::nparts]
    rng = random.Random(20260926 + part)
    t0 = time.time()
    logf = open("t3_local_part%d.log" % part, "w")
    def log(msg):
        line = "[%6.1fs] %s" % (time.time() - t0, msg)
        print(line, flush=True); logf.write(line + "\n"); logf.flush()
    out = {}; sage_sub = []
    for (r, s) in slopes:
        assert gcd(r, s) == 1 and abs(r) != abs(s)
        primes = sorted(set(BASE_PRIMES) | set(prime_factors(r - s)) | set(prime_factors(r + s)))
        st = run_local(r, s, primes, NA, -7, 6, rng, sage_sub, log)
        out["%d/%d" % (r, s)] = dict(local=st, real=real_check(r, s, rng))
        log("%d/%d real: %s" % (r, s, out["%d/%d" % (r, s)]["real"]))
    mism = sage_crosscheck(sage_sub)
    log("Sage-сверка подвыборки: случаев %d, расхождений %d" % (len(sage_sub), mism))
    tot = dict(tests=0, passed=0, viol=0, T_odd=0, L_odd=0, H_mismatch=0, real_bad=0)
    for k, v in out.items():
        for p, rec in v["local"].items():
            tot["tests"] += rec["tests"]; tot["passed"] += rec["passed"]
            tot["viol"] += rec.get("nviol", 0); tot["T_odd"] += rec["T_odd"]; tot["L_odd"] += rec["L_odd"]
            tot["H_mismatch"] += int(rec["H_realization_mismatch"])
        tot["real_bad"] += v["real"]["negative_cell"]
    tot["sage_sub"] = len(sage_sub); tot["sage_mismatch"] = mism; tot["time_s"] = round(time.time() - t0, 1)
    log("ИТОГ части %d: %s" % (part, tot))
    json.dump(dict(total=tot, per_slope=out), open("t3_local_part%d.json" % part, "w"), indent=1, default=str)
