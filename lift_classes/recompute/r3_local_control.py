# -*- coding: utf-8 -*-
"""r3: локальный контроль необходимости (частичный, НЕ доказательство).
Для набора наклонов r/s и простых p строим рациональные z (z = p^k * a/b, a,b взаимно просты с p),
точно (целочисленная арифметика) проверяем, что 8 произведений арифметической сетки f_ij=1+(ir+js)z
являются квадратами в Q_p, и для прошедших проверяем локальные утверждения теоремы Codex:
  (M)  матрица классов: противоположные клетки одного класса, [T]=класс клетки (1,1), [L]=класс клетки (1,-1),
       осевые клетки одного класса U;
  (P2) p=2: все 9 клеток — квадраты в Q_2, v_2(rz),v_2(sz)>=3;
  (P3) p=3: все 9 клеток — квадраты в Q_3, v_3(z)>=1;
  (Po) p нечётно: v_p(T) нечётно => p|(r-s);  v_p(L) нечётно => p|(r+s);
  (P34) p=3 mod 4: все 9 оценок чётны;
  (R)  вещественное место: 8 произведений >0  =>  все 9 клеток >0.
Затем (CRT) строим z, для которого 8 произведений — квадраты СРАЗУ во всех p из набора P и в R, и печатаем
локальные классы T, L (показывает, что локальные ограничения не пусты и что исключение 5,13 у 126/451
происходит только глобально, через согласование с Q_2, Q_3).
Код независим от файлов Codex. Все циклы ограничены явно."""
import random, json, time, sys
from math import gcd

random.seed(20260926)
t0 = time.time()

def vp(n, p):
    if n == 0: return 10 ** 9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v

def cls(num, den, p):
    """Класс ненулевого рационального num/den в Q_p^*/Q_p^*2 как битовая маска.
    p нечётно: бит0 = v mod 2, бит1 = [единичная часть — невычет].
    p=2: бит0 = v mod 2, бит1 = [u = 3 mod 4], бит2 = [u = +-3 mod 8]."""
    assert num != 0 and den != 0
    a, b = vp(num, p), vp(den, p)
    un, ud = num // p ** a, den // p ** b
    v = (a - b) & 1
    if p == 2:
        u = (un * ud) % 8          # ud нечётно, ud^-1 = ud mod 8
        return v | (((u - 1) // 2 % 2) << 1) | ((((u * u - 1) // 8) % 2) << 2)
    u = (un * ud) % p
    leg = pow(u, (p - 1) // 2, p)
    return v | ((0 if leg == 1 else 1) << 1)

IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]
LINES = [[(i, j) for j in IDX] for i in IDX] + [[(i, j) for i in IDX] for j in IDX] + \
        [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]]
T_CELLS = [(-1, 0), (1, 1), (0, -1)]
L_CELLS = [(-1, 0), (1, -1), (0, 1)]

def cell_nums(r, s, N, D):
    return {(i, j): D + (i * r + j * s) * N for (i, j) in CELLS}

def prod_num(cn, cl):
    x = 1
    for c in cl: x *= cn[c]
    return x

def eight_square_at(r, s, N, D, p):
    cn = cell_nums(r, s, N, D)
    if any(v == 0 for v in cn.values()): return None
    for ln in LINES:
        if cls(prod_num(cn, ln), D ** 3, p) != 0: return False
    return True

def check_claims(r, s, N, D, p):
    """Возвращает (список нарушений, нетривиальность) для z=N/D, прошедшего 8 условий в Q_p."""
    viol = []
    cn = cell_nums(r, s, N, D)
    c = {k: cls(v, D, p) for k, v in cn.items()}
    # самопроверка гомоморфизма: класс произведения = xor классов
    for ln in LINES:
        x = 0
        for k in ln: x ^= c[k]
        if x != 0: viol.append("xor-line")
    cT = cls(prod_num(cn, T_CELLS), D ** 3, p); cL = cls(prod_num(cn, L_CELLS), D ** 3, p)
    if c[(-1, -1)] != c[(1, 1)] or c[(-1, 1)] != c[(1, -1)] or c[(-1, 0)] != c[(1, 0)] or c[(0, -1)] != c[(0, 1)]:
        viol.append("M:opposite")
    if c[(1, 0)] != c[(0, 1)]: viol.append("M:axial")
    if cT != c[(1, 1)]: viol.append("M:T")
    if cL != c[(1, -1)]: viol.append("M:L")
    vT, vL = cT & 1, cL & 1
    if p == 2:
        if any(v != 0 for v in c.values()): viol.append("P2:cells")
        if vp(r * N, 2) - vp(D, 2) < 3 or vp(s * N, 2) - vp(D, 2) < 3: viol.append("P2:val")
    elif p == 3:
        if any(v != 0 for v in c.values()): viol.append("P3:cells")
        if vp(N, 3) - vp(D, 3) < 1: viol.append("P3:val")
    else:
        if vT and (r - s) % p != 0: viol.append("Po:T")
        if vL and (r + s) % p != 0: viol.append("Po:L")
        if p % 4 == 3 and any(v & 1 for v in c.values()): viol.append("P34")
        # уточнение Claude (не утверждение Codex): p=5 mod 8 и чётная v_p(r-s) => v_p(T) чётна; то же для L, r+s
        if p % 8 == 5 and vT and (r - s) % p == 0 and vp(r - s, p) % 2 == 0: viol.append("REF:T")
        if p % 8 == 5 and vL and (r + s) % p == 0 and vp(r + s, p) % 2 == 0: viol.append("REF:L")
    nontriv = any(v != 0 for v in c.values())
    return viol, nontriv, vT, vL

def primes_upto(n):
    return [q for q in range(2, n + 1) if all(q % d for d in range(2, int(q ** 0.5) + 1))]

def pf(n):
    n = abs(n); out = []; d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: out.append(n)
    return out

SLOPES = [(126, 451), (73, 362), (265, 298), (12, 85), (1, 4), (2, 3), (1, 6), (4, 9), (3, 10), (7, 10),
          (12, 17), (1, 12), (8, 21), (5, 13), (1, 3), (3, 5), (1, 24), (2, 27), (7, 24), (13, 24),
          (5, 12), (5, 8), (-126, 451), (-73, 362), (9, 16), (11, 14),
          (1, 26), (7, 32), (1, 170), (1, 290), (12, 13), (2, 23), (4, 29)]
BASEP = primes_upto(41)
NSAMP = int(sys.argv[1]) if len(sys.argv) > 1 else 120
KRANGE = range(-6, 7)

summary = {"tested": 0, "passed": 0, "passed_nontrivial": 0, "violations": 0, "viol_examples": []}
per_p = {}
odd_support = {}      # (slope,p) -> [T нечётна встречалась, L нечётна встречалась]
pass_store = {}       # (slope,p) -> список прошедших (N,D,vT,vL,nontriv), для CRT
for (r, s) in SLOPES:
    assert gcd(abs(r), s) == 1 and abs(r) != s
    P = sorted(set(BASEP) | set(pf(r)) | set(pf(s)) | set(pf(r - s)) | set(pf(r + s)))
    for p in P:
        key = (r, s, p)
        st = per_p.setdefault(p, [0, 0, 0])
        for k in KRANGE:
            for _ in range(NSAMP):
                while True:   # не более нескольких итераций: вероятность отказа < 1/2
                    a = random.randint(1, 10 ** 6) * random.choice((1, -1)); b = random.randint(1, 10 ** 6)
                    if a % p and b % p and gcd(a, b) == 1: break
                if k >= 0: N, D = a * p ** k, b
                else: N, D = a, b * p ** (-k)
                ok = eight_square_at(r, s, N, D, p)
                summary["tested"] += 1; st[0] += 1
                if not ok: continue
                summary["passed"] += 1; st[1] += 1
                viol, nontriv, vT, vL = check_claims(r, s, N, D, p)
                if nontriv: summary["passed_nontrivial"] += 1; st[2] += 1
                if viol:
                    summary["violations"] += 1
                    if len(summary["viol_examples"]) < 20:
                        summary["viol_examples"].append((r, s, p, N, D, viol))
                os_ = odd_support.setdefault(f"{r}/{s}@{p}", [0, 0])
                os_[0] += vT; os_[1] += vL
                lst = pass_store.setdefault(key, [])
                if len(lst) < 400: lst.append((N, D, vT, vL, nontriv))
    print(f"наклон {r}/{s}: простых {len(P)}; всего тестов {summary['tested']}, прошло {summary['passed']}, "
          f"нетрив {summary['passed_nontrivial']}, нарушений {summary['violations']}  t={time.time()-t0:.0f}s", flush=True)

print("по простым [тестов, прошло, нетрив. прошло]:")
print({p: v for p, v in sorted(per_p.items()) if p < 50})
odd_nonzero = {k: v for k, v in odd_support.items() if v[0] or v[1]}
print("где встретились нечётные v_p(T) / v_p(L) (число случаев):", odd_nonzero)

# ---- вещественное место
real_viol = 0; real_pass = 0; real_tested = 0
for (r, s) in SLOPES:
    for _ in range(4000):
        N = random.randint(-10 ** 6, 10 ** 6); D = random.randint(1, 10 ** 6) * (abs(r) + s) // random.choice((1, 2, 4))
        if N == 0 or D == 0: continue
        cn = cell_nums(r, s, N, D)
        if any(v == 0 for v in cn.values()): continue
        real_tested += 1
        if all(prod_num(cn, ln) > 0 for ln in LINES):
            real_pass += 1
            if any(v < 0 for v in cn.values()): real_viol += 1
print(f"R: тестов {real_tested}, 8 произведений >0: {real_pass}, из них с отрицательной клеткой: {real_viol}")
summary.update({"real_tested": real_tested, "real_pass": real_pass, "real_violations": real_viol})

# ---- CRT: одновременная локальная разрешимость на P ∪ {∞}
def crt_point(r, s, P, want, soft=()):
    """want: p -> фильтр для выбора локального z_p среди прошедших. Возвращает (N,D) или None."""
    choice = {}
    for p in P:
        allc = pass_store.get((r, s, p), [])
        cand = [x for x in allc if want.get(p, lambda x: True)(x)]
        if not cand and p in soft: cand = allc
        if not cand: return None, f"нет локального кандидата при p={p}"
        choice[p] = random.choice(cand)
    Dp = 1
    for p in P:
        N_, D_ = choice[p][0], choice[p][1]
        Dp *= p ** max(0, vp(D_, p) - vp(N_, p))
    prec = {p: vp(Dp, p) + 14 + (6 if p == 2 else 0) for p in P}
    M = 1
    for p in P: M *= p ** prec[p]
    Q = 4 * (abs(r) + s) * M + 1       # взаимно просто со всеми p из P, |z| < 1/(4(|r|+|s|))
    D = Dp * Q
    # N = z_p * D (mod p^prec)
    N = 0; mod = 1
    for p in P:
        zN, zD = choice[p][0], choice[p][1]
        m = p ** prec[p]
        bet = vp(zD, p); zD1 = zD // p ** bet
        dlt = vp(D, p); D1 = D // p ** dlt
        target = (zN * p ** (dlt - bet) * D1 * pow(zD1, -1, m)) % m
        # CRT слияние
        t = ((target - N) * pow(mod, -1, m)) % m
        N += mod * t; mod *= m
    if N > mod // 2: N -= mod
    g = gcd(N, D); N //= g; D //= g
    return (N, D), None

crt_results = []
plans = [
    ((126, 451), "L нечётна при 577", {577: lambda x: x[3] == 1}),
    ((126, 451), "T нечётна при 5 (локально)", {5: lambda x: x[2] == 1}),
    ((126, 451), "T нечётна при 13 (локально)", {13: lambda x: x[2] == 1}),
    ((73, 362), "L нечётна при 5 и 29", {5: lambda x: x[3] == 1, 29: lambda x: x[3] == 1}),
    ((73, 362), "T нечётна при 17", {17: lambda x: x[2] == 1}),
    ((12, 85), "T при 73, L при 97", {73: lambda x: x[2] == 1, 97: lambda x: x[3] == 1}),
    ((265, 298), "общий нетривиальный", {}),
]
for (r, s), label, want in plans:
    P = sorted(set(primes_upto(13)) | set(pf(r)) | set(pf(s)) | set(pf(r - s)) | set(pf(r + s)))
    # остальные простые: предпочитаем нетривиальные локальные точки
    soft = [p for p in P if p not in want]
    for p in soft:
        if p not in (2, 3):
            want[p] = (lambda x: x[4])
    best = None
    for attempt in range(30):
        pt, err = crt_point(r, s, P, want, soft)
        if pt is None: best = err; break
        N, D = pt
        okall = all(eight_square_at(r, s, N, D, p) for p in P)
        cn = cell_nums(r, s, N, D)
        pos = all(v * D > 0 for v in cn.values())
        if not (okall and pos): continue
        Tn, Ln = prod_num(cn, T_CELLS), prod_num(cn, L_CELLS)
        rec = {"slope": f"{r}/{s}", "label": label, "P": P, "z_num_digits": len(str(abs(N))), "z_den_digits": len(str(D)),
               "v_p(T)": {p: vp(Tn, p) - 3 * vp(D, p) for p in P}, "v_p(L)": {p: vp(Ln, p) - 3 * vp(D, p) for p in P},
               "T_square_Q2": cls(Tn, D ** 3, 2) == 0, "T_square_Q3": cls(Tn, D ** 3, 3) == 0,
               "L_square_Q2": cls(Ln, D ** 3, 2) == 0, "L_square_Q3": cls(Ln, D ** 3, 3) == 0,
               "T_pos": Tn * D > 0, "L_pos": Ln * D > 0}
        viol = []
        for p in P:
            v1, _, vT, vL = check_claims(r, s, N, D, p)
            viol += [f"{p}:{x}" for x in v1]
        rec["violations"] = viol
        best = rec; break
    crt_results.append(best)
    print("CRT", label, "->", best if isinstance(best, str) else
          {k: best[k] for k in ("slope", "z_num_digits", "z_den_digits", "T_square_Q2", "T_square_Q3",
                                "L_square_Q2", "L_square_Q3", "T_pos", "L_pos", "violations")})
    if isinstance(best, dict):
        oddT = [p for p, v in best["v_p(T)"].items() if v % 2]; oddL = [p for p, v in best["v_p(L)"].items() if v % 2]
        print("     нечётные v_p(T) при p из P:", oddT, "; нечётные v_p(L):", oddL)
        best["odd_T"] = oddT; best["odd_L"] = oddL

summary["crt"] = crt_results
summary["per_p"] = per_p
summary["odd_support"] = odd_nonzero
summary["time_s"] = round(time.time() - t0, 1)
print("ИТОГ:", {k: summary[k] for k in ("tested", "passed", "passed_nontrivial", "violations", "real_violations", "time_s")})
json.dump(summary, open("r3_local_control.json", "w"), indent=1, default=str)
