# -*- coding: utf-8 -*-
"""s11: универсальные простые p = 2, 3 (и 5, 7 при единичных клетках).

Факт (NOTE.md §7, доказано элементарно): у полного квадрата с центром a каждая линия через центр --
тройка квадратов (a-t, a, a+t), откуда t/a in 8Z_2 и t/a in 3Z_3 (коника u^2+w^2=2). Значит b/a, c/a
лежат в 8Z_2 и 3Z_3 при ЛЮБОМ наклоне: все девять клеток сравнимы с a (mod 8a) и (mod 3a).

Часть 1 (одна линия, p=2): класс delta_2(Q) канонической половинки (корни = 1 mod 4) по модулю delta_2(E_t[2])
  как функция s = t/x in 8Z_2. Вопрос: когда P = 2Q лежит в 4E_t(Q_2).
Часть 2 (полный квадрат, a = 1, p in {2,3,5,7}): все 2^8 знаковых узоров корней; классы 16 дельта-координат
  (x и x-t на 8 линиях); размерности, совместные соотношения; проверка правила
  delta_p(Q_l) = kappa_p(l) * delta_p(T(метки знаков)) и «транспонированного соотношения»
  (T^b_j - T^b_0)_{kappa(i)} = (T^c_i - T^c_0)_{kappa(j)},  i, j = +-1.
Запуск: python3 s11_universal_primes.py   (чистый Python)
"""
import random, itertools, json, sys
sys.path.insert(0, '.')
from grid import vp_int, legendre

random.seed(1111)
N = 240  # точность p-адики (по модулю p^N)
IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]
LINES = []  # (имя, кривая, (b,c)-коэф. t, (P1,P2,P3))
for j in IDX:
    LINES.append(("Eb_row%+d" % j, "b", (1, 0), ((-1, j), (0, j), (1, j))))
for i in IDX:
    LINES.append(("Ec_col%+d" % i, "c", (0, 1), ((i, -1), (i, 0), (i, 1))))
LINES.append(("Ebpc_diag", "b+c", (1, 1), ((-1, -1), (0, 0), (1, 1))))
LINES.append(("Ebmc_anti", "b-c", (1, -1), ((-1, 1), (0, 0), (1, -1))))


# ------------------------------------------------------------------ p-адика
def sqrt_padic(u, p, n=N):
    """корень из u (целое, u = 1 mod 8 при p=2, u = квадрат-единица при нечётном p), канонический:
    = 1 mod 4 при p = 2; при нечётном p -- тот, что = +sqrt(u mod p) с минимальным представителем."""
    M = p ** n
    if p == 2:
        assert u % 8 == 1
        x = 1
        for k in range(3, n):  # x^2 = u mod 2^k -> mod 2^(k+1)
            if (x * x - u) % (2 ** (k + 1)):
                x += 2 ** (k - 1)
        x %= M
        if x % 4 != 1:
            x = (-x) % M
        assert (x * x - u) % M == 0
        return x
    u0 = u % p
    assert legendre(u0, p) == 1
    x = next(y for y in range(1, p) if (y * y - u0) % p == 0)
    k = 1
    while k < n:
        k = min(2 * k, n)
        m = p ** k
        x = (x - (x * x - u) * pow(2 * x, -1, m)) % m
    return x % M


def pclass(z, p, n=N):
    """класс z (целое, представляющее p-адическое число mod p^n) в Q_p^*/Q_p^*2."""
    z %= p ** n
    assert z != 0
    v = vp_int(z, p)
    assert v < n - 5, "точность"
    u = z // p ** v
    if p == 2:
        return (v % 2, u % 8)
    return (v % 2, legendre(u, p))


def rclass(q_num, q_den, p):
    """класс рационального числа q_num/q_den в Q_p^*/Q_p^*2 (для констант)."""
    vn, vd = vp_int(abs(q_num), p), vp_int(abs(q_den), p)
    un, ud = q_num // p ** vn, q_den // p ** vd
    if p == 2:
        return ((vn - vd) % 2, (un * ud) % 8)
    return ((vn - vd) % 2, legendre(un * ud, p))


def cmul(x, y, p):
    if p == 2:
        return ((x[0] + y[0]) % 2, (x[1] * y[1]) % 8)
    return ((x[0] + y[0]) % 2, x[1] * y[1])


def class_bits(cl, p):
    if p == 2:
        return [cl[0], 1 if cl[1] in (3, 7) else 0, 1 if cl[1] in (5, 7) else 0]  # 5 и -1 -- базис
    return [cl[0], 1 if cl[1] == -1 else 0]


def torsion_classes(t_num, t_den, p):
    """delta_p(T) в координатах (x, x-t) для T in {O, T+, T-, T0}; метки (s1,s2): O=(0,0),T+=(1,0),T-=(0,1),T0=(1,1)."""
    cl = lambda n, d=1: rclass(n * t_den if False else n, d, p)
    one = rclass(1, 1, p)
    t = (t_num, t_den)
    c_t = rclass(t_num, t_den, p)
    c_mt = rclass(-t_num, t_den, p)
    c_m1 = rclass(-1, 1, p)
    c_2 = rclass(2, 1, p)
    c_m2t = rclass(-2 * t_num, t_den, p)
    return {(0, 0): (one, one), (1, 0): (c_t, c_2), (0, 1): (c_mt, c_m2t), (1, 1): (c_m1, c_mt)}


def rank_f2(rows):
    basis = []
    for x in rows:
        for bv in basis:
            x = min(x, x ^ bv)
        if x:
            basis.append(x)
    return len(basis)


def affine_dim(vecs):
    base = vecs[0]
    rows = []
    for v in vecs[1:]:
        m = 0
        for k, (x, y) in enumerate(zip(v, base)):
            if x != y:
                m |= 1 << k
        rows.append(m)
    return rank_f2(rows)


def rel_basis(vecs):
    """базис аффинных соотношений sum c_k v_k = const на всех векторах (над F_2)."""
    n = len(vecs[0])
    base = vecs[0]
    M = [[x ^ y for x, y in zip(v, base)] for v in vecs[1:]]
    piv, r = [], 0
    for col in range(n):
        pr = next((i for i in range(r, len(M)) if M[i][col]), None)
        if pr is None:
            continue
        M[r], M[pr] = M[pr], M[r]
        for i in range(len(M)):
            if i != r and M[i][col]:
                M[i] = [x ^ y for x, y in zip(M[i], M[r])]
        piv.append(col)
        r += 1
    free = [c for c in range(n) if c not in piv]
    out = []
    for f in free:
        cvec = [0] * n
        cvec[f] = 1
        for i, pc in enumerate(piv):
            if M[i][f]:
                cvec[pc] = 1
        const = sum(ci & bi for ci, bi in zip(cvec, base)) % 2
        out.append((cvec, const))
    return out


# ------------------------------------------------------------------ часть 1: одна линия при p = 2
def part1():
    res = {}
    for k in range(3, 12):
        for unit_mod in (1, 3, 5, 7):
            cnt, in4 = 0, 0
            for _ in range(40):
                # x = 1 (масштаб), t = 2^k * u, u = unit_mod mod 8
                u = random.randrange(1, 2 ** 20) * 8 + unit_mod
                t = 2 ** k * u
                al, ga, be = sqrt_padic(1 - t, 2), 1, sqrt_padic(1 + t, 2)
                A, B, C = al + ga, ga + be, al + be
                cq = (pclass(A * B, 2), pclass(A * C, 2))
                TC = torsion_classes(t, 1, 2)
                in4 += any(cq == tc for tc in TC.values())
                cnt += 1
            res["v2(t/x)=%d, unit=%d mod 8" % (k, unit_mod)] = {"samples": cnt, "P_in_4E(Q_2)": in4}
    return res


# ------------------------------------------------------------------ часть 2: полный квадрат
def labels_of(eps, line):
    P1, P2, P3 = line
    return ((eps[P1] != eps[P2]) * 1, (eps[P3] != eps[P2]) * 1)


def square_classes(p, b, c, eps, roots):
    out = {}
    for name, cv, coef, (P1, P2, P3) in LINES:
        r1, r2, r3 = (eps[P] * roots[P] for P in (P1, P2, P3))
        A, B, C = r1 + r2, r2 + r3, r1 + r3
        out[name] = (pclass(A * B, p), pclass(A * C, p))
    return out


def gen_bc(p):
    """случайные b, c с a = 1: b, c in 8Z_2 (p=2) или pZ_p (нечётное p); все клетки -- единицы."""
    lo = 3 if p == 2 else 1
    while True:
        vb, vc = random.randint(lo, lo + 3), random.randint(lo, lo + 3)
        ub = random.randrange(1, p ** 6)
        uc = random.randrange(1, p ** 6)
        if ub % p == 0 or uc % p == 0:
            continue
        b, c = random.choice((1, -1)) * p ** vb * ub, random.choice((1, -1)) * p ** vc * uc
        forms = [b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c]
        if 0 in forms:
            continue
        return b, c


ALL_EPS = []
for bits in itertools.product((1, -1), repeat=8):
    e = {(0, 0): 1}
    for P, s in zip([P for P in CELLS if P != (0, 0)], bits):
        e[P] = s
    ALL_EPS.append(e)


def part2(p, nbc=12):
    rec = {"p": p, "bc_samples": [], "rule_checks": 0, "rule_fail": 0, "transpose_checks": 0, "transpose_fail": 0,
           "injective_delta_T_all": True}
    allvecs = []
    kappa_table = {}
    for _ in range(nbc):
        b, c = gen_bc(p)
        roots = {P: sqrt_padic(1 + P[0] * b + P[1] * c, p) for P in CELLS}
        vecs = []
        for eps in ALL_EPS:
            cl = square_classes(p, b, c, eps, roots)
            # правило: delta_p(Q_l) = kappa(l) * delta_p(T(метки)), kappa(l) := класс при всех знаках +
            T_of = {}
            for name, cv, coef, L in LINES:
                t = coef[0] * b + coef[1] * c
                TC = torsion_classes(t, 1, p)
                if len(set(TC.values())) < 4:
                    rec["injective_delta_T_all"] = False
                lab = labels_of(eps, L)
                T_of[name] = lab
                if name not in kappa_table.get((b, c), {}):
                    kappa_table.setdefault((b, c), {})
                    base = square_classes(p, b, c, {P: 1 for P in CELLS}, roots)[name]
                    kappa_table[(b, c)][name] = base
                kap = kappa_table[(b, c)][name]
                pred = tuple(cmul(kap[k], TC[lab][k], p) for k in range(2))
                rec["rule_checks"] += 1
                rec["rule_fail"] += (pred != cl[name])
            # транспонированное соотношение на метках (i,j in {-1,+1}; kappa(-1) -> s1, kappa(+1) -> s2)
            kap_idx = {-1: 0, 1: 1}
            for i in (-1, 1):
                for j in (-1, 1):
                    lhs = (T_of["Eb_row%+d" % j][kap_idx[i]] + T_of["Eb_row+0"][kap_idx[i]]) % 2
                    rhs = (T_of["Ec_col%+d" % i][kap_idx[j]] + T_of["Ec_col+0"][kap_idx[j]]) % 2
                    rec["transpose_checks"] += 1
                    rec["transpose_fail"] += (lhs != rhs)
            v = []
            for name, cv, coef, L in LINES:
                for k in range(2):
                    v += class_bits(cl[name][k], p)
            vecs.append(v)
        allvecs += vecs
        # размерности для данного (b, c)
        nb = len(class_bits(((0, 1) if p != 2 else (0, 1)), p))
        per_line = 2 * nb
        def proj(vs, curves):
            idx = []
            for k, (name, cv, coef, L) in enumerate(LINES):
                if cv in curves:
                    idx += list(range(per_line * k, per_line * (k + 1)))
            return [[w[i] for i in idx] for w in vs]
        d_all = affine_dim(vecs)
        per = {cv: affine_dim(proj(vecs, {cv})) for cv in ("b", "c", "b+c", "b-c")}
        dbc = affine_dim(proj(vecs, {"b", "c"}))
        kap_nontriv = {name: (kappa_table[(b, c)][name] != (((0, 1), (0, 1)) if p != 2 else ((0, 1), (0, 1))))
                       for name, cv, coef, L in LINES}
        kap_in_T = {}
        for name, cv, coef, L in LINES:
            t = coef[0] * b + coef[1] * c
            TC = torsion_classes(t, 1, p)
            kap_in_T[name] = kappa_table[(b, c)][name] in TC.values()
        rec["bc_samples"].append({"b": b, "c": c, "v_p": [vp_int(abs(x), p) for x in (b, c, b + c, b - c)],
                                  "dim_all": d_all, "dim_per_curve": per, "dim_Eb_Ec": dbc,
                                  "joint_all": sum(per.values()) - d_all, "joint_Eb_Ec": per["b"] + per["c"] - dbc,
                                  "kappa_in_delta(E[2])": kap_in_T})
    rec["dim_union_all_bc"] = affine_dim(allvecs)
    return rec


if __name__ == "__main__":
    out = {}
    p1 = part1()
    out["part1_single_line_p2"] = p1
    print("Часть 1 (p=2, одна линия, x=1): P in 4E(Q_2)?")
    for k, v in p1.items():
        print("  ", k, v)
    out["part2"] = {}
    for p in (3, 5, 7, 2):
        r = part2(p)
        out["part2"][p] = r
        print("p=%d: правило kappa*delta(T(метки)): %d/%d нарушений; транспонирование: %d/%d нарушений; "
              "delta_p инъективно на E[2] у всех кривых: %s; dim объединения по всем (b,c): %d"
              % (p, r["rule_fail"], r["rule_checks"], r["transpose_fail"], r["transpose_checks"],
                 r["injective_delta_T_all"], r["dim_union_all_bc"]))
        for s in r["bc_samples"]:
            print("   v_p(b,c,b+c,b-c)=%s dim=%d per=%s dimEbEc=%d совм.(все)=%d совм.(Eb-Ec)=%d kappa in delta(E[2]): %s"
                  % (s["v_p"], s["dim_all"], s["dim_per_curve"], s["dim_Eb_Ec"], s["joint_all"], s["joint_Eb_Ec"],
                     "".join("1" if x else "0" for x in s["kappa_in_delta(E[2])"].values())))
    json.dump(out, open("s11_universal_primes.json", "w"), ensure_ascii=False, indent=1, default=str)
