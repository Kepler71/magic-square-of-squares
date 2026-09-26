# -*- coding: utf-8 -*-
"""s3: проверка явных локальных совместных соотношений (теорема B NOTE.md) на Q_p-точках полного квадрата (lambda = 0).

Типы простых (p нечётно, p не делит клетки):
  "b"   : p | b, p ∤ c, v_p(b) = e нечётно           -> (J_p): оценки E_b-классов задают вычеты E_c- и диагональных классов
  "c"   : симметрично (проверяется транспонированием сетки)
  "bc"  : p | b, p | c (общий простой, как 3)        -> всё задаётся знаковым узором rho(p), плакетные соотношения
  "b+c" : p | b+c, p ∤ bc                             -> связь вычетов средних половинок на E_b и E_c
Отрицательный контроль: строки E_b от одной точки (a,b,c), столбцы E_c от другой (a',b,c) с теми же b,c:
соотношения (J_p) должны нарушаться с положительной частотой.
Контроль теоремы C: a = t^2, v_p(t) = -K  => все 16 delta-координат -- квадраты в Q_p.
"""
import random, json, sys, itertools
sys.path.insert(0, '.')
from grid import IDX, lines_of, sums, padic_sqrt, local_class, legendre, vp_int

random.seed(925)
NPREC = 80


def chi(x, p):
    return legendre(x % p, p)


def make_roots(a, b, c, p, N=NPREC):
    r = {}
    for i in IDX:
        for j in IDX:
            cell = a + i * b + j * c
            s, ok = padic_sqrt(cell, p, N)
            if not ok:
                return None
            r[(i, j)] = s if random.random() < 0.5 else (-s) % p ** N
    return r


def classes(r, p, N=NPREC):
    """Для каждой линии: локальные классы трёх delta-координат и трёх сумм."""
    out = {}
    for name, tl, (al, ga, be) in lines_of(r):
        A, B, C = sums(al, ga, be)
        out[name] = {"coords": [local_class(A * B, p, N), local_class(A * C, p, N), local_class(B * C, p, N)],
                     "sums": [local_class(A, p, N), local_class(B, p, N), local_class(C, p, N)]}
    return out


def sgn_bit(r1, r2, p):
    """1, если r1 = -r2 mod p (обе -- единицы), 0 если r1 = r2 mod p; иначе None."""
    if (r1 - r2) % p == 0:
        return 0
    if (r1 + r2) % p == 0:
        return 1
    return None


def pw(ch, e):
    return ch if e % 2 else 1


# ------------------------------------------------------------------ тип "b"
def check_type_b(a, b, c, r, p):
    """Проверка формул (J_p) и формул для диагоналей. Возвращает список нарушений (пусто = всё верно)."""
    bad = []
    cl = classes(r, p)
    sig = {(i, j): sgn_bit(r[(i, j)], r[(0, j)], p) for i in (-1, 1) for j in IDX}
    assert all(v is not None for v in sig.values())
    # 1) E_b видит sigma: чётности оценок (x, x-b, x+b) строки j = (s_-1 + s_1, s_1, s_-1)
    for j in IDX:
        vpar = [cl["Eb_row%+d" % j]["coords"][k][0] for k in range(3)]
        pred = [(sig[(-1, j)] + sig[(1, j)]) % 2, sig[(1, j)], sig[(-1, j)]]
        if vpar != pred:
            bad.append(("Eb_val", j, vpar, pred))
    # восстановление sigma ТОЛЬКО по E_b-оценкам (как в формулировке теоремы)
    sig_from_Eb = {}
    for j in IDX:
        co = cl["Eb_row%+d" % j]["coords"]
        sig_from_Eb[(1, j)] = co[1][0]
        sig_from_Eb[(-1, j)] = co[2][0]
    cc, mc, m2c, m1 = chi(c, p), chi(-c, p), chi(-2 * c, p), chi(-1, p)
    # 2) (J_p): столбцы i = +-1 относительно столбца 0
    ref = cl["Ec_col+0"]["coords"]
    for i in (-1, 1):
        s1 = (sig_from_Eb[(i, -1)] + sig_from_Eb[(i, 0)]) % 2
        s2 = (sig_from_Eb[(i, 0)] + sig_from_Eb[(i, 1)]) % 2
        pred = [pw(cc, s1) * pw(mc, s2), pw(mc, s1) * pw(m2c, s1 + s2), pw(m1, s1) * pw(mc, s2) * pw(m2c, s1 + s2)]
        got = [cl["Ec_col%+d" % i]["coords"][k][1] * ref[k][1] for k in range(3)]
        vals = [cl["Ec_col%+d" % i]["coords"][k][0] for k in range(3)]
        if got != pred or any(vals):
            bad.append(("Jp", i, got, pred, vals))
    # 3) диагональ E_{b+c} против столбца 0 (координаты в том же порядке)
    s, s_ = sig_from_Eb[(-1, -1)], sig_from_Eb[(1, 1)]
    pred = [pw(cc, s) * pw(mc, s_), pw(mc, s) * pw(m2c, s + s_), pw(m1, s) * pw(mc, s_) * pw(m2c, s + s_)]
    got = [cl["Ebpc_diag"]["coords"][k][1] * ref[k][1] for k in range(3)]
    if got != pred:
        bad.append(("diag+", got, pred))
    # 4) антидиагональ E_{b-c}: координаты (x, x-t, x+t) против (x, x+c, x-c) столбца 0
    s, s_ = sig_from_Eb[(-1, 1)], sig_from_Eb[(1, -1)]
    pred = [pw(mc, s) * pw(cc, s_), pw(mc, s) * pw(m1, s_) * pw(m2c, s + s_), pw(mc, s_) * pw(m2c, s + s_)]
    got = [cl["Ebmc_anti"]["coords"][0][1] * ref[0][1],
           cl["Ebmc_anti"]["coords"][1][1] * ref[2][1],
           cl["Ebmc_anti"]["coords"][2][1] * ref[1][1]]
    if got != pred:
        bad.append(("diag-", got, pred))
    return bad, cl


def transpose(r):
    return {(i, j): r[(j, i)] for i in IDX for j in IDX}


# ------------------------------------------------------------------ тип "bc" (общий простой)
def predict_common(a, b, c, r, p):
    """Все локальные классы delta-координат из узора rho и единичных частей разностей клеток."""
    rho = {k: sgn_bit(r[k], r[(0, 0)], p) for k in r}
    assert all(v is not None for v in rho.values())
    cellval = {(i, j): (i, j) for i in IDX for j in IDX}
    pred = {}
    # суммы: (класс-оценка, символ) с общим множителем chi(r00), который в произведениях сокращается
    r00c = chi(r[(0, 0)], p)
    for name, tl, trip in lines_of(r):
        pass
    # перечислим линии с индексами клеток
    L = []
    for j in IDX:
        L.append(("Eb_row%+d" % j, [(-1, j), (0, j), (1, j)]))
    for i in IDX:
        L.append(("Ec_col%+d" % i, [(i, -1), (i, 0), (i, 1)]))
    L.append(("Ebpc_diag", [(-1, -1), (0, 0), (1, 1)]))
    L.append(("Ebmc_anti", [(-1, 1), (0, 0), (1, -1)]))

    def cell(k):
        return a + k[0] * b + k[1] * c

    def sum_class(k1, k2):
        # r1 + r2, r1 = (-1)^rho1 r00 mod p, r2 = (-1)^rho2 r00 mod p
        if rho[k1] == rho[k2]:
            return (0, chi(2, p) * (chi(-1, p) if rho[k1] else 1) * r00c)
        D = cell(k1) - cell(k2)  # r1^2 - r2^2 = D, r1 + r2 = D/(r1 - r2), r1 - r2 = 2 r1 mod p
        v = vp_int(abs(D), p)
        u = D // p ** v
        return (v % 2, chi(u, p) * chi(2, p) * (chi(-1, p) if rho[k1] else 1) * r00c)

    for name, (ka, kg, kb) in L:
        A_, B_, C_ = sum_class(ka, kg), sum_class(kg, kb), sum_class(ka, kb)
        mul = lambda X, Y: ((X[0] + Y[0]) % 2, X[1] * Y[1])
        pred[name] = [mul(A_, B_), mul(A_, C_), mul(B_, C_)]
    return pred, rho


# ------------------------------------------------------------------ тип "b+c"
def check_type_bpc(a, b, c, r, p):
    cl = classes(r, p)
    e3 = sgn_bit(r[(1, 0)], r[(0, -1)], p)
    e4 = sgn_bit(r[(-1, 0)], r[(0, 1)], p)
    assert e3 is not None and e4 is not None
    cc, mc, m1, c2 = chi(c, p), chi(-c, p), chi(-1, p), chi(2 * c, p)
    pred = [pw(mc, e4) * pw(cc, e3), pw(cc, e4) * pw(c2, e3 + e4), pw(cc, e3) * pw(m1, e4) * pw(c2, e3 + e4)]
    k0, kk0 = cl["Eb_row+0"]["coords"], cl["Ec_col+0"]["coords"]
    got = [k0[0][1] * kk0[0][1], k0[1][1] * kk0[2][1], k0[2][1] * kk0[1][1]]
    # множество допустимых значений без знания e3,e4
    allowed = set()
    for x3, x4 in itertools.product((0, 1), repeat=2):
        allowed.add((pw(mc, x4) * pw(cc, x3), pw(cc, x4) * pw(c2, x3 + x4), pw(cc, x3) * pw(m1, x4) * pw(c2, x3 + x4)))
    forced = (cc == 1 and m1 == 1 and chi(2, p) == 1)
    return got == pred, tuple(got) in allowed, forced, len(allowed), cl


def tuple_of(cl):
    return tuple(tuple(cl[nm]["coords"][k]) for nm in sorted(cl) for k in (0, 1))


# ------------------------------------------------------------------ основной прогон
out = {"type_b": {}, "type_c": {}, "type_bc": {}, "type_bpc": {}, "negative_control": {}, "theorem_C": {}}
PR = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 73, 89, 97, 101, 113]

for p in PR:
    # ---- тип b
    nb, nbad, images = 0, 0, {}
    while nb < 150:
        e = random.choice((1, 3))
        ub = random.randrange(1, p)
        b = p ** e * ub * random.choice((1, -1))
        c = random.randrange(1, 10 * p)
        if c % p == 0:
            continue
        r00 = random.randrange(1, 10 * p)
        if r00 % p == 0:
            continue
        a = r00 * r00
        if any((a + j * c) % p == 0 for j in IDX):
            continue
        r = make_roots(a, b, c, p)
        if r is None:
            continue
        bad, cl = check_type_b(a, b, c, r, p)
        nb += 1
        nbad += bool(bad)
        if nbad and bad and nbad < 3:
            print("TYPE b violation", p, a, b, c, bad)
    out["type_b"][p] = {"samples": nb, "violations": nbad}
    # ---- тип c (транспонирование)
    nc, ncbad = 0, 0
    while nc < 100:
        e = random.choice((1, 3))
        c = p ** e * random.randrange(1, p)
        b = random.randrange(1, 10 * p)
        if b % p == 0:
            continue
        r00 = random.randrange(1, 10 * p)
        if r00 % p == 0:
            continue
        a = r00 * r00
        if any((a + i * b) % p == 0 for i in IDX):
            continue
        r = make_roots(a, b, c, p)
        if r is None:
            continue
        bad, _ = check_type_b(a, c, b, transpose(r), p)
        nc += 1
        ncbad += bool(bad)
    out["type_c"][p] = {"samples": nc, "violations": ncbad}
    # ---- тип bc
    nbc, nbcbad, plaq_bad = 0, 0, 0
    while nbc < 100:
        eb, ec = random.choice((1, 2, 3)), random.choice((1, 2, 3))
        b = p ** eb * random.randrange(1, p) * random.choice((1, -1))
        c = p ** ec * random.randrange(1, p) * random.choice((1, -1))
        r00 = random.randrange(1, 10 * p)
        if r00 % p == 0 or (b + c) == 0 or (b - c) == 0:
            continue
        a = r00 * r00
        r = make_roots(a, b, c, p)
        if r is None:
            continue
        pred, rho = predict_common(a, b, c, r, p)
        cl = classes(r, p)
        nbc += 1
        ok = all(pred[nm][k] == tuple(cl[nm]["coords"][k]) for nm in pred for k in range(3))
        nbcbad += (not ok)
        # плакетные соотношения (только чётности оценок, при нечётных e_b, e_c)
        if eb % 2 and ec % 2:
            for i in (-1, 1):
                for j in (-1, 1):
                    lhs = (cl["Eb_row%+d" % j]["coords"][1 if i == 1 else 2][0] + cl["Eb_row+0"]["coords"][1 if i == 1 else 2][0]) % 2
                    rhs = (cl["Ec_col%+d" % i]["coords"][1 if j == 1 else 2][0] + cl["Ec_col+0"]["coords"][1 if j == 1 else 2][0]) % 2
                    plaq_bad += (lhs != rhs)
    out["type_bc"][p] = {"samples": nbc, "prediction_failures": nbcbad, "plaquette_failures": plaq_bad}
    # ---- тип b+c
    n4, n4exact, n4allowed, n4forced, sizes = 0, 0, 0, 0, set()
    tries = 0
    while n4 < 60 and tries < 200000:
        tries += 1
        b = random.randrange(1, 10 * p)
        if b % p == 0:
            continue
        e = random.choice((1, 3))
        c = -b + p ** e * random.randrange(1, p)
        if c % p == 0 or c == 0:
            continue
        r00 = random.randrange(1, 10 * p)
        if r00 % p == 0:
            continue
        a = r00 * r00
        cells = [a + i * b + j * c for i in IDX for j in IDX]
        if any(x % p == 0 for x in cells):
            continue
        r = make_roots(a, b, c, p)
        if r is None:
            continue
        ex, al, forced, nall, cl = check_type_bpc(a, b, c, r, p)
        n4 += 1
        n4exact += ex
        n4allowed += al
        n4forced += forced
        sizes.add(nall)
    out["type_bpc"][p] = {"samples": n4, "exact_formula_ok": n4exact, "in_allowed_set": n4allowed,
                          "forced_cases(p=1 mod 8, (c/p)=1)": n4forced, "allowed_set_sizes": sorted(sizes)}
    # ---- отрицательный контроль: строки от (a,b,c), столбцы от (a',b,c)
    nn, nviol = 0, 0
    while nn < 100:
        b = p * random.randrange(1, p)
        c = random.randrange(1, 10 * p)
        if c % p == 0:
            continue
        pts = []
        for _ in range(2):
            for _t in range(1000):
                r00 = random.randrange(1, 10 * p)
                if r00 % p == 0:
                    continue
                a = r00 * r00
                if any((a + j * c) % p == 0 for j in IDX):
                    continue
                r = make_roots(a, b, c, p)
                if r is not None:
                    pts.append(r)
                    break
        if len(pts) < 2:
            continue
        r1, r2 = pts
        mixed = dict(r1)
        # столбцы E_c берём от второй точки: подменяем r(i,j) для i = +-1 только в «столбцовой» роли нельзя --
        # поэтому считаем соотношение напрямую: sigma из строк r1, символы E_c из столбцов r2.
        cl1, cl2 = classes(r1, p), classes(r2, p)
        sig = {}
        for j in IDX:
            co = cl1["Eb_row%+d" % j]["coords"]
            sig[(1, j)], sig[(-1, j)] = co[1][0], co[2][0]
        cc, mc, m2c, m1 = chi(c, p), chi(-c, p), chi(-2 * c, p), chi(-1, p)
        ref = cl2["Ec_col+0"]["coords"]
        viol = False
        for i in (-1, 1):
            s1 = (sig[(i, -1)] + sig[(i, 0)]) % 2
            s2 = (sig[(i, 0)] + sig[(i, 1)]) % 2
            pred = [pw(cc, s1) * pw(mc, s2), pw(mc, s1) * pw(m2c, s1 + s2), pw(m1, s1) * pw(mc, s2) * pw(m2c, s1 + s2)]
            got = [cl2["Ec_col%+d" % i]["coords"][k][1] * ref[k][1] for k in range(3)]
            viol |= (got != pred)
        nn += 1
        nviol += viol
    out["negative_control"][p] = {"samples": nn, "violations": nviol}
    # ---- теорема C: a = t^2, v_p(t) = -K; проверяем в целых: умножим всё на p^(2K): a = 1, b,c -> b p^(2K), c p^(2K)
    nC, nCok = 0, 0
    for _ in range(30):
        b = random.randrange(1, 10 ** 6) * random.choice((1, -1))
        c = random.randrange(1, 10 ** 6) * random.choice((1, -1))
        if b in (0, c, -c, 2 * c, -2 * c) or c in (0, 2 * b, -2 * b):
            continue
        K = 3
        s = p ** (2 * K)
        # a = t^2, v_p(t) = -K; умножив всё на t^{-2} = p^{2K} (квадрат), получаем a = 1, b -> b p^{2K}, c -> c p^{2K}
        a_, b_, c_ = 1, b * s, c * s
        r = make_roots(a_, b_, c_, p)
        nC += 1
        if r is None:
            continue
        rr = {k: (v if (v - 1) % p == 0 else (-v) % p ** NPREC) for k, v in r.items()}  # корни ~ +1
        cl = classes(rr, p)
        allsq = all(cl[nm]["coords"][k] == (0, 1) for nm in cl for k in range(3))
        nCok += allsq
    out["theorem_C"][p] = {"samples": nC, "all_16_delta_coords_squares": nCok}
    print(p, {k: out[k][p] for k in out}, flush=True)

json.dump(out, open("s3_padic_local.json", "w"), indent=1)
