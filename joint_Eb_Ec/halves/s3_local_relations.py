# -*- coding: utf-8 -*-
"""s3: локальные совместные соотношения между образами половинок (теорема B NOTE.md) на Q_p-точках
полного квадрата (lambda = 0), p нечётно, все девять клеток -- p-адические единицы.

Модель (доказательство в NOTE.md, §3):
 (L1) линия l на E_t с p | t ("плохая"):  delta_p(Q_l) = delta_p(T(f)),  f = внутренние знаки
      (alpha = -gamma?, 0, beta = -gamma?) по модулю p;
 (L2) линии l, l' с p ∤ t, t' и клетками l' == клеткам l по модулю p (в том же или обратном порядке):
      delta_p(Q_l'') = delta_p(Q_l) * delta_p(T(f)),  l'' = l' (или l', развёрнутая, t'' = -t'),
      f = относительные знаки корней l'' и l по модулю p.
 T(f): (s1, s2) = (f1+f2, f2+f3) mod 2;  (0,0)->O, (1,0)->T+=(t,0), (0,1)->T-=(-t,0), (1,1)->T0=(0,0).
 delta(T0) = (-1,-t,t), delta(T+) = (t,2,2t), delta(T-) = (-t,-2t,2)  (координаты x, x-t, x+t).

Проверки:
 A. все предсказания (L1),(L2) на случайных Q_p-точках всех типов простых (знаки -- из настоящих корней);
 B. отрицательный контроль: строки (E_b) от одной точки, прочие линии от другой с теми же b, c;
    смешанный вектор классов сверяется с моделью при ВСЕХ 256 знаковых узорах rho -- «совместим ли хоть с одним»;
 C. (перенесено в s7, s7b, s10) полнота: при малых p эмпирика ограничена числом допустимых a mod p;
    модель L1+L2 НЕПОЛНА при p | b+-c (одно лишнее соотношение с E_{b-+c}), полная модель -- s10_pair_model.py.
Запуск: python3 s3_local_relations.py   (чистый Python)
"""
import random, json, sys, itertools, time
from fractions import Fraction
sys.path.insert(0, '.')
from grid import IDX, padic_sqrt, local_class, legendre, vp_int, rat_local_class

random.seed(20260926)
NPREC = 60
CELLS = [(i, j) for i in IDX for j in IDX]
LINES = []  # (имя, коэф. t при (b,c), (P1,P2,P3))
for j in IDX:
    LINES.append(("Eb_row%+d" % j, (1, 0), ((-1, j), (0, j), (1, j))))
for i in IDX:
    LINES.append(("Ec_col%+d" % i, (0, 1), ((i, -1), (i, 0), (i, 1))))
LINES.append(("Ebpc_diag", (1, 1), ((-1, -1), (0, 0), (1, 1))))
LINES.append(("Ebmc_anti", (1, -1), ((-1, 1), (0, 0), (1, -1))))
LNAMES = [l[0] for l in LINES]
UNRAM = [((0, 1), (0, s2), (0, s2)) for s2 in (1, -1)] + [((0, -1), (0, s2), (0, -s2)) for s2 in (1, -1)]


def cellval(a, b, c, P):
    return a + P[0] * b + P[1] * c


def tval(b, c, coef):
    return coef[0] * b + coef[1] * c


def mul(x, y):
    return ((x[0] + y[0]) % 2, x[1] * y[1])


def mul3(X, Y):
    return tuple(mul(X[k], Y[k]) for k in range(3))


_TC = {}


def torsion_classes(t, p):
    """Локальные классы delta_p(T), T in {O,T0,T+,T-} для E_t, t -- целое или рациональное."""
    key = (t, p)
    if key not in _TC:
        tt = Fraction(t)
        cl = lambda v: rat_local_class(Fraction(v), p)
        _TC[key] = {"O": ((0, 1), (0, 1), (0, 1)),
                    "T0": (cl(-1), cl(-tt), cl(tt)),
                    "T+": (cl(tt), cl(2), cl(2 * tt)),
                    "T-": (cl(-tt), cl(-2 * tt), cl(2))}
    return _TC[key]


def T_of_flip(f):
    s1, s2 = (f[0] + f[1]) % 2, (f[1] + f[2]) % 2
    return {(0, 0): "O", (1, 0): "T+", (0, 1): "T-", (1, 1): "T0"}[(s1, s2)]


def congruent(P, Q, b, c, p):
    return ((P[0] - Q[0]) * b + (P[1] - Q[1]) * c) % p == 0


def pred_list(b, c, p, sbit):
    """Предсказания модели. sbit(P, Q) -- относительный знак корней сравнимых клеток P, Q (0/1).
    Элемент: (вид, целевая линия, опорная линия или None, T, t опорной, развёрнута ли)."""
    preds = []
    for name, coef, (P1, P2, P3) in LINES:
        t = tval(b, c, coef)
        if t % p == 0:
            preds.append(("L1", name, None, T_of_flip((sbit(P1, P2), 0, sbit(P3, P2))), t, False))
    for (n1, c1, L1), (n2, c2, L2) in itertools.permutations(LINES, 2):
        t1, t2 = tval(b, c, c1), tval(b, c, c2)
        if t1 % p == 0 or t2 % p == 0:
            continue
        if all(congruent(P, Q, b, c, p) for P, Q in zip(L1, L2)):
            L2o, rev = L2, False
        elif all(congruent(P, Q, b, c, p) for P, Q in zip(L1, L2[::-1])):
            L2o, rev = L2[::-1], True
        else:
            continue
        f = tuple(sbit(P, Q) for P, Q in zip(L1, L2o))
        preds.append(("L2", n2, n1, T_of_flip(f), t1, rev))
    return preds


def consistent(cl, preds, p):
    for kind, tgt, ref, T, t, rev in preds:
        tc = torsion_classes(t, p)[T]
        g = cl[tgt]
        if kind == "L1":
            if g != tc:
                return False
        else:
            got = (g[0], g[2], g[1]) if rev else g
            if got != mul3(cl[ref], tc):
                return False
    return True


def sample_roots(p, a, b, c):
    r = {}
    for P in CELLS:
        v = cellval(a, b, c, P)
        if v % p == 0:
            return None
        s, ok = padic_sqrt(v, p, NPREC)
        if not ok:
            return None
        r[P] = s if random.random() < 0.5 else (-s) % p ** NPREC
    return r


def classes(r, p):
    out = {}
    for name, coef, (P1, P2, P3) in LINES:
        A, B, C = r[P1] + r[P2], r[P2] + r[P3], r[P1] + r[P3]
        out[name] = (local_class(A * B, p, NPREC), local_class(A * C, p, NPREC), local_class(B * C, p, NPREC))
    return out


def true_sbit(r, p):
    def s(P, Q):
        if (r[P] - r[Q]) % p == 0:
            return 0
        if (r[P] + r[Q]) % p == 0:
            return 1
        raise ValueError("клетки не сравнимы")
    return s


def rho_sbit(rho):
    return lambda P, Q: rho[P] ^ rho[Q]


ALL_RHO = []
for bits in itertools.product((0, 1), repeat=8):
    rho = {(0, 0): 0}
    for P, bt in zip([P for P in CELLS if P != (0, 0)], bits):
        rho[P] = bt
    ALL_RHO.append(rho)


def rnd_unit(p, K=4):
    while True:
        x = random.randrange(1, p ** K)
        if x % p:
            return x


TYPES = ["generic", "b", "c", "bc", "b+c", "b-c", "b+2c", "b-2c", "2b+c", "2b-c"]


def make_bc(p, typ):
    e = random.choice((1, 1, 2, 3))
    u = lambda: rnd_unit(p)
    if typ == "generic":
        return u(), u()
    if typ == "b":
        return p ** e * u(), u()
    if typ == "c":
        return u(), p ** e * u()
    if typ == "bc":
        return p ** e * u(), p ** random.choice((1, 2, 3)) * u()
    if typ == "b+c":
        b = u(); return b, -b + p ** e * u()
    if typ == "b-c":
        b = u(); return b, b + p ** e * u()
    if typ == "b+2c":
        c = u(); return -2 * c + p ** e * u(), c
    if typ == "b-2c":
        c = u(); return 2 * c + p ** e * u(), c
    if typ == "2b+c":
        b = u(); return b, -2 * b + p ** e * u()
    if typ == "2b-c":
        b = u(); return b, 2 * b + p ** e * u()


def forms(b, c):
    return {"b": b, "c": c, "b+c": b + c, "b-c": b - c, "b+2c": b + 2 * c, "b-2c": b - 2 * c,
            "2b+c": 2 * b + c, "2b-c": 2 * b - c}


def type_ok(p, typ, b, c):
    div = {k for k, v in forms(b, c).items() if v % p == 0}
    if typ == "generic":
        return not div
    if typ == "bc":
        return {"b", "c"} <= div
    return div == {typ}


def get_bc(p, typ):
    for _ in range(200):
        b, c = make_bc(p, typ)
        if type_ok(p, typ, b, c):
            return b, c
    return None


def find_point(p, b, c, tries=3000):
    for _ in range(tries):
        a = rnd_unit(p) ** 2
        r = sample_roots(p, a, b, c)
        if r is not None:
            return a, r
    return None


def find_point_type(p, typ, tries=20000):
    """Совместный выбор (a, b, c): новый (b, c) нужного типа на каждой попытке (для малых p)."""
    for _ in range(tries):
        bc = make_bc(p, typ)
        if not type_ok(p, typ, *bc):
            continue
        a = rnd_unit(p) ** 2
        r = sample_roots(p, a, *bc)
        if r is not None:
            return bc[0], bc[1], a, r
    return None


def bits_of(cl):
    v = []
    for nm in LNAMES:
        for k in (0, 1):
            v += [cl[nm][k][0], 1 if cl[nm][k][1] == -1 else 0]
    return v


def affine_rank(vecs):
    if not vecs:
        return -1
    base = vecs[0]
    rows = []
    for v in vecs[1:]:
        m = 0
        for k, (x, y) in enumerate(zip(v, base)):
            if x != y:
                m |= 1 << k
        rows.append(m)
    basis = []
    for x in rows:
        for bv in basis:
            x = min(x, x ^ bv)
        if x:
            basis.append(x)
    return len(basis)


def model_vectors(b, c, p, nsamp=3000):
    """Векторы, допустимые моделью: rho -- любой, опорные классы хороших линий -- любые неразветвлённые.
    Опорная линия группы -- первая в порядке LINES; прочие хорошие линии группы -- через L2."""
    goods = [nm for nm, coef, L in LINES if tval(b, c, coef) % p]
    vecs = []
    for _ in range(nsamp):
        rho = random.choice(ALL_RHO)
        preds = pred_list(b, c, p, rho_sbit(rho))
        cl = {}
        for kind, tgt, ref, T, t, rev in preds:
            if kind == "L1":
                cl[tgt] = torsion_classes(t, p)[T]
        # группы хороших линий: связность по L2
        for nm in goods:
            if nm in cl:
                continue
            cl[nm] = random.choice(UNRAM)
            changed = True
            while changed:
                changed = False
                for kind, tgt, ref, T, t, rev in preds:
                    if kind == "L2" and ref in cl and tgt not in cl:
                        val = mul3(cl[ref], torsion_classes(t, p)[T])
                        cl[tgt] = (val[0], val[2], val[1]) if rev else val
                        changed = True
        vecs.append(bits_of(cl))
    return vecs


if __name__ == "__main__":
    PR = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 61, 97]
    out = {"A_predictions": {}, "B_negative_control": {}, "C_completeness": [], "meta": {"NPREC": NPREC, "primes": PR}}
    t0 = time.time()
    for p in PR:
        stats = {}
        for typ in TYPES:
            ns, nchk, nbad, nskip, kinds = 0, 0, 0, 0, {}
            for _ in range(30):
                res = find_point_type(p, typ)
                if res is None:
                    nskip += 1
                    if nskip >= 2 and ns == 0:
                        break  # тип при этом p, по-видимому, неосуществим (все клетки -- ненулевые квадраты mod p)
                    continue
                b, c, a, r = res
                cl = classes(r, p)
                preds = pred_list(b, c, p, true_sbit(r, p))
                for pr in preds:
                    kinds[pr[0]] = kinds.get(pr[0], 0) + 1
                    if not consistent(cl, [pr], p):
                        nbad += 1
                        if nbad <= 3:
                            print("НАРУШЕНИЕ", p, typ, a, b, c, pr, flush=True)
                ns += 1
                nchk += len(preds)
            stats[typ] = {"points": ns, "checks": nchk, "violations": nbad, "skipped": nskip, "by_kind": kinds}
        out["A_predictions"][p] = stats
        # ---------------- B
        neg = {}
        for typ in ("b", "c", "bc", "b+c", "b-c"):
            nn, nincons, ntrue_ok = 0, 0, 0
            nmiss = 0
            for _ in range(25):
                res0 = find_point_type(p, typ)
                if res0 is None:
                    nmiss += 1
                    if nmiss >= 2 and nn == 0:
                        break
                    continue
                bc = res0[:2]
                r1, r2 = (res0[2], res0[3]), find_point(p, *bc)
                if r1 is None or r2 is None:
                    continue
                b, c = bc
                cl1, cl2 = classes(r1[1], p), classes(r2[1], p)
                mixed = {nm: (cl1[nm] if nm.startswith("Eb_row") else cl2[nm]) for nm in LNAMES}
                # позитивный контроль процедуры: настоящая точка совместима хотя бы с одним rho
                ntrue_ok += any(consistent(cl2, pred_list(b, c, p, rho_sbit(rho)), p) for rho in ALL_RHO)
                nincons += not any(consistent(mixed, pred_list(b, c, p, rho_sbit(rho)), p) for rho in ALL_RHO)
                nn += 1
            neg[typ] = {"pairs": nn, "true_point_consistent": ntrue_ok, "mixed_inconsistent": nincons}
        out["B_negative_control"][p] = neg
        print(p, "A:", {k: (v["points"], v["checks"], v["violations"]) for k, v in stats.items()},
              "B:", {k: (v["pairs"], v["true_point_consistent"], v["mixed_inconsistent"]) for k, v in neg.items()},
              "%.0fs" % (time.time() - t0), flush=True)
        json.dump(out, open("s3_local_relations.json", "w"), ensure_ascii=False, indent=1)
    # ---------------- C: полнота -- перенесена в s7 / s7b / s10 (при малых p эмпирика ограничена числом a mod p)
    for p in ():
        for typ in ("b", "c", "bc", "b+c", "b-c", "b+2c", "generic"):
            res0 = find_point_type(p, typ)
            if res0 is None:
                continue
            b, c = res0[:2]
            emp = []
            for _ in range(400):
                res = find_point(p, b, c, tries=500)
                if res is None:
                    continue
                a, r = res
                cl = classes(r, p)
                assert consistent(cl, pred_list(b, c, p, true_sbit(r, p)), p)
                emp.append(bits_of(cl))
            mod = model_vectors(b, c, p)
            both = affine_rank(emp + mod) if emp else -1
            rec = {"p": p, "type": typ, "b": b, "c": c, "empirical_points": len(emp),
                   "dim_empirical": affine_rank(emp), "dim_model": affine_rank(mod), "dim_union": both}
            out["C_completeness"].append(rec)
            print("C", rec, flush=True)
    json.dump(out, open("s3_local_relations.json", "w"), ensure_ascii=False, indent=1)
    print("всего %.0fs" % (time.time() - t0))
