# -*- coding: utf-8 -*-
"""s10: ПОЛНАЯ модель локальных классов половинок при нечётном p (все клетки -- p-адические единицы).

Клетки разбиваются на классы K сравнимости по модулю p (значение i*b + j*c mod p). Для каждого класса --
представитель корня w_K, остальные корни r_P = eps_P * w_K (mod p). Каждая сумма корней на линии:
  * P, Q в одном классе, eps_P = eps_Q:  r_P + r_Q = 2 eps_P w_K (единица);
  * P, Q в одном классе, eps_P != eps_Q: r_P + r_Q = (cell_P - cell_Q)/(r_P - r_Q), r_P - r_Q = 2 eps_P w_K
    (нормирование v_p(cell_P - cell_Q), единичная часть (cell_P - cell_Q)/p^v / (2 eps_P w_K));
  * разные классы: r_P + r_Q = eps_P (w_K + eps w_K'), eps = eps_P eps_Q;
    w_K - w_K' = (val_K - val_K')/(w_K + w_K')  ->  класс g_{KK'} := [w_K + w_K'] и константа.
Итог: класс каждой дельта-координаты = (нормирование) * (произведение символов «свободных функций»
[w_K], [g_{KK'}]) * (константа, зависящая от знаков eps аффинно над F_2).
Соотношения между классами = ядро этого отображения; они ДОКАЗАНЫ формулой (для всех точек данного типа).
Полнота (наблюдение): свободные функции на редуцированном многообразии независимы => эмпирическая оболочка
совпадает с модельной (сверка с Q_p-точками при больших p).
Проверки: (1) формула = фактический класс в каждой Q_p-точке; (2) размерности: модель s3 (L1+L2),
полная модель, эмпирика; (3) лишние соотношения полной модели сверх L1+L2 -- в читаемом виде.
Запуск: python3 s10_pair_model.py
"""
import json, random, time, itertools
from fractions import Fraction
import s3_local_relations as S
from s7_joint_relation_count import counts

S.NPREC = 14
random.seed(1010)
CELLS = S.CELLS


def cell_classes(b, c, p):
    cls = {}
    for P in CELLS:
        cls.setdefault((P[0] * b + P[1] * c) % p, []).append(P)
    reps = {}
    for k, L in cls.items():
        for P in L:
            reps[P] = (k, L[0])
    return cls, reps


def lc(x, p):
    return S.rat_local_class(Fraction(x), p)


class Model:
    """Символьная модель класса: словарь свободных функций (множество имён) + константный класс (v, chi)."""

    def __init__(self, b, c, p):
        self.b, self.c, self.p = b, c, p
        self.cls, self.rep = cell_classes(b, c, p)
        keys = sorted(self.cls)
        self.free = ["w%d" % k for k in keys] + ["g%d_%d" % (k1, k2) for k1, k2 in itertools.combinations(keys, 2)]

    def sum_class(self, P, Q, eps):
        """Класс r_P + r_Q: (множество свободных функций, константа (v, chi))."""
        p = self.p
        kP, RP = self.rep[P]
        kQ, RQ = self.rep[Q]
        eP, eQ = eps[P], eps[Q]
        cellP = P[0] * self.b + P[1] * self.c
        cellQ = Q[0] * self.b + Q[1] * self.c
        if kP == kQ:
            if eP == eQ:
                return {"w%d" % kP}, S.mul(lc(2, p), lc(eP, p))
            d = cellP - cellQ
            v = S.vp_int(abs(d), p)
            const = S.mul((v % 2, 1), lc(d // p ** v, p))
            const = S.mul(const, S.mul(lc(2, p), lc(eP, p)))  # деление на 2 eps_P w_K: класс тот же, что умножение
            return {"w%d" % kP}, const
        k1, k2 = sorted((kP, kQ))
        g = "g%d_%d" % (k1, k2)
        const = lc(eP, p)
        if eP * eQ == -1:
            # w_P - w_Q (в порядке P, Q) = (val_P - val_Q)/(w_P + w_Q)
            const = S.mul(const, lc((cellP - cellQ) % p or p, p) if (cellP - cellQ) % p else (0, 1))
        return {g}, const

    def coord_classes(self, eps):
        """Для каждой линии: [(free, const) для AB, AC]."""
        out = {}
        for name, coef, (P1, P2, P3) in S.LINES:
            A = self.sum_class(P1, P2, eps)
            B = self.sum_class(P2, P3, eps)
            C = self.sum_class(P1, P3, eps)
            AB = (A[0] ^ B[0], S.mul(A[1], B[1]))
            AC = (A[0] ^ C[0], S.mul(A[1], C[1]))
            out[name] = (AB, AC)
        return out

    def evaluate(self, eps, freevals):
        """freevals: имя -> символ (+1/-1). Возвращает классы (как в S.classes, но только 2 координаты)."""
        cc = self.coord_classes(eps)
        out = {}
        for nm, coords in cc.items():
            res = []
            for fr, const in coords:
                chi = const[1]
                for f in fr:
                    chi *= freevals[f]
                res.append((const[0], chi))
            out[nm] = tuple(res)
        return out


def actual_free(model, r, p):
    """Фактические значения свободных функций в Q_p-точке: w_K = корень представителя, g = w_K + w_K'."""
    vals = {}
    w = {}
    for k, L in model.cls.items():
        w[k] = r[L[0]]
        vals["w%d" % k] = S.legendre(w[k] % p, p)
    for k1, k2 in itertools.combinations(sorted(model.cls), 2):
        vals["g%d_%d" % (k1, k2)] = S.legendre((w[k1] + w[k2]) % p, p)
    eps = {}
    for P in CELLS:
        k, R = model.rep[P]
        eps[P] = 1 if (r[P] - r[R]) % p == 0 else -1
    return eps, vals


def bits2(cl2):
    v = []
    for nm in S.LNAMES:
        for k in (0, 1):
            v += [cl2[nm][k][0], 1 if cl2[nm][k][1] == -1 else 0]
    return v


out = []
t0 = time.time()
for p, typ, npts in [(13, "bc", 300), (1009, "b", 500), (1009, "c", 500), (1009, "bc", 300), (2003, "b+c", 800),
                     (2003, "b-c", 800), (4001, "b+2c", 800), (4001, "2b-c", 800), (50021, "generic", 800)]:
    res0 = S.find_point_type(p, typ, tries=400000)
    b, c = res0[:2]
    M = Model(b, c, p)
    emp, nform_ok, nform = [], 0, 0
    for _ in range(npts):
        res = S.find_point(p, b, c, tries=20000)
        if res is None:
            continue
        a, r = res
        cl = S.classes(r, p)
        eps, vals = actual_free(M, r, p)
        pred = M.evaluate(eps, vals)
        ok = all(pred[nm][k] == cl[nm][k] for nm in S.LNAMES for k in (0, 1))
        nform += 1
        nform_ok += ok
        emp.append(bits2(cl))
    # полная модель: все знаки eps (по классам, кроме представителей) и все значения свободных функций
    full = []
    nonrep = [P for P in CELLS if M.rep[P][1] != P]
    for _ in range(6000):
        eps = {P: 1 for P in CELLS}
        for P in nonrep:
            eps[P] = random.choice((1, -1))
        vals = {f: random.choice((1, -1)) for f in M.free}
        full.append(bits2(M.evaluate(eps, vals)))
    l12 = S.model_vectors(b, c, p, nsamp=4000)
    rec = {"p": p, "type": typ, "b": b, "c": c, "formula_checks": nform, "formula_ok": nform_ok,
           "dim_L1L2_model": S.affine_rank(l12), "dim_full_model": S.affine_rank(full),
           "dim_empirical": S.affine_rank(emp), "dim_full_union_emp": S.affine_rank(full + emp),
           "joint_full": counts(full), "joint_L1L2": counts(l12)}
    out.append(rec)
    print(p, typ, "формула:", nform_ok, "/", nform, "| dim L1L2:", rec["dim_L1L2_model"], "полная:", rec["dim_full_model"],
          "эмпирика:", rec["dim_empirical"], "полная+эмп:", rec["dim_full_union_emp"],
          "| совместных (все кривые): полная", rec["joint_full"]["joint_relations_all_curves"],
          "L1L2", rec["joint_L1L2"]["joint_relations_all_curves"],
          "| E_b-E_c: полная", rec["joint_full"]["joint_relations_Eb_Ec"], "L1L2", rec["joint_L1L2"]["joint_relations_Eb_Ec"],
          "%.0fs" % (time.time() - t0), flush=True)
json.dump(out, open("s10_pair_model.json", "w"), ensure_ascii=False, indent=1)
