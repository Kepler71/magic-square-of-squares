# -*- coding: utf-8 -*-
"""s4: строгое доказательство теоремы A (нет универсальных мультипликативных соотношений)
через нормирования вдоль граничных дивизоров многообразия X девяти квадратов.

X: r_ij^2 = a + i b + j c  (i, j in {-1,0,1}),  отображение X -> A^3 = Spec Q[a,b,c] конечно, степени 2^9.
Над общей точкой гиперплоскости H = {L = 0}, не лежащей на плоскостях клеток, все клетки -- единицы,
расширение K/Q(a,b,c) неразветвлено (Куммер с единичными подкоренными, char 0) =>
  (1) X гладко в общей точке каждой компоненты D над H, v_D(L) = 1;
  (2) компоненты D над H = выборы относительных знаков r_P = +-r_Q внутри групп совпадающих на H клеток
      (различные значения клеток на H независимы по модулю квадратов -- проверяется ниже);
  (3) (r_P + r_Q)(r_P - r_Q) = cell_P - cell_Q, поэтому
      v_D(r_P + r_Q) = 1, если P,Q совпадают на H и знак на D противоположный, иначе 0.
Над плоскостью клетки P: r_P^2 = cell_P, ветвление индекса 2, v_D(r_P) = 1, прочие функции -- единицы.
Если произведение функций = const * g^2 в K, то все его нормирования чётны. Полный ранг матрицы
нормирований mod 2 => соотношений нет. Константы имеют нулевые нормирования, поэтому учитываются автоматически.

Часть 2 (контроль правила (3), независимый расчёт в Sage): вдоль кривой (a0, b0, c0) + eps*(направление),
пересекающей H трансверсально в точке над F_p, корни -- степенные ряды по eps над GF(p);
нормирования произведений AB, AC сравниваются с комбинаторным правилом для каждой компоненты.
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s4_valuation_proof.py
"""
import itertools, json, random

IDX = (-1, 0, 1)
CELLS = [(i, j) for i in IDX for j in IDX]

LINES = []  # (имя, (P1,P2,P3)) с P1 = x - t, P2 = x, P3 = x + t
for j in IDX:
    LINES.append(("Eb_row%+d" % j, ((-1, j), (0, j), (1, j))))
for i in IDX:
    LINES.append(("Ec_col%+d" % i, ((i, -1), (i, 0), (i, 1))))
LINES.append(("Ebpc_diag", ((-1, -1), (0, 0), (1, 1))))
LINES.append(("Ebmc_anti", ((-1, 1), (0, 0), (1, -1))))

# гиперплоскости разностей: линейная форма (коэф. при b, c); клетка (i,j) даёт значение a + i b + j c
HYP = {"b": (1, 0), "c": (0, 1), "b+c": (1, 1), "b-c": (1, -1),
       "b+2c": (1, 2), "b-2c": (1, -2), "2b+c": (2, 1), "2b-c": (2, -1)}


def groups_on(h):
    """Группы клеток, значения которых совпадают на гиперплоскости h: L(b,c) = 0.
    Клетки P, Q совпадают на H <=> (iP-iQ, jP-jQ) пропорционален (коэф. при b, c) ядра.
    Ядро формы (lb, lc): направление (b, c) = (lc, -lb)."""
    lb, lc = HYP[h]
    db, dc = lc, -lb  # на H: (b, c) = s*(db, dc)
    val = {}
    for (i, j) in CELLS:
        k = i * db + j * dc  # значение клетки = a + s*k
        val.setdefault(k, []).append((i, j))
    return list(val.values())


def components(h):
    """Список компонент над h: словарь sign[P] in {+1,-1} относительно первой клетки группы."""
    gs = groups_on(h)
    free = [(g[0], P) for g in gs for P in g[1:]]
    comps = []
    for bits in itertools.product((1, -1), repeat=len(free)):
        s = {}
        for g in gs:
            s[g[0]] = 1
        for (ref, P), e in zip(free, bits):
            s[P] = e
        comps.append((s, gs))
    return comps


def same_group(P, Q, gs):
    return any(P in g and Q in g for g in gs)


def v_sum(P, Q, comp):
    s, gs = comp
    return 1 if (same_group(P, Q, gs) and s[P] * s[Q] == -1) else 0


FUNCS = []  # (имя, функция: (тип, данные) -> нормирование)
for name, (P1, P2, P3) in LINES:
    FUNCS.append((name + ":AB", ("prod", [(P1, P2), (P2, P3)])))
    FUNCS.append((name + ":AC", ("prod", [(P1, P2), (P1, P3)])))
for h in HYP:
    FUNCS.append(("form:" + h, ("form", h)))
for P in CELLS:
    FUNCS.append(("root:r%+d%+d" % P, ("root", P)))

# столбцы: компоненты над гиперплоскостями разностей + плоскости клеток
COLS = []
for h in HYP:
    for comp in components(h):
        COLS.append(("hyp", h, comp))
for P in CELLS:
    COLS.append(("cell", P, None))


def valuation(func, col):
    kind, data = func
    ctype, h, comp = col
    if ctype == "cell":
        return 1 if (kind == "root" and data == h) else 0
    if kind == "root":
        return 0
    if kind == "form":
        return 1 if data == h else 0
    return sum(v_sum(P, Q, comp) for (P, Q) in data)


def rank_f2(rows):
    basis = []
    for x in rows:
        for bv in basis:
            x = min(x, x ^ bv)
        if x:
            basis.append(x)
    return len(basis)


rows = []
for f in FUNCS:
    m = 0
    for k, col in enumerate(COLS):
        if valuation(f[1], col) % 2:
            m |= 1 << k
    rows.append(m)

out = {"n_functions": len(FUNCS), "n_divisors": len(COLS),
       "components_per_hyperplane": {h: len(components(h)) for h in HYP},
       "rank_all": rank_f2(rows),
       "rank_only_16_delta_coords": rank_f2(rows[:16]),
       "rank_Eb_Ec_12_delta_coords": rank_f2(rows[:12])}

# ---------------- независимость значений клеток на H по модулю квадратов (условие (2)):
# на H значения клеток -- различные линейные формы от (a, s); попарно непропорциональны => неприводимы и
# независимы в Q(a,s)^*/квадраты (разные простые элементы кольца многочленов).
indep = {}
for h in HYP:
    ks = sorted(set(i * HYP[h][1] - j * HYP[h][0] for (i, j) in CELLS))
    indep[h] = len(ks)  # число различных форм a + s*k
out["distinct_cell_values_on_H"] = indep

# ---------------- часть 2: контроль правила нормирований степенными рядами над GF(p) (Sage)
try:
    from sage.all import GF, PowerSeriesRing
    p = 10007
    F = GF(p)
    R = PowerSeriesRing(F, 'e', default_prec=12)
    e = R.gen()
    random.seed(20260926)
    checks, fails = 0, 0
    for h in ["b", "c", "b+c", "b-c"]:
        lb, lc = HYP[h]
        db, dc = lc, -lb
        for trial in range(6):
            # точка на H: b0 = s*db, c0 = s*dc, все клетки -- ненулевые квадраты в F_p
            while True:
                a0, s0 = F.random_element(), F.random_element()
                b0, c0 = s0 * db, s0 * dc
                vals = {P: a0 + P[0] * b0 + P[1] * c0 for P in CELLS}
                if all(v != 0 and v.is_square() for v in vals.values()) and s0 != 0:
                    break
            # трансверсальное направление: (b, c) = (b0, c0) + e*(nb, nc), L(nb, nc) != 0
            nb, nc = (1, 0) if lb != 0 else (0, 1)
            comps = components(h)
            for comp in random.sample(comps, min(8, len(comps))):
                s, gs = comp
                roots = {}
                for g in gs:
                    base = vals[g[0]].sqrt()
                    for P in g:
                        cell = vals[P] + e * (P[0] * nb + P[1] * nc)
                        # sqrt(cell) = base*sqrt(1 + e*.../base^2) как ряд
                        u = cell / vals[g[0]]
                        roots[P] = s[P] * base * u.sqrt()
                for name, (P1, P2, P3) in LINES:
                    A = roots[P1] + roots[P2]
                    B = roots[P2] + roots[P3]
                    C = roots[P1] + roots[P3]
                    for tag, val_series, pairs in (("AB", A * B, [(P1, P2), (P2, P3)]),
                                                   ("AC", A * C, [(P1, P2), (P1, P3)])):
                        pred = sum(v_sum(P, Q, comp) for (P, Q) in pairs)
                        got = val_series.valuation()
                        checks += 1
                        fails += (got != pred)
    out["power_series_check"] = {"p": p, "checks": checks, "failures": fails}
except ImportError:
    out["power_series_check"] = "sage недоступен"

print(json.dumps(out, ensure_ascii=False, indent=1))
json.dump(out, open("s4_valuation_proof.json", "w"), ensure_ascii=False, indent=1)
