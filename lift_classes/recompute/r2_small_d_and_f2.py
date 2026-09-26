# -*- coding: utf-8 -*-
"""r2: (а) наименьший нетривиальный положительный бесквадратный d = 1 mod 24 со всеми простыми = 1 mod 4;
(б) конечные леммы доказательства, проверяемые полным перебором:
   - решение 8 уравнений классов над F_2 (9 клеток, центр включён) -> пространство размерности 2 вида
       S U S+U / U 0 U / S+U U S ; [T]=S, [L]=S+U; веса ненулевых векторов чётностей = {4,6,6};
   - ни один ненулевой допустимый вектор не лежит на <=3 клетках, лежащих на одной прямой (в т.ч. любых 3 клетках);
   - в F_5 и F_7 нет трёх различных ненулевых квадратов в АП (замечание Codex о бесконечном цикле s3).
Всё — чистый Python, все циклы с явной границей."""
import json, itertools

def fac(n):
    out = {}; p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: out[n] = out.get(n, 0) + 1
    return out

res = {}
BOUND = 100000
adm, sqf24 = [], []
for d in range(2, BOUND + 1):
    if d % 24 != 1: continue
    f = fac(d)
    if all(e == 1 for e in f.values()):
        sqf24.append(d)
        if all(p % 4 == 1 for p in f):
            adm.append(d)
res["first_admissible"] = adm[:15]
res["min_admissible"] = adm[0]
res["first_squarefree_1mod24_without_mod4_condition"] = sqf24[:10]
res["numbers_1mod24_below_73"] = [d for d in range(1, 73) if d % 24 == 1]
print("первые допустимые d:", adm[:15])
print("минимум:", adm[0])
print("бесквадратные d=1 mod 24 без условия mod 4:", sqf24[:10])
print("числа = 1 mod 24 меньше 73:", res["numbers_1mod24_below_73"], "; 25=5^2, 49=7^2")

# ---- F_2: клетки (i,j), i,j in {-1,0,1}; линии арифметической сетки
cells = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
idx = {c: k for k, c in enumerate(cells)}
lines = []
for i in (-1, 0, 1): lines.append([(i, j) for j in (-1, 0, 1)])   # фикс. i
for j in (-1, 0, 1): lines.append([(i, j) for i in (-1, 0, 1)])   # фикс. j
lines.append([(-1, -1), (0, 0), (1, 1)]); lines.append([(-1, 1), (0, 0), (1, -1)])
T_cells = [(-1, 0), (1, 1), (0, -1)]   # (1-rz)(1+(r+s)z)(1-sz)
L_cells = [(-1, 0), (1, -1), (0, 1)]   # (1-rz)(1+(r-s)z)(1+sz)
sols = []
for bits in range(2 ** 9):
    x = [(bits >> k) & 1 for k in range(9)]
    if all(sum(x[idx[c]] for c in ln) % 2 == 0 for ln in lines):
        sols.append(x)
print("решений системы 8 линий над F_2 (центр свободен):", len(sols))
ok_shape = True; weights = []; T_eq_S = True; L_eq_SU = True
for x in sols:
    S = x[idx[(1, 1)]]; U = x[idx[(1, 0)]]
    expect = {(-1, -1): S, (1, 1): S, (-1, 1): S ^ U, (1, -1): S ^ U,
              (-1, 0): U, (1, 0): U, (0, -1): U, (0, 1): U, (0, 0): 0}
    ok_shape &= all(x[idx[c]] == v for c, v in expect.items())
    T_eq_S &= (sum(x[idx[c]] for c in T_cells) % 2 == S)
    L_eq_SU &= (sum(x[idx[c]] for c in L_cells) % 2 == (S ^ U))
    if any(x): weights.append(sum(x))
print("форма матрицы S/U подтверждена:", ok_shape, "; [T]=S:", T_eq_S, "; [L]=S+U:", L_eq_SU)
print("веса ненулевых решений:", sorted(weights))
# носитель на <=3 клетках
supp3 = [x for x in sols if any(x) and sum(x) <= 3]
print("ненулевых решений с носителем <=3 клеток:", len(supp3))
res.update({"f2_solutions": len(sols), "f2_shape_ok": ok_shape, "T_eq_S": T_eq_S, "L_eq_SU": L_eq_SU,
            "nonzero_weights": sorted(weights), "nonzero_support_le3": len(supp3)})
# в F_2^k (k=2,3 — локальные группы классов) линейная алгебра покомпонентная: те же выводы.

# ---- АП из трёх различных ненулевых квадратов в F_5, F_7, F_11, F_13
ap = {}
for p in (5, 7, 11, 13):
    Q = {x * x % p for x in range(1, p)}
    found = [(a, c) for a in range(p) for c in range(1, p)
             if {(a - c) % p, a % p, (a + c) % p} <= Q and len({(a - c) % p, a, (a + c) % p}) == 3]
    ap[p] = len(found)
print("число АП (a-c,a,a+c) из различных ненулевых квадратов:", ap)
res["ap_of_three_distinct_nonzero_squares"] = ap
json.dump(res, open("r2_small_d_and_f2.json", "w"), indent=1)
