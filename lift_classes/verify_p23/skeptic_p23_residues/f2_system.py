"""Шаг 1 (линейная алгебра над F2), независимая проверка.
Клетки f_ij = 1 + i*b + j*c, индекс (i+1)*3+(j+1); центр = 4 (класс 0).
Проверяем: пространство решений 8 уравнений (строки/столбцы/диагонали) на 8 нецентральных
неизвестных = {M(S,U)}; веса ненулевых векторов; [T]=S, [L]=S+U; какие подмножества линий
дают то же пространство (т.е. сколько произведений реально нужно)."""
from itertools import product, combinations
import json

def idx(i, j): return (i+1)*3 + (j+1)
LINES = {}
for i in (-1, 0, 1): LINES[f"row{i:+d}"] = [idx(i, j) for j in (-1, 0, 1)]
for j in (-1, 0, 1): LINES[f"col{j:+d}"] = [idx(i, j) for i in (-1, 0, 1)]
LINES["diag"] = [idx(-1, -1), 4, idx(1, 1)]
LINES["anti"] = [idx(-1, 1), 4, idx(1, -1)]
NONC = [k for k in range(9) if k != 4]

def solutions(line_names):
    sols = []
    for bits in product((0, 1), repeat=8):   # 256 векторов, явная граница
        x = [0]*9
        for k, bt in zip(NONC, bits): x[k] = bt
        if all(sum(x[k] for k in LINES[n]) % 2 == 0 for n in line_names):
            sols.append(tuple(x))
    return sols

def M(S, U):
    x = [0]*9
    x[idx(-1, -1)] = S; x[idx(-1, 0)] = U; x[idx(-1, 1)] = S ^ U
    x[idx(0, -1)] = U; x[idx(0, 1)] = U
    x[idx(1, -1)] = S ^ U; x[idx(1, 0)] = U; x[idx(1, 1)] = S
    return tuple(x)

out = {}
all8 = list(LINES)
sol = solutions(all8)
MS = sorted({M(S, U) for S in (0, 1) for U in (0, 1)})
out["n_solutions_8lines"] = len(sol)
out["solutions_equal_M(S,U)"] = sorted(sol) == MS
out["weights_nonzero"] = sorted(sum(v) for v in sol if any(v))
# [T] = [f(-1,0)]+[f(1,1)]+[f(0,-1)],  [L] = [f(-1,0)]+[f(1,-1)]+[f(0,1)]
okTL = True
for S, U in product((0, 1), repeat=2):
    x = M(S, U)
    t = (x[idx(-1, 0)] + x[idx(1, 1)] + x[idx(0, -1)]) % 2
    l = (x[idx(-1, 0)] + x[idx(1, -1)] + x[idx(0, 1)]) % 2
    okTL &= (t == S) and (l == (S ^ U))
    # угол (-1,-1) несёт S, угол (-1,1) несёт S+U
    okTL &= x[idx(-1, -1)] == S and x[idx(-1, 1)] == (S ^ U)
out["T_is_S_and_L_is_S+U"] = okTL
# нет ненулевого решения с носителем <=3 клеток (используется в Шаге 2)
out["min_weight_nonzero"] = min(sum(v) for v in sol if any(v))
# минимальные подмножества линий, дающие то же пространство решений
same = {}
for k in range(1, 9):
    same[k] = [list(c) for c in combinations(all8, k) if sorted(solutions(c)) == MS]
out["n_subsets_with_same_space_by_size"] = {k: len(v) for k, v in same.items()}
out["minimal_size"] = min(k for k, v in same.items() if v)
out["example_minimal_subsets"] = same[out["minimal_size"]][:6]
print(json.dumps(out, indent=1, ensure_ascii=False))
json.dump(out, open(__file__.replace(".py", ".json"), "w"), indent=1, ensure_ascii=False)
