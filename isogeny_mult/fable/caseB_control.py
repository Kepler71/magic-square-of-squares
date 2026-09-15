# Fable 15.09. Контроль решателя случая B: подсаживаем решение. Для пары (i2,i3) и кручения t берём точку A₀ = n·G + t′,
# вычисляем z₂(A₀) и подменяем мёбиусову матрицу M[i3] на M′ так, чтобы z₃′(A₀ − t) = z₂(A₀) (сдвиг константы);
# решатель обязан вернуть z₂(A₀) среди кандидатов.
from sage.all import *
import sys
sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable')
from demjanenko2 import *

r, s = 204, 247
S = slope_cells(r, s)
Ts = [cells_of(Sg, r, s) for Sg in CLASSES[BIG[1]]]
com = Common(Ts, verbose=False)
i2, i3 = 2, 0
G = com.gens[0]
found = 0; total = 0
for n in (1, 2, 3):
    for ti, t in enumerate(com.tors):
        A0 = n * G + com.tors[(ti + 1) % len(com.tors)]
        z0 = com.z_of_point(i2, A0)
        # z₃′(A) := z₃(A) + δ, δ = z0 − z₃(A0 − t); M: z ↦ x_min, поэтому M′ = M ∘ (z ↦ z − δ) = [[M00, M01 − δ M00],[M10, M11 − δ M10]]
        z3 = com.z_of_point(i3, A0 - t)
        if z3 is infinity or z0 is infinity: continue
        delta = z0 - z3
        M = com.M[i3]; Mp = matrix(QQ, [[M[0, 0], M[0, 1] - delta * M[0, 0]], [M[1, 0], M[1, 1] - delta * M[1, 0]]])
        saved = com.M[i3]; com.M[i3] = Mp
        cands, triv = case_B_candidates(com, i2, i3, verbose=False)
        com.M[i3] = saved
        total += 1; ok = z0 in cands; found += ok
        print(f'n = {n}, t = #{ti}: подсаженное z найдено: {ok} (кандидатов {len(cands)})')
print(f'найдено {found} из {total}')
