"""Полнота GPU-сита: для случайных строк i полный CPU-перебор j > i тем же модульным фильтром (numpy, независимо от ядра)
и сравнение с множеством выживших GPU в этих строках. Плюс все строки i < 50 и последние 50."""
import sys, numpy as np
B = int(sys.argv[1]); nrows = int(sys.argv[2]) if len(sys.argv) > 2 else 200
lines = open(f'gen_{B}_meta.txt').read().split('\n'); n, k = map(int, lines[0].split()); PR = list(map(int, lines[1].split()))
mw = [list(map(int, lines[2 + l].split())) for l in range(k)]
if len(mw[0]) == 256:   # 128-битные маски
    bits = [[(row[2*x] | (row[2*x+1] << 64)) for x in range(128)] for row in mw]
else:
    bits = [row for row in mw]
res = np.fromfile(f'gen_{B}.bin', dtype=np.uint8).reshape(n, k)
pairs = np.fromfile(f'pairs_{B}.bin', dtype=np.uint32).reshape(-1, 2)
rng = np.random.default_rng(12345)
rows = sorted(set(list(range(50)) + list(range(n - 51, n - 1)) + list(rng.choice(n - 1, nrows, replace=False))))
# таблица допустимости allowed[l][x][y] как bool
allow = [np.array([[(bits[l][x] >> y) & 1 for y in range(PR[l] + 1)] for x in range(PR[l] + 1)], dtype=bool) for l in range(k)]
gpu = {}
sel = np.isin(pairs[:, 0], np.array(rows, dtype=np.uint32))
for i, j in pairs[sel]: gpu.setdefault(int(i), set()).add(int(j))
bad = 0; tot = 0
for i in rows:
    J = np.arange(i + 1, n)
    ok = np.ones(len(J), dtype=bool)
    for l in range(k):
        ok &= allow[l][res[i, l], res[J, l]]
    cpu = set((J[ok]).tolist()); tot += len(J)
    if cpu != gpu.get(i, set()): bad += 1; print("MISMATCH row", i, len(cpu), len(gpu.get(i, set())))
print(f"B={B}: rows checked {len(rows)}, pairs checked on CPU {tot}, mismatching rows {bad}")
