# Эвристика BKLPR (Bhargava–Kane–Lenstra–Poonen–Rains): Ш[2^oo] ~ coker_tors(A), A — случайная кососимметрическая над Z_2
# нечётного размера (корранг 1 => ранг 1). Для ранга r: Prob_r(G) ∝ #G^{1-r} / #Sp(G) = #G^{-(r-1)} * Prob_1(G).
# Наблюдение: dim Sel_2 (без кручения) = 5 = r + dim Ш[2]. Считаем P(dim Ш[2] = 5 - r | r) для r = 1, 3, 5.
import random
random.seed(int(1))
n, K, N = 15, 12, 60000          # размер, 2-адическая точность 2^K, число выборок
M = 2^K
counts = {}; wsum = {1: 0.0, 3: 0.0, 5: 0.0}; wd = {1: {}, 3: {}, 5: {}}
for it in range(N):
    A = matrix(ZZ, n, n)
    for i in range(n):
        for j in range(i+1, n):
            v = random.randrange(M); A[i, j] = v; A[j, i] = -v
    d = pari(A).matsnf()
    vals = [ZZ(x).valuation(2) if x != 0 else K for x in d]
    vals = [min(v, K) for v in vals]
    zeros = sum(1 for v in vals if v >= K)          # корранг над Q_2 (приближённо: инвариант делится на 2^K)
    if zeros != 1: continue                          # ранг-1 модель: корранг ровно 1
    tors = [v for v in vals if 0 < v < K]
    dimS = len(tors); logG = sum(tors)
    for r in (1, 3, 5):
        w = 2.0^(-(r - 1)*logG)
        wsum[r] += w; wd[r][dimS] = wd[r].get(dimS, 0.0) + w
for r in (1, 3, 5):
    tot = wsum[r]
    dist = {k_: v/tot for k_, v in sorted(wd[r].items())}
    print(f"rank {r}: P(dim Ш[2] = d) ≈ " + ", ".join(f"{k_}: {v:.2e}" for k_, v in dist.items() if v > 0))
like = {r: wd[r].get(5 - r, 0.0)/wsum[r] for r in (1, 3, 5)}
print("likelihood of observed dim Sel = 5:", {r: f"{v:.2e}" for r, v in like.items()})
