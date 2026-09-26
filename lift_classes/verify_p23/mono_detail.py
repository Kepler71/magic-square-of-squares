# Детализация нарушений монотонности (какой метод, какой переход).
import p23lib as L
res = {}
for p, plan in ((2, [(k, m) for k in range(1, 6) for m in range(0, 3)]),
                (3, [(k, m) for k in range(1, 4) for m in range(0, 2)])):
    for (k, m) in plan:
        K = k + m; K1 = K + 1
        if p ** (2 * K1) > 3 * 10 ** 5:
            continue
        PK = p ** K
        par = {}
        for B in range(PK):
            for C in range(PK):
                cls = L.ball_cells(B, C, m, p, K)
                par[(B, C)] = (L.classify_naive(cls, p), L.classify_pairs(cls, p))
        for B1 in range(p ** K1):
            for C1 in range(p ** K1):
                if m >= 1 and B1 % p == 0 and C1 % p == 0:
                    continue
                pa = par[(B1 % PK, C1 % PK)]
                cls = L.ball_cells(B1, C1, m, p, K1)
                ch = (L.classify_naive(cls, p), L.classify_pairs(cls, p))
                for meth, a, b in zip(('naive', 'pairs'), pa, ch):
                    if a in ('EXCL', 'PASS') and a != b:
                        key = (p, k, m, meth, a, b, L.minval_ball(B1 % PK, C1 % PK, m, p, K))
                        res.setdefault(key, []).append((B1 % PK, C1 % PK))
for key, v in sorted(res.items()):
    print(key, len(v), sorted(set(v))[:8])
