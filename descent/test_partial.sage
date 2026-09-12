# Контроли для ek_descent_partial.sage (Claude, 2026-09-12).
# 1) E0/Q: y^2 = x(x^2+ax+b), одна рациональная точка порядка 2; k = Q(sqrt d).
#    rank E0(k) = rank E0(Q) + rank E0^{(d)}(Q)  (точно, PARI ellrank, требуем lo == hi).
load('/home/kep/magicKube/descent/ek_descent_partial.sage')
import sys, os
if os.environ.get('PROOF'): proof.number_field(True)   # контроль без GRH

def true_rank_Q(a, b):
    E = EllipticCurve(QQ, [0, a, 0, b, 0])
    r = pari(E).ellrank()
    lo, hi = ZZ(r[0]), ZZ(r[1])
    assert lo == hi, f"ellrank не точен для {(a,b)}: {lo}..{hi}"
    return lo

def true_rank_k(a, b, d):
    # твист на d: y^2 = x(x^2 + a d x + b d^2)
    return true_rank_Q(a, b) + true_rank_Q(a*d, b*d^2)

CASES = [
    # (a, b, d)
    (1, -3, 5), (1, -3, -1), (1, -3, 2), (1, 1, 5), (1, 1, -7),
    (-2, -4, 3), (3, 1, 5), (3, 1, -3), (-5, 5, 2), (2, -7, 13),
    (1, -5, -1), (5, 5, 6), (7, 5, 17), (-3, -11, 2), (4, -1, 5), (6, 7, -1),
    (1, -6, 5), (-1, -6, 2), (9, 2, 3), (-4, 2, 6),
]

def main():
    import os
    seed = int(os.environ.get('SEED', '1'))
    random.seed(seed)
    rows = []
    for (a, b, d) in CASES:
        k = QuadraticField(d, 'r')
        Dl = k(a^2 - 4*b)
        if Dl.is_square():
            print(f"пропуск (a,b,d)={(a,b,d)}: Delta квадрат в k"); continue
        t0 = time.time()
        tr = true_rank_k(a, b, d)
        dim, bd, C = partial_descent_bound(k, a, b, verbose=True)
        rows.append((a, b, d, tr, dim, bd, time.time()-t0))
        print(f"### a={a} b={b} d={d}: истинный ранг {tr}, dim Sel^2 = {dim}, граница {bd}"
              f" {'OK' if bd >= tr else 'ПРОТИВОРЕЧИЕ!!!'}{' (точно)' if bd == tr else ''}  ({time.time()-t0:.0f}s)")
        assert bd >= tr, "ГРАНИЦА НИЖЕ ИСТИННОГО РАНГА — ошибка в реализации"
    print()
    print("| a | b | d | rank E(k) | dim Sel^2 | граница | точна |")
    print("|---|---|---|---|---|---|---|")
    for (a, b, d, tr, dim, bd, t) in rows:
        print(f"| {a} | {b} | {d} | {tr} | {dim} | {bd} | {'да' if bd == tr else 'нет'} |")

main()
