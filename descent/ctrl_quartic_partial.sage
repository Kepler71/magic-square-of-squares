# Контроль (а): кривые y^2 = x(x^2+ax+b) над k = Q(sqrt d) с ОДНОЙ точкой порядка 2,
# истинный ранг известен точно (PARI ellrank, lo == hi):  rank E(k) = rank E0(Q) + rank E0^{(d)}(Q).
# Проверяем: rank <= dim Sel^2 - 1 - rank CTP  (и насколько граница точна).
# Запуск:  SEED=1 NEED=8 sage ctrl_quartic_partial.sage
import os, sys, time, traceback
load('/home/kep/magicKube/descent/ek_descent_partial.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
proof.number_field(False)


def true_rank_Q(a, b):
    E = EllipticCurve(QQ, [0, a, 0, b, 0])
    if E.discriminant() == 0:
        return None
    r = pari(E).ellrank()
    lo, hi = ZZ(r[0]), ZZ(r[1])
    return lo if lo == hi else None


def true_rank_k(a, b, d):
    r1 = true_rank_Q(a, b)
    if r1 is None:
        return None
    r2 = true_rank_Q(a * d, b * d ^ 2)
    if r2 is None:
        return None
    return r1 + r2


def run_case(a, b, d, min_excess=2, verbose=False):
    k = QuadraticField(d, 'r')
    if k(a ^ 2 - 4 * b).is_square():
        return None
    tr = true_rank_k(a, b, d)
    if tr is None:
        return None
    dim, bd, PD = partial_descent_bound(k, a, b, verbose=verbose)
    exc = dim - 1 - tr
    if exc < min_excess:
        return ('skip', a, b, d, tr, dim, exc)
    E = EllipticCurve(k, [0, k(a), 0, k(b), 0])
    I, J = IJ_of_curve(E)
    F = FisherCTP(k, I, J)
    dl = deltas_partial(F, E, PD, PD.Sel)
    t0 = time.time()
    M, quart = ctp_matrix(F, dl, verbose=False, with_diag=True)
    n = M.nrows()
    sym = (M == M.transpose()); dg0 = all(M[i, i] == 0 for i in range(n))
    rk = M.rank()
    return ('case', a, b, d, tr, dim, exc, rk, dim - 1 - rk, sym, dg0, round(time.time() - t0))


CASES = [(1, -3, 5), (1, 1, 5), (-2, -4, 3), (3, 1, -3), (-5, 5, 2), (2, -7, 13), (5, 5, 6),
         (7, 5, 17), (-3, -11, 2), (6, 7, -1), (9, 2, 3), (-4, 2, 6), (1, -3, -1), (1, -3, 2)]

if __name__ == '__main__' or True:
    random.seed(int(os.environ.get('SEED', '1')))
    need = int(os.environ.get('NEED', '8'))
    scan = os.environ.get('SCAN', '1') == '1'
    rows = []
    tried = set()
    cands = list(CASES)
    if scan:
        extra = []
        for d in [2, 3, 5, 6, 7, 10, 11, 13, 14, 15, 17, -1, -2, -3, -7]:
            for a in range(-25, 26):
                for b in range(-25, 26):
                    if b == 0 or a * a - 4 * b == 0:
                        continue
                    if gcd(gcd(a, b), 1) != 1:
                        continue
                    extra.append((a, b, d))
        random.shuffle(extra)
        cands += extra
    found = 0
    for (a, b, d) in cands:
        if (a, b, d) in tried:
            continue
        tried.add((a, b, d))
        try:
            res = run_case(a, b, d)
        except Exception as ex:
            print(f"ERR a={a} b={b} d={d}: {type(ex).__name__}: {str(ex)[:150]}")
            continue
        if res is None:
            continue
        if res[0] == 'skip':
            print("skip", res[1:])
            continue
        print("CASE", res[1:])
        rows.append(res)
        found += 1
        if found >= need:
            break
    print()
    print("| a | b | d | rank E(k) | dim Sel^2 | rank CTP | граница | точна | симм | диаг0 |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        _, a, b, d, tr, dim, exc, rk, bnd, sym, dg0, tt = r
        print(f"| {a} | {b} | {d} | {tr} | {dim} | {rk} | {bnd} | {'да' if bnd == tr else 'нет'} | {sym} | {dg0} |")
        assert bnd >= tr, "ГРАНИЦА НИЖЕ ИСТИННОГО РАНГА"
