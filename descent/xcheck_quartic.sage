# Контроль (б): сверка спаривания Фишера (ctp_quartic.sage) с методом Касселса (ctp.sage)
# на кривых с ПОЛНЫМ 2-кручением над квадратичным полем. Матрицы должны совпасть поэлементно.
# Запуск:  SEED=1 CASES=... sage xcheck_quartic.sage
import os, sys, time, traceback
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')

CASES = [(10, [0, 18, -11]), (7, [0, 11, -12]), (2, [0, 30, -13]), (13, [0, 31, -25]),
         (5, [0, 7, -32]), (3, [0, 21, -16])]

def run(d, rts, seed=1, do_cassels=True, with_diag=True):
    random.seed(int(seed))
    k = QuadraticField(d, 'w')
    e = [k(t) for t in rts]
    E = EllipticCurve(k, [0, -sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], -prod(e)])
    Cv = Curve3(k, rts)
    Sel = Cv.selmer_full()
    n = Sel.dimension()
    print(f"=== d={d} roots={rts}: dim Sel^2 = {n}")
    I, J = IJ_of_curve(E)
    F = FisherCTP(k, I, J)
    t0 = time.time()
    dl = deltas_full2(F, E, rts, Cv, Sel)
    Mq, quart = ctp_matrix(F, dl, verbose=False, with_diag=with_diag)
    tq = time.time() - t0
    print(f"  Фишер:  rank = {Mq.rank()}, симметрия {Mq == Mq.transpose()}, "
          f"нулевая диагональ {all(Mq[i,i] == 0 for i in range(n))}  ({tq:.0f}s)")
    print(Mq)
    if do_cassels:
        t0 = time.time()
        els = []
        for w in Sel.basis():
            a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n))
            b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
            els.append((k(a), k(b)))
        ct = CTP(Cv); ct.randomize_conic = True; ct.cassels_form = True
        Mc = matrix(GF(2), n, n)
        for i in range(n):
            vals, _ = ct.pair(els[i], els, reps=1)
            for j in range(n):
                Mc[i, j] = vals[j]
        print(f"  Касселс: rank = {Mc.rank()}  ({time.time()-t0:.0f}s)")
        print(Mc)
        same = (Mq == Mc)
        print(f"  *** МАТРИЦЫ СОВПАЛИ: {same}; ранги {Mq.rank()} vs {Mc.rank()}")
        return Mq, Mc, same
    return Mq, None, None

if __name__ == '__main__' or True:
    seed = int(os.environ.get('SEED', '1'))
    sel = os.environ.get('CASES')
    cases = CASES if not sel else [CASES[int(t)] for t in sel.split(',')]
    ok = []
    for d, rts in cases:
        try:
            _, _, same = run(d, rts, seed=seed)
            ok.append((d, rts, same))
        except Exception as ex:
            print(f"!! d={d} {rts}: {type(ex).__name__}: {ex}")
            traceback.print_exc()
            ok.append((d, rts, 'ERR'))
    print()
    for t in ok:
        print(t)
