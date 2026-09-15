# Перепись совместного решета: все наклоны r/s, SMIN ≤ s ≤ SMAX, фиксированные параметры (8 множителей ранга 1–2, N=12, ℓ ≤ 500).
# python3 census_sieve.py SMIN SMAX [nf=8] [lmax=500]
from jsieve import *
from run_slope import degenerate_points, build_factors
import sys, json, os
from math import gcd


def work(args):
    r, s, nf, lmax = args
    t0 = time.time()
    try:
        cells = slope_cells(r, s); degz = degenerate_points(r, s)
        facs, info = build_factors(r, s)
        sel = [i for i, f in enumerate(facs) if f.rank_proved and 1 <= f.rank <= 2]
        sel.sort(key=lambda i: (facs[i].rank, len(facs[i].T), i)); sel = sel[:nf]
        F = [facs[i] for i in sel]
        nr0 = sum(1 for f in facs if f.rank_proved and f.rank == 0)
        if len(F) < nf: return dict(slope=f'{r}/{s}', err='мало множителей', nfac=len(F), rank0=nr0)
        JS = JointSieve(r, s, cells, F, verbose=False)
        p0 = [f.decompose(f.point(0, 1)) for f in F]   # класс точки z=0 (квадрат из единиц) по образующим
        known = []
        for z in degz:
            roots = {lam: (1 + lam * z).sqrt() for lam in cells}; known += JS.sign_classes(z, roots)
        N0 = {i: lcm(f.tinv) for i, f in enumerate(F)}
        primes = JS.good_primes(3, lmax, range(len(F)))
        lifts = [(i, 2) for i in range(len(F))] + [(i, 3) for i in range(len(F))]
        cand, idx, cols, N, log = JS.run(list(range(len(F))), N0, primes, lifts=lifts)
        S = set(map(tuple, cand.tolist())); kn = set(JS.class_tuple(cl, idx, N) for _, cl in known)
        return dict(slope=f'{r}/{s}', degz=[str(z) for z in degz], rank0=nr0, nfac_all=len(facs), used=[f.name for f in F], ranks=[f.rank for f in F],
                    tors=[f.tinv for f in F], p0=[[int(x) for x in c] for c in p0], sat_index=[int(f.sat_index) for f in F], N=[int(N[i]) for i in range(len(F))], nprimes=len(primes), survivors=len(S), known=len(kn),
                    missing=len(kn - S), extra=len(S - kn), extra_list=[list(map(int, t)) for t in sorted(S - kn)][:50],
                    maxcand=max(max(st[1] for st in L['stats']) for L in log if L['stats']), time=time.time() - t0)
    except Exception as e:
        return dict(slope=f'{r}/{s}', err=str(e)[:200], time=time.time() - t0)


if __name__ == '__main__':
    import multiprocessing as mp
    SMIN, SMAX = int(sys.argv[1]), int(sys.argv[2])
    nf = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    lmax = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    todo = [(r, s, nf, lmax) for s in range(SMIN, SMAX + 1) for r in range(1, s) if gcd(r, s) == 1 and 2 * r != s]
    fn = f'census_{SMIN}_{SMAX}.jsonl'; done = set()
    if os.path.exists(fn): done = {json.loads(l)['slope'] for l in open(fn)}
    todo = [t for t in todo if f'{t[0]}/{t[1]}' not in done]
    print('к расчёту', len(todo), flush=True)
    out = open(fn, 'a'); n = 0; nextra = 0
    with mp.get_context('fork').Pool(8) as pool:
        for res in pool.imap_unordered(work, todo):
            out.write(json.dumps(res) + '\n'); out.flush(); n += 1
            if res.get('extra'): nextra += 1; print('НЕВЫРОЖДЕННЫЕ КЛАССЫ:', res['slope'], res['extra'], res['extra_list'][:3], flush=True)
            if 'err' in res: print('ошибка', res['slope'], res['err'], flush=True)
            if n % 25 == 0: print(n, '/', len(todo), 'с невырожденными', nextra, flush=True)
    print('готово', n, 'с невырожденными', nextra, flush=True)
