# Совместное решето для наклона r/s со всеми восемью клетками.
# python3 run_slope.py r s [nfactors=8] [lmax=500] [maxrank=2] [--norank0]
from jsieve import *
import sys, json


def degenerate_points(r, s):
    """Рациональные точки X_k над z=0 всегда; проверяем z = −1/λ (клетка 0, остальные семь должны быть квадратами)."""
    cells = slope_cells(r, s); out = [QQ(0)]
    for lam in cells:
        z = QQ(-1) / lam
        if all((1 + mu * z) >= 0 and (1 + mu * z).is_square() for mu in cells): out.append(z)
    return out


def build_factors(r, s, verbose=False, maxrank=2, norank0=False):
    cells = slope_cells(r, s)
    facs = []; info = []
    for k in (3, 4):
        for T in itertools.combinations(sorted(cells), k):
            try:
                f = Factor(T)
            except Exception as e:
                info.append(dict(T=list(map(int, T)), err=str(e)[:80])); continue
            info.append(dict(T=list(map(int, T)), lo=f.lo, hi=f.hi, rank=f.rank, proved=f.rank_proved, tors=f.tinv, sat_index=int(f.sat_index)))
            facs.append(f)
    return facs, info


def main():
    r, s = int(sys.argv[1]), int(sys.argv[2])
    nf = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    lmax = int(sys.argv[4]) if len(sys.argv) > 4 else 500
    maxrank = int(sys.argv[5]) if len(sys.argv) > 5 else 2
    norank0 = '--norank0' in sys.argv
    t0 = time.time()
    cells = slope_cells(r, s)
    degz = degenerate_points(r, s)
    print(f'наклон {r}/{s}: клетки {cells}; вырожденные рациональные z на X_k: {degz}', flush=True)
    facs, info = build_factors(r, s)
    ranks = [f.rank for f in facs]
    print(f'множителей {len(facs)} (из 126), время {time.time()-t0:.1f} с; ранги: ' + ' '.join(f'{ranks.count(k)}×{k}' for k in sorted(set(ranks))) +
          f'; недоказанных {sum(1 for f in facs if not f.rank_proved)}; насыщение индексы>1: {sum(1 for f in facs if f.sat_index>1)}', flush=True)
    sel = [i for i, f in enumerate(facs) if f.rank_proved and f.rank <= maxrank and (f.rank >= 1 or not norank0)]
    sel.sort(key=lambda i: (facs[i].rank, len(facs[i].T), i))
    sel = sel[:nf]
    F = [facs[i] for i in sel]
    print('используем: ' + '; '.join(f'{f.name} (r={f.rank}, T={f.tinv})' for f in F), flush=True)
    JS = JointSieve(r, s, cells, F)
    known = []
    for z in degz:
        roots = {lam: (1 + lam * z).sqrt() for lam in cells}
        known += JS.sign_classes(z, roots)
    N0 = {i: lcm(f.tinv) for i, f in enumerate(F)}
    primes = JS.good_primes(3, lmax, range(len(F)))
    lifts = [(i, 2) for i in range(len(F))] + [(i, 3) for i in range(len(F))] + ([(i, 2) for i in range(len(F))] if "--n24" in sys.argv or "--big" in sys.argv else [])
    if "--big" in sys.argv:   # N = 4·2·3·2·2·5 = 480 (кручение (2,2)) или 960 (при (2,4))
        lifts += [(i, 2) for i in range(len(F))] + [(i, 5) for i in range(len(F))]
    JS.verbose = False
    cand, idx, cols, N, log = JS.run(list(range(len(F))), N0, primes, lifts=lifts, known=None)
    for st in log:
        print(f"  {st['step']}: выживших {st['survivors']} (N={[st['N'][i] for i in range(len(F))]})", flush=True)
    S = set(map(tuple, cand.tolist()))
    kn = set(JS.class_tuple(cl, idx, N) for _, cl in known)
    missing = kn - S
    extra = S - kn
    print(f'ИТОГ {r}/{s}: выживших классов {len(S)}, классов вырожденных точек {len(kn)}, потеряно известных {len(missing)}, НЕВЫРОЖДЕННЫХ выживших {len(extra)}; N={[N[i] for i in range(len(F))]}; простых {len(primes)} (≤{lmax}); время {time.time()-t0:.1f} с', flush=True)
    if extra and len(extra) <= 40:
        for t in sorted(extra): print('   невырожденный класс:', t)
    res = dict(slope=f'{r}/{s}', cells=cells, degenerate_z=[str(z) for z in degz], factors=[dict(T=list(map(int, f.T)), rank=f.rank, tors=f.tinv, sat_index=int(f.sat_index), gens=[[str(c) for c in P] for P in f.gens]) for f in F],
               ranks_all=info, N=[N[i] for i in range(len(F))], primes=primes, survivors=len(S), known=len(kn), missing=len(missing), extra=len(extra),
               extra_list=[list(t) for t in sorted(extra)][:200], log=log, time=time.time() - t0)
    def conv(o):
        if isinstance(o, dict): return {str(k): conv(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)): return [conv(v) for v in o]
        if isinstance(o, (bool, str, float)) or o is None: return o
        try: return int(o)
        except Exception: return str(o)
    with open(f'res_{r}_{s}.json', 'w') as fh: json.dump(conv(res), fh)


if __name__ == '__main__':
    main()
