"""Recount: slopes closed by the 12 pairs (odd criterion / exact Selmer), grid n<=400, G1 pairs."""
import sys, math, time, json
from multiprocessing import Pool
sys.path.insert(0, '/home/kep/magicKube/criterion_proof/fable_review')
from fr_selmer import *

CACHE = {}


def pair_info(mn):
    m, n = mn
    S1, S2 = selmer_sets(m, n)
    be = int(round(math.log2(len(S1) * len(S2)))) - 2
    O1, O2 = odd_sets(m, n)
    bo = int(round(math.log2(len(O1) * len(O2)))) - 2
    incl = set(S1) <= set(O1) and set(S2) <= set(O2)
    return (m, n, be, bo, incl, len(S1), len(S2), UNDETERMINED[0])


def main():
    t0 = time.time()
    which = sys.argv[1]
    if which == 'slopes':
        smax = int(sys.argv[2])
        sl = list(slopes(2, smax))
        allpairs = sorted(set(p for r, s in sl for p in twelve_pairs(r, s)))
        print('slopes', len(sl), 'distinct pairs', len(allpairs), flush=True)
        with Pool(20) as pool:
            res = pool.map(pair_info, allpairs, chunksize=50)
        info = {(m, n): (be, bo, incl, s1, s2, und) for m, n, be, bo, incl, s1, s2, und in res}
        print('inclusion violations:', sum(not v[2] for v in info.values()), 'undetermined:', sum(v[5] for v in info.values()), flush=True)
        out = {}
        for lo, hi in [(2, 48), (49, 100), (101, 200)]:
            if lo > smax:
                continue
            tot = 0; c_odd = 0; c_ex = 0; c_odd_pyth_only = 0; c_ex_pyth_only = 0; open_ex = []; open_odd = []
            for r, s in sl:
                if not (lo <= s <= min(hi, smax)):
                    continue
                tot += 1
                prs = twelve_pairs(r, s)
                odd_np = any(info[p][1] == 0 and not pythagorean(*p) for p in prs)
                odd_any = any(info[p][1] == 0 for p in prs)
                ex_np = any(info[p][0] == 0 and not pythagorean(*p) for p in prs)
                ex_any = any(info[p][0] == 0 for p in prs)
                c_odd += odd_np; c_ex += ex_np
                c_odd_pyth_only += (odd_any and not odd_np); c_ex_pyth_only += (ex_any and not ex_np)
                if not ex_np:
                    open_ex.append((r, s))
                if not odd_np:
                    open_odd.append((r, s))
            print(f'denominators {lo}-{min(hi,smax)}: slopes {tot}; closed odd (non-pyth pair) {c_odd} (+{c_odd_pyth_only} only via pythagorean pair); '
                  f'closed exact (non-pyth) {c_ex} (+{c_ex_pyth_only} only via pythagorean pair)', flush=True)
            out[f'{lo}-{hi}'] = dict(total=tot, closed_odd=c_odd, closed_odd_pyth_only=c_odd_pyth_only, closed_exact=c_ex,
                                     closed_exact_pyth_only=c_ex_pyth_only, open_exact=open_ex, open_odd=open_odd)
        json.dump(out, open(f'/home/kep/magicKube/criterion_proof/fable_review/coverage_{smax}.json', 'w'))
    elif which == 'grid':
        N = int(sys.argv[2])
        allpairs = [(m, n) for n in range(2, N + 1) for m in range(1, n) if math.gcd(m, n) == 1]
        print('grid pairs', len(allpairs), flush=True)
        with Pool(20) as pool:
            res = pool.map(pair_info, allpairs, chunksize=100)
        print('odd bound 0:', sum(r[3] == 0 for r in res), '; exact bound 0:', sum(r[2] == 0 for r in res),
              '; inclusion violations:', sum(not r[4] for r in res), '; undetermined:', sum(r[7] for r in res), flush=True)
        json.dump([list(r[:7]) for r in res], open(f'/home/kep/magicKube/criterion_proof/fable_review/grid_{N}.json', 'w'))
    elif which == 'g1':
        N = int(sys.argv[2])
        base = [(m, n) for n in range(2, N + 1) for m in range(1, n) if math.gcd(m, n) == 1]
        sq = [(m * m, n * n) for m, n in base]
        print('G1 pairs', len(base), flush=True)
        with Pool(20) as pool:
            res = pool.map(pair_info, sq, chunksize=20)
        print('G1 odd bound 0:', sum(r[3] == 0 for r in res), '; exact bound 0:', sum(r[2] == 0 for r in res),
              '; inclusion violations:', sum(not r[4] for r in res), '; undetermined:', sum(r[7] for r in res), flush=True)
        json.dump([list(r[:7]) for r in res], open(f'/home/kep/magicKube/criterion_proof/fable_review/g1_{N}.json', 'w'))
    print(f'time {time.time()-t0:.0f}s')


if __name__ == '__main__':
    main()
