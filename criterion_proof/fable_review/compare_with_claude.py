"""Phase C: element-wise comparison of Claude's selmer_exact / odd_criterion with own fr_selmer."""
import sys, math, time
from multiprocessing import Pool
sys.path.insert(0, '/home/kep/magicKube/criterion_proof/fable_review')
sys.path.insert(0, '/home/kep/magicKube/criterion_proof')
import fr_selmer
import selmer_exact, odd_criterion

def cmp(mn):
    m, n = mn
    S1, S2 = fr_selmer.selmer_sets(m, n)
    C1, C2 = selmer_exact.selmer(m, n)
    O1, O2 = fr_selmer.odd_sets(m, n)
    a, b, _ = odd_criterion.odd_bound(m, n)
    return (m, n, set(S1) == set(C1) and set(S2) == set(C2), (len(O1), len(O2)) == (a, b))

if __name__ == '__main__':
    t0 = time.time()
    grid = [(m, n) for n in range(2, 401) for m in range(1, n) if math.gcd(m, n) == 1]
    sl = list(fr_selmer.slopes(2, 200))
    sp = sorted(set(p for r, s in sl for p in fr_selmer.twelve_pairs(r, s)))
    g1 = [(m*m, n*n) for n in range(2, 201) for m in range(1, n) if math.gcd(m, n) == 1]
    allp = sorted(set(grid) | set(sp) | set(g1))
    with Pool(20) as pool:
        res = pool.map(cmp, allp, chunksize=100)
    print('pairs compared:', len(res), '; Selmer sets differ:', [(m, n) for m, n, ok, _ in res if not ok][:10], sum(not ok for *_, ok, _ in res),
          '; odd counts differ:', sum(not o for *_, o in res), f'; {time.time()-t0:.0f}s')
