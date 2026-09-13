"""Compare own Selmer bound with PARI ellrank on random pairs; check odd-criterion inclusions."""
import random, math, sys, time, json
sys.path.insert(0, '/home/kep/magicKube/criterion_proof/fable_review')
from fr_selmer import *
import cypari2
pari = cypari2.Pari()
pari.allocatemem(2 * 10**9)
pari.default("parisizemax", 4 * 10**9)

random.seed(1309)
pairs = set()
while len(pairs) < 300:
    m = random.randint(1, 400); n = random.randint(1, 400)
    if m != n and math.gcd(m, n) == 1:
        pairs.add((min(m, n), max(m, n)))
pairs = sorted(pairs)

t0 = time.time()
rows = []
viol_lower = 0; viol_upper = 0; incl_viol = 0; odd_lt_exact = 0
for (m, n) in pairs:
    S1, S2 = selmer_sets(m, n)
    be = int(round(math.log2(len(S1) * len(S2)))) - 2
    O1, O2 = odd_sets(m, n)
    bo = int(round(math.log2(len(O1) * len(O2)))) - 2
    if not set(S1) <= set(O1) or not set(S2) <= set(O2):
        incl_viol += 1
        print('INCLUSION VIOLATED', m, n, S1, O1, S2, O2)
    if bo < be:
        odd_lt_exact += 1
    E = pari.ellinit([0, m*m + n*n, 0, m*m*n*n, 0])
    try:
        rk = pari.ellrank(E)
        r1, r2 = int(rk[0]), int(rk[1])
    except Exception as ex:
        r1, r2 = -1, -1
        print('PARI FAIL', m, n, ex)
    if r1 >= 0 and be < r1:
        viol_lower += 1; print('BOUND BELOW PARI LOWER', m, n, be, r1)
    if r2 >= 0 and be < r2:
        viol_upper += 1; print('ISOGENY BOUND BELOW FULL 2-DESCENT BOUND', m, n, be, r2)
    rows.append(dict(m=m, n=n, sel_phi=len(S1), sel_phihat=len(S2), exact_bound=be, odd_bound=bo, pari_lower=r1, pari_upper=r2))
print(f'pairs {len(pairs)}, time {time.time()-t0:.0f}s')
print('bound below PARI lower:', viol_lower, '; isogeny bound below PARI upper:', viol_upper)
print('inclusion Sel ⊆ S_odd violated:', incl_viol, '; odd bound < exact bound:', odd_lt_exact)
print('exact bound 0:', sum(r['exact_bound'] == 0 for r in rows), '; odd bound 0:', sum(r['odd_bound'] == 0 for r in rows),
      '; PARI rank 0 (r1=r2=0):', sum(r['pari_lower'] == 0 and r['pari_upper'] == 0 for r in rows),
      '; PARI undecided:', sum(r['pari_lower'] != r['pari_upper'] for r in rows))
print('exact bound == PARI upper:', sum(r['exact_bound'] == r['pari_upper'] for r in rows),
      '; exact bound == PARI lower:', sum(r['exact_bound'] == r['pari_lower'] for r in rows))
print('undetermined depth events:', UNDETERMINED[0])
json.dump(rows, open('/home/kep/magicKube/criterion_proof/fable_review/pari_compare.json', 'w'), indent=0)
