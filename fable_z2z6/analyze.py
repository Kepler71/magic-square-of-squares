import json, glob, sys
from collections import Counter
BMAX = sys.argv[1] if len(sys.argv)>1 else '60'
recs = []
for f in sorted(glob.glob(f'grid/shard_{BMAX}_*.jsonl')):
    for line in open(f): recs.append(json.loads(line))
print("pairs:", len(recs))
err3 = [r for r in recs if 'isog3_error' in r]; print("isog3 errors:", len(err3), [ (r['a'],r['b'],r['isog3_error'][:80]) for r in err3][:5])
errp = [r for r in recs if 'pari_error' in r]; print("pari errors:", len(errp))
viol = 0; c = Counter()
for r in recs:
    if 'isog3' not in r or 'pari' not in r: continue
    b3 = r['isog3']['bound']; b2 = min(x['bound'] for x in r['isog2']); b2s = [x['bound'] for x in r['isog2']]
    lo, hi = r['pari']
    if b3 < lo or b2 < lo: viol += 1; print("VIOLATION", r['a'], r['b'], b3, b2s, r['pari'])
    c['pari_hi0'] += (hi == 0)
    c['pari_lo>0'] += (lo > 0)
    c['isog3_0'] += (b3 == 0)
    c['isog2_0'] += (b2 == 0)
    c['isog2_root0_0'] += (b2s[0]==0); c['isog2_root1_0'] += (b2s[1]==0); c['isog2_root2_0'] += (b2s[2]==0)
    c['any_0'] += (b3 == 0 or b2 == 0)
    c['isog3_0_not_isog2'] += (b3 == 0 and b2 > 0)
    c['isog2_0_not_isog3'] += (b2 == 0 and b3 > 0)
    c['dim_phi>0'] += (r['isog3']['dim_phi'] > 0)
    c['tors12'] += (r['tors']==12)
    c['sha3_gap'] += (b3 > hi)
    c['sha2_gap'] += (b2 > hi)
    c['parity_mismatch_3_2'] += ((b3 - b2) % 2 == 1)
print("violations:", viol); 
for k,v in sorted(c.items()): print(f"  {k}: {v}")
