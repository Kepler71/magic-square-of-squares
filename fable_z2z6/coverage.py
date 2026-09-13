# Покрытие наклонов r/s (s<=200) теоремой о 3-изогении и точными 2-изогенными оценками; сравнение со старым покрытием Z/2xZ/4.
import json, glob, sys
from math import gcd
sys.path.insert(0, '/home/kep/magicKube/fable_z2z6')
from explicit_formula import explicit_dims
recs = {}
for f in sorted(glob.glob('grid/shard_200_*.jsonl')):
    for line in open(f):
        r = json.loads(line); recs[(r['a'], r['b'])] = r
print("pairs loaded:", len(recs))
old = json.load(open('/home/kep/magicKube/criterion_proof/coverage_exact.json'))
old_open = set()
for k in old: old_open |= set(old[k]['open'])
def info(a, b):
    r = recs[(a, b)]
    b3 = r['isog3']['bound'] if 'isog3' in r else 99
    b2 = min(x['bound'] for x in r['isog2'])
    ex = explicit_dims(a, b)[2]
    assert ex == b3, (a, b, ex, b3)
    return b3, b2
ranges = [(2, 48), (49, 100), (101, 200)]
res = {}
newly = {}
for lo, hi in ranges:
    tot = c3 = c2 = cany = cold = cnew = c3not2 = c2not3 = 0
    opens = []
    for s in range(lo, hi+1):
        for r in range(1, s):
            if gcd(r, s) != 1 or 2*r == s: continue
            tot += 1
            pairs = {(min(r, s-r), max(r, s-r)), (r, s)}
            b3s = []; b2s = []
            for (a, b) in pairs:
                b3, b2 = info(a, b); b3s.append(b3); b2s.append(b2)
            k3 = min(b3s) == 0; k2 = min(b2s) == 0
            c3 += k3; c2 += k2; cany += (k3 or k2)
            c3not2 += (k3 and not k2); c2not3 += (k2 and not k3)
            sl = f"{r}/{s}"
            isold = sl not in old_open
            cold += isold
            if (k3 or k2) and not isold: cnew += 1; newly.setdefault((lo,hi), []).append((sl, min(b3s), min(b2s)))
            if not (k3 or k2) and not isold: opens.append(sl)
    res[(lo,hi)] = dict(total=tot, isog3=c3, isog2=c2, any=cany, old=cold, new_not_old=cnew, isog3_not_isog2=c3not2, isog2_not_isog3=c2not3, still_open_both=len(opens), open_both=opens)
    print(f"s in [{lo},{hi}]: total {tot}; 3-isog {c3}; 2-isog(best root) {c2}; any {cany}; old Z/2xZ/4 {cold}; NEW (closed now, open before) {cnew}; 3not2 {c3not2}; 2not3 {c2not3}; open after both families {len(opens)}")
json.dump({str(k): v for k, v in res.items()} | {"newly": {str(k): v for k, v in newly.items()}}, open('coverage_result.json', 'w'), indent=1, ensure_ascii=False)
for k, v in newly.items(): print(k, "newly closed:", v)
# трудные наклоны
for sl in ["11/142", "48/163", "79/110", "19/60"]:
    r, s = map(int, sl.split('/'))
    pairs = {(min(r, s-r), max(r, s-r)), (r, s)}
    print(sl, [(p, info(*p)) for p in pairs])
