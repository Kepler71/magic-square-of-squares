# сетка 1<=a<b<=60, gcd=1: 3-изогения, 2-изогения (3 корня), PARI ellrank. Шарды.
import sys, json, os, time
sys.path.insert(0, '/home/kep/magicKube/fable_z2z6')
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
from two_isog import two_isog_bounds
shard, nsh = int(sys.argv[1]), int(sys.argv[2]); NOPARI = (len(sys.argv) > 4 and sys.argv[4] == "nopari")
BMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 60
pairs = [(a,b) for b in range(2, BMAX+1) for a in range(1, b) if gcd(a,b) == 1]
mine = pairs[shard::nsh]
os.makedirs('/home/kep/magicKube/fable_z2z6/grid', exist_ok=True)
out = open(f'/home/kep/magicKube/fable_z2z6/grid/shard_{BMAX}_{shard}.jsonl', 'w')
t0 = time.time()
for idx, (a,b) in enumerate(mine):
    rec = dict(a=a, b=b)
    try:
        r3 = selmer_3isog(a, b)
        rec['isog3'] = r3
    except Exception as e:
        rec['isog3_error'] = str(e)
    r2 = two_isog_bounds(a, b)
    rec['isog2'] = [dict(bound=x['bound'], nS=len(x['S']), nT=len(x['T']), S=x['S'], T=x['T']) for x in r2]
    M = b**3*(2*a+b); N = a**3*(a+2*b)
    E = EllipticCurve(QQ, [0, M+N, 0, M*N, 0])
    if not NOPARI:
        try:
            rk = pari(E).ellrank()
            rec['pari'] = [int(rk[0]), int(rk[1])]
        except Exception as e:
            rec['pari_error'] = str(e)
    rec['tors'] = int(E.torsion_order())
    out.write(json.dumps(rec) + '\n'); out.flush()
    if idx % 10 == 0:
        print(f"shard {shard}: {idx+1}/{len(mine)} ({a},{b}) t={time.time()-t0:.0f}s", flush=True)
print(f"shard {shard} done, t={time.time()-t0:.0f}s")
