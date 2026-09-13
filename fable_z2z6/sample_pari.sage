import json, glob, random
random.seed(int(2026))
recs = []
for f in sorted(glob.glob('/home/kep/magicKube/fable_z2z6/grid/shard_200_*.jsonl')):
    for line in open(f):
        r = json.loads(line)
        if r['b'] > 60: recs.append(r)
sample = random.sample(recs, 300)
viol = 0; n = 0; hi0 = 0; both0 = 0
out = open('/home/kep/magicKube/fable_z2z6/sample_pari.jsonl', 'w')
for r in sample:
    a, b = r['a'], r['b']
    M = b**3*(2*a+b); N = a**3*(a+2*b)
    E = EllipticCurve(QQ, [0, M+N, 0, M*N, 0])
    try:
        alarm(300); rk = pari(E).ellrank(); cancel_alarm()
    except Exception as e:
        cancel_alarm(); out.write(json.dumps(dict(a=a,b=b,error=str(e)[:80]))+'\n'); continue
    lo, hi = int(rk[0]), int(rk[1])
    b3 = r['isog3']['bound']; b2 = min(x['bound'] for x in r['isog2'])
    n += 1
    if b3 < lo or b2 < lo: viol += 1; print("VIOLATION", a, b, b3, b2, lo, hi)
    hi0 += (hi == 0); both0 += (min(b3, b2) == 0)
    out.write(json.dumps(dict(a=a,b=b,b3=b3,b2=b2,pari=[lo,hi]))+'\n'); out.flush()
print(f"sample {n}: violations {viol}; PARI hi=0: {hi0}; ours min=0: {both0}")
