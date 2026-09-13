# Контроль: для пар ранга 0 (PARI hi=0), b<=40: разрыв (3-изогенная оценка - 0) > 0  <=>  3 | Sha_an(E)*Sha_an(E')
import json, glob, sys
load('/home/kep/magicKube/fable_z2z6/three_isog.sage')
recs = []
for f in sorted(glob.glob('/home/kep/magicKube/fable_z2z6/grid/shard_60_*.jsonl')):
    for line in open(f): recs.append(json.loads(line))
out = open('/home/kep/magicKube/fable_z2z6/sha_control.jsonl', 'w')
n = 0; ok = 0; bad = []
for r in recs:
    if r['b'] > 40 or r['pari'][1] != 0: continue
    a, b = r['a'], r['b']
    E, Ep, phi, al, be = curve_data(a, b)
    Em = E.minimal_model(); Epm = Ep.minimal_model()
    try:
        alarm(120)
        sE = Em.sha().an(); sEp = Epm.sha().an()
        cancel_alarm()
    except Exception as e:
        cancel_alarm(); out.write(json.dumps(dict(a=a,b=b,error=str(e)[:100]))+'\n'); out.flush(); continue
    gap = r['isog3']['bound']
    pred = (gap > 0)
    obs = ((ZZ(sE)*ZZ(sEp)) % 3 == 0)
    n += 1; ok += (pred == obs)
    if pred != obs: bad.append((a,b,gap,sE,sEp))
    out.write(json.dumps(dict(a=a,b=b,gap=gap,shaE=int(sE),shaEp=int(sEp),ok=(pred==obs)))+'\n'); out.flush()
print("checked", n, "consistent", ok, "inconsistent", bad)
