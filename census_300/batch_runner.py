# добивает перепись пакетами в отдельных процессах; упавший пакет дробится до одиночных; одиночный segfault записывается как crash
import sys, json, subprocess, os
from math import gcd
SMIN,SMAX=int(sys.argv[1]),int(sys.argv[2]); fn=f'c{SMIN}_{SMAX}.jsonl'
allsl=[f"{r}/{s}" for s in range(SMIN,SMAX+1) for r in range(1,s) if gcd(r,s)==1]
done={json.loads(l)['slope'] for l in open(fn)} if os.path.exists(fn) else set()
todo=[t for t in allsl if t not in done]; print('осталось',len(todo),flush=True)
one="""import sys,json
sys.path.insert(0,'/home/kep/magicKube/census_300')
from c300 import work
for sl in sys.argv[1:]:
    print(json.dumps(work(sl)),flush=True)
"""
def run(batch):
    p=subprocess.run([sys.executable,'-c',one]+batch,capture_output=True,text=True,timeout=3600)
    got=[json.loads(l) for l in p.stdout.splitlines() if l.startswith('{')]
    with open(fn,'a') as f:
        for g in got: f.write(json.dumps(g)+"\n")
    rest=batch[len(got):]
    if rest:
        if len(rest)==1:
            with open(fn,'a') as f: f.write(json.dumps(dict(slope=rest[0],closed=False,crash=True))+"\n")
            print('crash',rest[0],flush=True)
        else:
            run(rest[:1]); run(rest[1:]) if len(rest)>1 else None
import concurrent.futures as cf
B=[todo[i:i+25] for i in range(0,len(todo),25)]
with cf.ThreadPoolExecutor(10) as ex:
    list(ex.map(run,B))
print('готово',flush=True)
