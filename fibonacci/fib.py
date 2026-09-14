# Claude, 14.09: наклоны Fₙ/Fₙ₊₁ за пределами переписи. Для каждого — тот же поиск, что в переписи (c300.work:
# род 1 из 3–4 клеток, шесть клеток; ellrank), в отдельном процессе с тайм-аутом; затем L(E,1) для кандидатов.
import sys, json, subprocess, concurrent.futures as cf
F=[1,1]
while len(F)<70: F.append(F[-1]+F[-2])
NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 45
items=[(n,f"{F[n]}/{F[n+1]}") for n in range(3,NMAX)]
one="""import sys,json
sys.path.insert(0,'/home/kep/magicKube/census_300')
from c300 import work
print(json.dumps(work(sys.argv[1])),flush=True)"""
def run(it):
    n,sl=it
    try:
        p=subprocess.run([sys.executable,'-c',one,sl],capture_output=True,text=True,timeout=1500)
        got=[json.loads(l) for l in p.stdout.splitlines() if l.startswith('{')]
        rec=got[0] if got else dict(slope=sl,closed=False,crash=p.stderr[-200:])
    except subprocess.TimeoutExpired:
        rec=dict(slope=sl,closed=False,timeout=True)
    rec['n']=n; return rec
out=open('fib.jsonl','w')
with cf.ThreadPoolExecutor(8) as ex:
    for rec in ex.map(run,items):
        out.write(json.dumps(rec)+"\n"); out.flush()
        tag='ЗАКРЫТ '+str(rec.get('by')) if rec.get('closed') else ('ALERT' if 'ALERT' in rec else ('тайм-аут' if rec.get('timeout') else 'не закрыт'))
        print(rec['n'],rec['slope'],tag,flush=True)
print('готово',flush=True)
