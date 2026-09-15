# Fable 15.09. Устойчивый прогон rank_classes: подпроцессы по кускам (PARI ellrank иногда падает — cysignals crash);
# упавший кусок пересчитывается по одному наклону; наклон, роняющий PARI, записывается с пометкой crash.
# python3 rank_driver.py SMIN SMAX [nproc] [chunk]   → ranks_SMIN_SMAX.jsonl (+ ranks_crash.txt)
import sys, os, json, subprocess, time
from math import gcd
smin, smax = int(sys.argv[1]), int(sys.argv[2]); nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
chunk = int(sys.argv[4]) if len(sys.argv) > 4 else 40
D = '/home/kep/magicKube/isogeny_mult/fable'
fn = f'{D}/ranks_{smin}_{smax}.jsonl'
done = set()
for f in os.listdir(D):
    if f.startswith('ranks_') and f.endswith('.jsonl'):
        for l in open(f'{D}/{f}'):
            try: done.add(json.loads(l)['slope'])
            except Exception: pass
todo = [f'{r}/{s}' for s in range(smin, smax + 1) for r in range(1, s) if gcd(r, s) == 1 and 2 * r != s and f'{r}/{s}' not in done]
print('к счёту:', len(todo), flush=True)
WORKER = f'''
import sys, json
sys.path.insert(0, "{D}")
from rank_classes import work
for sl in sys.argv[2:]:
    r, s = map(int, sl.split("/")); res = work((r, s))
    with open(sys.argv[1], "a") as f: f.write(json.dumps(res) + "\\n")
'''
open(f'{D}/_rank_worker.py', 'w').write(WORKER)
queue = [todo[i:i + chunk] for i in range(0, len(todo), chunk)]
running = []; t0 = time.time(); ndone = 0; crashes = []
def launch(lst, single):
    p = subprocess.Popen([sys.executable, f'{D}/_rank_worker.py', fn] + lst, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    running.append((p, lst, single, time.time()))
while queue or running:
    while queue and len(running) < nproc:
        lst = queue.pop(0); launch(lst, len(lst) == 1)
    time.sleep(0.5)
    for item in running[:]:
        p, lst, single, ts = item
        rc = p.poll()
        if rc is None:
            if time.time() - ts > 900: p.kill(); rc = -9
            else: continue
        running.remove(item)
        if rc == 0: ndone += len(lst)
        else:
            if single: crashes.append(lst[0]); open(f'{D}/ranks_crash.txt', 'a').write(lst[0] + f' rc={rc}\n')
            else:
                # пересчитать по одному те, что не записаны
                got = set()
                for l in open(fn) if os.path.exists(fn) else []:
                    try: got.add(json.loads(l)['slope'])
                    except Exception: pass
                for sl in lst:
                    if sl not in got: queue.insert(0, [sl])
                    else: ndone += 1
        if ndone and ndone % 400 < len(lst): print(ndone, f'{time.time()-t0:.0f} с', flush=True)
print('готово', ndone, 'падений', len(crashes), f'{time.time()-t0:.0f} с', flush=True)
