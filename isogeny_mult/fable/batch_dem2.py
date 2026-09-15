# Fable 15.09. Пакетный прогон demjanenko2 по списку наклонов: для каждого — классы с доказанным рангом == RANK
# (по переписи ranks_*.jsonl / all54), первые NCL таких классов; отчёт в batch_<name>.jsonl.
from sage.all import *
import sys, json, os, time
sys.path.insert(0, '/home/kep/magicKube/isogeny_mult/fable')
from demjanenko2 import run_class, CLASSES, BIG, cells_of, slope_cells
D = '/home/kep/magicKube/isogeny_mult/fable'
lst = [l.strip() for l in open(sys.argv[1]) if l.strip()]; name = sys.argv[2]; RANK = int(sys.argv[3]); NCL = int(sys.argv[4]) if len(sys.argv) > 4 else 1
rows = {}
for f in sorted(os.listdir(D)):
    if f.startswith('ranks_') and f.endswith('.jsonl'):
        for l in open(f'{D}/{f}'):
            d = json.loads(l); rows[d['slope']] = d['classes']
fn = f'{D}/batch_{name}.jsonl'
done = set(json.loads(l)['slope'] for l in open(fn)) if os.path.exists(fn) else set()
for sl in lst:
    if sl in done: continue
    r, s = map(int, sl.split('/')); S = slope_cells(r, s); t0 = time.time()
    cl = rows.get(sl)
    if cl is None: print(sl, 'нет в переписи рангов'); continue
    cis = [i for i, (m, T, lo, hi) in enumerate(cl) if lo is not None and lo == hi == RANK][:NCL]
    reps = []
    for ci in cis:
        Ts = [cells_of(Sg, r, s) for Sg in CLASSES[BIG[ci]]]
        try:
            rep = run_class(S, Ts, verbose=False)
        except Exception as e:
            rep = dict(ok=False, why='ошибка: ' + repr(e)[:150])
        rep['class'] = ci; reps.append(rep)
    closed = any(o.get('ok') and not o['nondeg_solutions'] for o in reps)
    sols = [o['nondeg_solutions'] for o in reps if o.get('ok')]
    out = dict(slope=sl, rank=RANK, classes=cis, closed=closed, reps=reps, time=time.time() - t0)
    with open(fn, 'a') as f: f.write(json.dumps(out, default=str) + '\n')
    print(sl, 'классы', cis, 'ЗАКРЫТ:', closed, 'решения:', sols, [(o.get('H0'), o.get('ncand'), o.get('why')) for o in reps], f'{time.time()-t0:.1f} с', flush=True)
