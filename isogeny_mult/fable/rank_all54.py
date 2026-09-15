# Fable 15.09. Ранги всех 54 классов (представителей) для списка наклонов (файл, по строке r/s) → all54_<name>.jsonl
from sage.all import *
import sys, json, subprocess, os, time
from classes import CLASSES, curve_of
D = '/home/kep/magicKube/isogeny_mult/fable'
ORDER = sorted(CLASSES, key=lambda c: (-len(CLASSES[c]), c))


def work(sl):
    r, s = map(int, sl.split('/')); out = []
    for c in ORDER:
        S = CLASSES[c][0]
        try:
            E, L, c0, T = curve_of(S, r, s)
            rk = E.pari_curve().ellrank(); out.append([len(CLASSES[c]), T, int(rk[0]), int(rk[1])])
        except Exception as e:
            out.append([len(CLASSES[c]), None, None, str(e)[:50]])
    return dict(slope=sl, classes=out)


if __name__ == '__main__':
    if sys.argv[1] == '--worker':
        fn = sys.argv[2]
        for sl in sys.argv[3:]:
            res = work(sl)
            with open(fn, 'a') as f: f.write(json.dumps(res) + '\n')
        sys.exit(0)
    lst = [l.strip() for l in open(sys.argv[1]) if l.strip()]; name = sys.argv[2]; nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    fn = f'{D}/all54_{name}.jsonl'
    done = set(json.loads(l)['slope'] for l in open(fn)) if os.path.exists(fn) else set()
    todo = [x for x in lst if x not in done]; print('к счёту', len(todo), flush=True)
    chunk = 5; queue = [todo[i:i + chunk] for i in range(0, len(todo), chunk)]; running = []; t0 = time.time()
    while queue or running:
        while queue and len(running) < nproc:
            lst_ = queue.pop(0)
            running.append((subprocess.Popen([sys.executable, f'{D}/rank_all54.py', '--worker', fn] + lst_, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL), lst_, time.time()))
        time.sleep(0.5)
        for it in running[:]:
            p, lst_, ts = it; rc = p.poll()
            if rc is None and time.time() - ts > 1200: p.kill(); rc = -9
            if rc is None: continue
            running.remove(it)
            if rc != 0 and len(lst_) > 1:
                got = set(json.loads(l)['slope'] for l in open(fn)) if os.path.exists(fn) else set()
                for sl in lst_:
                    if sl not in got: queue.insert(0, [sl])
            elif rc != 0: open(f'{D}/all54_crash.txt', 'a').write(lst_[0] + '\n')
    print('готово', f'{time.time()-t0:.0f} с', flush=True)
