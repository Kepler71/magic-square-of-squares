# Fable 15.09. Ранги (PARI ellrank: нижняя/верхняя) представителей 11 больших классов для наклонов переписи.
# python3 rank_classes.py SMIN SMAX [nproc]  → ranks_SMIN_SMAX.jsonl
from sage.all import *
import sys, json, time, os
import multiprocessing as mp
from classes import CLASSES, BIG, curve_of


def work(sl):
    r, s = sl
    out = []
    for c in BIG:
        S = CLASSES[c][0]
        try:
            E, L, c0, T = curve_of(S, r, s)
            if E.discriminant() == 0: out.append([len(CLASSES[c]), T, None, None]); continue
            rk = E.pari_curve().ellrank()
            out.append([len(CLASSES[c]), T, int(rk[0]), int(rk[1])])
        except Exception as e:
            out.append([len(CLASSES[c]), None, None, str(e)[:60]])
    return dict(slope=f'{r}/{s}', classes=out)


if __name__ == '__main__':
    smin, smax = int(sys.argv[1]), int(sys.argv[2]); nproc = int(sys.argv[3]) if len(sys.argv) > 3 else 8
    slopes = [(r, s) for s in range(smin, smax + 1) for r in range(1, s) if gcd(r, s) == 1 and 2 * r != s]
    fn = f'/home/kep/magicKube/isogeny_mult/fable/ranks_{smin}_{smax}.jsonl'
    done = set()
    if os.path.exists(fn):
        done = set(json.loads(l)['slope'] for l in open(fn))
    slopes = [x for x in slopes if f'{x[0]}/{x[1]}' not in done]
    print('наклонов к счёту:', len(slopes), flush=True)
    t0 = time.time(); n = 0
    with open(fn, 'a') as f, mp.get_context('fork').Pool(nproc) as pool:
        for res in pool.imap_unordered(work, slopes, chunksize=4):
            f.write(json.dumps(res) + '\n'); f.flush(); n += 1
            if n % 200 == 0: print(n, f'{time.time()-t0:.0f} с', flush=True)
    print('готово', n, f'{time.time()-t0:.0f} с', flush=True)
