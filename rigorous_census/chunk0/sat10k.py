# Claude, 26.09.2026. Усиление sat100.py: PARI ellsaturation(E,[P0],10^4) — p-насыщение для всех p <= 10^4
# (алгебраически, через редукции по модулю простых). Если P0 = mQ + t, m > 1, то у m есть простой p | m;
# p <= 10^4 исключено здесь, а p > 10^4 требует hhat(Q) = hhat(P0)/m^2 <= hhat(P0)/10007^2.
# Сравнение с lambda_Sage (численная нижняя граница CPS для hhat на E(Q) \ Tors) — остающаяся зависимость (ПО).
from sage.all import *
import json, glob, time, multiprocessing as mp
base = '/home/kep/magicKube/rigorous_census/chunk0'
PB = 10**4
PNEXT = next_prime(PB)


def one(arg):
    sl, o = arg
    E = EllipticCurve([QQ(a) for a in o['Emin']])
    P0 = E([QQ(v) for v in o['P0']])
    t = time.time()
    V = pari.ellsaturation(pari(E), [pari([P0[0], P0[1]])], PB)
    Q = E(list(map(QQ, V[0])))
    same = (Q - P0).order() != oo or (Q + P0).order() != oo
    need_hi = o['hP0_hi'] / PNEXT**2          # верхний конец строгой вилки hhat(P0), делённый на p^2
    return dict(slope=sl, T=o['T'], sat_p_le_1e4=bool(same), lam_sage=o['lambda_sage'],
                hQ_would_be_le=float(need_hi), lam_over_need=float(o['lambda_sage'] / need_hi), sec=float(round(time.time() - t, 2)))


if __name__ == '__main__':
    t0 = time.time()
    tasks = []
    for fn in sorted(glob.glob(f'{base}/rig_*_*.json')):
        d = json.load(open(fn))
        tasks += [(d['slope'], o) for o in d['results'] if o.get('closed_rigorous')]
    print('кривых', len(tasks), flush=True)
    out = []
    with mp.get_context('fork').Pool(3) as pool:
        for i, v in enumerate(pool.imap_unordered(one, tasks)):
            out.append(v)
            print(f"[{i+1}/{len(tasks)} {time.time()-t0:.0f}s] {v['slope']} {v['T']} p<=1e4 насыщена: {v['sat_p_le_1e4']} "
                  f"λ_Sage/нужное = {v['lam_over_need']:.3g}", flush=True)
    json.dump(out, open(f'{base}/sat10k.json', 'w'), ensure_ascii=False, indent=1)
    print('итог: насыщено при p<=1e4:', sum(v['sat_p_le_1e4'] for v in out), '/', len(out),
          ' min λ_Sage/нужное = %.3g' % min(v['lam_over_need'] for v in out), ' %.0fs' % (time.time() - t0), flush=True)
