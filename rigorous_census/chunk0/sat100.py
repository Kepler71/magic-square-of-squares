# Claude, 26.09.2026. Доп. проверка образующих: PARI ellsaturation(E,[P0],100) — p-насыщение для всех p <= 100
# (алгебраическое, по редукциям), независимо от численной оценки индекса. Если P0 = mQ + t, m > 1, то у m есть
# простой делитель p; p <= 100 исключено здесь, p > 100 требует hhat(Q) <= hhat(P0)/101^2.
from sage.all import *
import json, glob, time
base = '/home/kep/magicKube/rigorous_census/chunk0'
out = {}; t0 = time.time()
for fn in sorted(glob.glob(f'{base}/rig_*_*.json')):
    d = json.load(open(fn))
    for o in d['results']:
        if not o.get('closed_rigorous'):
            continue
        E = EllipticCurve([QQ(a) for a in o['Emin']])
        P0 = E([QQ(v) for v in o['P0']])
        V = pari.ellsaturation(pari(E), [pari([P0[0], P0[1]])], 100)
        Q = E(list(map(QQ, V[0])))
        same = (Q - P0).order() != oo or (Q + P0).order() != oo
        need = float(P0.height()) / 101**2
        out.setdefault(d['slope'], []).append(dict(T=o['T'], sat100_same=bool(same), lam_sage=o['lambda_sage'],
                                                   hQ_would_be_le=need, lam_over_need=o['lambda_sage'] / need))
        print(d['slope'], o['T'], 'p<=100 насыщена:', same, 'λ_Sage/нужное = %.0f' % (o['lambda_sage'] / need), flush=True)
json.dump(out, open(f'{base}/sat100.json', 'w'), ensure_ascii=False, indent=1)
allv = [x for v in out.values() for x in v]
print('итог: насыщено при p<=100:', sum(x['sat100_same'] for x in allv), '/', len(allv),
      ' min λ_Sage/нужное = %.0f' % min(x['lam_over_need'] for x in allv), ' %.0fs' % (time.time() - t0), flush=True)
