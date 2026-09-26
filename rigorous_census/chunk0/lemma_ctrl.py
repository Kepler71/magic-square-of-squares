# Claude, 26.09.2026. Контроль леммы Безу на одном шаге удвоения — в целых числах, без логарифмов.
# Для x = m/n (НОД 1), H(x) = max(|m|,|n|): лемма утверждает  H(x)^4 <= K * H(x(2Q))  и  H(x(2Q)) <= L * H(x)^4.
# Точки: Q = n P0 + t, 1 <= n <= NMAX, t in Tors (x(2Q) — через групповой закон Sage, т.е. независимо от F/G).
# Также: насколько лемма близка к равенству (запас в натуральных логарифмах, численно — только для информации).
from sage.all import *
import json, glob, time
base = '/home/kep/magicKube/rigorous_census/chunk0'
NMAX = 15
RR2 = RealField(200)


def H(q):
    q = QQ(q); return max(abs(q.numerator()), q.denominator())


t0 = time.time(); tot = 0; fail = 0; slack_lo = []; slack_hi = []
for fn in sorted(glob.glob(f'{base}/rig_*_*.json')):
    d = json.load(open(fn))
    for o in d['results']:
        if not o.get('closed_rigorous'):
            continue
        E = EllipticCurve([QQ(a) for a in o['Emin']])
        P0 = E([QQ(v) for v in o['P0']])
        K = ZZ(o['cert']['K']); L = ZZ(o['cert']['L'])
        tors = E.torsion_points()
        Q = E(0); smin_lo = None; smin_hi = None
        for n in range(1, NMAX + 1):
            Q = Q + P0
            for t in tors:
                Qt = Q + t
                h1 = H(Qt[0]); h2 = H((2 * Qt)[0])
                tot += 1
                if not (h1**4 <= K * h2 and h2 <= L * h1**4):
                    fail += 1; print('НАРУШЕНИЕ', d['slope'], o['T'], n, t, flush=True)
                s_lo = float(RR2(K * h2).log() - 4 * RR2(h1).log()); s_hi = float(RR2(L).log() + 4 * RR2(h1).log() - RR2(h2).log())
                smin_lo = s_lo if smin_lo is None else min(smin_lo, s_lo)
                smin_hi = s_hi if smin_hi is None else min(smin_hi, s_hi)
        slack_lo.append(smin_lo / float(RR2(K).log())); slack_hi.append(smin_hi / float(RR2(L).log()))
    print(f"[{time.time()-t0:.0f}s] {d['slope']}: точек всего {tot}, нарушений {fail}", flush=True)
print(f'итог: {tot} проверок, нарушений {fail}; min запас нижней стороны / log K = {min(slack_lo):.3f}, '
      f'верхней / log L = {min(slack_hi):.3f}', flush=True)
