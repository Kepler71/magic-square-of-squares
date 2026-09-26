# Claude, 26.09.2026. Попытка сквозного положительного контроля: наборы S'' = ±T с ИЗВЕСТНЫМ решением z* ≠ 0
# (все шесть клеток 1 ± λz* — квадраты). λ = a² − 1 при a² + b² = 2 ⇒ 1 + λ = a², 1 − λ = b² (z = 1),
# после приведения к целым λ' = λ·D, z* = 1/D. Ищем T с рангом E_T = 1 (PARI ellrank [1,1]) и прогоняем
# check_factor с клетками S'' — метод обязан выдать z* как «невырожденное решение».
import sys, os, json, time, itertools
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import rigdem
from rigdem import *

rigdem.LOGF = open(os.path.join(here, 'ctrl_positive.log'), 'a')
MAXT = int(sys.argv[1]) if len(sys.argv) > 1 else 400
TLIM = float(sys.argv[2]) if len(sys.argv) > 2 else 600
t0 = time.time()
lams = set()
for num in range(-12, 13):
    for den in range(1, 9):
        if gcd(num, den) != 1:
            continue
        t = QQ(num) / den
        a = (t ** 2 - 2 * t - 1) / (1 + t ** 2)
        lam = a ** 2 - 1
        if lam != 0 and abs(lam) != 1:
            lams.add(lam)
lams = sorted(lams, key=lambda q: (max(abs(q.numerator()), q.denominator()), q))
log('ctrl_positive: λ с конического семейства: %d; берём первые 14 по высоте' % len(lams))
lams = lams[:14]
out = []; tried = 0; rank1 = 0
for trip in itertools.combinations(lams, 3):
    if time.time() - t0 > TLIM or tried >= MAXT:
        break
    if len(set(abs(x) for x in trip)) < 3:
        continue
    D = lcm([x.denominator() for x in trip])
    T = sorted(ZZ(x * D) for x in trip)
    g = gcd(T)
    # z* = 1/D; делить T на g нельзя без изменения z*: z* = g/D при T/g
    Tg = [x // g for x in T]; zs = QQ(g) / D
    S2 = sorted(set(Tg + [-x for x in Tg]))
    assert full_ok(zs, S2) and full_ok(-zs, S2)
    tried += 1
    md = model(Tg)
    try:
        er = pari(md['E1'].minimal_model()).ellrank()
    except Exception:
        continue
    rk = [int(er[0]), int(er[1])]
    if rk != [1, 1]:
        continue
    rank1 += 1
    rigdem.cells = lambda r, s, S2=S2: S2
    res = check_factor(0, 0, Tg, ctrl_brute=0)
    found = [v for tag, v in res.get('nondegenerate', [])]
    ok = str(zs) in found or str(-zs) in found
    log('  T=%s z*=%s: ранг 1; метод нашёл %s; z* найден: %s; M0=%s' % (Tg, zs, found, ok, res.get('M0')))
    out.append(dict(T=[int(x) for x in Tg], zstar=str(zs), found=found, detected=ok, M0=res.get('M0'),
                    tags=res.get('nondegenerate')))
log('ctrl_positive: перебрано %d троек, ранг 1 у %d, за %.0f с' % (tried, rank1, time.time() - t0))
json.dump(dict(tried=tried, rank1=rank1, cases=out), open(os.path.join(here, 'ctrl_positive.json'), 'w'), ensure_ascii=False, indent=1)
