# Claude, 26.09: строгое закрытие 21 наклона s ≤ 200, у которых единственное основание — один прогон Magma.
# Для наклона: все 3-клеточные T ⊂ {±r,±s,±(s−r),±(s+r)}, T ≠ −T, ранг множителя 1 (PARI ellrank [1,1]) → rigdem.check_factor
# (строгая B из сертификата Безу для удвоения). Останавливаемся на первом закрывающем T; сохраняем до 3 закрывающих.
import sys, json, itertools, time
sys.path.insert(0, '/home/kep/magicKube/rigorous_census/magma21')
import rigdem
from sage.all import QQ, pari
slope = sys.argv[1]; r, s = map(int, slope.split('/'))
vals = [r, s - r, s, s + r]; C = sorted(set(vals + [-v for v in vals]))
cands = []
for T in itertools.combinations(C, 3):
    if sorted(-t for t in T) == sorted(T): continue
    cands.append(list(T))
out = dict(slope=slope, tried=[], closed=False)
nclosed = 0; t0 = time.time()
for T in cands:
    if time.time() - t0 > 3000 or nclosed >= 3: break
    try:
        res = rigdem.check_factor(r, s, T)
    except Exception as ex:
        res = dict(T=T, error=str(ex)[:300])
    short = {k: res.get(k) for k in ('ellrank', 'M0', 'B_rig', 'nondegenerate', 'closed', 'skip', 'error')}
    short['T'] = T
    out['tried'].append(short)
    if res.get('nondegenerate'): out['ALERT'] = res['nondegenerate']
    if res.get('closed'): nclosed += 1; out['closed'] = True
    print(slope, T, {k: short[k] for k in ('ellrank', 'M0', 'closed', 'skip')}, flush=True)
out['n_closing'] = nclosed; out['sec'] = round(time.time() - t0)
json.dump(out, open('m21_%d_%d.json' % (r, s), 'w'), ensure_ascii=False, indent=1, default=str)
print('ИТОГ', slope, 'closed =', out['closed'], 'закрывающих T:', nclosed, 'ALERT:', out.get('ALERT'), flush=True)
