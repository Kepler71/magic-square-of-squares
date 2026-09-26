# Claude, 26.09.2026 (сводка STATUS). Своя сверка census_500_master.json: полнота множества наклонов,
# классы основных оснований, наклоны с единственным основанием «Magma», покрытие 46 трудных наклонов строгой B,
# текущее состояние rig_*.json во всех трёх чанках (мастер собран в 20:20, чанк 0 перезаписан в 20:40).
import json, glob, collections
from math import gcd
R = '/home/kep/magicKube/rigorous_census'
d = json.load(open(f'{R}/census_500_master.json'))['slopes']
exp = {f'{r}/{s}' for s in range(3, 501) for r in range(1, s) if gcd(r, s) == 1}
print('ожидается', len(exp), 'в мастере', len(d), 'множества равны:', set(d) == exp)
print('основные основания по классам:', dict(collections.Counter(v['m'] for v in d.values())))
print('по методам:', dict(collections.Counter(v['m'] + ':' + v['b'] for v in d.values())))
print('основное с сомнением:', [k for k, v in d.items() if v.get('doubt')])
print('основное класса Е (численно):', [k for k, v in d.items() if v['m'] == 'Е'])
only_magma = [k for k, v in d.items() if v['b'] == 'g2_magma'
              and not any(a['m'] in 'АБД' or a['b'] == 'dem_rank1_strictB' for a in v.get('alt', []))]
print('единственное основание — Magma род 2:', len(only_magma), only_magma)
cur = {}
for ch in ('chunk0', 'chunk1', 'chunk2'):
    for fn in glob.glob(f'{R}/{ch}/rig_*_*.json'):
        x = json.load(open(fn))
        cur[x['slope']] = (ch, x.get('closed_rigorous', x.get('closed_rigorously')), x.get('ALERT', x.get('alert')))
print('rig_*.json сейчас:', len(cur), 'closed=True и без ALERT:', sum(1 for v in cur.values() if v[1] is True and not v[2]))
hard = [k for k, v in d.items() if v['m'] in 'ВГД' and int(k.split('/')[1]) > 200]
print('s>200 с основным В/Г/Д:', len(hard), '; из них покрыты строгой B (сейчас):', sum(1 for k in hard if k in cur and cur[k][1] is True))
print('s>200 с основным В/Г/Д без строгой B:', [k for k in hard if k not in cur])
