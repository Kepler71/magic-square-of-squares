# Проверяющий: свой анализ (ранг <w_u>, лемма, классы D, шаги Z') на шести клетках Саллоуса. Использует функции my_configs.py.
import sys, json
src = open('/home/kep/magicKube/joint_Eb_Ec/center4/check/my_configs.py').read()
head = src.split('t0 = time.time()')[0]
sys.argv = ['x', '1000']
exec(head)
res = {}
for a, tr in ((1105, (264, 576, 1104)), (4420, (1056, 2304, 4416))):
    rels = relations(tr)
    o = analyse(a, tr, rels)
    o['maxr'] = str(o['maxr'])
    res[a] = o
    print(a, tr, 'rels', rels, 'rank', o['rank'], 'lemma', o['lemma_ok'], 'D', o['D_ok'], 'steps', o['steps_ok'])
    for c in o['classes']: print('   ', c)
# прогноз: какие объединения простых 5,13,17 допустимы как класс при NL=4, a=1105
from itertools import combinations
for r in (1, 2, 3):
    for S in combinations((5, 13, 17), r):
        Pi = 1
        for p in S: Pi *= p
        print('класс', S, 'Pi^2 =', Pi*Pi, '<= 4a=4420 ?', Pi*Pi <= 4420)
json.dump(res, open('my_sallows.json', 'w'), indent=1, default=str)
