from sage.all import *
import json
from kolyvagin import curve
R=[json.loads(l) for l in open('kolyvagin.jsonl')]
sqf=0; nsq=[]
for r in R:
    if r['closed']:
        b=r['by']; E,lc=curve(*b['abc'],b['kind']); N=E.minimal_model().conductor()
        if N.is_squarefree(): sqf+=1
        else: nsq.append((r['slope'],b['abc'],b['kind'],str(N.factor()),b['L_ratio']))
print('кондуктор свободен от квадратов:',sqf,'; нет:',len(nsq))
for x in nsq: print(x)
print('ошибки:',[(r['slope'],t['abc'],t['kind'],t['err']) for r in R for t in r['tried'] if 'err' in t])
