# ранги E+ и E- для всех троек клеток трёх открытых наклонов (кандидаты биэллиптического квадратичного Шаботи: ранги 1 и 1)
from sage.all import *
import itertools, json, sys
sys.path.insert(0,'/home/kep/magicKube/census_300')
from c300 import e_pm
out=[]
for sl in sys.argv[1:] or ['265/298','73/362','126/451']:
    r,s=map(int,sl.split('/')); vals=sorted(set([r,s-r,s,s+r]))
    for a,b,c in itertools.combinations(vals,3):
        row=dict(slope=sl,abc=[a,b,c])
        for kind in '+-':
            E,lc=e_pm(a,b,c,kind)
            lo,hi=[int(t) for t in pari(E).ellrank(4)[:2]]
            row[kind]=[lo,hi,int(E.root_number())]
        row['QC_candidate']= row['+'][:2]==[1,1] and row['-'][:2]==[1,1]
        out.append(row); print(row,flush=True)
json.dump(out,open('ranks_open.json' if len(sys.argv)<2 else 'ranks_'+sys.argv[1].replace('/','_')+'.json','w'),indent=1)
