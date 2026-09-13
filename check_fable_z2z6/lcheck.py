# Claude, 14.09: независимая проверка теоремы 1 Fable. Формулу Fable беру только как «заявление»;
# проверяю аналитикой: граница 0 ⇒ L(E,1)≠0 (L_ratio); L(E,1)=0 ⇒ граница ≥1. Плюс (50,59).
import sys, json
sys.path.insert(0,'/home/kep/magicKube/fable_z2z6')
from explicit_formula import explicit_dims
from sage.all import *
import multiprocessing as mp
from math import gcd
import random as rnd
def one(ab):
    a,b=ab; M=b**3*(2*a+b); N=a**3*(a+2*b)
    E=EllipticCurve([0,M+N,0,M*N,0]).minimal_model()
    bound=explicit_dims(a,b)[2]
    Lr=E.lseries().L_ratio()
    return (a,b,int(bound),str(Lr), int(E.root_number()))
if __name__=='__main__':
    rnd.seed(14092026)
    allp=[(a,b) for b in range(2,201) for a in range(1,b) if gcd(a,b)==1]
    z=[p for p in allp if explicit_dims(*p)[2]==0]
    nz=[p for p in allp if explicit_dims(*p)[2]>0]
    sample=[(50,59)]+rnd.sample(z,150)+rnd.sample(nz,150)
    viol=[]; rows=[]
    with mp.get_context('fork').Pool(12) as pool:
        for a,b,bound,Lr,w in pool.imap_unordered(one,sample):
            rows.append((a,b,bound,Lr,w))
            if bound==0 and Lr=='0': viol.append((a,b,'граница 0, но L(E,1)=0'))
            if bound<0: viol.append((a,b,'отрицательная граница'))
            if w==-1 and bound%2==0 and bound==0: viol.append((a,b,'граница 0 при корневом числе −1'))
    print('выборка',len(rows),'; граница 0:',sum(1 for r in rows if r[2]==0),'; из них L≠0:',sum(1 for r in rows if r[2]==0 and r[3]!='0'))
    print('граница>0 и L≠0 (ранг 0, но теорема не видит):',sum(1 for r in rows if r[2]>0 and r[3]!='0'))
    print('(50,59):',[r for r in rows if r[:2]==(50,59)])
    print('нарушений:',viol)
    json.dump(dict(rows=rows,viol=viol),open('lcheck.json','w'))
