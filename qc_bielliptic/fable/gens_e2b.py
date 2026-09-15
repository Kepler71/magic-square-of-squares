# точка бесконечного порядка на E2min через изогенный класс: mwrank на каждой кривой класса (быстро), перенос по изогениям
import sys,time,json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f); Em=E2.minimal_model(); print('E2min',Em.ainvs(),flush=True)
cl=Em.isogeny_class(); curves=list(cl.curves); M=cl.isogenies(); i0=curves.index(Em); print('класс изогений: кривых',len(curves),'матрица степеней',cl.matrix(),flush=True)
# BFS путей к Em
import collections
paths={i0:[]}; dq=collections.deque([i0])
while dq:
    i=dq.popleft()
    for j in range(len(curves)):
        if j not in paths and M[j][i]!=0:
            paths[j]=[M[j][i]]+paths[i]; dq.append(j)
found=None
for j,Ej in enumerate(curves):
    t0=time.time()
    try:
        alarm(240); g=list(Ej.gens(proof=False)); cancel_alarm()
    except Exception as e:
        cancel_alarm(); g=[]
    print(' кривая',j,Ej.ainvs()[-2:],'gens',g,round(time.time()-t0,1),'с',flush=True)
    if g:
        P=g[0]
        for phi in paths[j]: P=phi(P)
        assert P.curve()==Em and P.order()==oo
        print('ТОЧКА на E2min:',P,'высота',float(P.height()),flush=True); found=P; break
out=json.load(open('gens_open.json'))
if found is not None:
    out[f'{a},{b},{c}']['E2']=[str(found)]; out[f'{a},{b},{c}']['E2_how']='isogeny class + mwrank'
    json.dump(out,open('gens_open.json','w'),indent=1)
