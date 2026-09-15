# вторая попытка для E2: PARI ellrank(effort=4) на каждой кривой изогенного класса (alarm 280 с на кривую), перенос по изогениям
import sys,time,json,collections; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f); Em=E2.minimal_model()
cl=Em.isogeny_class(); curves=list(cl.curves); M=cl.isogenies(); i0=curves.index(Em)
paths={i0:[]}; dq=collections.deque([i0])
while dq:
    i=dq.popleft()
    for j in range(len(curves)):
        if j not in paths and M[j][i]!=0: paths[j]=[M[j][i]]+paths[i]; dq.append(j)
found=None
for j,Ej in enumerate(curves):
    t0=time.time(); pts=[]
    try:
        alarm(280); r=pari(Ej).ellrank(4); cancel_alarm(); pts=[Ej(list(P)) for P in r[3]]
    except Exception as e:
        cancel_alarm(); print(' кривая',j,'прервано/ошибка',str(e)[:80],flush=True)
    pts=[P for P in pts if P.order()==oo]
    print(' кривая',j,'ellrank(4) точки',len(pts),round(time.time()-t0,1),'с',flush=True)
    if pts:
        P=pts[0]
        for phi in paths[j]: P=phi(P)
        assert P.curve()==Em and P.order()==oo
        print('ТОЧКА на E2min:',P,'высота',float(P.height()),flush=True); found=P; break
if found is not None:
    out=json.load(open('gens_open.json')); out[f'{a},{b},{c}']['E2']=[str(found)]; out[f'{a},{b},{c}']['E2_how']='isogeny class + pari ellrank(4)'
    json.dump(out,open('gens_open.json','w'),indent=1)
