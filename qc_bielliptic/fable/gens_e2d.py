# третья попытка для E2: 2-накрытия PARI (ell2cover) + поиск точек на квартиках hyperellratpoints с растущей границей
import sys,time,json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f); Em=E2.minimal_model(); print('E2min',Em.ainvs(),flush=True)
t0=time.time(); cov=pari(Em).ell2cover(); print('2-накрытий:',len(cov),round(time.time()-t0,1),'с',flush=True)
found=None
for B in (10**4,10**5,10**6,10**7):
    for i in range(len(cov)):
        q,m=cov[i][0],cov[i][1]; t0=time.time()
        try:
            alarm(900); pts=pari.hyperellratpoints(q,B); cancel_alarm()
        except Exception as e:
            cancel_alarm(); print(' B=',B,'накрытие',i,'прервано',str(e)[:60],flush=True); continue
        print(' B=',B,'накрытие',i,'точек',len(pts),round(time.time()-t0,1),'с',flush=True)
        for P in pts:
            try:
                Q=pari.substvec(m,['x','y'],[P[0],P[1]]) if False else m.subst('x',P[0]).subst('y',P[1])
                PQ=Em(list(Q))
                if PQ.order()==oo: found=PQ; print('ТОЧКА на E2min:',PQ,'высота',float(PQ.height()),flush=True); break
            except Exception as e: print('  перенос не удался',str(e)[:80],flush=True)
        if found is not None: break
    if found is not None: break
if found is not None:
    out=json.load(open('gens_open.json')); out[f'{a},{b},{c}']['E2']=[str(found)]; out[f'{a},{b},{c}']['E2_how']='ell2cover + hyperellratpoints'
    json.dump(out,open('gens_open.json','w'),indent=1)
