# генераторы E1=E^-, E2~E^+ (модель кода Bianchi) для пяти троек (1,1) трёх открытых наклонов;
# сначала mwrank (Sage gens(proof=False)), если пусто — PARI ellrank с большим усилием; точки переносятся на минимальную модель
import sys,time,json; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
out={}
def find(Em):
    t0=time.time()
    try:
        alarm(900); g=list(Em.gens(proof=False)); cancel_alarm()
    except Exception as e:
        cancel_alarm(); g=[]
    if g: return g,'mwrank',time.time()-t0
    for eff in (4,8,12):
        t0=time.time()
        try:
            alarm(1200); r=pari(Em).ellrank(eff); cancel_alarm()
        except Exception as e:
            cancel_alarm(); continue
        pts=[Em(list(P)) for P in r[2]] if len(r[2]) else []
        pts=[P for P in pts if P.order()==oo]
        if pts:
            s=Em.saturation(pts); return list(s[0]),f'ellrank({eff}) sat index {s[1]}',time.time()-t0
    return [],'нет',0
for a,b,c in [(33,265,298),(33,298,563),(265,298,563),(73,289,362),(126,451,577)]:
    f=fpoly(a,b,c); E1,E2=E12(f); row={}
    for E,name in ((E1,'E1'),(E2,'E2')):
        Em=E.minimal_model(); g,how,dt=find(Em)
        row[name]=[str(P) for P in g]; row[name+'_how']=how
        print((a,b,c),name,g,how,'высота',[float(P.height()) for P in g],round(dt,1),'с',flush=True)
    good=[int(p) for p in primes(5,60) if E1.has_good_reduction(p) and E2.has_good_reduction(p) and E1.is_ordinary(p) and E2.is_ordinary(p)]
    row['good_ordinary']=good; print((a,b,c),'p:',good,flush=True)
    out[f'{a},{b},{c}']=row
    json.dump(out,open('gens_open.json','w'),indent=1)
