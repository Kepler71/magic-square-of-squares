import sys,time; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f)
for E,name in ((E1,'E1'),(E2,'E2')):
    Em=E.minimal_model(); t0=time.time()
    r=pari(Em).ellrank(3); print(name,'ellrank(3)',r[:2],'точки',r[2],round(time.time()-t0,1),'с',flush=True)
    t0=time.time()
    try:
        g=Em.gens(proof=False); print(name,'gens (mwrank)',g,round(time.time()-t0,1),'с',flush=True)
    except Exception as e: print(name,'gens fail',e,flush=True)
