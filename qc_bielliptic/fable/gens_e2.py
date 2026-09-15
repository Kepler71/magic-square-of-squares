import sys,time; sys.path.insert(0,'/home/kep/magicKube/qc_bielliptic/fable')
from models import *
a,b,c=map(int,sys.argv[1:4]); f=fpoly(a,b,c); E1,E2=E12(f); Em=E2.minimal_model(); print(Em.ainvs(),flush=True)
t0=time.time(); r=pari(Em).ellrank(10); print('ellrank(10)',r[:2],'pts',r[2],round(time.time()-t0,1),'с',flush=True)
t0=time.time()
try:
    g=Em.gens(proof=False,descent_second_limit=16,verbose=False); print('mwrank sl16',g,[float(P.height()) for P in g],round(time.time()-t0,1),'с',flush=True)
except Exception as e: print('mwrank fail',e)
