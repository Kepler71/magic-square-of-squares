from jsieve import *
import traceback
r,s=333,442
cells=slope_cells(r,s)
bad=[]; odd=[]
for k in (3,4):
    for T in itertools.combinations(sorted(cells),k):
        try: f=Factor(T)
        except Exception as e:
            bad.append((T,type(e).__name__,str(e)[:60],traceback.format_exc().splitlines()[-3][:120])); continue
        if not f.rank_proved: odd.append((T,f.lo,f.hi,f.rank))
print('не построились',len(bad)); 
for b in bad[:5]: print(b)
print('недоказанные (T, lo, hi, найдено точек):',odd)
# P_T(0): кручение или образующая?
for (r,s) in [(204,247),(333,442)]:
    cells=slope_cells(r,s); cnt={'tors':0,'gen':0,'mult':0}; ex=[]
    for k in (3,4):
        for T in itertools.combinations(sorted(cells),k):
            try: f=Factor(T)
            except Exception: continue
            if not f.rank_proved: continue
            P=f.point(0,1); cls=f.decompose(P); n=cls[:f.rank]
            if all(x==0 for x in n): cnt['tors']+=1
            elif f.rank==1 and abs(n[0])==1: cnt['gen']+=1
            elif f.rank==1: cnt['mult']+=1; ex.append((f.name,f.rank,cls))
            else: cnt['r>=2']=cnt.get('r>=2',0)+1
    print(r,s,'P_T(0):',cnt,'примеры кратных:',ex[:6])
