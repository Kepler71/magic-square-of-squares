import sys, json, random
sys.path.insert(0,'.')
load('both_squares.sage') if False else None
exec(preparse(open('both_squares.sage').read().split("restored=[]")[0]))
def ok_general(m,n):
    if not proves_rank0(m,n): return False
    if sq(m) and sq(n) and sq(isqrt(m)**2+isqrt(n)**2): return False   # пифагоров случай — отдельно
    return True
cov={}; closed_list=[]
for lo,hi in ((2,48),(49,100),(101,200)):
    tot=c=0
    for s in range(lo,hi+1):
        for r in range(1,s):
            if gcd(r,s)!=1 or 2*r==s: continue
            tot+=1; k=QQ(r)/s
            rows=[(nm,cf(k),red(*pr(r,s))) for nm,cf,pr in CLASSES]
            good=[(nm,cs,mn) for nm,cs,mn in rows if mn[0]>0 and mn[0]!=mn[1] and ok_general(*mn)]
            pyth=[(nm,cs,mn) for nm,cs,mn in rows if mn[0]>0 and mn[0]!=mn[1] and proves_rank0(*mn) and not ok_general(*mn)]
            if good: c+=1; closed_list.append((r,s,good[0]))
            elif pyth: c+=1; closed_list.append((r,s,pyth[0]))    # пифагоровы проверены поточечно (7/32)
    cov[f"{lo}-{hi}"]=(tot,c); print(f"знаменатели {lo}..{hi}: закрыто доказанно {c} из {tot} ({100.0*c/tot:.1f}%)")
# Контроль всей теоремы поточечно: 400 случайных закрытых наклонов, реальные точки кручения через Sage
random.seed(int(17)); sample=random.sample(closed_list,int(400)); bad=[]; tors_seen={}
for r,s,(nm,cs,mn) in sample:
    k=QQ(r)/s
    E,ps=class_points([QQ(c) for c in cs])
    t=E.torsion_subgroup().invariants(); tors_seen[t]=tors_seen.get(t,0)+1
    # ранг 0 здесь не пересчитываем: он доказан критерием; проверяем шаги (в)-(д)
    if any(p is not None and is_nondeg_square_set(k,p) for p in ps): bad.append((r,s,nm))
print(f"контроль на 400 закрытых наклонах: кручение {tors_seen}; точек, дающих невырожденный квадрат: {bad if bad else 0}")
