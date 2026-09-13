import json
crit=json.load(open("unify_claude/crit_g1.json"))
def EB(m,n):  return EllipticCurve(QQ,[0,-((m^2-n^2)^2+(m^2+n^2)^2),0,(m^2-n^2)^2*(m^2+n^2)^2,0])
r0=0; caught=0; false=[]
for key,pure in crit.items():
    m,n=map(int,key.split(','))
    r=pari(EB(m,n).minimal_model().a_invariants()).ellinit().ellrank(); lo,hi=int(r[0]),int(r[1])
    if hi==0: r0+=1
    if pure and hi==0: caught+=1
    if pure and hi!=0: false.append((m,n,lo,hi))
print(f"ранг B = 0 (доказан PARI): {r0} пар из {len(crit)}; критерий ловит {caught}; ложных 'закрыто': {false}")
