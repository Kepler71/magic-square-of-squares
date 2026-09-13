import json
R.<kk> = QQ[]
coef = [R(1), -(1+kk), kk, -(1-kk), (1-kk), -kk, (1+kk), R(-1)]
good = [o['idx'] for o in json.load(open("p0_torsion.json")) if o['torsion_all']]
def ainv(S, kv):
    cs=[QQ(coef[i](kv)) for i in S]
    if len(set(cs))<len(cs): return None
    if len(S)==3:
        Pp.<p>=QQ[]; f=prod(1+c*p for c in cs); A,B,C=f[3],f[2],f[1]; return [0,B,0,A*C,A^2]
    Px.<x>=QQ[]; c1=cs[0]; g=c1*prod((1-c/c1)*x+c for c in cs[1:]); A,B,C,D=g[3],g[2],g[1],g[0]
    return [0,B,0,A*C,A^2*D]
DMAX=int(150)
walls=[]; walls_cubic=[]; total=0
for b in range(2,DMAX+1):
    for a in range(1,b):
        if gcd(a,b)!=1: continue
        kv=QQ(a)/b; total+=1
        s=[]
        for S in good:
            ai=ainv(S,kv)
            if ai is None: s.append(None); continue
            E=EllipticCurve(QQ,ai).minimal_model()
            s.append(int(pari(E.a_invariants()).ellinit().ellrootno()))
        cub=[x for S,x in zip(good,s) if len(S)==3 and x is not None]
        allc=[x for x in s if x is not None]
        if cub and all(x==-1 for x in cub): walls_cubic.append(str(kv))
        if allc and all(x==-1 for x in allc): walls.append(str(kv))
print(f"наклонов со знаменателем <= {DMAX}: {total}")
print(f"знак -1 у всех 12 троек: {len(walls_cubic)}  (доля {float(len(walls_cubic))/float(total):.3f})")
print(f"знак -1 у всех 22 кривых: {len(walls)}  (доля {float(len(walls))/float(total):.3f})")
print("примеры стены (все 22):", walls[:15])
json.dump(dict(total=int(total), walls_all22=walls, walls_cubic12=walls_cubic), open("rootno_wall.json","w"))
