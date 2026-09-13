import json, itertools
R.<kk> = QQ[]
coef = [R(1), -(1+kk), kk, -(1-kk), (1-kk), -kk, (1+kk), R(-1)]
good = [o['idx'] for o in json.load(open("p0_torsion.json")) if o['torsion_all']]
def curve(S, kv):
    cs=[QQ(coef[i](kv)) for i in S]
    if len(S)==3:
        Pp.<p>=QQ[]; f=prod(1+c*p for c in cs); A,B,C=f[3],f[2],f[1]; return EllipticCurve(QQ,[0,B,0,A*C,A^2])
    Px.<x>=QQ[]; c1=cs[0]; g=c1*prod((1-c/c1)*x+c for c in cs[1:]); A,B,C,D=g[3],g[2],g[1],g[0]
    return EllipticCurve(QQ,[0,B,0,A*C,A^2*D])
for ks in ['19/60','41/60','38/63','29/65','33/97','80/97','1/48','40/43']:
    kv=QQ(ks); ivs=[]; errs=0
    for S in good:
        try:
            r=pari(curve(S,kv).minimal_model().a_invariants()).ellinit().ellrank(); ivs.append((int(r[0]),int(r[1])))
        except Exception as e: errs+=1
    zero=sum(1 for lo,hi in ivs if hi==0); proven_pos=sum(1 for lo,hi in ivs if lo>=1)
    print(f"k={ks:6} кривых {len(ivs)}, ошибок {errs}: ранг 0 доказан у {zero}, ранг>=1 доказан у {proven_pos}, не решено {len(ivs)-zero-proven_pos}")
