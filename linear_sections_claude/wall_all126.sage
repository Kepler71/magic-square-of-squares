import itertools
R.<kk> = QQ[]
coef = [R(1), -(1+kk), kk, -(1-kk), (1-kk), -kk, (1+kk), R(-1)]
def curve(S, kv):
    cs=[QQ(coef[i](kv)) for i in S]
    if len(S)==3:
        Pp.<p>=QQ[]; f=prod(1+c*p for c in cs); A,B,C=f[3],f[2],f[1]; return EllipticCurve(QQ,[0,B,0,A*C,A^2])
    Px.<x>=QQ[]; c1=cs[0]; g=c1*prod((1-c/c1)*x+c for c in cs[1:]); A,B,C,D=g[3],g[2],g[1],g[0]
    return EllipticCurve(QQ,[0,B,0,A*C,A^2*D])
for ks in ['80/97','73/105','79/110','43/145']:
    kv=QQ(ks); plus=[]
    for size in (3,4):
        for S in itertools.combinations(range(8),size):
            E=curve(S,kv).minimal_model()
            if int(pari(E.a_invariants()).ellinit().ellrootno())==1: plus.append(S)
    print(f"k={ks}: кривых со знаком +1 среди всех 126: {len(plus)}")
    zero=[]; errs=[]; intervals={}
    for S in plus:
        E=curve(S,kv).minimal_model()
        try:
            r=pari(E.a_invariants()).ellinit().ellrank()
            lo,hi=int(r[0]),int(r[1])
            intervals[(lo,hi)]=intervals.get((lo,hi),0)+1
        except Exception as e:
            errs.append(repr(e)[:120]); lo,hi=None,None
        if hi==0: zero.append(S)
    print(f"      ошибок ellrank: {len(errs)} {errs[:2]}; интервалы ранга: {intervals}")
    print(f"      из них с доказанным rank=0 (PARI ellrank [0,0]): {len(zero)}  {[list(s) for s in zero[:6]]}")
