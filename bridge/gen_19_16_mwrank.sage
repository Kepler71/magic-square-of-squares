import time
m,n=19,16; s=QQ(m^2+n^2)/2; b=s*m^2*n^2
r=[ZZ(-4*b),ZZ(-4*s*m^4),ZZ(-4*s*n^4)]
E=EllipticCurve(QQ,[0,-(r[0]+r[1]+r[2]),0,r[0]*r[1]+r[0]*r[2]+r[1]*r[2],-r[0]*r[1]*r[2]]).minimal_model()
cl=E.isogeny_class()
for lim in [14,16,18,20,22]:
    for idx,Ei in enumerate(cl.curves):
        t0=time.time()
        try:
            G=Ei.gens(proof=False, descent_second_limit=lim, algorithm='mwrank_shell')
            nt=[P for P in G if P.order()==oo]
            print("second_limit=%d curve %d : nontorsion found = %s   [%.0fs]"%(lim,idx,nt,time.time()-t0))
        except Exception as ex:
            print("second_limit=%d curve %d : FAILED %s [%.0fs]"%(lim,idx,str(ex)[:80],time.time()-t0))
        sys.stdout.flush()
