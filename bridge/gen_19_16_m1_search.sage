import time
t0=time.time()
m=19;n=16
s=QQ(m^2+n^2)/2; b=s*m^2*n^2
R.<X>=QQ[]
f=(X+b)*(X+s*m^4)*(X+s*n^4)
E=EllipticCurve([0,f[2],0,f[1],f[0]])
Em=E.minimal_model()
print("Em",Em); import sys; sys.stdout.flush()

curves = [("Em", Em)]
for i,phi in enumerate(Em.isogenies_prime_degree(2)):
    curves.append(("iso%d"%i, phi.codomain().minimal_model()))
# second layer: 2-isogenies of the isogenous curves
seen = [c[1] for c in curves]
extra=[]
for nm,C in list(curves[1:]):
    for j,psi in enumerate(C.isogenies_prime_degree(2)):
        D=psi.codomain().minimal_model()
        if not any(D.is_isomorphic(S) for S in seen):
            seen.append(D); extra.append((nm+"_%d"%j, D))
curves += extra
print("curves in 2-isogeny class explored:", [c[0] for c in curves])
sys.stdout.flush()

for h in [12,16,20,24]:
    print("=== point_search height bound", h, " t=%.1f"%(time.time()-t0))
    sys.stdout.flush()
    for nm,C in curves:
        try:
            pts=C.point_search(h, rank_bound=1)
        except Exception as e:
            print("   ",nm,"ERR",e); sys.stdout.flush(); continue
        nz=[P for P in pts if P.order()==oo]
        print("   ",nm,"-> found",len(pts),"gens-of-search; infinite order:",nz, " t=%.1f"%(time.time()-t0))
        sys.stdout.flush()
        if nz:
            print("   !!!! INFINITE ORDER POINT ON",nm,":",nz)
            sys.stdout.flush()
