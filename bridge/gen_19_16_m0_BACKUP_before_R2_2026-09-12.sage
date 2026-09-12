#!/usr/bin/env sage
# =====================================================================
#  gen_19_16_m0.sage   --  (m,n) = (19,16), G1 family
#  Approach 1 of 3: mwrank/eclib quartic search with LARGE search bound,
#  plus direct ratpoints search over the whole 2-isogeny class.
#  Goal: find a generator G of E(Q) (rank 1) and compute all eight
#        delta-classes of E(Q)/2E(Q).
#  All arithmetic on found points is EXACT (rational), never numerical.
# =====================================================================
import sys, time

def sqclass(q):
    q = QQ(q); assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

# ---------------- the curve ----------------
m, n = 19, 16
s = QQ(m^2+n^2)/2                 # 617/2
b = s*m^2*n^2                     # 28510336
e = [-b, -s*m^4, -s*n^4]          # e1,e2,e3  in the order fixed by the project
E4 = [4*x for x in e]             # integral model  x = 4X, y = 8V ; 4 is a SQUARE
R = PolynomialRing(QQ,'x'); x = R.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
c = f.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
target = (1, sqclass(s), sqclass(s))

print("="*70)
print("(m,n) = (%d,%d)   s = %s   b = %s"%(m,n,s,b))
print("e1,e2,e3          =", e)
print("4e (integral model)=", E4)
print("E  =", E)
print("Emin =", E.minimal_model())
print("torsion =", E.torsion_subgroup().invariants())
print("REQUIRED class delta = (1,s,s) =", target)
print("="*70)

def delta(P):
    if P.is_zero(): return (1,1,1)
    X = P[0]/P[2]
    out=[]
    for i in range(3):
        d = X - E4[i]
        if d == 0:
            j,k = [t for t in range(3) if t!=i]
            d = (E4[i]-E4[j])*(E4[i]-E4[k])
        out.append(sqclass(d))
    return tuple(out)

def full_report(G):
    """G: a point of infinite order on E. Exact checks (a)-(g)."""
    print()
    print("#"*70); print("GENERATOR FOUND"); print("#"*70)
    X = G[0]/G[2]; Y = G[1]/G[2]
    print("G = (%s , %s)"%(X,Y))
    # (a) EXACT substitution into the Weierstrass equation
    lhs = Y^2
    rhs = (X-E4[0])*(X-E4[1])*(X-E4[2])
    print("(a) EXACT check  Y^2 - prod(X-4e_i) =", lhs-rhs, "  -> on curve:", lhs==rhs)
    assert lhs == rhs
    # (b) infinite order
    print("(b) order of G =", G.order(), "  canonical height =", G.height())
    assert G.order() == Infinity
    # saturation
    try:
        sat = E.saturation([G])
        print("    saturation:", sat)
        G = sat[0][0]
        print("    saturated generator G =", G)
    except Exception as ex:
        print("    saturation failed:", ex)
    # (c) all eight classes
    T = E.torsion_points()
    cls = {}
    for t in T:
        cls[str(t)] = delta(t)
        cls["G+"+str(t)] = delta(G+t)
    print("(c) the eight delta classes:")
    for k,v in cls.items(): print("     %-28s -> %s"%(k,v))
    vals = list(cls.values())
    print("    number of DISTINCT classes =", len(set(vals)), " (need 8 = 2^(r+2), r=1)")
    # (d) required class present?
    present = target in set(vals)
    print("(d) required class", target, "PRESENT:", present)
    return G, set(vals), present

# ---------------- searches ----------------
t_start = time.time()
found = None

print()
print(">>> stage 1: direct ratpoints search on E and its 2-isogenous curves")
cl = E.minimal_model().isogeny_class()
for hb in [10, 12, 14, 16]:
    for i,C in enumerate(cl.curves):
        t0=time.time()
        pts = C.point_search(hb, rank_bound=1)
        print("   curve %d  height_bound=%d  -> %d pts   (%.1fs)"%(i,hb,len(pts),time.time()-t0))
        sys.stdout.flush()
        if pts:
            print("      points:", pts)
            found = (i,C,pts); break
    if found: break

print()
print(">>> stage 2: mwrank/eclib two_descent with large second_limit / search bound")
Emin = E.minimal_model()
from sage.libs.eclib.interface import mwrank_EllipticCurve
mw = mwrank_EllipticCurve(list(Emin.ainvs()))
for lim1, lim2, sl in [(10,10,8),(12,12,10),(14,14,12),(16,16,14),(18,18,16),(20,20,18)]:
    t0=time.time()
    try:
        mw.two_descent(verbose=False, selmer_only=False,
                       first_limit=lim1, second_limit=lim2,
                       n_aux=-1, second_descent=True)
        g = mw.gens()
        print("   two_descent(first=%d,second=%d): rank=%s gens=%s   (%.1fs)"%(
              lim1,lim2,mw.rank(),g,time.time()-t0))
        sys.stdout.flush()
        if g:
            found = ('mwrank', Emin, [Emin(p) for p in g]); break
    except Exception as ex:
        print("   two_descent(first=%d,second=%d) FAILED: %s  (%.1fs)"%(lim1,lim2,ex,time.time()-t0))
        sys.stdout.flush()

print()
print("total search time %.1f s"%(time.time()-t_start))
if found:
    tag, C, pts = found
    # move the point back to E
    P = pts[0]
    if C != E:
        try:
            iso = C.isogeny_class().curves  # not used
        except Exception: pass
    print("RAW FOUND on", C.ainvs(), ":", pts)
else:
    print("NO point of infinite order found by stages 1-2.")
