#!/usr/bin/env sage
# =====================================================================
#  gen_19_16_m1_audit.sage
#  AUDIT of the result of gen_19_16_m1.sage.
#  The whole exclusion argument rests on ONE external input:
#        rank E(Q) <= 1        (2-Selmer upper bound)
#  Everything else (the 8 classes) is now explicit.  So here we
#   1. recompute the 2-Selmer rank by several independent routes,
#   2. decide whether the TARGET class is even in Sel^2 (local
#      solvability of its 2-covering torsor, computed by hand),
#   3. sanity-check the whole delta picture.
# =====================================================================
import sys, time
t00 = time.time()
def flush(*a):
    print(*a); sys.stdout.flush()

def sqclass(q):
    q = QQ(q); assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

m, n = 19, 16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
e = [-b, -s*m^4, -s*n^4]
E4 = [4*x for x in e]
R = PolynomialRing(QQ,'x'); x = R.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
c = f.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
Emin = E.minimal_model()
target = (1, sqclass(s), sqclass(s))
flush("E    =", E.ainvs())
flush("Emin =", Emin.ainvs())
flush("target class =", target)

# ---------------------------------------------------------------
flush("\n=== 1. RANK UPPER BOUND, several routes ===")

t0=time.time()
try:
    rb = Emin.rank_bounds()
    flush("Sage Emin.rank_bounds()          = %s   (%.1fs)"%(rb, time.time()-t0))
except Exception as ex:
    flush("rank_bounds failed:", ex)

t0=time.time()
try:
    sr = Emin.selmer_rank()
    flush("Sage Emin.selmer_rank()          = %s   (%.1fs)"%(sr, time.time()-t0))
    flush("   => rank <= selmer_rank - rank(E(Q)[2]) = %s - 2 = %s"%(sr, sr-2))
except Exception as ex:
    flush("selmer_rank failed:", ex)

t0=time.time()
try:
    rk = Emin.rank(only_use_mwrank=True, proof=True)
    flush("Sage Emin.rank(proof=True)       = %s   (%.1fs)"%(rk, time.time()-t0))
except Exception as ex:
    flush("rank(proof=True) failed: %s   (%.1fs)"%(ex, time.time()-t0))

t0=time.time()
try:
    pr = pari(Emin).ellrank()
    flush("PARI ellrank(Emin)               = %s   (%.1fs)"%(pr, time.time()-t0))
except Exception as ex:
    flush("PARI ellrank failed:", ex)

t0=time.time()
try:
    ar = Emin.analytic_rank()
    flush("analytic rank (via L-series)     = %s   (%.1fs)"%(ar, time.time()-t0))
except Exception as ex:
    flush("analytic_rank failed:", ex)

# mwrank with the 'certain' flag
from sage.libs.eclib.interface import mwrank_EllipticCurve
for lim in [12, 16]:
    t0=time.time()
    mw = mwrank_EllipticCurve(list(Emin.ainvs()))
    mw.two_descent(verbose=False, selmer_only=False,
                   first_limit=lim, second_limit=lim, n_aux=-1, second_descent=True)
    flush("mwrank lim=%d : rank=%s  certain=%s  selmer_rank=%s  #gens=%d  (%.1fs)"%(
        lim, mw.rank(), mw.certain(), mw.selmer_rank(), len(mw.gens()), time.time()-t0))

# ---------------------------------------------------------------
flush("\n=== 2. IS THE TARGET CLASS IN Sel^2 ? ===")
# The 2-covering attached to (d1,d2,d3), d1 d2 d3 in Q*^2, for
#   y^2 = (x-a1)(x-a2)(x-a3) :
#     d1 z1^2 - d2 z2^2 = a2 - a1
#     d1 z1^2 - d3 z3^2 = a3 - a1
# (x = a1 + d1 z1^2).  Membership in Sel^2 <=> this has points over
# every Q_p and over R.
a = E4
d = [QQ(t) for t in target]
flush("roots a1,a2,a3 =", a)
flush("(d1,d2,d3) =", d, "  product square:", (d[0]*d[1]*d[2]).is_square())

A = a[1]-a[0]
B = a[2]-a[0]
flush("a2-a1 = %s = %s"%(A, factor(A)))
flush("a3-a1 = %s = %s"%(B, factor(B)))

# real solvability
flush("\n-- real place --")
flush("   d1>0:", d[0]>0, " need d1 z1^2 - d2 z2^2 = A =", A, " and d1 z1^2 - d3 z3^2 = B =", B)

# p-adic solvability by brute force lifting (Hensel / exhaustive mod p^k)
def solvable_padic(d, A, B, p, prec):
    """search for (z1,z2,z3) in Z_p^3 (up to scaling) solving the pair,
    by exhaustive search modulo p^prec on suitably scaled variables."""
    mod = p^prec
    Zm = Zmod(mod)
    # allow z_i = p^{-k} u_i ; clear denominators: multiply the system by p^{2k}
    for k in range(0, 8):
        sc = p^(2*k)
        found = False
        rng = range(mod)
        # z1 ranges over units*p^j ; do a plain search, prec small
        for z1 in rng:
            v1 = Zm(d[0]*z1^2 - sc*A)
            # need d2 z2^2 = v1
            ok2 = any(Zm(d[1]*z2^2) == v1 for z2 in rng)
            if not ok2: continue
            v2 = Zm(d[0]*z1^2 - sc*B)
            ok3 = any(Zm(d[2]*z3^2) == v2 for z3 in rng)
            if ok3:
                found = True; break
        if found:
            return True, k
    return False, None

# do it properly with Sage's own machinery instead: use the 2-covering as a
# genus-one quartic / conic intersection and test local solvability via
# Qp points on the associated quartic  y^2 = d1*(d1 z^2 ... ) -- simpler:
# the class (d1,d2,d3) is in the image of the local descent map at p iff
# the torsor has a Qp point.  We use the quartic model:
#   from  d1 z1^2 - d2 z2^2 = A  and  d1 z1^2 - d3 z3^2 = B
# set z1 = 1 (scaling) -> conic intersection; instead use the standard
# quartic:  d1 * Y^2 = d1^2 U^4 + ... ; Sage provides this via
# EllipticCurve.descend_via_2_isogeny? Not general.  We test directly with
# the intersection of two conics in P^3 using Sage's local solubility for
# quartics:  eliminate z1.
flush("\n-- p-adic places: test the pair of conics for Qp-points --")
# Homogeneous: d1 z1^2 - d2 z2^2 - A w^2 = 0 ,  d1 z1^2 - d3 z3^2 - B w^2 = 0
# in P^3 with coordinates (z1:z2:z3:w).  Equivalent quartic (w=1 patch
# eliminated): from the first,  z2^2 = (d1 z1^2 - A w^2)/d2 ; substitute in
# nothing -- we instead do the standard trick: the torsor is isomorphic to
# the quartic  y^2 = g(u) with
#     g(u) = d1*d2*(d1 u^2 - A)*(d1 u^2 - B)/(d2*d3)  ... (up to squares)
# Cleanest: the curve  y^2 = d1*(d1 u^2 - A)*(d1 u^2 - B) * d1  is wrong.
# Use: z1 = u w.  Then need BOTH  (d1 u^2 - A)/d2  and  (d1 u^2 - B)/d3  to
# be squares.  Product must be a square:  (d1u^2-A)(d1u^2-B)/(d2 d3) = square
# Since d1 d2 d3 = square,  d2 d3 ~ d1.  So the quartic is
#     y^2 = d1 (d1 u^2 - A)(d1 u^2 - B)              (*)
# and the torsor maps 2:1 onto (*).  ELS of the torsor implies ELS of (*),
# and for the purpose of an OBSTRUCTION a failure of (*) is decisive.
K.<u> = QQ[]
quart = d[0]*(d[0]*u^2 - A)*(d[0]*u^2 - B)
flush("quartic (*) : y^2 =", quart)
flush("   factored  :", factor(quart))

from sage.schemes.hyperelliptic_curves.constructor import HyperellipticCurve
bad = [2,3,5,7,19,617]
flush("   testing local solubility of (*) at", bad, "and R")
for p in bad:
    try:
        Cq = HyperellipticCurve(quart)
        ok = Cq.has_rational_point(algorithm='local', bound=p) if False else None
    except Exception as ex:
        ok = None
    flush("     p=%s -> (skipped, see explicit test below)"%p)

# explicit: does (d1 u^2 - A)/d2 and (d1 u^2 - B)/d3 admit a common Qp solution?
def is_square_Qp(t, p, prec=40):
    t = QQ(t)
    if t == 0: return True
    v = t.valuation(p)
    if v % 2: return False
    tu = t / p^v
    if p == 2:
        return Zmod(8)(tu.numerator()*tu.denominator().inverse_mod(8)) == 1
    num = tu.numerator(); den = tu.denominator()
    return kronecker(num*den, p) == 1

def torsor_has_Qp_point(d, A, B, p, N=None):
    """search u in Q_p (approximated by rationals u = c/p^j and c/1) for which
    (d1 u^2 - A)/d2 and (d1 u^2 - B)/d3 are both squares in Qp."""
    if N is None: N = p^6 if p < 20 else 4000
    cand = []
    for j in range(-6, 7):
        for cnum in range(0, min(N, 5000)):
            cand.append(QQ(cnum)*p^j)
            if cnum: cand.append(-QQ(cnum)*p^j)
    for uu in cand:
        t1 = (d[0]*uu^2 - A)/d[1]
        t2 = (d[0]*uu^2 - B)/d[2]
        if t1 == 0 or t2 == 0:
            # still a point provided the other is a square
            if (t1 == 0 and is_square_Qp(t2,p)) or (t2 == 0 and is_square_Qp(t1,p)):
                return True, uu
            continue
        if is_square_Qp(t1, p) and is_square_Qp(t2, p):
            return True, uu
    # also the point at infinity: u -> oo  means z1/w -> oo, i.e. w=0:
    # d1 z1^2 = d2 z2^2 = d3 z3^2  -> need d1 d2 and d1 d3 squares in Qp
    if is_square_Qp(d[0]*d[1], p) and is_square_Qp(d[0]*d[2], p):
        return True, "infinity"
    return False, None

flush("\n-- explicit search for Qp points on the target torsor --")
for p in [2,3,5,7,19,617]:
    t0 = time.time()
    ok, wit = torsor_has_Qp_point(d, A, B, p)
    flush("   p=%-4s solvable: %-5s witness u=%s   (%.1fs)"%(p, ok, wit, time.time()-t0))
# real place
okR = None
# need d1 u^2 - A > 0 wrt sign of d2 etc.
flush("   real: d=(%s,%s,%s), A=%s, B=%s"%(d[0],d[1],d[2],A,B))
found_real = False
for uu in [QQ(0), QQ(1), QQ(10)^6, QQ(1)/1000, QQ(-1)]:
    t1 = (d[0]*uu^2 - A)/d[1]; t2 = (d[0]*uu^2 - B)/d[2]
    if t1 >= 0 and t2 >= 0:
        found_real = True; flush("   real solvable, u=%s"%uu); break
if not found_real:
    flush("   real: no u found among the probes (needs a proper interval argument)")

flush("\nTOTAL %.1fs"%(time.time()-t00))
