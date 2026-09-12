#!/usr/bin/env sage
# =====================================================================
#  gen_19_16_m1.sage   --  (m,n) = (19,16), G1 family
#  APPROACH 2 of 3: go through the 2-isogeny class (full rational
#  2-torsion => three 2-isogenies, and a second layer), find a point of
#  infinite order on an ISOGENOUS curve where the height is smaller,
#  then PULL IT BACK to E by the dual isogeny.
#
#  Everything about the final point is checked EXACTLY (rational
#  arithmetic), never numerically.
#
#  Output: logs_19_16_m1/stageA.log
# =====================================================================
import sys, time
t00 = time.time()

def flush(*a):
    print(*a); sys.stdout.flush()

def sqclass(q):
    q = QQ(q)
    assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

# ---------------- the curve ----------------
m, n = 19, 16
s = QQ(m^2+n^2)/2                 # 617/2
b = s*m^2*n^2                     # 28510336
e = [-b, -s*m^4, -s*n^4]          # e1,e2,e3 in the project's fixed order
E4 = [4*x for x in e]             # x = 4X, y = 8V; the factor 4 is a SQUARE
R = PolynomialRing(QQ,'x'); x = R.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
c = f.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
target = (1, sqclass(s), sqclass(s))

flush("="*70)
flush("(m,n) = (%d,%d)   s = %s   b = %s"%(m,n,s,b))
flush("e1,e2,e3           =", e)
flush("4e (integral model)=", E4)
flush("E      =", E)
flush("Emin   =", E.minimal_model())
flush("disc   =", factor(E.discriminant()))
flush("cond   =", factor(E.conductor()))
flush("torsion=", E.torsion_subgroup().invariants())
flush("REQUIRED class delta = (1,s,s) =", target)
flush("="*70)

# ---------------- delta map on the E4 model ----------------
def delta(P):
    if P.is_zero():
        return (1,1,1)
    X = P[0]/P[2]
    out = []
    for i in range(3):
        d = X - E4[i]
        if d == 0:
            j,k = [t for t in range(3) if t != i]
            d = (E4[i]-E4[j])*(E4[i]-E4[k])
        out.append(sqclass(d))
    return tuple(out)

# ---------------- build the 2-isogeny tree with pullback maps -------
flush("\n>>> building the 2-isogeny tree (pullback maps to E kept exactly)")
tree = []          # (name, curve C, function pull: point on C -> point on E)

def ident(P): return P
tree.append(("E", E, ident))

lvl1 = []
for i, phi in enumerate(E.isogenies_prime_degree(2)):
    C = phi.codomain()
    dual = phi.dual()                    # C -> E
    def mk(dual=dual):
        return lambda P: dual(P)
    pull = mk()
    nm = "iso%d"%i
    tree.append((nm, C, pull))
    lvl1.append((nm, C, phi, pull))

for nm, C, phi, pull in lvl1:
    for j, psi in enumerate(C.isogenies_prime_degree(2)):
        D = psi.codomain()
        if D.is_isomorphic(E):
            continue
        if any(D.is_isomorphic(T[1]) for T in tree):
            continue
        dual2 = psi.dual()               # D -> C
        def mk2(dual2=dual2, pull=pull):
            return lambda P: pull(dual2(P))
        tree.append((nm+"_%d"%j, D, mk2()))

flush("curves in the explored 2-isogeny class:")
for nm, C, _ in tree:
    Cm = C.minimal_model()
    flush("   %-10s %s   min=%s" % (nm, C.ainvs(), Cm.ainvs()))

# ---------------- search on each curve via mwrank two_descent -------
from sage.libs.eclib.interface import mwrank_EllipticCurve

flush("\n>>> mwrank two_descent on every curve of the tree (minimal models)")
found = None
for lim in [8, 10, 12, 14, 16, 18]:
    for nm, C, pull in tree:
        Cm = C.minimal_model()
        iso_back = Cm.isomorphism_to(C)          # Cm -> C
        t0 = time.time()
        try:
            mw = mwrank_EllipticCurve(list(Cm.ainvs()))
            mw.two_descent(verbose=False, selmer_only=False,
                           first_limit=lim, second_limit=lim,
                           n_aux=-1, second_descent=True)
            g = mw.gens()
        except Exception as ex:
            flush("   lim=%2d  %-10s FAILED %s  (%.1fs)"%(lim,nm,ex,time.time()-t0))
            continue
        flush("   lim=%2d  %-10s rank=%s  #gens=%d  (%.1fs)"%(
              lim, nm, mw.rank(), len(g), time.time()-t0))
        if g:
            for pt in g:
                Q = Cm(pt)
                if Q.order() != Infinity:
                    continue
                flush("      >>> nontorsion on %s (min model): %s"%(nm, Q))
                Q1 = iso_back(Q)                  # onto the working model C
                P  = pull(Q1)                     # onto E by the dual isogeny
                flush("      >>> pulled back to E: %s"%(P,))
                if P.order() == Infinity:
                    found = (nm, Cm, Q, P)
                    break
                else:
                    flush("      (pullback is torsion -- discarded)")
            if found:
                break
    if found:
        break

flush("\ntotal search time %.1f s"%(time.time()-t00))

if not found:
    flush("NO point of infinite order pulled back to E.")
    sys.exit(0)

nm, Cm, Q, P = found
flush("\n" + "#"*70)
flush("RAW POINT OF INFINITE ORDER ON E, obtained from curve %s"%nm)
flush("#"*70)
flush("P = %s"%(P,))

# ---------------- saturation ----------------
flush("\n>>> saturating on E")
sat, index, reg = E.saturation([P])
flush("saturation index = %s   regulator = %s"%(index, reg))
G = sat[0]
flush("G (saturated generator) = %s"%(G,))
save((E4, E.ainvs(), G.xy()), "/home/kep/magicKube/bridge/gen_19_16_m1_G.sobj")

# ---------------- (a) EXACT check on the curve ----------------
flush("\n(a) EXACT verification of the Weierstrass equation (rational arithmetic)")
Xg = QQ(G[0]/G[2]); Yg = QQ(G[1]/G[2])
lhs = Yg^2
rhs = (Xg-E4[0])*(Xg-E4[1])*(Xg-E4[2])
flush("    X = %s"%Xg)
flush("    Y = %s"%Yg)
flush("    Y^2 - prod(X - 4e_i) = %s"%(lhs-rhs))
flush("    ON CURVE (exact):", lhs == rhs)
assert lhs == rhs
# and the ORIGINAL model  V^2 = (Xo+b)(Xo+s m^4)(Xo+s n^4),  Xo = X/4, V = Y/8
Xo = Xg/4; Vo = Yg/8
lhs0 = Vo^2
rhs0 = (Xo+b)*(Xo+s*m^4)*(Xo+s*n^4)
flush("    original model: Xo = X/4 = %s"%Xo)
flush("    original model: V  = Y/8 = %s"%Vo)
flush("    V^2 - (Xo+b)(Xo+s m^4)(Xo+s n^4) = %s"%(lhs0-rhs0))
flush("    ON ORIGINAL CURVE (exact):", lhs0 == rhs0)
assert lhs0 == rhs0

# ---------------- (b) infinite order ----------------
flush("\n(b) order and height")
flush("    order(G) = %s"%G.order())
flush("    canonical height h(G) = %s"%G.height())
assert G.order() == Infinity

# ---------------- (c) all eight classes ----------------
flush("\n(c) the eight delta classes of E(Q)/2E(Q)  (r=1 => 2^(1+2)=8)")
T = E.torsion_points()
flush("    torsion points: %s"%(T,))
rows = []
for t in T:
    rows.append((str(t), delta(t)))
for t in T:
    rows.append(("G+"+str(t), delta(G+t)))
for k,v in rows:
    flush("     %-46s -> %s"%(k,v))
vals = [v for _,v in rows]
dis = set(vals)
flush("    number of DISTINCT classes = %d   (need 8)"%len(dis))
flush("    pairwise distinct:", len(dis) == len(vals))
flush("    the 8 classes: %s"%(sorted(dis),))

# group check: the set must be a subgroup of (Q*/Q*^2)^3
def mul(u,v): return tuple(sqclass(QQ(a)*QQ(bb)) for a,bb in zip(u,v))
closed = all(mul(u,v) in dis for u in dis for v in dis)
flush("    closed under multiplication (is a subgroup):", closed)

# ---------------- (d) required class present? ----------------
flush("\n(d) REQUIRED class (1,s,s) = %s"%(target,))
flush("    PRESENT in the image of E(Q)/2E(Q):", target in dis)
flush("\nTOTAL time %.1f s"%(time.time()-t00))
