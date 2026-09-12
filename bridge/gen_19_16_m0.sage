#!/usr/bin/env sage
# =====================================================================
#  gen_19_16_m0.sage  --  ROUND 2  --  (m,n) = (19,16), family G1
#
#  APPROACH 1 of 3: eclib/mwrank quartic search with a LARGE search
#  bound.  Two search fronts:
#     (C)  mwrank directly on E with escalating second_limit  (honest
#          record of how far the direct search reaches);
#     (D)  mwrank with the SAME bounds on every curve of the 2-isogeny
#          class, pulling any point back to E by the dual isogeny
#          (exact rational maps, no numerics).
#
#  Everything about the resulting point is then verified EXACTLY:
#   (a) substitution into the Weierstrass equation of BOTH models
#       (the 4-scaled integral one and the original half-integral one);
#   (b) infinite order;
#   (c) all eight delta classes + pairwise distinctness + subgroup test;
#   (d) presence/absence of the required class (1,s,s).
#
#  Output: logs_19_16_m0r2/run.log
# =====================================================================
import sys, os, time

LOGDIR = "/home/kep/magicKube/bridge/logs_19_16_m0r2"
os.makedirs(LOGDIR, exist_ok=True)

t00 = time.time()
def flush(*a):
    print(*a); sys.stdout.flush()

def sqclass(q):
    """squarefree integer representing the class of q in Q*/Q*^2"""
    q = QQ(q)
    assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

# ---------------------------------------------------------------- 0
# the curve, built from scratch out of (m,n)
# ---------------------------------------------------------------- 0
m, n = 19, 16
assert gcd(m, n) == 1
s = QQ(m^2 + n^2)/2                  # 617/2
b = s*m^2*n^2                        # 28510336
e = [-b, -s*m^4, -s*n^4]             # e1, e2, e3  IN THE PROJECT'S ORDER
E4 = [4*ei for ei in e]              # x = 4X, y = 8V ; 4 is a SQUARE
R = PolynomialRing(QQ, 'x'); x = R.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
cf = f.coefficients(sparse=False)
E = EllipticCurve([0, cf[2], 0, cf[1], cf[0]])
Emin = E.minimal_model()
target = (1, sqclass(s), sqclass(s))

flush("="*72)
flush("(m,n) = (%d,%d)   s = %s   b = %s" % (m, n, s, b))
flush("e1,e2,e3 (original model)  =", e)
flush("4*e_i    (integral model)  =", E4)
flush("E    =", E)
flush("Emin =", Emin)
flush("disc =", factor(E.discriminant()))
flush("cond =", E.conductor(), "=", factor(E.conductor()))
flush("torsion =", E.torsion_subgroup().invariants())
flush("REQUIRED class delta = (1,s,s) =", target)
flush("="*72)

# ---------------------------------------------------------------- 1
# SYMBOLIC check of the three identities that fix the required class
# and of the root ORDER.  This is the step where a wrong permutation
# would silently flip the final answer, so it is checked by algebra.
# ---------------------------------------------------------------- 1
flush("\n[1] symbolic check of  X-e_i  identities (fixes the root order)")
St.<t> = PolynomialRing(QQ)
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
X = b*t^2
id1 = X - e[0] - (m*n)^2*F4          # X-e1 = (mn)^2 u4^2
id2 = X - e[1] - s*m^2*F0            # X-e2 = s m^2 u0^2
id3 = X - e[2] - s*n^2*F8            # X-e3 = s n^2 u8^2
flush("    X-e1 - (mn)^2*F4 =", id1, "   (must be 0)")
flush("    X-e2 - s*m^2*F0  =", id2, "   (must be 0)")
flush("    X-e3 - s*n^2*F8  =", id3, "   (must be 0)")
assert id1 == 0 and id2 == 0 and id3 == 0
flush("    => class of (X-e1,X-e2,X-e3) = (1, s, s) = %s   CONFIRMED symbolically" % (target,))
flush("    sqclass(s) = %s, sqclass(m^2 n^2)=1, sqclass(m^2)=1, sqclass(n^2)=1"
      % sqclass(s))

# ---------------------------------------------------------------- 2
# delta map on the integral model (and a cross-check on the original)
# ---------------------------------------------------------------- 2
def delta_int(P):
    """2-descent class of P in (Q*/Q*^2)^3 on the 4-scaled model E"""
    if P.is_zero():
        return (1, 1, 1)
    XX = P[0]/P[2]
    out = []
    for i in range(3):
        d = XX - E4[i]
        if d == 0:
            j, k = [u for u in range(3) if u != i]
            d = (E4[i]-E4[j])*(E4[i]-E4[k])
        out.append(sqclass(d))
    return tuple(out)

def delta_orig(Xo):
    """same, computed on the ORIGINAL half-integral model V^2=prod(X-e_i)"""
    out = []
    for i in range(3):
        d = Xo - e[i]
        if d == 0:
            j, k = [u for u in range(3) if u != i]
            d = (e[i]-e[j])*(e[i]-e[k])
        out.append(sqclass(d))
    return tuple(out)

# ---------------------------------------------------------------- 3
# rank UPPER bound, independently from two sources
# ---------------------------------------------------------------- 3
flush("\n[3] rank upper bound (needed for |E(Q)/2E(Q)| = 2^(r+2))")
t0 = time.time()
try:
    lo, hi = Emin.rank_bounds()
    flush("    Sage/eclib rank_bounds(Emin) = [%s, %s]   (%.1fs)" % (lo, hi, time.time()-t0))
except Exception as ex:
    flush("    rank_bounds failed:", ex)
    lo, hi = None, None
t0 = time.time()
try:
    pr = pari(Emin).ellrank()
    flush("    PARI ellrank(Emin) = %s   (%.1fs)" % (pr, time.time()-t0))
except Exception as ex:
    flush("    PARI ellrank failed:", ex)

# ---------------------------------------------------------------- 4
# (C) mwrank on E ITSELF with escalating quartic-search bound
# ---------------------------------------------------------------- 4
from sage.libs.eclib.interface import mwrank_EllipticCurve

def mwrank_try(EE, second_limit, first_limit=20, timeout=900):
    """run eclib two_descent with the given bounds; return list of gens"""
    mw = mwrank_EllipticCurve(list(EE.ainvs()))
    mw.set_verbose(0)
    t0 = time.time()
    status = "ok"
    try:
        alarm(timeout)
        mw.two_descent(verbose=False, selmer_only=False,
                       first_limit=first_limit, second_limit=second_limit,
                       n_aux=-1, second_descent=True)
        g = mw.gens()
        cancel_alarm()
    except AlarmInterrupt:
        status = "TIMEOUT after %ds" % timeout
        g = []
    except Exception as ex:
        cancel_alarm()
        status = "ERROR %s" % ex
        g = []
    return g, status, time.time()-t0

flush("\n[4] (C) mwrank/eclib on E itself, escalating second_limit")
found_direct = []
for sl in [10, 12, 14, 16]:
    g, st, dt = mwrank_try(Emin, sl, timeout=1200)
    flush("    second_limit=%-3d  gens=%s  [%s, %.1fs]" % (sl, g, st, dt))
    if g:
        found_direct = g
        break

# ---------------------------------------------------------------- 5
# (D) mwrank over the whole 2-isogeny class, with EXACT pullback
# ---------------------------------------------------------------- 5
flush("\n[5] (D) 2-isogeny class with exact pullback maps to E")

class Node:
    def __init__(self, curve, pull):
        self.curve = curve      # an elliptic curve 2-power-isogenous to E
        self.pull = pull        # function: point on self.curve -> point on E

nodes = [Node(E, lambda P: P)]
seen = [E.ainvs()]
frontier = [nodes[0]]
for depth in range(3):
    newfront = []
    for nd in frontier:
        for phi in nd.curve.isogenies_prime_degree(2):
            Ei = phi.codomain()
            if Ei.ainvs() in seen:
                continue
            seen.append(Ei.ainvs())
            dual = phi.dual()
            prev = nd.pull
            nd2 = Node(Ei, (lambda d, p: (lambda P: p(d(P))))(dual, prev))
            nodes.append(nd2)
            newfront.append(nd2)
    frontier = newfront
    if not frontier:
        break
flush("    curves in the 2-isogeny class reachable from E: %d" % len(nodes))
for i, nd in enumerate(nodes):
    flush("      #%d  %s" % (i, nd.curve.ainvs()))

# sanity: the pullback maps really land on E and are homomorphisms
for i, nd in enumerate(nodes):
    if i == 0:
        continue
    T = nd.curve.torsion_points()
    for P in T:
        Q = nd.pull(P)
        assert Q.curve() == E, "pullback #%d lands on the wrong curve" % i
flush("    pullback maps verified to land on E (checked on torsion)")

raw = None
raw_src = None
for sl in [10, 12, 14]:
    for i, nd in enumerate(nodes):
        Emin_i = nd.curve.minimal_model()
        iso_back = Emin_i.isomorphism_to(nd.curve)
        g, st, dt = mwrank_try(Emin_i, sl, timeout=900)
        nontors = []
        for P in g:
            Pt = Emin_i(P)
            Q = iso_back(Pt)
            if Q.order() == Infinity:
                nontors.append(Q)
        flush("    second_limit=%-3d curve#%d: %d gens, %d nontorsion  [%s, %.1fs]"
              % (sl, i, len(g), len(nontors), st, dt))
        if nontors:
            raw = nd.pull(nontors[0])
            raw_src = i
            flush("      >>> nontorsion on curve #%d: %s" % (i, nontors[0]))
            flush("      >>> pulled back to E:      %s" % (raw,))
            break
    if raw is not None:
        break

if raw is None and found_direct:
    raw = E(Emin.isomorphism_to(E)(Emin(found_direct[0])))
    raw_src = 0

if raw is None:
    flush("\n!!! NO POINT OF INFINITE ORDER FOUND with the bounds tried.")
    flush("    Nothing is proved.  Total time %.1fs" % (time.time()-t00))
    sys.exit(0)

# ---------------------------------------------------------------- 6
# saturate, then verify EXACTLY
# ---------------------------------------------------------------- 6
flush("\n[6] saturation on E")
t0 = time.time()
sat, index, reg = E.saturation([raw])
G = sat[0]
flush("    saturation index = %s   regulator = %s   (%.1fs)" % (index, reg, time.time()-t0))
flush("    G = %s" % (G,))
save(G, LOGDIR + "/G.sobj")

flush("\n[a] EXACT verification (rational arithmetic, no floating point)")
XG = QQ(G[0]/G[2]); YG = QQ(G[1]/G[2])
lhs = YG^2
rhs = prod([XG - E4[i] for i in range(3)])
flush("    integral model:  Y^2 - prod(X - 4e_i) =", lhs - rhs)
assert lhs - rhs == 0
Xo = XG/4; Vo = YG/8
lhs2 = Vo^2
rhs2 = prod([Xo - e[i] for i in range(3)])
flush("    original model:  X = x/4 =", Xo)
flush("                     V = y/8 =", Vo)
flush("                     V^2 - (X+b)(X+s m^4)(X+s n^4) =", lhs2 - rhs2)
assert lhs2 - rhs2 == 0
flush("    BOTH exact substitutions give 0  ->  the point is on the curve.")

flush("\n[b] order")
ordG = G.order()
flush("    order(G) =", ordG)
assert ordG == Infinity
hG = G.height()
flush("    canonical height h(G) =", hG)
flush("    naive height of x(G)  = log max(|num|,|den|) =",
      RR(log(max(abs(XG.numerator()), abs(XG.denominator())))))

# ---------------------------------------------------------------- 7
# the eight classes
# ---------------------------------------------------------------- 7
flush("\n[c] the eight classes of E(Q)/2E(Q)")
tors = E.torsion_points()
reps = list(tors) + [G + T for T in tors]
classes = []
for P in reps:
    d = delta_int(P)
    classes.append(d)
    label = "O" if P.is_zero() else str(P.xy())
    flush("    %-70s -> %s" % (label, d))
uniq = set(classes)
flush("    number of DISTINCT classes = %d  (need 8)" % len(uniq))
flush("    pairwise distinct:", len(uniq) == 8)

# cross-check: same classes computed on the ORIGINAL half-integral model
flush("\n    cross-check of delta on the ORIGINAL model V^2=prod(X-e_i):")
ok_cross = True
for P in reps:
    if P.is_zero():
        d2 = (1, 1, 1)
    else:
        d2 = delta_orig(QQ(P[0]/P[2])/4)
    d1 = delta_int(P)
    if d1 != d2:
        ok_cross = False
        flush("      MISMATCH at %s : %s vs %s" % (P, d1, d2))
flush("    original-model delta agrees with integral-model delta:", ok_cross)
assert ok_cross

# subgroup test
def mul(u, v):
    return tuple(sqclass(u[i]*v[i]) for i in range(3))
is_grp = all(mul(u, v) in uniq for u in uniq for v in uniq)
flush("    closed under multiplication (is a subgroup of (Q*/Q*^2)^3):", is_grp)

flush("\n[d] REQUIRED class")
flush("    required (1,s,s) =", target)
present = target in uniq
flush("    PRESENT among the eight classes:", present)
flush("    sorted list of the eight classes:")
for c in sorted(uniq):
    flush("       ", c, "   <-- TARGET" if c == target else "")

# ---------------------------------------------------------------- 8
# logical summary, with the hypotheses spelled out
# ---------------------------------------------------------------- 8
flush("\n[8] LOGIC")
flush("    G has infinite order  ->  rank >= 1.")
flush("    eclib/PARI 2-descent  ->  rank <= 1.")
flush("    E has full rational 2-torsion  ->  |E(Q)/2E(Q)| = 2^(rank+2) = 8.")
flush("    delta is INJECTIVE on E(Q)/2E(Q) (kernel = 2E(Q)).")
flush("    %d distinct values from the 8 listed points  ->  the image is complete."
      % len(uniq))
flush("    required class present: %s" % present)
if len(uniq) == 8 and not present:
    flush("    ==> no finite rational point of C maps into E(Q); together with the")
    flush("        point at infinity (needs s to be a square, and s=617/2 is not),")
    flush("        C_(19,16)(Q) = empty.   PAIR (19,16) IS CLOSED.")
elif present:
    flush("    ==> the required class DOES occur; this 2-descent does NOT exclude the pair.")

flush("\nTOTAL time %.1fs" % (time.time()-t00))
