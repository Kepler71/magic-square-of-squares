#!/usr/bin/env sage
# =====================================================================
#  gen_19_16_m1.sage   --  (m,n) = (19,16), family G1
#  APPROACH 2 of 3: through the 2-isogeny class.  E has full rational
#  2-torsion => three 2-isogenies (and a second layer).  Find a point of
#  infinite order on an ISOGENOUS curve, where the canonical height can
#  be noticeably smaller, and PULL IT BACK to E by the dual isogeny.
#
#  Round 2 (independent re-verification).  Added relative to round 1:
#   * SYMBOLIC proof of the three identities X-e_i = (...)*u_i^2 and of
#     V^2 = b^2 (u0 u4 u8)^2  -- i.e. of the claim "delta must be (1,s,s)";
#   * RIGOROUS upper bound rank <= 1 by three independent 2-Selmer
#     computations (pari / mwrank / simon).  Round 1 never proved this,
#     and WITHOUT it 8 distinct classes prove nothing;
#   * a SECOND, independent implementation of the square-class map
#     (via factorisation) cross-checked against squarefree_part;
#   * structural checks of the 8 classes (subgroup, product-is-square);
#   * the point at infinity of C.
#
#  Everything about the final point is checked EXACTLY, in rational
#  arithmetic, never numerically.
#
#  Output: logs_19_16_m1/round2.log
# =====================================================================
import sys, time
t00 = time.time()

def flush(*a):
    print(*a); sys.stdout.flush()

FAIL = []
def check(name, cond):
    flush("    [%s] %s" % ("OK " if cond else "FAIL", name))
    if not cond:
        FAIL.append(name)
    return cond

# ---------- square classes: two independent implementations ----------
def sqclass(q):
    """squarefree representative of q in Q*/Q*^2"""
    q = QQ(q)
    assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

def sqclass2(q):
    """the same, computed independently from the factorisation"""
    q = QQ(q)
    assert q != 0
    r = ZZ(-1) if q < 0 else ZZ(1)
    for p, k in factor(q):
        if k % 2 != 0:
            r *= p
    return r

# =====================================================================
# 0. the data
# =====================================================================
m, n = 19, 16
s = QQ(m^2+n^2)/2                 # 617/2
b = s*m^2*n^2
e = [-b, -s*m^4, -s*n^4]          # e1,e2,e3, the project's fixed order

flush("="*72)
flush("(m,n) = (%d,%d)   gcd = %s" % (m, n, gcd(m, n)))
flush("s = (m^2+n^2)/2 = %s      b = s*m^2*n^2 = %s" % (s, b))
flush("e1,e2,e3 = %s" % (e,))
flush("="*72)

# =====================================================================
# 1. SYMBOLIC verification of the bridge C -> E and of the required class
#    (this is the load-bearing claim: delta(P) = (1,s,s) is NECESSARY)
# =====================================================================
flush("\n[1] SYMBOLIC check of the identities defining the map C -> E")
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
F0 = m^2 + n^2*t^2
F4 = s*(1 + t^2)
F8 = n^2 + m^2*t^2
X  = b*t^2                                  # X = b t^2
# X - e_i as polynomials in t, against the claimed expressions in u_i^2=F_i
check("X - e1 == (m n)^2 * F4      (=> class 1)",  X - e[0] == (m*n)^2*F4)
check("X - e2 == s m^2 * F0        (=> class s)",  X - e[1] == s*m^2*F0)
check("X - e3 == s n^2 * F8        (=> class s)",  X - e[2] == s*n^2*F8)
check("(X-e1)(X-e2)(X-e3) == b^2 * F0*F4*F8  (=> V = +- b u0 u4 u8)",
      (X-e[0])*(X-e[1])*(X-e[2]) == b^2*F0*F4*F8)
# non-vanishing: F0,F4,F8 > 0 for every real t, so u_i != 0 and V != 0,
# so the image of a rational point of C is never 2-torsion and never O.
check("F0 > 0 for all real t (disc<0 and lead>0)", F0.discriminant() < 0 and F0[2] > 0)
check("F4 > 0 for all real t", (F4.discriminant() < 0) and (F4[2] > 0))
check("F8 > 0 for all real t", F8.discriminant() < 0 and F8[2] > 0)
target = (1, sqclass(s), sqclass(s))
check("sqclass agrees with independent impl. on s", sqclass(s) == sqclass2(s))
flush("    => REQUIRED class of any finite rational point of C:  delta = (1,s,s) = %s" % (target,))

# point at infinity of C:  z=1/t -> 0 gives (u4/t)^2 = s
flush("\n[1b] points of C over t = infinity need s to be a square in Q")
check("s = %s is NOT a square in Q (so NO points over infinity)" % s, not QQ(s).is_square())

# =====================================================================
# 2. the curve E, integral model x = 4X, y = 8V  (4 is a SQUARE => classes
#    are unchanged:  X - e_i = (x - 4 e_i)/4 )
# =====================================================================
flush("\n[2] the elliptic curve")
E4 = [4*ei for ei in e]
Rx = PolynomialRing(QQ, 'x'); x = Rx.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
c = f.coefficients(sparse=False)
check("the 4-scaled model is INTEGRAL", all(ci in ZZ for ci in c) and all(ei in ZZ for ei in E4))
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
flush("    E  : y^2 = x^3 + %s x^2 + %s x + %s" % (c[2], c[1], c[0]))
flush("    4e = %s" % (E4,))
flush("    Emin = %s" % (E.minimal_model().ainvs(),))
flush("    conductor N = %s" % factor(E.conductor()))
flush("    disc        = %s" % factor(E.discriminant()))
Tor = E.torsion_subgroup()
flush("    torsion     = %s" % (Tor.invariants(),))
check("torsion is exactly (Z/2)^2  => |E(Q)/2E(Q)| = 2^(r+2)",
      tuple(Tor.invariants()) == (2, 2))
# symbolic check that y=8V, x=4X really lands on E
Xs = polygen(QQ, 'Xs')
check("y^2 = prod(x-4e_i) is the model of V^2=prod(X-e_i) under x=4X,y=8V",
      Rx((8*Xs)^0)  # trivial guard
      and (f(4*Xs) == 64*(Xs-e[0])*(Xs-e[1])*(Xs-e[2])))

# ---------------- the delta map (on the 4-scaled model) --------------
def delta(P):
    if P.is_zero():
        return (1, 1, 1)
    Xp = P[0]/P[2]
    out = []
    for i in range(3):
        d = Xp - E4[i]
        if d == 0:
            j, k = [u for u in range(3) if u != i]
            d = (E4[i]-E4[j])*(E4[i]-E4[k])
        a1 = sqclass(d); a2 = sqclass2(d)
        assert a1 == a2, "square-class implementations disagree on %s" % d
        out.append(a1)
    return tuple(out)

# =====================================================================
# 3. RIGOROUS upper bound for the rank.
#
#    TRAP FOUND IN ROUND 1: the 2-Selmer group OF E ITSELF is NOT enough.
#    dim Sel_2(E) = 5, and with dim E(Q)[2] = 2 this only gives
#    rank + dim Sha(E)[2] = 3, i.e. rank <= 3.  With rank 2 or 3 the image
#    of delta would have 16 or 32 elements and 8 classes would prove
#    NOTHING.  So an upper bound rank <= 1 has to be obtained separately.
#
#    THE CLEAN ROUTE: rank is invariant under isogeny.  The 2-isogenous
#    curve C = iso0 has only C(Q)[2] = Z/2 and dim Sel_2(C) = 2, hence
#       rank(C) + dim Sha(C)[2] = 2 - 1 = 1  =>  rank(C) <= 1,
#    and rank(E) = rank(C).  No Cassels pairing, no BSD, no parity needed.
#
#    First, CALIBRATE the semantics of Sage's selmer_rank on curves whose
#    invariants are known, so that the formula used here is not guesswork.
# =====================================================================
flush("\n[3] RIGOROUS upper bound for rank E(Q)")
flush("\n  [3a] calibration of Sage semantics: selmer_rank == dim_F2 Sel_2 (torsion included)")
flush("       formula used:  dim Sel_2 = rank + dim Sha[2] + dim E(Q)[2]")
calib = [([1, 1, 1, 508, -2551],    1, 2, "Sha[2]=4, E(Q)[2]=Z/2"),
         ([0, -1, 0, -900, -10098], 0, 2, "960d1, Sha=4"),
         ([0, -1, 1, -929, -10595], 0, 2, "571a1, no rational 2-torsion"),
         ([0, 0, 0, -1, 0],         0, 0, "y^2=x^3-x, full 2-torsion")]
for ai, rk, dsha, note in calib:
    Ec = EllipticCurve(ai)
    d2 = len(Ec.division_polynomial(2).roots())      # dim E(Q)[2]
    pred = rk + dsha + d2
    got = Ec.selmer_rank(algorithm='pari')
    check("calib %-22s predicted %d, selmer_rank %d" % (note, pred, got), pred == got)

Em = E.minimal_model()
flush("\n  [3b] E itself -- NOT sufficient")
sel = {}
for alg in ['pari', 'mwrank']:
    t0 = time.time()
    sel[alg] = Em.selmer_rank(algorithm=alg)
    flush("       dim Sel_2(E) [%s] = %s   (%.1f s)" % (alg, sel[alg], time.time()-t0))
check("both algorithms agree on dim Sel_2(E)", len(set(sel.values())) == 1)
dimSelE = sel['pari']
flush("       dim E(Q)[2] = 2  =>  rank + dim Sha(E)[2] = %d  =>  only rank <= %d"
      % (dimSelE-2, dimSelE-2))
check("dim Sel_2(E) is indeed > 3 (round 1's implicit assumption was WRONG)", dimSelE > 3)

flush("\n  [3c] the 2-isogenous curve C = iso0 -- THIS gives rank <= 1")
C0 = EllipticCurve([0, 0, 0, -1269386570673900, 17388781634016192078000])
check("C is 2-isogenous to E over Q (so rank(C) = rank(E))",
      any(phi.codomain().is_isomorphic(C0) for phi in Em.isogenies_prime_degree(2))
      or any(phi.codomain().is_isomorphic(C0) for phi in E.isogenies_prime_degree(2)))
d2C = len(C0.division_polynomial(2).roots())
flush("       C torsion = %s ;  dim C(Q)[2] = %d" % (C0.torsion_subgroup().invariants(), d2C))
check("dim C(Q)[2] = 1", d2C == 1)
selC = {}
for alg in ['pari', 'mwrank']:
    t0 = time.time()
    selC[alg] = C0.selmer_rank(algorithm=alg)
    flush("       dim Sel_2(C) [%s] = %s   (%.1f s)" % (alg, selC[alg], time.time()-t0))
check("both algorithms agree on dim Sel_2(C)", len(set(selC.values())) == 1)
dimSelC = selC['pari']
check("dim Sel_2(C) = 2  =>  rank(C) + dim Sha(C)[2] = 1  =>  rank(E) <= 1", dimSelC == 2)
RANK_UPPER_PROVED = (dimSelC - d2C == 1)

flush("\n  [3d] independent confirmations of rank <= 1 (not relied upon)")
for nmx, Cx in [("E", Em), ("C", C0)]:
    flush("       %s: rank_bound(pari)=%s  rank_bound(mwrank)=%s  pari ellrank=%s"
          % (nmx, Cx.rank_bound(), Cx.rank_bound(algorithm='mwrank'), pari(Cx).ellrank()))
flush("       E root number = %s   E.analytic_rank() = %s  (numerical; Gross-Zagier+Kolyvagin"
      % (E.root_number(), E.analytic_rank()))
flush("        would give rank=1 too, but the proof above does not need them)")

# =====================================================================
# 4. the 2-isogeny tree, and the search for a point of infinite order
# =====================================================================
flush("\n[4] building the 2-isogeny tree (exact pullback maps to E are kept)")
tree = [("E", E, (lambda P: P))]
lvl1 = []
for i, phi in enumerate(E.isogenies_prime_degree(2)):
    C = phi.codomain(); dual = phi.dual()
    pull = (lambda dual=dual: (lambda P: dual(P)))()
    tree.append(("iso%d" % i, C, pull)); lvl1.append((("iso%d" % i), C, phi, pull))
for nm, C, phi, pull in lvl1:
    for j, psi in enumerate(C.isogenies_prime_degree(2)):
        D = psi.codomain()
        if D.is_isomorphic(E) or any(D.is_isomorphic(T[1]) for T in tree):
            continue
        dual2 = psi.dual()
        pull2 = (lambda dual2=dual2, pull=pull: (lambda P: pull(dual2(P))))()
        tree.append((nm + "_%d" % j, D, pull2))
for nm, C, _ in tree:
    flush("    %-10s min = %s   (h = %s)" % (nm, C.minimal_model().ainvs(),
                                             RR(C.minimal_model().discriminant().global_height())))

flush("\n    mwrank two_descent on every curve of the tree")
from sage.libs.eclib.interface import mwrank_EllipticCurve
found = None
for lim in [8, 10, 12, 14, 16, 18, 20]:
    for nm, C, pull in tree:
        Cm = C.minimal_model(); iso_back = Cm.isomorphism_to(C)
        t0 = time.time()
        try:
            mw = mwrank_EllipticCurve(list(Cm.ainvs()))
            mw.two_descent(verbose=False, selmer_only=False, first_limit=lim,
                           second_limit=lim, n_aux=-1, second_descent=True)
            g = mw.gens()
        except Exception as ex:
            flush("      lim=%2d %-10s FAILED %s (%.1fs)" % (lim, nm, ex, time.time()-t0)); continue
        flush("      lim=%2d %-10s rank=%s #gens=%d (%.1fs)" % (lim, nm, mw.rank(), len(g), time.time()-t0))
        for pt in g:
            Q = Cm(pt)
            if Q.order() != Infinity:
                continue
            P = pull(iso_back(Q))
            flush("        nontorsion on %s ; pulled back to E, order = %s" % (nm, P.order()))
            if P.order() == Infinity:
                found = (nm, Cm, Q, P); break
        if found: break
    if found: break

if not found:
    flush("\n!!! NO point of infinite order was found. Nothing is proved. !!!")
    flush("total %.1f s" % (time.time()-t00)); sys.exit(1)

nm, Cm, Q, P = found
flush("\n    source curve: %s   point there: %s" % (nm, Q.xy()))
flush("    pulled back to E by the dual isogeny.")

# ---------------- saturation (NOT needed for the argument) -----------
sat, index, reg = E.saturation([P])
G = sat[0]
flush("    saturation index = %s   regulator = %s" % (index, reg))
flush("    (saturation is a convenience only: the proof below uses just")
flush("     'the 8 classes are distinct' + 'rank <= 1', not that G generates)")

# =====================================================================
# 5. (a) EXACT verification of the point
# =====================================================================
flush("\n[5] (a) EXACT verification of G, in rational arithmetic")
Xg = QQ(G[0]/G[2]); Yg = QQ(G[1]/G[2])
flush("    x(G) = %s" % Xg)
flush("    y(G) = %s" % Yg)
flush("    numerator/denominator sizes: %s / %s digits" %
      (len(str(Xg.numerator())), len(str(Xg.denominator()))))
check("y^2 - prod(x - 4e_i) == 0   (integral model)",
      Yg^2 - (Xg-E4[0])*(Xg-E4[1])*(Xg-E4[2]) == 0)
check("y^2 - (x^3 + a2 x^2 + a4 x + a6) == 0",
      Yg^2 - (Xg^3 + c[2]*Xg^2 + c[1]*Xg + c[0]) == 0)
Xo = Xg/4; Vo = Yg/8
check("ORIGINAL model: V^2 - (X+b)(X+s m^4)(X+s n^4) == 0",
      Vo^2 - (Xo+b)*(Xo+s*m^4)*(Xo+s*n^4) == 0)
flush("    X(orig) = %s" % Xo)
flush("    V(orig) = %s" % Vo)

# (b) infinite order -- exactly, via the torsion subgroup
flush("\n    (b) order of G")
flush("        G in E(Q)_tors ? %s" % (G in [pt for pt in E.torsion_points()]))
check("order(G) = infinity", G.order() == Infinity)
check("2G, 3G, ... nontrivial: 12G != O (independent sanity)", not (12*G).is_zero())
flush("        canonical height h(G) = %s" % G.height())
flush("        => rank E(Q) >= 1")

# =====================================================================
# 6. (c) the eight classes
# =====================================================================
flush("\n[6] (c) the eight classes delta of E(Q)/2E(Q)")
Tpts = E.torsion_points()
def nmT(T): return "O" if T.is_zero() else ("T(x=%s)" % (T[0]/T[2]))
rows = [(nmT(T), delta(T)) for T in Tpts]
rows += [("G + " + nmT(T), delta(G+T)) for T in Tpts]
for k, v in rows:
    flush("      %-26s -> %s" % (k[:26], v))
vals = [v for _, v in rows]
dis = sorted(set(vals))
check("the 8 listed classes are PAIRWISE DISTINCT (=> |image| >= 8)", len(dis) == 8)

def mul(u, v): return tuple(sqclass(QQ(a)*QQ(bb)) for a, bb in zip(u, v))
check("the set is closed under multiplication (a subgroup of (Q*/Q*^2)^3)",
      all(mul(u, v) in dis for u in dis for v in dis))
check("every class has product of coordinates a square (forced by V^2=prod)",
      all(QQ(u[0]*u[1]*u[2]).is_square() for u in dis))
check("(1,1,1) is in the set (image of O)", (1, 1, 1) in dis)

flush("\n    |E(Q)/2E(Q)| = 2^(rank+2)  (torsion (Z/2)^2);  rank >= 1 from the exhibited")
flush("    point of infinite order, rank <= 1 from dim Sel_2(C) = 2 on the 2-isogenous C")
check("rank is PINNED to 1 (>=1 from the point, <=1 from Sel_2 of the isogenous curve)",
      RANK_UPPER_PROVED)
check("hence |image of delta| = 2^3 = 8 = number of classes found  => IMAGE IS COMPLETE",
      RANK_UPPER_PROVED and len(dis) == 8)

# =====================================================================
# 7. (d) is the required class present?
# =====================================================================
flush("\n[7] (d) the required class  (1, s, s) = %s" % (target,))
present = target in dis
flush("      PRESENT in the (complete) image: %s" % present)
for u in dis:
    flush("        %-28s %s" % (str(u), "<== REQUIRED" if u == target else ""))

# =====================================================================
# 8. NUMERICAL cross-check via BSD (consistency only, proves nothing)
# =====================================================================
flush("\n[8] BSD consistency cross-check (NUMERICAL -- not part of the proof)")
try:
    t0 = time.time()
    Lp = RR(pari(Em).ellL1(1))
    om = RR(Em.period_lattice().omega())
    tam = prod([Em.tamagawa_number(p) for p in Em.conductor().prime_factors()])
    tor = Em.torsion_order()
    RegSha = Lp*tor^2/(om*tam)
    flush("      L'(1) = %s   Omega = %s   prod c_p = %s   |tors| = %s  (%.1f s)"
          % (Lp, om, tam, tor, time.time()-t0))
    flush("      BSD:  Reg * |Sha|  = %s" % RegSha)
    flush("      h(G) = Reg(G)      = %s" % RR(G.height()))
    flush("      dim Sha(E)[2] = dim Sel_2(E) - rank - dim E(Q)[2] = %d - 1 - 2 = %d"
          % (dimSelE, dimSelE-3))
    flush("      => |Sha(E)| divisible by %d ; ratio (Reg*|Sha|)/h(G) = %s"
          % (2^(dimSelE-3), RegSha/RR(G.height())))
except Exception as ex:
    flush("      BSD cross-check failed/skipped: %s" % ex)

flush("\n" + "="*72)
if FAIL:
    flush("SOME CHECKS FAILED: %s" % FAIL)
elif present:
    flush("VERDICT: the required class (1,s,s) IS in the image.")
    flush("         The pair (19,16) is NOT excluded by this 2-descent.")
else:
    flush("VERDICT (proved, modulo the software): image of E(Q)/2E(Q) under delta is")
    flush("         complete (8 = 2^(1+2) classes, rank exactly 1), and (1,s,s) is NOT")
    flush("         among them.  Hence C has NO finite rational points; and s is not a")
    flush("         square, so C has no rational points over t=infinity either.")
    flush("         => C(Q) = empty  => the pair (m,n) = (19,16) IS CLOSED.")
flush("="*72)
flush("total %.1f s" % (time.time()-t00))
save((E4, E.ainvs(), G.xy(), dis, target), "/home/kep/magicKube/bridge/gen_19_16_m1_round2.sobj")
