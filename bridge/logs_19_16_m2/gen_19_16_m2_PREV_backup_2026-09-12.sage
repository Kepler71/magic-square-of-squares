# gen_19_16_m2.sage
# G1, pair (m,n) = (19,16): find a generator of E(Q) (rank 1) and compute the
# full image of the 2-descent map delta, to decide whether the class (1,s,s)
# required by C_(19,16) is in that image.
#
# Method 3 of 3: Heegner point (PARI ellheegner), then exact verification.
#
# Run:  sage gen_19_16_m2.sage
# Optional: pass a precomputed generator via the file GEN_FILE below.

import time, os, json

# ----------------------------------------------------------------------
# 0.  The family data and the curve
# ----------------------------------------------------------------------
m, n = 19, 16
s = QQ(m^2 + n^2) / 2                 # 617/2
b = s * m^2 * n^2                     # 28510336

print("=" * 70)
print("G1 pair (m,n) = (%d,%d)   s = %s   b = %s" % (m, n, s, b))

# E : V^2 = (X+b)(X + s m^4)(X + s n^4),  roots e1=-b, e2=-s m^4, e3=-s n^4.
# s is half-integral, so clear denominators by (X,V) -> (x/4, y/8):
#     y^2 = (x + 4b)(x + 4 s m^4)(x + 4 s n^4)
# The scale factor on each (X - e_i) is 4, a square, so square classes are
# unchanged.
r = [-4*b, -4*s*m^4, -4*s*n^4]        # roots of the scaled model, in order e1,e2,e3
assert all(_ in ZZ for _ in r)
r = [ZZ(_) for _ in r]
print("scaled roots (order e1,e2,e3):", r)

Rx.<xv> = QQ[]
fscaled = (xv - r[0]) * (xv - r[1]) * (xv - r[2])
Escaled = EllipticCurve([0, fscaled[2], 0, fscaled[1], fscaled[0]])
E = Escaled.minimal_model()
iso = Escaled.isomorphism_to(E)
u_, r_, s_, t_ = iso.u, iso.r, iso.s, iso.t
print("E (minimal) :", E)
print("isomorphism (u,r,s,t) =", (u_, r_, s_, t_))
assert u_ == 1 and s_ == 0 and t_ == 0, "unexpected isomorphism shape"

# roots carried over to the minimal model, in the SAME order e1,e2,e3
R = [ZZ(ri - r_) for ri in r]
print("roots of E in order (e1,e2,e3):", R)
for Ri in R:
    assert E.is_on_curve(Ri, 0)
assert sorted(R) == sorted(E.two_division_polynomial().roots(multiplicities=False))

# ----------------------------------------------------------------------
# 1.  The class required by C, re-derived from scratch (not copied)
# ----------------------------------------------------------------------
# C : u0^2 = m^2 + n^2 t^2,  u4^2 = s(1+t^2),  u8^2 = n^2 + m^2 t^2
# X = b t^2 gives, as polynomial identities in t:
#   X - e1 = b(1+t^2)        = (mn)^2 * [ s(1+t^2) ]        = (mn)^2 u4^2
#   X - e2 = s m^2 (m^2+n^2 t^2) = s m^2 * u0^2
#   X - e3 = s n^2 (n^2+m^2 t^2) = s n^2 * u8^2
Rt.<tv> = QQ[]
X_of_t = b * tv^2
F0 = m^2 + n^2*tv^2
F4 = s * (1 + tv^2)
F8 = n^2 + m^2*tv^2
assert X_of_t - (-b)      == (m*n)^2 * F4, "identity X-e1 failed"
assert X_of_t - (-s*m^4)  == s*m^2  * F0, "identity X-e2 failed"
assert X_of_t - (-s*n^4)  == s*n^2  * F8, "identity X-e3 failed"
print("[proved] the three identities X-e_i hold identically in t")

def sqclass(q):
    """squarefree integer representing the class of q in Q*/Q*^2"""
    q = QQ(q)
    assert q != 0
    return ZZ(q.numerator() * q.denominator()).squarefree_part()

delta_star = (sqclass(1), sqclass(s), sqclass(s))
print("required class delta* = (1, s, s) =", delta_star)

# points of C over t = infinity:  with z=1/t, (u4 z)^2 -> s at z=0
print("s is a square in Q? ", QQ(s).is_square(),
      " -> points of C over t=infinity exist?", QQ(s).is_square())

# ----------------------------------------------------------------------
# 2.  Torsion and the delta map
# ----------------------------------------------------------------------
T = E.torsion_subgroup()
print("torsion structure:", T.invariants())
assert T.order() == 4

def delta(P):
    """delta(P) in (Q*/Q*^2)^3 w.r.t. the ordered roots R = [e1,e2,e3]."""
    if P == E(0):
        return (1, 1, 1)
    x = P[0]
    out = []
    for i in range(3):
        if x == R[i]:
            j, k = [u for u in range(3) if u != i]
            out.append(sqclass((R[i] - R[j]) * (R[i] - R[k])))
        else:
            out.append(sqclass(x - R[i]))
    return tuple(out)

# ----------------------------------------------------------------------
# 3.  Rank facts (each computed here, with its own status label)
# ----------------------------------------------------------------------
print("-" * 70)
print("selmer_rank  (dim_F2 Sel^2)      =", E.selmer_rank())
print("rank_bound   (mwrank, 2nd desc.) =", E.rank_bound())
print("root number                      =", E.root_number())
print("analytic rank (PARI)             =", pari(E).ellanalyticrank())

# ----------------------------------------------------------------------
# 4.  The generator
# ----------------------------------------------------------------------
GEN_FILE = "/home/kep/magicKube/bridge/gen_19_16_m2_point.json"

G = None
if os.path.exists(GEN_FILE):
    d = json.load(open(GEN_FILE))
    G = E(QQ(d["x"]), QQ(d["y"]))
    print("generator loaded from", GEN_FILE)
else:
    print("computing Heegner point (PARI ellheegner) ...")
    t0 = time.time()
    pari.allocatemem(2*10^9, 24*10^9, silent=True)
    pari.set_real_precision(300)
    H = pari(E).ellheegner()
    print("  ellheegner time: %.1f s" % (time.time() - t0))
    G = E(QQ(H[0]), QQ(H[1]))
    json.dump({"x": str(G[0]), "y": str(G[1])}, open(GEN_FILE, "w"))

# ---- (a) EXACT verification, rational arithmetic, no floating point ----
xg, yg = QQ(G[0]), QQ(G[1])
a1, a2, a3, a4, a6 = [QQ(c) for c in E.a_invariants()]
lhs = yg^2 + a1*xg*yg + a3*yg
rhs = xg^3 + a2*xg^2 + a4*xg + a6
print("-" * 70)
print("EXACT check  y^2 + a1 x y + a3 y - (x^3+a2 x^2+a4 x+a6) =", lhs - rhs)
assert lhs == rhs, "point is NOT on the curve"
print("[proved] G lies on E, verified by exact rational arithmetic")
print("x(G) numerator digits:", len(str(abs(xg.numerator()))))

# ---- (b) infinite order ----
print("order of G (Sage) :", G.order())
assert G.order() == Infinity or G.order() == oo
h = G.height()
print("canonical height  :", h)
print("[proved] G has infinite order (torsion subgroup is %s and G is not in it)"
      % (T.invariants(),))

# ---- saturation: is G a generator, or a multiple of one? ----
print("-" * 70)
try:
    t0 = time.time()
    Gsat, idx, _ = E.saturation([G])
    print("saturation index :", idx, "  (%.1f s)" % (time.time() - t0))
    G = Gsat[0]
    print("saturated G      : x =", G[0])
    print("canonical height :", G.height())
except Exception as e:
    print("saturation failed:", e)

# ----------------------------------------------------------------------
# 5.  The eight classes
# ----------------------------------------------------------------------
T1 = E(R[0], 0)
T2 = E(R[1], 0)
T3 = E(R[2], 0)
assert T1 + T2 == T3

reps = [("O",        E(0)),
        ("T1",       T1),
        ("T2",       T2),
        ("T1+T2=T3", T3),
        ("G",        G),
        ("G+T1",     G + T1),
        ("G+T2",     G + T2),
        ("G+T1+T2",  G + T3)]

print("-" * 70)
classes = []
for name, P in reps:
    d = delta(P)
    classes.append(d)
    print("%-10s delta = %s" % (name, d))

# homomorphism cross-check: delta(G+T) must equal delta(G)*delta(T) mod squares
def mul(d1, d2):
    return tuple(sqclass(d1[i]*d2[i]) for i in range(3))
for j, (name, P) in enumerate(reps[4:]):
    expected = mul(classes[4], classes[j])
    assert classes[4 + j] == expected, ("homomorphism check failed", name)
print("[checked] delta is a homomorphism on the eight representatives")

distinct = len(set(classes))
print("-" * 70)
print("number of DISTINCT classes :", distinct, "of 8")
print("required class", delta_star, "present?", delta_star in set(classes))

print("=" * 70)
if distinct == 8:
    print("image of E(Q)/2E(Q) has order 8 = 2^(rank+2) with rank=1 -> image is COMPLETE")
    if delta_star in set(classes):
        print("VERDICT: required class IS present -> pair (19,16) is NOT excluded by this test")
    else:
        print("VERDICT: required class is ABSENT -> C_(19,16) has no finite rational point;")
        print("         together with 'no points over t=infinity' this gives C_(19,16)(Q) = empty")
else:
    print("VERDICT: image incomplete (%d/8) -> no conclusion" % distinct)
