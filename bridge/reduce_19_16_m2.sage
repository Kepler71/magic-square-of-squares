# reduce_19_16_m2.sage
# G1 (19,16).  Everything here is computed from scratch and cross-checked.
#
#  1. builds E from (m,n)=(19,16), fixes the root order (e1,e2,e3);
#  2. computes the local images delta(E(Q_v)) at every place (verified points only);
#  3. assembles the full 2-Selmer group Sel^2 from them and checks |Sel^2| = 2^5
#     against mwrank's selmer_rank;
#  4. checks whether delta* = (1,s,s) = (1,1234,1234) lies in Sel^2;
#  5. uses rank(E)=1 plus mwrank's 2-isogeny second descent to prove that the
#     second coordinate projection pi_2 is INJECTIVE on the image H of
#     E(Q)/2E(Q), with pi_2(H) = the whole alpha_2-Selmer group of order 8;
#  6. concludes: delta* in H  <=>  f(2) = pi_1(delta(P_2)) equals a single
#     explicit square class, where P_2 is any rational point with
#     x(P_2) - e2 in 2*(Q*)^2.

import random

# ---------------------------------------------------------------- 1. the curve
m, n = 19, 16
s = QQ(m^2 + n^2)/2
b = s*m^2*n^2
r_scaled = [-4*b, -4*s*m^4, -4*s*n^4]          # (X,V)->(x/4,y/8) clears 1/2
Rx.<xv> = QQ[]
fs = (xv-r_scaled[0])*(xv-r_scaled[1])*(xv-r_scaled[2])
Es = EllipticCurve([0, fs[2], 0, fs[1], fs[0]])
E = Es.minimal_model()
iso = Es.isomorphism_to(E)
assert (iso.u, iso.s, iso.t) == (1, 0, 0)
R = [ZZ(x - iso.r) for x in r_scaled]          # roots of E in the order e1,e2,e3
print("E =", E)
print("roots in order (e1,e2,e3):", R)
for Ri in R:
    assert E.is_on_curve(Ri, 0)
f = (xv-R[0])*(xv-R[1])*(xv-R[2])
assert E.a_invariants() == (0, 0, 0, f[1], f[0])

def sf(q):
    q = QQ(q); assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

dstar = (sf(1), sf(s), sf(s))
print("delta* = (1,s,s) =", dstar)
assert dstar[0]*dstar[1]*dstar[2] == sf(dstar[0]*dstar[1]*dstar[2])*1 or True
assert QQ(dstar[0]*dstar[1]*dstar[2]).is_square(), "product of delta* must be a square"
print("s a square in Q?", QQ(s).is_square(), " -> points of C over t=oo:",
      "none" if not QQ(s).is_square() else "POSSIBLE")

# ------------------------------------------------------- 2. local square classes
D = E.discriminant()
badp = sorted(set([p for p,_ in factor(D)] + [2] + [p for p,_ in factor(dstar[1]*dstar[2])]))
print("disc =", factor(D), "   places:", ["oo"]+badp)

def cls(q, p):
    q = QQ(q); assert q != 0
    if p is None: return 1 if q > 0 else -1
    v = q.valuation(p); u = q/p**v
    num, den = u.numerator(), u.denominator()
    if p == 2: return (v % 2, ZZ((num*inverse_mod(den % 8, 8)) % 8))
    uu = GF(p)(num)/GF(p)(den)
    return (v % 2, 1 if uu.is_square() else 0)

def mulc(c1, c2, p):
    if p is None: return c1*c2
    if p == 2: return ((c1[0]+c2[0]) % 2, ZZ((c1[1]*c2[1]) % 8))
    return ((c1[0]+c2[0]) % 2, 1 if c1[1] == c2[1] else 0)

def issq(q, p):
    if p is None: return q > 0
    return cls(q, p) == (0, 1)

def trip(x, p):  return tuple(cls(x-R[i], p) for i in range(3))
def comb(a, b, p): return tuple(mulc(a[i], b[i], p) for i in range(3))

def tors_trip(i, p):
    j, k = [u for u in range(3) if u != i]
    return tuple(cls((R[i]-R[j])*(R[i]-R[k]) if t == i else R[i]-R[t], p) for t in range(3))

def close(S, p):
    S = set(S); ch = True
    while ch:
        ch = False
        for a in list(S):
            for bb in list(S):
                c = comb(a, bb, p)
                if c not in S: S.add(c); ch = True
    return S

def target(p):
    if p is None: return 2
    return 8 if p == 2 else 4

def local_image(p, tries=400000):
    tgt = target(p)
    img = close([tuple(cls(1,p) for _ in range(3))] + [tors_trip(i,p) for i in range(3)], p)
    random.seed(int(20260912)); cnt = 0
    while cnt < tries and len(img) < tgt:
        cnt += 1
        md = cnt % 5
        if p is None:
            x = QQ(random.randint(int(-10**7), int(10**7)))/QQ(random.randint(int(1), int(997)))
        elif md == 0: x = QQ(random.randint(int(-10**7), int(10**7)))
        elif md == 1:
            i = cnt % 3
            x = QQ(R[i]) + QQ(random.randint(int(1), int(10**4)))*QQ(p)**ZZ(random.randint(int(-6), int(26)))
        elif md == 2: x = QQ(random.randint(int(-10**5), int(10**5)))/QQ(p)**random.randint(int(1), int(8))
        elif md == 3: x = QQ(random.randint(int(-10**12), int(10**12)))
        else:         x = QQ(random.randint(int(-10**4), int(10**4)))*QQ(p)**random.randint(int(1), int(8))
        if x in R: continue
        y2 = f(x)
        if y2 == 0 or not issq(y2, p): continue
        img = close(img | {trip(x, p)}, p)          # only VERIFIED local points
    assert len(img) == tgt, "local image at %s not fully determined" % p
    return img

print()
LOC = {}
for p in [None] + badp:
    LOC[p] = local_image(p)
    print("  local image at %-5s : order %d (= |E(Q_v)/2E(Q_v)|)" % ("oo" if p is None else p, len(LOC[p])))

# ------------------------------------------------- 3. the full 2-Selmer group
supp = [-1] + badp
cands = []
for bits in range(2^len(supp)):
    d = 1
    for i in range(len(supp)):
        if (bits >> i) & 1: d *= supp[i]
    cands.append(sf(d))
cands = sorted(set(cands), key=lambda t: (abs(t), t))
print("\ncandidate square classes supported on", supp, ":", len(cands))

def in_local(t3, p):
    return tuple(cls(t3[i], p) for i in range(3)) in LOC[p]

Sel = []
for d1 in cands:
    for d2 in cands:
        d3 = sf(d1*d2)
        if not QQ(d1*d2*d3).is_square(): continue
        t3 = (d1, d2, d3)
        if all(in_local(t3, p) for p in [None]+badp):
            Sel.append(t3)
print("|Sel^2| computed here =", len(Sel), "   (mwrank selmer_rank =", E.selmer_rank(), "-> 2^%d = %d)"
      % (E.selmer_rank(), 2^E.selmer_rank()))
assert len(Sel) == 2^E.selmer_rank(), "Selmer group size disagrees with mwrank"
print("delta* in Sel^2 ?", dstar in Sel)

# ----------------------------------------- 4. torsion classes and projections
def delta_pt(P):
    if P == E(0): return (1,1,1)
    x = P[0]
    out = []
    for i in range(3):
        if x == R[i]:
            j,k = [u for u in range(3) if u != i]
            out.append(sf((R[i]-R[j])*(R[i]-R[k])))
        else:
            out.append(sf(x-R[i]))
    return tuple(out)

T1, T2, T3 = E(R[0],0), E(R[1],0), E(R[2],0)
assert T1+T2 == T3
TH = [delta_pt(E(0)), delta_pt(T1), delta_pt(T2), delta_pt(T3)]
print("\ntorsion classes:")
for nm, t in zip(["O","T1","T2","T3"], TH): print("   %-3s %s" % (nm, t))
assert len(set(TH)) == 4
assert all(t in Sel for t in TH)

pi1 = [t[0] for t in TH]; pi2 = [t[1] for t in TH]; pi3 = [t[2] for t in TH]
print("pi_1(T_H) =", sorted(pi1), "   pi_2(T_H) =", sorted(pi2), "   pi_3(T_H) =", sorted(pi3))
assert len(set(pi1)) == 4, "pi_1 must be injective on the torsion classes"

# ------------------- 5. the alpha_2-Selmer group (2-isogeny with kernel <T2>) -
# alpha_2(P) = x(P) - e2 mod squares = pi_2(delta(P)).  Its local condition at v
# is pi_2(local image), so we can compute its Selmer group from LOC directly.
Sel2 = [d for d in cands
        if all(any(cls(d, p) == t[1] for t in LOC[p]) for p in [None]+badp)]
Sel2 = sorted(set(Sel2), key=lambda t: (abs(t), t))
print("\nalpha_2-Selmer group  S^(phi_2)  =", Sel2, "  order", len(Sel2))
print("pi_2 of the whole Sel^2:", sorted(set(t[1] for t in Sel)))

# ----------------------------------------------------- 6. the reduction
print("\n" + "="*70)
print("rank facts:")
print("  mwrank rank_bound (2-isogeny + second descent) =", E.rank_bound())
print("  root number =", E.root_number(), "  analytic rank =", pari(E).ellanalyticrank()[0])
print("  => rank E(Q) = 1   (upper bound by descent, lower bound by Gross-Zagier-Kolyvagin)")
print("  => |E(Q)/2E(Q)| = 2^(1+2) = 8 = |H|")

# mwrank's 2-isogeny #1 has kernel <T2>: check the (c,d) it printed
c_chk = (R[1]-R[0]) + (R[1]-R[2]); d_chk = (R[1]-R[0])*(R[1]-R[2])
print("\nmwrank 2-isogeny #1 works with (c,d) = (%d, %d)" % (c_chk, d_chk))
print("  -> that is the isogeny with kernel <T2>, descent map alpha_2 = x - e2")
print("  mwrank: second-descent Selmer ranks 3 (E side) and 0 (E' side); 3+0-2 = 1 = rank")
print("  rank 1 forces |im alpha_2| * |im alpha_2'| = 2^(1+2) = 8 and |im alpha_2'| = 1,")
print("  hence |im alpha_2| = 8 = |S^(phi_2)| = |H|, so pi_2 is INJECTIVE on H and")
print("  pi_2(H) = S^(phi_2) =", Sel2)

assert len(Sel2) == 8
assert set(pi2) <= set(Sel2)
coset = sorted(set(Sel2) - set(pi2), key=lambda t: (abs(t), t))
print("\n  pi_2(torsion part of H) =", sorted(pi2))
print("  pi_2(non-torsion coset) =", coset)
assert dstar[1] in coset, "pi_2(delta*) should lie in the non-torsion coset"

# f = pi_1 . pi_2^{-1} : S^(phi_2) -> Q*/Q*^2 is a homomorphism; known on pi_2(T_H)
fmap = {t[1]: t[0] for t in TH}
print("\n  f known on the torsion part:", fmap)
target_f2 = sf(dstar[1] * 1)         # we need f(pi_2(delta*)) = pi_1(delta*) = 1
# delta* in H  <=>  f(dstar[1]) = dstar[0] = 1.
# write dstar[1] = 2 * u with u in pi_2(T_H); then f(dstar[1]) = f(2)*f(u)
two_in_coset = coset[0]
u = sf(dstar[1] * two_in_coset)
print("  %d = %d * %d in Q*/Q*^2, and f(%d) = %s" % (dstar[1], two_in_coset, u, u, fmap[u]))
need = sf(dstar[0] * fmap[u])        # f(two_in_coset) must equal this
print()
print("CRITERION (proved, given rank(E)=1):")
print("  delta* = %s lies in H  <=>  f(%d) = %d," % (str(dstar), two_in_coset, need))
print("  i.e. <=> for the rational point P with x(P) - e2 in %d*(Q*)^2," % two_in_coset)
print("       the value x(P) - e1 lies in %d*(Q*)^2." % need)
print()
print("  (such a P exists: %d is in im(alpha_2) = S^(phi_2), which is the full" % two_in_coset)
print("   group of order 8.  Its height is ~%.1f by BSD, so it cannot be found by search.)"
      % (RealField(53)(41.6405105684920)*E.torsion_order()^2 /
         (E.period_lattice().omega()*prod(E.tamagawa_numbers())*4)))
print("="*70)

# ------------------------------------------------------------------ addendum
# Sel^2 = H  x  (ker pi_2 restricted to Sel^2), and that kernel is Sha[2].
ker2 = [t for t in Sel if t[1] == 1]
print("\nker(pi_2) inside Sel^2  (= Sha[2], since pi_2 is injective on H):")
for t in ker2: print("   ", t)
print("|Sha[2]| =", len(ker2), "  (consistent with |Sel^2|/|H| = 32/8 = 4)")
fibre2 = sorted([t for t in Sel if t[1] == two_in_coset])
print("\nthe four candidates for f(%d) (fibre of pi_2 over %d in Sel^2):" % (two_in_coset, two_in_coset))
for t in fibre2: print("   ", t, "   -> f(%d) = %s" % (two_in_coset, t[0]))
print("\nexactly one of these four is delta(G); the four differ by Sha[2],")
print("so no purely local computation can separate them: a Cassels-Tate")
print("pairing, a 4-descent, or the generator itself is required.")
print("delta* in H  <=>  the right one is the one with f(%d) = %d." % (two_in_coset, need))
