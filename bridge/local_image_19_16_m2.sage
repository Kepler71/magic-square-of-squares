# local_image_19_16_m2.sage
# Independent audit for G1 (19,16): is delta* = (1,1234,1234) in the 2-Selmer
# group of E?  We compute the LOCAL image delta(E(Q_v)) at each relevant place
# using ONLY verified local points, and stop when the subgroup reaches its
# theoretical order  |E(Q_v)/2E(Q_v)|  (2 at oo, 4 at odd p, 8 at p=2).
# Since delta is injective on E(Q_v)/2E(Q_v), reaching that order proves the
# subgroup found IS the whole local image, so membership of delta* is decided.

import random

E = EllipticCurve([0,0,0,-1613212693425900, 7222574158576522830000])
R = [4534950, -42239820, 37704870]        # roots in the order e1, e2, e3
for Ri in R:
    assert E.is_on_curve(Ri, 0)
Rs = sorted(R)
dstar = (1, 1234, 1234)

Rx.<X> = QQ[]
f = (X - R[0]) * (X - R[1]) * (X - R[2])
# f must be the 2-division cubic of E (a1=a2=a3=0 model)
assert E.a_invariants() == (0, 0, 0, f[1], f[0]) and f[2] == 0

D = E.discriminant()
badp = sorted(set([p for p, _ in factor(D)] + [2] + [p for p, _ in factor(1234)]))
print("disc =", factor(D))
print("places:", ["oo"] + badp)
print("v_2 of roots:", [ZZ(x).valuation(2) for x in R])
print()

# ---------------- square classes in Q_v^*/(Q_v^*)^2 ----------------
def cls(q, p):
    q = QQ(q)
    assert q != 0
    if p is None:
        return 1 if q > 0 else -1
    v = q.valuation(p)
    u = q / p**v
    num, den = u.numerator(), u.denominator()
    if p == 2:
        return (v % 2, ZZ((num * inverse_mod(den % 8, 8)) % 8))
    uu = GF(p)(num) / GF(p)(den)
    return (v % 2, 1 if uu.is_square() else 0)

def mulc(c1, c2, p):
    if p is None:
        return c1 * c2
    if p == 2:
        return ((c1[0] + c2[0]) % 2, ZZ((c1[1] * c2[1]) % 8))
    return ((c1[0] + c2[0]) % 2, 1 if c1[1] == c2[1] else 0)

def issq(q, p):
    """is the nonzero rational q a square in Q_v?"""
    if p is None:
        return q > 0
    return cls(q, p) == (0, 1)

def combine(a, b, p):
    return tuple(mulc(a[i], b[i], p) for i in range(3))

def trip(x, p):
    """delta of a non-2-torsion point with abscissa x (assumed a local point)"""
    return tuple(cls(x - R[i], p) for i in range(3))

def tors_trip(i, p):
    j, k = [u for u in range(3) if u != i]
    out = []
    for t in range(3):
        val = (R[i] - R[j]) * (R[i] - R[k]) if t == i else R[i] - R[t]
        out.append(cls(val, p))
    return tuple(out)

def close(S, p):
    S = set(S)
    changed = True
    while changed:
        changed = False
        for a in list(S):
            for b in list(S):
                c = combine(a, b, p)
                if c not in S:
                    S.add(c); changed = True
    return S

def target_order(p):
    if p is None: return 2      # disc > 0: E(R) has 2 components, E(R)/2E(R) = Z/2
    if p == 2:    return 8      # |E(Q2)[2]| * |2|_2^{-1} = 4*2
    return 4                    # p odd: |E(Qp)[2]| * 1

# ---------------- build each local image from verified points ----------------
def local_image(p, tries=400000, verbose=False):
    tgt = target_order(p)
    img = set()
    img.add(tuple(cls(1, p) for _ in range(3)))       # delta(O)
    for i in range(3):
        img.add(tors_trip(i, p))                       # delta(T_i): always legitimate
    img = close(img, p)
    witnesses = {}
    random.seed(int(20260912))
    cnt = 0
    while cnt < tries and len(img) < tgt:
        cnt += 1
        if p is None:
            x = QQ(random.randint(int(-10**7), int(10**7))) / QQ(random.randint(int(1), int(1000)))
        else:
            mode = cnt % 5
            if mode == 0:
                x = QQ(random.randint(int(-10**7), int(10**7)))
            elif mode == 1:
                # p-adically close to a root: x = R_i + p^j * c
                i = cnt % 3
                j = ZZ(random.randint(int(-6), int(26)))
                c = ZZ(random.randint(int(1), int(10**4)))
                x = QQ(R[i]) + QQ(c) * QQ(p)**j
            elif mode == 2:
                x = QQ(random.randint(int(-10**5), int(10**5))) / QQ(p)**random.randint(int(1), int(8))
            elif mode == 3:
                x = QQ(random.randint(int(-10**12), int(10**12)))
            else:
                x = QQ(random.randint(int(-10**4), int(10**4))) * QQ(p)**random.randint(int(1), int(8))
        if x in R:
            continue
        y2 = f(x)
        if y2 == 0:
            continue
        if not issq(y2, p):
            continue
        # VERIFIED: (x, sqrt(y2)) is a genuine Q_v-point of E
        t = trip(x, p)
        if t not in img:
            witnesses[t] = x
        img = close(img | {t}, p)
    return img, cnt, tgt, witnesses

print("%-6s %-8s %-8s %-9s %s" % ("place", "|image|", "target", "samples", "delta* in local image?"))
res = {}
for p in [None] + badp:
    img, used, tgt, wit = local_image(p)
    ds = tuple(cls(d, p) for d in dstar)
    inimg = ds in img
    res[p] = (len(img), tgt, inimg, wit)
    lab = "oo" if p is None else str(p)
    mark = "" if len(img) == tgt else "  <<< ORDER NOT REACHED, INCONCLUSIVE"
    print("%-6s %-8d %-8d %-9d %s%s" % (lab, len(img), tgt, used, inimg, mark))

print()
for p in [None] + badp:
    lab = "oo" if p is None else str(p)
    wit = res[p][3]
    if wit:
        print("  witnesses at v=%s:" % lab)
        for t, x in wit.items():
            print("     x =", x, " -> ", t)

print()
determined = all(res[p][0] == res[p][1] for p in res)
allin = all(res[p][2] for p in res)
print("every local image fully determined (order matches theory):", determined)
print("delta* lies in EVERY local image  =>  delta* in Sel^2(E):", allin)
if determined and not allin:
    bad = [("oo" if p is None else p) for p in res if not res[p][2]]
    print()
    print("delta* FAILS to be locally in the image at:", bad)
    print("=> delta* is NOT in Sel^2, so certainly not in the image of E(Q)/2E(Q);")
    print("   C_(19,16) has no finite rational point, by a purely LOCAL obstruction.")
