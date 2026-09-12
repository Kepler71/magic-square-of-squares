# -*- coding: utf-8 -*-
# Explicit full 2-Selmer group of E : y^2 = (x-r1)(x-r2)(x-r3) as triples of square
# classes, via the local images  L_v = delta(E(Q_v))  in  (Q_v*/Q_v*^2)^3.
# |L_v| = |E(Q_v)[2]| / |2|_v  = 4 (p odd), 8 (p=2), 2 (v=infinity).

import itertools, sys

def loc_class(a, p):
    """square class of a nonzero rational in Q_p*/Q_p*^2, as a hashable tag."""
    a = QQ(a)
    v = a.valuation(p)
    u = a / QQ(p)**v
    t = ZZ(u.numerator()*u.denominator())
    if p == 2:
        return (v % 2, t % 8)
    return (v % 2, kronecker(t, p))

def loc_triple(c, p):
    return tuple(loc_class(x, p) for x in c)

def real_triple(c):
    return tuple(1 if x > 0 else -1 for x in c)

def sqfree(x):
    x = QQ(x)
    return ZZ(x.numerator()*x.denominator()).squarefree_part()

def delta_of(x, r):
    return tuple(sqfree(x - ri) for ri in r)

def torsion_deltas(r):
    out = [(ZZ(1), ZZ(1), ZZ(1))]
    for i in range(3):
        j, k = [t for t in range(3) if t != i]
        c = [0, 0, 0]
        c[i] = sqfree((r[i]-r[j])*(r[i]-r[k]))
        c[j] = sqfree(r[i]-r[j])
        c[k] = sqfree(r[i]-r[k])
        out.append(tuple(c))
    return out

def grp_close(elts, mulf):
    S = set(elts)
    changed = True
    while changed:
        changed = False
        for a in list(S):
            for b in list(S):
                c = mulf(a, b)
                if c not in S:
                    S.add(c); changed = True
    return S

def local_image_p(r, p, want, nsample=40000):
    """L_p as a set of tags (triples of local square classes)."""
    mul = lambda a, b: tuple((( (x[0]+y[0]) % 2, (x[1]*y[1]) % (8 if p==2 else 10**9) if p==2
                               else x[1]*y[1]) ) for x, y in zip(a, b))
    def mul2(a, b):
        out = []
        for x, y in zip(a, b):
            if p == 2:
                out.append(((x[0]+y[0]) % 2, (x[1]*y[1]) % 8))
            else:
                out.append(((x[0]+y[0]) % 2, x[1]*y[1]))
        return tuple(out)
    L = set([loc_triple(t, p) for t in torsion_deltas(r)])
    L = grp_close(L, mul2)
    if len(L) >= want:
        return L, True
    # sample points of E(Q_p): x in Q_p with f(x) a square
    f = lambda x: (x-r[0])*(x-r[1])*(x-r[2])
    def is_sq(a):
        if a == 0: return False
        v = a.valuation(p)
        if v % 2: return False
        u = a/QQ(p)**v; t = ZZ(u.numerator()*u.denominator())
        return (t % 8 == 1) if p == 2 else kronecker(t, p) == 1
    cnt = 0
    for k in range(-6, 9):
        for a in range(-300, 301):
            if a == 0: continue
            x = QQ(a) * QQ(p)**k
            cnt += 1
            if cnt > nsample: break
            fx = f(x)
            if fx != 0 and is_sq(fx):
                tg = loc_triple(delta_of(x, r), p)
                if tg not in L:
                    L = grp_close(L | set([tg]), mul2)
                    if len(L) >= want:
                        return L, True
    return L, (len(L) >= want)

def local_image_R(r):
    L = set([real_triple(t) for t in torsion_deltas(r)])
    L = grp_close(L, lambda a, b: tuple(x*y for x, y in zip(a, b)))
    f = lambda x: (x-r[0])*(x-r[1])*(x-r[2])
    for x in [RR(q) for q in srange(-2*10**9, 2*10**9, 10**7)]:
        if f(x) > 0:
            tg = tuple(1 if (x-ri) > 0 else -1 for ri in r)
            if tg not in L:
                L = grp_close(L | set([tg]), lambda a, b: tuple(u*v for u, v in zip(a, b)))
    return L

# ------------------------------------------------------------------ curve
m, n = 19, 16
s = QQ(m**2+n**2)/2; b = s*m**2*n**2
r = [ZZ(-4*b), ZZ(-4*s*m**4), ZZ(-4*s*n**4)]
print("r =", r)
E = EllipticCurve(QQ, [0, -(r[0]+r[1]+r[2]), 0,
                       r[0]*r[1]+r[0]*r[2]+r[1]*r[2], -r[0]*r[1]*r[2]])
bad = sorted(set([q for q, _ in factor(2*E.discriminant())]))
print("bad primes:", bad)
target = (ZZ(1), sqfree(s), sqfree(s))
print("target class:", target)

Lv = {}
for p in bad:
    want = 8 if p == 2 else 4
    L, ok = local_image_p(r, p, want)
    Lv[p] = L
    print("  L_%d : |L| = %d (expected %d)  %s" % (p, len(L), want, "OK" if ok else "*** SHORT ***"))
    if not ok:
        raise RuntimeError("local image at %s incomplete" % p)
LR = local_image_R(r)
print("  L_R : |L| = %d (expected 2)" % len(LR))

# ------------------------------------------------------------- Selmer group
prim = bad[:]                      # support of the square classes
gens = [ZZ(-1)] + [ZZ(q) for q in prim]
allcls = []
for sub in Subsets(range(len(gens))):
    v = prod([gens[i] for i in sub], ZZ(1))
    allcls.append(ZZ(v).squarefree_part())
allcls = sorted(set(allcls))
print("  #square classes supported on {-1}u%s = %d" % (prim, len(allcls)))

sel = []
for c1 in allcls:
    for c2 in allcls:
        c3 = sqfree(c1*c2)
        c = (c1, c2, c3)
        ok = True
        if real_triple(c) not in LR:
            continue
        for p in bad:
            if loc_triple(c, p) not in Lv[p]:
                ok = False; break
        if ok:
            sel.append(c)
print("\n|Sel_2(E/Q)| =", len(sel), " = 2^%s" % ZZ(len(sel)).exact_log(2),
      "  (Sage selmer_rank says %s)" % E.selmer_rank())
tors = torsion_deltas(r)
print("torsion classes in Sel_2 ?", all(t in sel for t in tors))
print("TARGET (1,%s,%s) in Sel_2 ?  %s" % (target[1], target[2], target in sel))

# the four elements with prescribed second coordinate
for c2v in [2]:
    cand = [c for c in sel if c[1] == c2v]
    print("\nelements of Sel_2 with c2 = %s : %s" % (c2v, cand))
tgt_norm = [c for c in sel if c in [tuple(sqfree(a*t) for a, t in zip(target, tt)) for tt in tors]]
print("\ntarget * torsion (the 4 classes that would make (19,16) survive):")
for tt in tors:
    cc = tuple(sqfree(a*t) for a, t in zip(target, tt))
    print("   ", cc, "  in Sel_2 ?", cc in sel)
