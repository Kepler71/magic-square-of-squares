# Fable review: (1) list of l for the 1/l theorem, PARI ranks; (2) for random slopes, each of the 13
# cell classes: Jacobian (PARI ellfromeqn) isomorphic over Q to claimed E_{m,n}? and obvious degenerate
# points on the cell curve count == torsion order of E_{m,n}.
import random, json, sys
random.seed(int(913))

def Emn(m, n):
    return EllipticCurve(QQ, [0, m*m + n*n, 0, m*m*n*n, 0])

# ---------- (1) ----------
def cond(l):
    if l % 8 not in (3, 5):
        return False
    return all(q % 4 == 3 for q, e in factor(l*l - 1) if q != 2)

ls = [l for l in prime_range(3, 3001) if cond(l)]
print("ell list (<=3000):", ls, "count", len(ls))
claimed = [3, 5, 13, 37, 43, 197, 277, 283, 397, 557, 643, 683, 757, 797, 827, 907, 947, 997, 1237, 1453,
           1693, 1867, 1933, 1987, 2267, 2357, 2477, 2557, 2917]
print("matches claimed list:", ls == claimed)
bad = []
for l in ls:
    E = Emn(1, l)
    rk = pari(E).ellrank()
    if not (int(rk[0]) == 0 and int(rk[1]) == 0) or E.torsion_order() != 8:
        bad.append((l, rk, E.torsion_order()))
print("ell with PARI rank not proven 0 or torsion != 8:", bad)

# ---------- (2) ----------
R.<p> = PolynomialRing(QQ)

def cell_classes(k):
    # 5 cubic triples {a,b,a+b} and 8 quartic classes, with the claimed pair (as (m,n) integers after scaling)
    r, s = k.numerator(), k.denominator()
    cls = []
    cls.append(([1, k, 1+k], (r, s)))
    cls.append(([k, 1-k, 1], (r, s-r)))
    cls.append(([1, -(1-k), k], (s, s-r)))
    cls.append(([1, -(1+k), -k], (s, r+s)))
    cls.append(([1+k, -k, 1], (r, r+s)))
    cls.append(([1, -1, k, -k], (s-r, r+s)))
    cls.append(([1+k, -(1+k), 1-k, -(1-k)], (r, s)))
    cls.append(([k, -k, 1-k, -(1-k)], (s, abs(2*r-s))))
    cls.append(([k, -k, 1+k, -(1+k)], (s, 2*r+s)))
    cls.append(([1, -1, 1-k, -(1-k)], (r, 2*s-r)))
    cls.append(([1, -1, 1+k, -(1+k)], (r, 2*s+r)))
    cls.append(([k, 1+k, 1-k, -k], (s, 2*r)))
    cls.append(([1, 1+k, -1, -(1-k)], (r, 2*s)))
    return cls

allc = lambda k: [1, -1, k, -k, 1+k, -(1+k), 1-k, -(1-k)]

def jacobian(cs):
    f = prod(1 + c*p for c in cs)
    # PARI ellfromeqn needs polynomial in x,y: y^2 - f(x)
    x = pari('x'); y = pari('y')
    fp = pari(f)(x)
    ell = pari.ellfromeqn(y**2 - fp)
    return EllipticCurve(QQ, [QQ(a) for a in ell])

def obvious_points(cs, k):
    f = prod(1 + c*p for c in cs)
    pts = []          # (p, kind)
    pts.append((0, 'p=0'))
    pts.append((0, 'p=0'))
    for c in cs:
        pts.append((-1/c, 'root'))
    lead = prod(cs)
    if len(cs) == 3:
        pts += [(oo, 'inf')]
    elif lead.is_square():
        pts += [(oo, 'inf'), (oo, 'inf')]
    for c in allc(k):
        v = f(-2/c)
        if v != 0 and v.is_square():
            pts += [(-2/c, 'p=-2/%s' % c)] * 2
    return pts

def degenerate(pv, k):
    if pv == 0 or pv == oo:
        return True
    cells = [1 + c*pv for c in allc(k)]
    # not all squares  <=> degenerate for our purpose; also flag the reason
    return not all(QQ(c).is_square() for c in cells)

slopes_all = [(r, s) for s in range(3, 201) for r in range(1, s) if gcd(r, s) == 1 and 2*r != s]
sample = random.sample(slopes_all, int(120))
sample += [(1, l) for l in ls if l <= 200] + [(l-1, l) for l in ls if l <= 200]
iso_fail = []; count_fail = []; nondeg = []; checked = 0
for r, s in sample:
    k = QQ(r)/QQ(s)
    for cs, (m, n) in cell_classes(k):
        g = gcd(m, n); m, n = m//g, n//g
        if m == n or m == 0 or n == 0:
            continue
        E = Emn(m, n)
        J = jacobian(cs)
        checked += 1
        if not J.is_isomorphic(E):
            iso_fail.append((r, s, cs, (m, n)))
            continue
        pts = obvious_points(cs, k)
        T = E.torsion_order()
        if len(pts) != T:
            count_fail.append((r, s, cs, (m, n), len(pts), T))
        for pv, kind in pts:
            if not degenerate(pv, k):
                nondeg.append((r, s, cs, pv))
print("classes checked:", checked, "slopes:", len(sample))
print("isomorphism failures:", iso_fail[:10], len(iso_fail))
print("obvious-point count != torsion order:", count_fail[:20], len(count_fail))
from collections import Counter
print("by class index / torsion:", Counter((cell_classes(QQ(1)/QQ(3)).__len__(), t) for *_, t in count_fail))
print("non-degenerate obvious points:", nondeg[:10], len(nondeg))

# ---------- (3) closure check tied to own Selmer bound ----------
sys.path.insert(0, '/home/kep/magicKube/criterion_proof/fable_review')
import fr_selmer
viol = []; closed_classes = 0; slopes_closed = set(); extra_rank = []
for r, s in sample:
    k = QQ(r)/QQ(s)
    for cs, (m, n) in cell_classes(k):
        g = gcd(m, n); m, n = int(m//g), int(n//g)
        if m == n or m == 0 or n == 0:
            continue
        if m > n: m, n = n, m
        if fr_selmer.exact_bound(m, n) != 0:
            continue
        E = Emn(m, n); T = E.torsion_order()
        pts = obvious_points(cs, k)
        if len(pts) != T or any(not degenerate(pv, k) for pv, kind in pts):
            viol.append((r, s, cs, (m, n), len(pts), T))
        else:
            closed_classes += 1; slopes_closed.add((r, s))
print("closed classes (own bound 0, count == torsion, all degenerate):", closed_classes, "; slopes closed:", len(slopes_closed), "of", len(sample))
print("violations:", viol)
for r, s, cs, (m, n), c, T in count_fail:
    rk = pari(Emn(m, n)).ellrank()
    extra_rank.append((m, n, int(rk[0]), int(rk[1]), fr_selmer.exact_bound(int(m), int(n))))
print("pairs with extra points: (m,n,PARI lower,upper,own bound):", sorted(set(extra_rank)))
print("any extra-point pair with own bound 0:", any(b == 0 for *_, b in extra_rank))
