#!/usr/bin/env sage
# Independent computation of S^(2)(E/Q) for (m,n)=(19,16) as triples of square classes,
# by brute force over all locally-admissible triples.  No reliance on eclib/PARI internals.
import itertools

m,n = 19,16
s = QQ(m^2+n^2)/2
b = s*m^2*n^2
e = [-b, -s*m^4, -s*n^4]
E4 = [ZZ(4*x) for x in e]
R = PolynomialRing(QQ,'x'); x=R.gen()
f = (x-E4[0])*(x-E4[1])*(x-E4[2])
c = f.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]])
print("E =",E); print("roots:",E4)
print("disc =", factor(E.discriminant()))
S = [2,3,5,7,19,617]                      # bad primes
places = ['inf'] + S

def sq_p(q,p):
    """canonical label of q in Q_p^*/(Q_p^*)^2"""
    q = QQ(q); assert q!=0
    a = q.valuation(p); u = q/p^a
    if p==2:
        num=u.numerator(); den=u.denominator()
        uu = (num*inverse_mod(den,8)) % 8
        return (a%2, uu)
    else:
        num=u.numerator(); den=u.denominator()
        uu = (num*inverse_mod(den,p)) % p
        return (a%2, legendre_symbol(uu,p))

def sq_inf(q): return 1 if q>0 else -1

def loc(dv, v):
    if v=='inf': return tuple(sq_inf(t) for t in dv)
    return tuple(sq_p(t,v) for t in dv)

def delta_pt(X):
    """delta triple of a point with x-coordinate X (X not a root)"""
    return [X-E4[0], X-E4[1], X-E4[2]]

def delta_tors(i):
    out=[]
    for k in range(3):
        if k==i:
            j,l=[t for t in range(3) if t!=i]
            out.append((E4[i]-E4[j])*(E4[i]-E4[l]))
        else:
            out.append(E4[i]-E4[k])
    return out

tors_triples = [[1,1,1]] + [delta_tors(i) for i in range(3)]
print(); print("torsion delta triples (as rationals):")
for t in tors_triples: print("   ", t, "->", tuple(ZZ(u).squarefree_part() for u in t))

# ---------- local images ----------
expected = {'inf':2, 2:8}
for p in S:
    if p!=2: expected[p]=4
localimg = {}
for v in places:
    img = set(loc(t,v) for t in tors_triples)
    # sample extra points of E(Q_v)
    tries=0
    bound = 400000
    while len(img) < expected[v] and tries < bound:
        tries += 1
        if v=='inf':
            X = QQ(randint(-3*10^8, 3*10^8))
        else:
            k = randint(-6,6)
            X = QQ(randint(-10^9,10^9))/ (v^max(0,-k)) * v^max(0,k) if False else QQ(randint(-10^9,10^9))*v^k
        if X in E4: continue
        val = (X-E4[0])*(X-E4[1])*(X-E4[2])
        if val==0: continue
        if v=='inf':
            ok = (val>0)
        else:
            ok = (sq_p(val,v)==(0,1) if v!=2 else sq_p(val,2)==(0,1))
        if ok:
            img.add(loc(delta_pt(X),v))
    localimg[v]=img
    print("place %-4s : local image size %d (expected %d)  [%d samples]"%(str(v),len(img),expected[v],tries))

ok_all = all(len(localimg[v])==expected[v] for v in places)
print("all local images complete:", ok_all)

# ---------- enumerate Selmer ----------
sf = [-1]+S
def sqfree_units():
    out=[]
    for bits in itertools.product([0,1],repeat=len(sf)):
        val=1
        for bi,pi in zip(bits,sf):
            if bi: val*=pi
        out.append(ZZ(val))
    return out
U = sqfree_units()
print("candidate d-values:",len(U))
sel=[]
for d1 in U:
    for d2 in U:
        d3 = ZZ(d1*d2).squarefree_part()
        trip=[d1,d2,d3]
        good=True
        for v in places:
            if loc(trip,v) not in localimg[v]: good=False;break
        if good: sel.append(tuple(trip))
print("#S^(2) =", len(sel), "  2-rank =", log(len(sel),2))
target=(1, ZZ(s.numerator()*s.denominator()).squarefree_part(), ZZ(s.numerator()*s.denominator()).squarefree_part())
print("target", target, "in S^(2):", target in set(sel))
print()
print("elements of S^(2) of the shape (1,beta,beta):")
for t in sel:
    if t[0]==1 and t[1]==t[2]: print("   ",t)
print()
print("all 32 Selmer elements:")
for t in sorted(sel): print("   ",t)
