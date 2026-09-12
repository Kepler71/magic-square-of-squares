from itertools import permutations
from fractions import Fraction

# cells 0..8 = a b c d e f g h i
# magic space basis: e-vector, X-vector (corner pair a,i step), Y-vector (corner pair c,g step)
E = [1,1,1,1,1,1,1,1,1]
X = [1,-1,0,-1,0,1,0,1,-1]
Y = [0,-1,1,1,0,-1,-1,1,0]
BASIS=[E,X,Y]

LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
def is_magic(v):
    s=sum(v[i] for i in LINES[0])
    return all(sum(v[i] for i in L)==s for L in LINES)

# sanity: basis vectors are magic
assert all(is_magic(b) for b in BASIS), "basis not magic"

def apply(sigma, v):
    # new[sigma[i]] = v[i]
    w=[0]*9
    for i in range(9): w[sigma[i]]=v[i]
    return w

def in_span(w):
    # solve w = c0*E + c1*X + c2*Y over Q by least squares / gaussian
    import itertools
    # build augmented matrix 9x4
    rows=[[Fraction(E[i]),Fraction(X[i]),Fraction(Y[i]),Fraction(w[i])] for i in range(9)]
    # gaussian elimination
    piv=0; pivcols=[]
    for c in range(3):
        r=None
        for rr in range(piv,9):
            if rows[rr][c]!=0: r=rr;break
        if r is None: continue
        rows[piv],rows[r]=rows[r],rows[piv]
        pv=rows[piv][c]
        rows[piv]=[x/pv for x in rows[piv]]
        for rr in range(9):
            if rr!=piv and rows[rr][c]!=0:
                f=rows[rr][c]
                rows[rr]=[a-f*b for a,b in zip(rows[rr],rows[piv])]
        pivcols.append(c); piv+=1
    for rr in range(piv,9):
        if rows[rr][3]!=0 and all(rows[rr][c]==0 for c in range(3)): return False
    return True

good=[]
for sigma in permutations(range(9)):
    if all(in_span(apply(sigma,b)) for b in BASIS):
        good.append(sigma)
print("permutations of the 9 cells that map EVERY magic square to a magic square:", len(good))

PAIRS={'corner(a,i)':(0,8),'corner(c,g)':(2,6),'edge(b,h)':(1,7),'edge(d,f)':(3,5)}
def pairtype(p):
    p=tuple(sorted(p))
    for k,v in PAIRS.items():
        if tuple(sorted(v))==p: return k
    return None

crossings=0
for sigma in good:
    img={}
    for name,(u,v) in PAIRS.items():
        t=pairtype((sigma[u],sigma[v]))
        assert t is not None, (sigma,name)
        img[name]=t
        if ('corner' in name) != ('corner' in t): crossings+=1
print("images of pairs: any edge<->corner crossing?", crossings>0, " (count of crossing images:",crossings,")")

# print the group action table compactly
seen=set()
for sigma in good:
    key=tuple(pairtype((sigma[u],sigma[v])) for (u,v) in PAIRS.values())
    seen.add(key)
print("distinct induced actions on the 4 pairs:", len(seen))
for k in sorted(seen): print("   ", k)
# does sigma fix the centre?
print("all fix the centre cell e:", all(s[4]==4 for s in good))
print("all map {corners}->{corners}:", all(sorted(s[i] for i in (0,2,6,8))==[0,2,6,8] for s in good))
