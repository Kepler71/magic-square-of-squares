print("="*78); print("1. G1: nine cells, magic, pair types")
R.<m,n,t> = QQ[]
s = (m^2+n^2)/2
F0 = m^2+n^2*t^2; F4 = s*(1+t^2); F8 = n^2+m^2*t^2
L  = F4-2*m*n*t;  U  = F4+2*m*n*t
M = [[F0,(m*t+n)^2,L],[(m*t-n)^2,F4,(m+n*t)^2],[U,(m-n*t)^2,F8]]
names = ['a','b','c','d','e','f','g','h','i']
flat = [M[0][0],M[0][1],M[0][2],M[1][0],M[1][1],M[1][2],M[2][0],M[2][1],M[2][2]]
LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
sums=set(sum(flat[i] for i in Ln) for Ln in LINES)
print("  all 8 line sums equal 3*centre :", len(sums)==1 and sums.pop()==3*F4)
for (u,v),lab in [((0,8),'(a,i) CORNER'),((2,6),'(c,g) CORNER'),((1,7),'(b,h) EDGE'),((3,5),'(d,f) EDGE')]:
    auto = flat[u].is_square() and flat[v].is_square()
    print(f"  pair {names[u]}{names[v]} {lab:14s} sum-2e ok = {flat[u]+flat[v]-2*F4==0}   both cells automatic squares = {auto}")
print("  => the curve C: u0^2=F0,u4^2=F4,u8^2=F8 is cells (a,e,i) = MAIN DIAGONAL")
print("     F0 = cell a (corner), F4 = cell e (centre), F8 = cell i (corner); free/unconstrained: c=L, g=U (the OTHER corner pair)")

print("="*78); print("2. structural identities: corner = mean of two edge cells; edge = corner+corner-centre")
e_,X_,Y_ = var('e_ X_ Y_')
A_=[e_+X_, e_-X_-Y_, e_+Y_, e_-X_+Y_, e_, e_+X_-Y_, e_-Y_, e_+X_+Y_, e_-X_]  # a..i
print("  a == (h+f)/2 :", bool((A_[0]-(A_[7]+A_[5])/2).simplify_full()==0))
print("  i == (b+d)/2 :", bool((A_[8]-(A_[1]+A_[3])/2).simplify_full()==0))
print("  c == (h+d)/2 :", bool((A_[2]-(A_[7]+A_[3])/2).simplify_full()==0))
print("  g == (b+f)/2 :", bool((A_[6]-(A_[1]+A_[5])/2).simplify_full()==0))
print("  h == a+c-e   :", bool((A_[7]-(A_[0]+A_[2]-A_[4])).simplify_full()==0))
print("  b == i+g-e   :", bool((A_[1]-(A_[8]+A_[6]-A_[4])).simplify_full()==0))
print("  (no corner is a mean of two corners, no edge is a mean of two edges - checked below)")
cells_idx={'corner':[0,2,6,8],'edge':[1,3,5,7]}
from itertools import combinations
for typ in ('corner','edge'):
    bad=[]
    for u,v in combinations(cells_idx[typ],2):
        for w in range(9):
            if w in (u,v): continue
            if bool((A_[w]-(A_[u]+A_[v])/2).simplify_full()==0): bad.append((names[u],names[v],names[w]))
    print(f"  means of two {typ} cells that equal another cell:", bad)

print("="*78); print("3. sections family (edge-fixed) vs corner-fixed analogue -- symbolic")
Rt.<T> = QQ[]
var('b0 h0 n0')
def build(kind, b0,h0,n0):
    # centre 1 (scale n0^2 later); fixed pair cells bS^2,hS^2 ; second pair of the SAME type parametrised by T
    bS,hS = b0/n0, h0/n0
    Bf = (1+2*T-T^2)/(1+T^2); Hf=(1-2*T-T^2)/(1+T^2)
    if kind=='edge':      # fixed EDGE pair (b,h); other EDGE pair (d,f) = Bf^2,Hf^2 ; free cells = 4 CORNERS
        free = [(hS^2+Hf^2)/2,(hS^2+Bf^2)/2,(bS^2+Hf^2)/2,(bS^2+Bf^2)/2]
    else:                 # fixed CORNER pair (a,i); other CORNER pair (c,g)=Bf^2,Hf^2 ; free cells = 4 EDGES
        free = [bS^2+Bf^2-1, hS^2+Hf^2-1, hS^2+Bf^2-1, bS^2+Hf^2-1]
    return [f.simplify_full() for f in free]
for kind in ('edge','corner'):
    fr = build(kind,b0,h0,n0)
    print(f"  --- {kind}-fixed: four free cells x n0^2 (1+T^2)^2 ---")
    for f in fr:
        g = (f*n0^2*(1+T^2)^2).simplify_full().expand()
        print("     ", g.factor() if False else g)
