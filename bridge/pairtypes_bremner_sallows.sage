print("="*78); print("4. Bremner-Sallows example: which pairs are complete, and of which type")
sq=[373^2,289^2,565^2,360721,425^2,23^2,205^2,527^2,222121]
names=['a','b','c','d','e','f','g','h','i']
LINES=[(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
print("  magic:", len(set(sum(sq[i] for i in L) for L in LINES))==1, " sum =", sum(sq[i] for i in LINES[0]), "= 3*e:", sum(sq[i] for i in LINES[0])==3*sq[4])
print("  squares:", {names[i]:ZZ(sq[i]).is_square() for i in range(9)})
for (u,v),lab in [((0,8),'CORNER (a,i)'),((2,6),'CORNER (c,g)'),((1,7),'EDGE (b,h)'),((3,5),'EDGE (d,f)')]:
    comp = ZZ(sq[u]).is_square() and ZZ(sq[v]).is_square()
    print(f"   pair {names[u]}{names[v]} {lab:14s}: sum=2e {sq[u]+sq[v]==2*sq[4]}, COMPLETE={comp}, step={abs(sq[u]-sq[4])}")
print("  => two complete pairs: one EDGE (b,h) and one CORNER (c,g)  ==> MIXED type = G3, not G1 and not G2")
print("  ratios of the complete pairs (cell roots : centre root):")
print("     EDGE   (b,h) = (289 : 527 : 425) -> primitive", tuple(x//gcd([289,527,425]) for x in (289,527,425)))
print("     CORNER (c,g) = (565 : 205 : 425) -> primitive", tuple(x//gcd([565,205,425]) for x in (565,205,425)))
print("  check 2n^2=P^2+Q^2 :", 289^2+527^2==2*425^2, 565^2+205^2==2*425^2)

print("="*78); print("5. steps of the four pairs, and the relation among any THREE complete pairs")
R.<Xv,Yv>=QQ[]
steps={'corner(a,i)':Xv,'corner(c,g)':Yv,'edge(b,h)':Xv+Yv,'edge(d,f)':Xv-Yv}
from itertools import combinations
for miss in steps:
    three=[k for k in steps if k!=miss]
    v=[steps[k] for k in three]
    rel=None
    for i in range(3):
        j,k=[x for x in range(3) if x!=i]
        for sg in (1,-1):
            if v[i]==v[j]+sg*v[k]: rel=f"{three[i]} = {three[j]} {'+' if sg>0 else '-'} {three[k]}"
        for sg in (1,-1):
            if 2*v[i]==v[j]+sg*v[k]: rel=(rel or "")+f" | 2*{three[i]} = {three[j]} {'+' if sg>0 else '-'} {three[k]}"
    print(f"  missing {miss:14s} -> relation among the other three: {rel}")

print("="*78); print("6. sections: edge-fixed (proved) vs corner-fixed (untouched) -- constants, field k, degeneration")
secs=[(17,7,13),(7,1,5),(23,7,17),(31,17,25),(41,1,29),(47,23,37),(49,31,41),(73,17,53),(71,49,61),(89,23,65),(79,47,65)]
def sqf(x):
    x=ZZ(x); return sign(x)*prod(p^(e%2) for p,e in x.abs().factor())
print(f"  {'section':14s} {'A_edge':>7s} {'C_edge':>7s} {'k_edge':>10s} | {'A_cor':>6s} {'C_cor':>6s} {'k_corner':>10s}  AC_cor square?")
for (b0,h0,n0) in secs:
    A=(h0^2+n0^2)/2; C=(b0^2+n0^2)/2
    dk=sqf((C-A)*(C+A))
    Ac=h0^2; Cc=b0^2; dkc=sqf((Cc-Ac)*(Cc+Ac))
    print(f"  ({b0},{h0},{n0})".ljust(16)+f"{A:7d} {C:7d}  Q(sqrt{dk:4d}) | {Ac:6d} {Cc:6d}  Q(sqrt{dkc:5d})   {(Ac*Cc).is_square()}   (A*C edge square? {(A*C).is_square()})")
print("  note k_edge = Q(sqrt(b0^2-h0^2)) and k_corner = Q(sqrt(2(b0^2-h0^2))) -- always different (ratio 2):")
print("   ", all(sqf((( (b0^2+n0^2)/2 - (h0^2+n0^2)/2 ))*(( (b0^2+n0^2)/2 + (h0^2+n0^2)/2 ))) == sqf(b0^2-h0^2) and sqf((b0^2-h0^2)*(b0^2+h0^2))==sqf(2*(b0^2-h0^2)) for (b0,h0,n0) in secs))
