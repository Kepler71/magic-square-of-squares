# Независимая разведка: модель, I,J, phi, классы 2-накрытий PARI. (11,4)
from sage.all import *
m, n = 11, 4
s = QQ(137)/2
b = s*m^2*n^2
e = [-b, -s*m^4, -s*n^4]
print("b =", b, " e =", e)
# целочисленная депрессированная модель: r_i = 4 e_i + shift, sum r = 0
sh = -sum(4*ei for ei in e)/3
print("shift =", sh, " integer:", sh in ZZ)
r = [4*ei + sh for ei in e]
print("roots r =", r, "sum =", sum(r))
R.<X> = QQ[]
f = prod(X - ri for ri in r)
A = f[1]; B = f[0]
print("A,B =", A, B)
M = EllipticCurve([0,0,0,A,B])
print("M =", M, " disc =", factor(M.discriminant()))
print("minimal? ", M.is_minimal(), M.minimal_model())
I = -48*A; J = -1728*B
print("I,J =", I, J)
print("c4,c6 =", M.c4(), M.c6(), " I==c4:", I==M.c4(), " J==2c6:", J==2*M.c6())
phi = [-12*ri for ri in r]
print("phi =", phi)
print("check phi cubic:", [ph^3 - 3*I*ph + J for ph in phi])
print("root diffs:", [(r[0]-r[1]).factor(), (r[0]-r[2]).factor(), (r[1]-r[2]).factor()])
# 2-накрытия PARI
cov = pari(M).ell2cover()
print("ell2cover count:", len(cov))
def quartI(g):
    a,bb,c,d,ee = g
    return 12*a*ee-3*bb*d+c*c
def quartJ(g):
    a,bb,c,d,ee = g
    return 72*a*c*ee-27*a*d*d-27*bb*bb*ee+9*bb*c*d-2*c^3
def zinv(g, ph):
    a,bb,c,d,ee = g
    return (4*a*ph + 3*bb*bb - 8*a*c)/3
def sqcl(q):
    q = QQ(q)
    if q == 0: return 0
    sgn = 1 if q>0 else -1
    num = q.numerator()*q.denominator()
    return sgn*prod([p for p,k in factor(abs(num)) if k%2==1])
Rx.<x> = QQ[]
for i,c in enumerate(cov):
    q = Rx(c[0])
    g = [QQ(q[4-j]) for j in range(5)]
    Ig, Jg = quartI(g), quartJ(g)
    sc = QQ(J*Ig/(Jg*I))
    print("cover",i,"raw g =", g, " I,J ratio:", Ig/I, Jg/J, " scale:", sc, " square:", sc.is_square())
    g = [sc*a for a in g]
    assert quartI(g)==I and quartJ(g)==J, (quartI(g), quartJ(g))
    zs = [zinv(g, ph) for ph in phi]
    print("   g =", g)
    print("   z =", zs, " classes =", [sqcl(z) for z in zs])
