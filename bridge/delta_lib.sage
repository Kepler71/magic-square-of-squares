# delta-class machinery for the G1 family, exact rational arithmetic
def sqclass(q):
    """squarefree integer representative of q in Q*/Q*^2; q must be nonzero rational"""
    q = QQ(q)
    assert q != 0
    return ZZ(q.numerator()*q.denominator()).squarefree_part()

def make_E(m,n):
    s = QQ(m^2+n^2)/2
    b = s*m^2*n^2
    e = [-b, -s*m^4, -s*n^4]          # e1,e2,e3 in the stated order
    E4 = [4*x for x in e]             # integral model x=4X, y=8V ; 4 is a square
    R = PolynomialRing(QQ,'x'); x = R.gen()
    f = (x-E4[0])*(x-E4[1])*(x-E4[2])
    c = f.coefficients(sparse=False)
    E = EllipticCurve([0,c[2],0,c[1],c[0]])
    return E, s, b, e, E4

def delta(P, E4):
    """delta class of a point P on the model with roots E4 = [4e1,4e2,4e3]."""
    if P.is_zero():
        return (1,1,1)
    X = P[0]/P[2]; Y = P[1]/P[2]
    out = []
    for i in range(3):
        d = X - E4[i]
        if d == 0:
            j,k = [t for t in range(3) if t != i]
            d = (E4[i]-E4[j])*(E4[i]-E4[k])
        out.append(sqclass(d))
    return tuple(out)
