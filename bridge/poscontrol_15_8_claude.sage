# ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ: та же машинерия на (15,8), где точка C ИЗВЕСТНА (t=1).
def sqclass(x):
    x = QQ(x)
    z = x.numerator() * x.denominator()
    sgn = 1 if z > 0 else -1
    out = Integer(1)
    for p, e in Integer(abs(z)).factor():
        if e % 2:
            out *= p
    return Integer(sgn) * out

def build(M, N):
    S = QQ(M**2 + N**2) / 2
    B = S * M**2 * N**2
    E1 = QQ(-B); E2 = QQ(-S * M**4); E3 = QQ(-S * N**4)
    Rx = PolynomialRing(QQ, 'X'); X = Rx.gen()
    cub = (X - E1) * (X - E2) * (X - E3)
    co = cub.coefficients(sparse=False)
    E = EllipticCurve([0, co[2], 0, co[1], co[0]])
    return S, B, (E1, E2, E3), E

def delta(E, es, P):
    if P == E(0):
        return (Integer(1), Integer(1), Integer(1))
    x = P.xy()[0]
    out = []
    for i in range(3):
        d = x - es[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (es[i] - es[j]) * (es[i] - es[k])
        out.append(sqclass(d))
    return tuple(out)

for (M, N) in [(15, 8), (19, 5)]:
    M = Integer(M); N = Integer(N)
    S, B, es, E = build(M, N)
    tgt = (sqclass((M * N)**2), sqclass(S * M**2), sqclass(S * N**2))
    print("\n(m,n)=(%s,%s): s=%s b=%s  цель=%s" % (M, N, S, B, tgt))
    hits = []
    for q in range(1, 60):
        for p in range(0, 60):
            if gcd(p, q) != 1:
                continue
            A = M**2 * q**2 + N**2 * p**2
            C8 = N**2 * q**2 + M**2 * p**2
            F4n = S * (q**2 + p**2)
            if Integer(A).is_square() and Integer(C8).is_square() and QQ(F4n).is_square():
                hits.append(QQ(p) / q)
    print("  рациональные t с тремя квадратами (|p|,|q|<60):", hits)
    for tv in hits:
        X = B * tv**2
        u0 = QQ(M**2 + N**2 * tv**2).sqrt()
        u4 = QQ(S * (1 + tv**2)).sqrt()
        u8 = QQ(N**2 + M**2 * tv**2).sqrt()
        V = B * u0 * u4 * u8
        P = E(X, V)
        d = delta(E, es, P)
        print("    t=%s : u0=%s u4=%s u8=%s" % (tv, u0, u4, u8))
        print("           P=(%s, %s) лежит на E: %s" % (X, V, P in E))
        print("           delta(P)=%s ; равен цели: %s" % (d, d == tgt))
