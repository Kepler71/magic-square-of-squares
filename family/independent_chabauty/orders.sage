load('/home/kep/magicKube/family/independent_chabauty/setup.sage')
import sys
b, h, n, D, u0 = sys.argv[1:6]
S = Section(ZZ(b), ZZ(h), ZZ(n), ZZ(D)); E = S.E; G = S.point_from_u(QQ(u0))
for p in primes(5, 120):
    if kronecker(S.dk, p) != 1: continue
    if E.discriminant().norm() % p == 0: continue
    Fp = GF(p); r = Fp(S.dk).sqrt()
    ords = []
    for sgn in (1, -1):
        red = lambda a: Fp(S.k(a).list()[0]) + sgn*r*Fp(S.k(a).list()[1])
        Eb = EllipticCurve(Fp, [red(a) for a in E.ainvs()])
        ords.append((Eb.cardinality(), Eb([red(G[0]), red(G[1])]).order()))
    print(f"p={p}: |E(F_P)|={[o[0] for o in ords]}, ord(G)={[o[1] for o in ords]}, N={lcm([o[1] for o in ords])}")
