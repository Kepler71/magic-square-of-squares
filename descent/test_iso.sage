load('isogeny_descent.sage')
R2 = PolynomialRing(QQ, 't'); t = R2.gen()
def section_pair(A, C, D, matching):
    rootpoly = {'0': t, 'oo': None, '1': t - 1, '-1': t + 1, 'C/A': t - QQ(C)/A, 'A/C': t - QQ(A)/C}
    Fs = []
    for (u, v) in matching:
        fu = rootpoly[u]; fv = rootpoly[v]
        Fs.append(fu if fv is None else (fv if fu is None else fu*fv))
    f0 = D*t*(t^2 - 1)*(A*t - C)*(C*t - A)
    lead = f0 // prod(Fs); Fs[0] = lead*Fs[0]
    assert prod(Fs) == f0
    return Fs
tests = [("s=1/5 (C0,K0)  [rank 1; Magma C0 3, K0 1]", 109, 229, 65, [('0','oo'), ('1','-1'), ('C/A','A/C')]),
         ("(23,7,17) (C0,K10) [rank 1; Magma K10 1]", 169, 409, 34, [('0','C/A'), ('oo','-1'), ('1','A/C')]),
         ("(23,7,17) (C0,K0)  [Magma K0 5]", 169, 409, 34, [('0','oo'), ('1','-1'), ('C/A','A/C')])]
for name, A, C, D, m in tests:
    print("===", name)
    try:
        RichelotPair(section_pair(A, C, D, m)).run()
    except Exception as ex:
        import traceback; traceback.print_exc()
