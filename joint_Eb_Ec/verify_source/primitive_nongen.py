# Численная проверка теоремы о примитивных делителях (чётные n > 2) на НЕ-образующих точках:
# G, G+T (T in E[2]\O), 3G, а на кривых ранга 2 -- G1+G2, G1-G2, 2G1+G2. Модель y^2 = x^3 - N^2 x, N бесквадратно.
# Примитивность проверяется без факторизации: из B_n удаляются все простые, делящие какой-либо B_k, k < n.
import sys, json, time
from sage.all import EllipticCurve, QQ, ZZ, gcd, is_squarefree, pari
NMAX = int(sys.argv[1]); NMAXIDX = int(sys.argv[2])
def Bseq(P, nmax):
    B = [None]; Q = P
    for k in range(1, nmax + 1):
        if k > 1: Q = Q + P
        B.append(ZZ(Q[0].denominator()))
    return B
def has_prim(B, n):
    r = B[n]
    for k in range(1, n):
        g = gcd(r, B[k])
        while g > 1:
            r //= g; g = gcd(r, g)
        if r == 1: return False
    return r > 1
fails = []; checked = 0; t0 = time.time(); curves = 0
for N in range(5, NMAX + 1):
    if not is_squarefree(N): continue
    E = EllipticCurve(QQ, [-N*N, 0])
    try:
        gens = E.gens(proof=False)
    except Exception as ex:
        print("N", N, "gens fail", ex); continue
    if len(gens) == 0: continue
    curves += 1
    T = [E(0, 0), E(N, 0), E(-N, 0)]
    pts = []
    for i, G in enumerate(gens):
        pts += [("G%d" % i, G)] + [("G%d+T%d" % (i, j), G + t) for j, t in enumerate(T)] + [("3G%d" % i, 3*G)]
    if len(gens) >= 2:
        G1, G2 = gens[0], gens[1]
        pts += [("G1+G2", G1 + G2), ("G1-G2", G1 - G2), ("2G1+G2", 2*G1 + G2)]
    for name, P in pts:
        nmax = NMAXIDX if not name.startswith(("3G", "2G1")) else max(12, NMAXIDX // 2)
        B = Bseq(P, nmax)
        for n in range(4, nmax + 1, 2):
            checked += 1
            if not has_prim(B, n):
                fails.append((N, name, str(P), n)); print("НЕТ ПРИМИТИВНОГО:", N, name, P, n, flush=True)
print("кривых ранга>=1:", curves, "проверено пар (точка,n):", checked, "нарушений:", len(fails), "%.0fs" % (time.time()-t0))
print("PARI", pari.version())
