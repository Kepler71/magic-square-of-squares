# Сверка: та же матрица в форме Касселса [Cas98, Lemma 7.4] — три конуса, без нормировки и без перестановки.
import sys, time
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
random.seed(int(sys.argv[1]) if len(sys.argv) > 1 else int(21))
k.<r> = QuadraticField(165)
rts = [2939651, 163724*r, -2939651]
Cv = Curve3(k, rts); Sel = Cv.selmer_full()
els = []
for w in Sel.basis():
    a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n)); b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
    els.append((k(a), k(b)))
n = len(els)
for mode in (False, True):
    ct = CTP(Cv); ct.randomize_conic = True; ct.cassels_form = mode
    M = matrix(GF(2), n, n); t0 = time.time()
    for i in range(n):
        vals, nP = ct.pair(els[i], els, reps=1)
        for j in range(n): M[i, j] = vals[j]
    print(f"{'Cassels' if mode else 'наша'} форма: rank {M.rank()}, симметрия {M == M.transpose()}, диагональ {all(M[i,i]==0 for i in range(n))} ({time.time()-t0:.0f}s)")
    print(M)
    if mode: M2 = M
    else: M1 = M
print("матрицы совпадают:", M1 == M2)
