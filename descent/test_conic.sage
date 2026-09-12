load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
random.seed(int(2))
k.<r> = QuadraticField(165)
Cv = Curve3(k, [2939651, 163724*r, -2939651])
ct = CTP(Cv)
Sel = Cv.selmer_full()
import time
for w in Sel.basis()[:3]:
    a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n)); b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
    t0 = time.time()
    L, eps3 = ct.tangent_forms(k(a), k(b))
    print("tangent forms ok", {j: len(str(L[j])) for j in L}, f"{time.time()-t0:.1f}s")
