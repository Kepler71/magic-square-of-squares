# CTP на Sel^2(E1/k) для сечения семейства: rank E1(k) <= dim Sel^2 - 2 - rank CTP.  Аргументы: A C D seed
import sys, time
A_, C_, D_, seed = [ZZ(v) for v in sys.argv[1:5]]
sys.argv = ['x']
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
random.seed(int(seed))
t0 = time.time()
k, rts = E1_roots(A_, C_, D_)
Cv = Curve3(k, rts)
Sel = Cv.selmer_full()
els = []
for w in Sel.basis():
    a = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n)); b = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
    els.append((k(a), k(b)))
e1, e2, e3 = Cv.e
tors = [Cv.kappa(e1), Cv.kappa(e2)]
ct = CTP(Cv); ct.randomize_conic = True
n = len(els); M = matrix(GF(2), n, n); tk = []
for i in range(n):
    vals, nP = ct.pair(els[i], els + tors, reps=2)
    for j in range(n): M[i, j] = vals[j]
    tk.append(vals[n:])
    print(f"row {i}: {vals[:n]} torsion {vals[n:]} places {nP} ({time.time()-t0:.0f}s)", flush=True)
print(f"A={A_} C={C_} D={D_} k={k.defining_polynomial()} roots={rts}")
print(M)
print(f"dim Sel2 = {n}; symmetric {M == M.transpose()}; zero diagonal {all(M[i,i] == 0 for i in range(n))}; torsion in kernel {all(all(v == 0 for v in r) for r in tk)}")
print(f"rank CTP = {M.rank()}  =>  rank E1(k) <= {n - 2 - M.rank()}   ({time.time()-t0:.0f}s, seed {seed})")
