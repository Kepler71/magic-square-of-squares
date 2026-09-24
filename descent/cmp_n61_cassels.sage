# Поэлементная сверка: матрица Фишера (с приёмом, сохранена) против матрицы Касселса (ctp.sage).
import os, time
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
random.seed(int(os.environ.get('SEED', '1')))
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
A, C, D = (3061, 4381, 671)
k, rts = E1_roots(A, C, D)
Cv = Curve3(k, rts); Sel = Cv.selmer_full(); n = Sel.dimension()
els = []
for w in Sel.basis():
    a_ = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n))
    b_ = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
    els.append((k(a_), k(b_)))
ct = CTP(Cv); ct.randomize_conic = True; ct.cassels_form = True
t0 = time.time()
Mc = matrix(GF(2), n, n)
for i in range(n):
    vals, _ = ct.pair(els[i], els, reps=1)
    for j in range(n):
        Mc[i, j] = vals[j]
print(f"Касселс: rank = {Mc.rank()}  ({time.time()-t0:.0f}s)"); print(Mc)
Mq = load('/home/kep/magicKube/descent/Mq_n61_trick.sobj')
print("Фишер (с приёмом):"); print(Mq)
print(f"*** МАТРИЦЫ СОВПАЛИ ПОЭЛЕМЕНТНО: {Mq == Mc};  ранги: Фишер {Mq.rank()}, Касселс {Mc.rank()}")
