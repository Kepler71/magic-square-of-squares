# Контроль (б2): прямая сверка Фишера с ctp.sage на E1 сечения (71,49,61) = n61 (полное 2-кручение).
# Известно (ctp.sage): dim Sel^2 = 7, rank CTP = 4  =>  rank E1(k) <= 1.
# Запуск: SEED=1 SEC=n61 sage xcheck_n61.sage
import os, time, traceback
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())

SECTIONS = {'s15': (109, 229, 65), 'n17': (169, 409, 34), 'n61': (3061, 4381, 671),
            'n79': (3217, 5233, 455), 'n89': (2377, 6073, 910)}

name = os.environ.get('SEC', 'n61')
random.seed(int(os.environ.get('SEED', '1')))
A, C, D = SECTIONS[name]
k, rts = E1_roots(A, C, D)
e = [k(t) for t in rts]
E = EllipticCurve(k, [0, -sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], -prod(e)])
Cv = Curve3(k, rts)
Sel = Cv.selmer_full()
n = Sel.dimension()
print(f"=== {name}: A={A} C={C} D={D}, k = {k.defining_polynomial()}, roots = {e}")
print(f"    dim Sel^2(E1) = {n}")
I, J = IJ_of_curve(E)
F = FisherCTP(k, I, J)
t0 = time.time()
dl = deltas_full2(F, E, rts, Cv, Sel)
Mq, quart = ctp_matrix(F, dl, verbose=True, with_diag=True)
print(f"Фишер: rank CTP = {Mq.rank()}, симметрия {Mq == Mq.transpose()}, "
      f"диагональ 0: {all(Mq[i,i] == 0 for i in range(n))}  ({time.time()-t0:.0f}s)")
print(Mq)
print(f"=> rank E1(k) <= {n - 2 - Mq.rank()}")

if os.environ.get('CASSELS', '1') == '1':
    t0 = time.time()
    els = []
    for w in Sel.basis():
        a_ = prod(Cv.gens[i]^int(w[i]) for i in range(Cv.n))
        b_ = prod(Cv.gens[i]^int(w[Cv.n + i]) for i in range(Cv.n))
        els.append((k(a_), k(b_)))
    ct = CTP(Cv); ct.randomize_conic = True; ct.cassels_form = True
    Mc = matrix(GF(2), n, n)
    for i in range(n):
        vals, _ = ct.pair(els[i], els, reps=1)
        for j in range(n):
            Mc[i, j] = vals[j]
    print(f"Касселс (ctp.sage): rank CTP = {Mc.rank()}  ({time.time()-t0:.0f}s)")
    print(Mc)
    print(f"*** МАТРИЦЫ СОВПАЛИ: {Mq == Mc}")
