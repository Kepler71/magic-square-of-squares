# Сверка Фишер против Касселса на n61 (dim Sel^2 = 7) — с приёмом Фишера для выбора точки.
# Раньше сходились 2 строки из 7 (упиралось в факторизацию 46-значных чисел).
import os, time
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ctp.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
load('/home/kep/magicKube/descent/ctp_trick.sage')
random.seed(int(os.environ.get('SEED', '1')))
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
A, C, D = (3061, 4381, 671)      # n61 = сечение (71,49,61)
k, rts = E1_roots(A, C, D)
e = [k(t) for t in rts]
E = EllipticCurve(k, [0, -sum(e), 0, e[0]*e[1] + e[0]*e[2] + e[1]*e[2], -prod(e)])
Cv = Curve3(k, rts); Sel = Cv.selmer_full(); n = Sel.dimension()
print(f"n61: k = {k.defining_polynomial()}, dim Sel^2(E1) = {n}", flush=True)
I, J = IJ_of_curve(E)
F = FisherCTP(k, I, J); F.conic_timeout = int(os.environ.get('CTO', '180'))
F.fact_cache_path = '/home/kep/magicKube/descent/fact_n61_trick.sobj'
install_trick(F, bound=int(os.environ.get('BND', '10')), ncheck=int(os.environ.get('NCHK', '1')), verbose=True)
t0 = time.time()
dl = deltas_full2(F, E, rts, Cv, Sel)
Mq, quart, sym_ok = ctp_matrix(F, dl, verbose=True, with_diag=True,
                               cache='/home/kep/magicKube/descent/quart_n61_trick.sobj')
print(f"Фишер (с приёмом): rank CTP = {Mq.rank()}, симметрия {Mq == Mq.transpose()}, "
      f"диагональ 0: {all(Mq[i,i] == 0 for i in range(n))}  ({time.time()-t0:.0f}s)", flush=True)
print(Mq)
print(f"=> rank E1(k) <= {n - 2 - Mq.rank()}", flush=True)
save(Mq, '/home/kep/magicKube/descent/Mq_n61_trick.sobj')
