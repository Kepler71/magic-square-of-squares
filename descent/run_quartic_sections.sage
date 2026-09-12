# Боевой счёт: CTP по Фишеру на Sel^2(E'_i) для 2-изогенных кривых E'_i = E1/<(e_i,0)> сечений семейства.
# rank E1(k) = rank E'_i(k) <= dim Sel^2(E'_i) - 1 - rank CTP.
# Запуск:  SEC=n89 WHICH=0,1,2 SEED=1 sage run_quartic_sections.sage
import os, time, traceback
load('/home/kep/magicKube/descent/ek_descent_partial.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
proof.number_field(False)

SECTIONS = {'s15': (109, 229, 65), 'n17': (169, 409, 34), 'n61': (3061, 4381, 671),
            'n79': (3217, 5233, 455), 'n89': (2377, 6073, 910)}

name = os.environ.get('SEC', 'n89')
which = [int(t) for t in os.environ.get('WHICH', '0,1,2').split(',')]
random.seed(int(os.environ.get('SEED', '1')))
A, C, D = SECTIONS[name]
k, rts = E1_roots(A, C, D)
e = [k(t) for t in rts]
print(f"=== {name}: A={A} C={C} D={D}, k = {k.defining_polynomial()}, roots = {e}")
best = None
for i in which:
    ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
    a = 2*ei - ej - ek_; b = (ei - ej)*(ei - ek_)
    ap, bp = -2*a, a^2 - 4*b            # E'_i: Y^2 = X(X^2 + ap X + bp)
    print(f"-- i={i+1}: ядро <(e{i+1},0)>, E'_i: a' = {ap}, b' = {bp}")
    t0 = time.time()
    try:
        dim, bd, PD = partial_descent_bound(k, ap, bp, verbose=True)
        Ep = EllipticCurve(k, [0, k(ap), 0, k(bp), 0])
        I, J = IJ_of_curve(Ep)
        F = FisherCTP(k, I, J); F.conic_verbose = True; F.conic_timeout = int(os.environ.get("CTO", "180"))
        F.base_S = sorted(set(list(PD.S) + [P for q in primes(2, 60) for P in k.primes_above(q)]), key=lambda P: (P.norm(), str(P)))
        dl = deltas_partial(F, Ep, PD, PD.Sel)
        M, quart = ctp_matrix(F, dl, verbose=True, with_diag=True)
        n = M.nrows()
        sym = (M == M.transpose()); dg0 = all(M[t, t] == 0 for t in range(n))
        rk = M.rank()
        print(M)
        print(f"** {name} i={i+1}: dim Sel^2(E') = {dim}, rank CTP = {rk} (симм {sym}, диаг0 {dg0})"
              f"  =>  rank E1(k) <= {dim - 1 - rk}   ({time.time()-t0:.0f}s)")
        bnd = dim - 1 - rk
        best = bnd if best is None else min(best, bnd)
    except Exception as ex:
        print(f"!! {name} i={i+1}: {type(ex).__name__}: {ex}   ({time.time()-t0:.0f}s)")
        traceback.print_exc()
print(f"=== ИТОГ {name}: rank E1(k) <= {best}")
