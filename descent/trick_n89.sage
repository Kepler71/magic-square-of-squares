# Приём Фишера на боевых квартиках n89 (E'_1, k = Q(sqrt 462), dim Sel^2 = 8).
# Эталон — матрица из CTP_QUARTIC_2026-09-12.md §4: M[1,2]=M[1,6]=M[1,7]=1, остальные 0.
import os, time
load('/home/kep/magicKube/descent/ek_descent_partial.sage')
load('/home/kep/magicKube/descent/ctp_quartic.sage')
load('/home/kep/magicKube/descent/ctp_trick.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
proof.number_field(False)
A, C, D = (2377, 6073, 910)
k, rts = E1_roots(A, C, D)
e = [k(t) for t in rts]
i = 0
ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
a = 2*ei - ej - ek_; b = (ei - ej)*(ei - ek_)
ap, bp = -2*a, a^2 - 4*b
Ep = EllipticCurve(k, [0, k(ap), 0, k(bp), 0])
I, J = IJ_of_curve(Ep)
F = FisherCTP(k, I, J); F.conic_timeout = 180
F.fact_cache_path = '/home/kep/magicKube/descent/fact_n89_1.sobj'
F.base_S = sorted(set([P for q in primes(2, 60) for P in k.primes_above(q)]), key=lambda P: (P.norm(), str(P)))
quart = load('/home/kep/magicKube/descent/quart_n89_1.sobj')
print('квартик в кэше:', len(quart), flush=True)
ref = {(1,2): 1, (1,6): 1, (1,7): 1, (0,1): 0, (2,3): 0, (0,2): 0}
pairs = [tuple(int(t) for t in p.split(',')) for p in os.environ.get('PAIRS', '1,2;1,6;2,3').split(';')]
NCHK = int(os.environ.get('NCHK', '2'))
for (i1, j1) in pairs:
    g1 = quart[(i1,)]; g2 = quart[(j1,)]; g3 = quart[(min(i1,j1), max(i1,j1))]
    t0 = time.time()
    try:
        v, used = pair_trick(F, g1, g2, g3, bound=int(os.environ.get('BND','10')), verbose=True, ncheck=NCHK)
        vr, _ = pair_trick(F, g2, g1, g3, bound=int(os.environ.get('BND','10')), verbose=False, ncheck=1)
        dt = time.time() - t0
        r = ref.get((min(i1,j1), max(i1,j1)))
        sym = 'симметрия OK' if v == vr else f'НЕСИММЕТРИЧНО (обратно {vr})'
        print(f"<g{i1},g{j1}> = {v}   эталон {r}   {'СОВПАЛО' if r is None or v == r else 'РАСХОЖДЕНИЕ'}   {sym}   {dt:.0f}s", flush=True)
    except Exception as ex:
        import traceback; traceback.print_exc()
        print(f"<g{i1},g{j1}>: ОШИБКА {type(ex).__name__}: {str(ex)[:200]}", flush=True)
