# Контроль 2: сверка спуска над алгеброй (ek_descent_partial) с изогенным спуском (ek_descent) — Claude, 2026-09-12.
# phi: E -> E' = E/<(0,0)>, dual phi^: E' -> E, [2]_{E'} = phi o phi^.
# Точная последовательность  Sel^{phi^}(E') -> Sel^2(E') -> Sel^{phi}(E),  ядро первой стрелки = E(k)[phi]/phi^(E'(k)[2]) = Z/2
# (E'(k)[2] = {O,(0,0)} и phi^((0,0)) = O, так как Delta' = 16b не квадрат) =>
#     dim Sel^2(E')  <=  dim Sel^{phi^}(E') + dim Sel^{phi}(E) - 1.
# Обе границы ранга должны быть согласованы: dim Sel^2(E') - 1  vs  Sel^phi + Sel^phi^ - 2.
load('/home/kep/magicKube/descent/ek_descent.sage')
load('/home/kep/magicKube/descent/ek_descent_partial.sage')      # переопределяет LocalSq на обобщённый (совместим)
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
import os, traceback

SECTIONS = {'s15': (109, 229, 65), 'n17': (169, 409, 34), 'n61': (3061, 4381, 671),
            'n79': (3217, 5233, 455), 'n89': (2377, 6073, 910)}

name = os.environ.get('SEC', 's15')
which = [int(t) for t in os.environ.get('WHICH', '0,1,2').split(',')]
random.seed(int(os.environ.get('SEED', '1')))
A, C, D = SECTIONS[name]
k, rts = E1_roots(A, C, D)
e = [k(t) for t in rts]
Cv = Curve3(k, rts)
Sf = Cv.selmer_full()
print(f"=== {name}: full 2-descent на E1: dim Sel^2 = {Sf.dimension()} -> rank <= {Sf.dimension()-2}")
for i in which:
    ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
    a = 2*ei - ej - ek_; b = (ei - ej)*(ei - ek_)
    t0 = time.time()
    s1, s2 = Cv.isogeny_descent_dual(i)
    print(f"-- i={i+1}: Sel^phi(E) = {s1}, Sel^phi^(E') = {s2} -> rank <= {s1+s2-2}  ({time.time()-t0:.0f}s)")
    # контроль изогении в Sage
    E = EllipticCurve(k, [0, a, 0, b, 0]); Ep = EllipticCurve(k, [0, -2*a, 0, a^2 - 4*b, 0])
    ok_iso = E.isogeny(E(0, 0)).codomain().is_isomorphic(Ep)
    try:
        dim, bd, Cp = partial_descent_bound(k, -2*a, a^2 - 4*b, verbose=False)
        print(f"   dim Sel^2(E') = {dim} -> rank <= {bd};  изогения E -> E' проверена: {ok_iso};"
              f"  неравенство dim Sel^2(E') <= {s1}+{s2}-1 = {s1+s2-1}: {'OK' if dim <= s1+s2-1 else 'НАРУШЕНО!!!'}"
              f"  ({time.time()-t0:.0f}s)")
        assert ok_iso and dim <= s1 + s2 - 1
    except Exception as ex:
        print(f"   !! {type(ex).__name__}: {ex}")
        traceback.print_exc()
