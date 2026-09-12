# Спуск над этальной алгеброй для 2-изогенных кривых E' = E1/<(e_i,0)> сечений семейства (Claude, 2026-09-12).
# E1: y^2 = prod(x - e_i) над k (полное 2-кручение). Сдвиг x -> x + e_i: y^2 = x(x^2 + a x + b),
#     a = 2e_i - e_j - e_k, b = (e_i - e_j)(e_i - e_k).
# E' = E1/<(0,0)>: Y^2 = X(X^2 - 2a X + (a^2 - 4b)); её Delta' = 16 b — не квадрат (проверено), значит применим
# спуск над L = k x k(sqrt b). rank E'(k) = rank E1(k) (2-изогения).
load('/home/kep/magicKube/descent/ek_descent_partial.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
import os, traceback

SECTIONS = {'s15': (109, 229, 65), 'n17': (169, 409, 34), 'n61': (3061, 4381, 671),
            'n79': (3217, 5233, 455), 'n89': (2377, 6073, 910)}

def run(name, which=(0, 1, 2), seed=1):
    random.seed(seed)
    A, C, D = SECTIONS[name]
    k, rts = E1_roots(A, C, D)
    e = [k(t) for t in rts]
    print(f"=== {name}: A={A} C={C} D={D}, k = {k.defining_polynomial()}, roots = {e}")
    best = None
    for i in which:
        ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
        a = 2*ei - ej - ek_; b = (ei - ej)*(ei - ek_)
        print(f"-- ядро <(e{i+1},0)>: a = {a}, b = {b}, b квадрат: {b.is_square()}")
        t0 = time.time()
        try:
            dim, bd, C_ = partial_descent_bound(k, -2*a, a^2 - 4*b, verbose=True)
            print(f"** {name} i={i+1}: dim Sel^2(E') = {dim} -> rank E1(k) <= {bd}   ({time.time()-t0:.0f}s)")
            best = bd if best is None else min(best, bd)
        except Exception as ex:
            print(f"!! {name} i={i+1}: {type(ex).__name__}: {ex}   ({time.time()-t0:.0f}s)")
            traceback.print_exc()
    print(f"=== ИТОГ {name}: rank E1(k) <= {best}")
    return best

name = os.environ.get('SEC', 's15')
which = [int(t) for t in os.environ.get('WHICH', '0,1,2').split(',')]
run(name, which, seed=int(os.environ.get('SEED', '1')))
