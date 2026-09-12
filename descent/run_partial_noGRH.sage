load('/home/kep/magicKube/descent/ek_descent_partial.sage')
exec(open('/home/kep/magicKube/descent/e1roots.py').read())
import os
proof.number_field(True)          # БЕЗ GRH
random.seed(int(11))
SECTIONS = {'n79': (3217, 5233, 455), 'n89': (2377, 6073, 910)}
name = os.environ['SEC']
A, C, D = SECTIONS[name]
k, rts = E1_roots(A, C, D); e = [k(t) for t in rts]
for i in range(3):
    ei = e[i]; ej, ek_ = [e[t] for t in range(3) if t != i]
    a = 2*ei - ej - ek_; b = (ei - ej)*(ei - ek_)
    t0 = time.time()
    dim, bd, Cp = partial_descent_bound(k, -2*a, a^2 - 4*b, verbose=False)
    print(f"** {name} i={i+1} (proof=True): dim Sel^2(E') = {dim} -> rank <= {bd}  ({time.time()-t0:.0f}s)")
