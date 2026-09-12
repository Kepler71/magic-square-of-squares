# Инволюции P^1, переставляющие 6 корней f1 парами (без неподвижных корней) -> бисэллиптичность C1.
import functools, itertools
print = functools.partial(print, flush=True)
K.<i> = QuadraticField(-1)
R.<t> = K[]
f1 = 65*(t^2 - 1)*(t^2 + 1)*(32161*t^2 - 49922*t + 32161)
roots = [r for r, _ in f1.roots()]
assert len(roots) == 6
def matchings(lst):
    if not lst: yield []; return
    a = lst[0]
    for j in range(1, len(lst)):
        rest = lst[1:j] + lst[j+1:]
        for m in matchings(rest):
            yield [(a, lst[j])] + m
found = []
for m in matchings(roots):
    (a, b), (c, d), (e, f) = m
    # M(z) = (al z + be)/(ga z - al): M(a) = b  <=>  al(a+b) + be - ga a b = 0
    A = matrix(K, [[a + b, 1, -a*b], [c + d, 1, -c*d]])
    ker = A.right_kernel().basis()
    if len(ker) != 1: continue
    al, be, ga = ker[0]
    if al^2 + be*ga == 0: continue                       # вырожденное преобразование
    Mf = lambda z: (al*z + be)/(ga*z - al)
    if Mf(e) == f:
        # нормируем коэффициенты и смотрим поле определения
        v = vector(K, [al, be, ga]); piv = next(x for x in v if x != 0); v = v/piv
        rat = all(x in QQ for x in v)
        found.append((m, v, rat))
print("involutions permuting roots in 3 pairs:", len(found))
for m, v, rat in found:
    print("   (alpha:beta:gamma) =", v, " defined over Q:", rat)
