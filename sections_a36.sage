# Сечения с a = 36 (среди них все минимальные): подсистема над Q, минимальные простые.
import functools, time
print = functools.partial(print, flush=True)
Rs.<s> = QQ[]
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
print("s^6 coeff of X^3+A4X+A6 at a=36:", 36^3 + A4[4]*36 + A6[6])
B.<b,c,d,e,f,g> = PolynomialRing(QQ, order='degrevlex')
Bt.<t> = B[]
X = 36*t^2 + b*t + c; Y = d*t^3 + e*t^2 + f*t + g
J = B.ideal((Y^2 - (X^3 + Bt(A4.change_ring(B))(t)*X + Bt(A6.change_ring(B))(t))).coefficients())
t0 = time.time()
print("dim", J.dimension(), " length", J.vector_space_dimension(), f"({time.time()-t0:.1f}s)")
t0 = time.time()
mp = J.minimal_associated_primes()
print("minimal primes:", len(mp), f"({time.time()-t0:.1f}s)")
for Pj in sorted(mp, key=lambda Z: Z.vector_space_dimension()):
    print("---- orbit size", Pj.vector_space_dimension())
    for gen in Pj.groebner_basis():
        print("     ", gen)
