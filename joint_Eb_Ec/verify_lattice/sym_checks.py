# Символьные проверки звеньев §3 UNIFIED_THEOREMS (Codex): кручение, удвоение, деление пополам, чётность v_l(x).
from sage.all import *
R = PolynomialRing(QQ, 'x,N'); x, N = R.gens()
E_gen = None
# 1) psi_3 для y^2 = x^3 + A x + B при A=-N^2, B=0
A, B = -N**2, 0
psi3 = 3*x**4 + 6*A*x**2 + 12*B*x - A**2
print("psi3 =", psi3, "| совпадает с 3x^4-6N^2x^2-N^4:", psi3 == 3*x**4 - 6*N**2*x**2 - N**4)
# корни в t = x^2/N^2: 3t^2 - 6t - 1 = 0
T = PolynomialRing(QQ,'t').gen()
f = 3*T**2 - 6*T - 1
print("3t^2-6t-1 неприводим над Q:", f.is_irreducible(), "корни:", f.roots(QQbar))
# 2) x(2P) = (x^2+N^2)^2 / (4 x (x^2-N^2)); проверка через общую формулу удвоения
Nv = 7
E = EllipticCurve([0,0,0,-Nv**2,0])
for P in [E.lift_x(QQ(25)) if E.is_x_coord(25) else None]:
    pass
# символьно: lambda = (3x^2+A)/(2y), x2 = lambda^2 - 2x, y^2 = x^3 + A x
num = (3*x**2 - N**2)**2 - 8*x*(x**3 - N**2*x)
print("x(2P) числитель (x^2+N^2)^2 :", num == (x**2 + N**2)**2)
# x(2P) = +-N  <=> (x^2 -+ 2Nx - N^2)^2 = 0
for s in [1,-1]:
    lhs = (x**2+N**2)**2 - s*N*4*(x**3 - N**2*x)
    rhs = (x**2 - s*2*N*x - N**2)**2
    print("x(2P)=%+dN <=> квадрат:" % s, lhs == rhs)
g = PolynomialRing(QQ,'u').gen(); print("u^2-2u-1 неприводим:", (g**2-2*g-1).is_irreducible())
# 3) деление пополам: alpha^2=x-b, gamma^2=x, beta^2=x+b; Q=((a+g)(g+be), (a+g)(g+be)(a+be)); x(2Q)=x
S = PolynomialRing(QQ, 'al,ga,be'); al, ga, be = S.gens()
xx = ga**2; bb = ga**2 - al**2
rel = [be**2 - (xx + bb)]
I = S.ideal(rel)
u = (al+ga)*(ga+be); v = (al+ga)*(ga+be)*(al+be)
on_curve = (v**2 - (u**3 - bb**2*u))
print("Q на E_b (по модулю beta^2=x+b):", I.reduce(on_curve) == 0)
# x(2Q) = ((u^2+b^2)^2)/(4 v^2) ; проверим (u^2+b^2)^2 - 4 x v^2 == 0 mod I
print("x(2Q) = x:", I.reduce((u**2+bb**2)**2 - 4*xx*v**2) == 0)
# 4) чётность отрицательной нормировки x на y^2=x^3-N^2x и знаменатели -- квадраты: численно на точках
for Nv in [5,6,7,14,34,210]:
    E = EllipticCurve([0,0,0,-Nv**2,0])
    gens = E.gens()
    ok = True
    for G in gens:
        for k in range(1,9):
            d = (k*G)[0].denominator()
            if not d.is_square(): ok = False
    print("N=%d: знаменатели x(kG), k<=8, квадраты:" % Nv, ok, " torsion:", E.torsion_order(), " rank:", len(gens))
