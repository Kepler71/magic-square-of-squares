# Построение H0(k,u) по рецепту HANDOFF.md §3 и контроль заявленных свойств.
# Результат пишется в h0.txt (строка многочлена) для последующей работы в Sage.
from sympy import *

x, y, z, t, k, u = symbols('x y z t k u')
A0 = -(x**2+2*x*z-1)*(x**2*z-2*x-z)
A1 = -(x**2*y-x**2*z-2*x*y*z-2*x-y+z)*(x**2*y*z+x**2+2*x*y-2*x*z-y*z-1)
A2 = (x+y)*(x*y-1)
A3 = -(x**2*y+x**2*z-2*x*y*z+2*x-y-z)*(x**2*y*z-x**2+2*x*y+2*x*z-y*z+1)
F = expand((1+x**2)*(1+z**2)*A2 - (1+y**2)*A0 - A1)
G = expand(A3 - (1+y**2)*A0 + A1)
P = Poly(expand(F.subs({x: t, z: k*t})), y)
print("deg_y F on ray:", P.degree())
Au = expand(P.coeff_monomial(y**2)/t).subs(t**2, u)
Cu = expand(P.coeff_monomial(1)/t).subs(t**2, u)
Bu = P.coeff_monomial(y).subs(t**2, u)
for name, e in (("Au", Au), ("Bu", Bu), ("Cu", Cu)):
    assert not expand(e).has(t), f"{name} still contains t: substitution t^2->u incomplete"
N = (2*k**2-3*k)*u**2 + (-2*k**2+18*k-2)*u + (2-3*k)
D = (6*k**2-k)*u**2 + (-6*k**2+6*k-6)*u + (6-k)
H = expand(u*(Au*N+Cu*D)**2 - Bu**2*N*D)
H0 = expand(cancel(H/u))
H0 = expand(H0/gcd_list(Poly(H0, k, u).coeffs()))

PH = Poly(H0, k, u)
print("deg_k, deg_u:", PH.degree(k), PH.degree(u), " total:", PH.total_degree())

print("\n[1] involution H0(k,u) = k^8 u^8 H0(1/k,1/u):",
      expand(k**8*u**8*H0.subs({k: 1/k, u: 1/u}, simultaneous=True) - H0) == 0)
print("    N:", expand(k**2*u**2*N.subs({k: 1/k, u: 1/u}, simultaneous=True) - N) == 0,
      " D:", expand(k**2*u**2*D.subs({k: 1/k, u: 1/u}, simultaneous=True) - D) == 0)

print("\n[2] H0(1,u) =", factor(H0.subs(k, 1)))
print("    claimed 4(u-3)^2(u+1)^4(3u-1)^2:",
      expand(H0.subs(k, 1) - 4*(u-3)**2*(u+1)**4*(3*u-1)**2) == 0)
print("[3] H0(0,u) =", factor(H0.subs(k, 0)))
print("    claimed -4(u-1)^2(3u-1)(u^2-22u+25):",
      expand(H0.subs(k, 0) + 4*(u-1)**2*(3*u-1)*(u**2-22*u+25)) == 0)

lc = Poly(H0, u).LC()
print("\n[4] leading coeff in u:", factor(lc))

a, b = symbols('a b')
L = Poly(expand(H0.subs({k: 1+a, u: -1+b}, simultaneous=True)), a, b)
m = min(sum(mon) for mon in L.monoms())
cone = sum(c*a**i*b**j for (i, j), c in L.terms() if i+j == m)
print("\n[5] multiplicity at (1,-1):", m, " tangent cone:", factor(cone))

print("\n[6] factor over Q (may take a while)...")
fl = factor_list(H0)
print("    content:", fl[0], " factors (deg_k,deg_u,mult):",
      [(Poly(f, k, u).degree(k), Poly(f, k, u).degree(u), e) for f, e in fl[1]])

with open("h0.txt", "w") as fh:
    fh.write(str(H0) + "\n")
print("\nwritten h0.txt, terms:", len(PH.terms()))
