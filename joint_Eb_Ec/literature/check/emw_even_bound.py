# Утв. 6: пересчёт финального неравенства EMW (arXiv:math/0409540v2, доказательство Thm 2.2, чётный случай):
#   n^2 (3/4 − ρ(n)) <= 4 (0.621 η(n) + 1.216 ω(n) + 9.0776),  ρ = Σ_{p|n} 1/p^2, η = 2 Σ_{p|n} log p, ω = #простых.
# Утверждение EMW: отсюда n <= 11, т.е. Z_e <= 10 => в шаге 5 хватает при 2M >= 12, M >= 6.
from math import log
from sympy import primefactors
surv = []
for n in range(2, 400, 2):
    ps = primefactors(n)
    rho = sum(1/p**2 for p in ps); eta = 2*sum(log(p) for p in ps); om = len(ps)
    lhs = n*n*(0.75 - rho); rhs = 4*(0.621*eta + 1.216*om + 9.0776)
    if lhs <= rhs: surv.append(n)
    if n in (10, 12): print(f"n={n}: LHS={lhs:.3f}, RHS={rhs:.3f}")
print("чётные n, не исключённые неравенством:", surv)
