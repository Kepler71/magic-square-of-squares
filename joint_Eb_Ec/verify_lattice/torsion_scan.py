# Кручение E_N: y^2 = x^3 - N^2 x для бесквадратных N <= NMAX: ожидается порядок 4 (= E_N[2]).
from sage.all import *
import sys
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
bad = []; cnt = 0
for Nv in range(1, NMAX+1):
    if not Integer(Nv).is_squarefree(): continue
    E = EllipticCurve([0,0,0,-Nv**2,0])
    t = E.torsion_subgroup().invariants()
    cnt += 1
    if tuple(t) != (2,2): bad.append((Nv, t))
print("бесквадратных N <=", NMAX, ":", cnt, "; с кручением != (Z/2)^2:", bad)
