# Является ли модель y^2 = x^3 - N^2 x глобально минимальной (важно для совпадения B_n с «минимальной моделью» VY/Verzobio)?
from sage.all import *
bad = []
for Nv in range(1, 2001):
    if not Integer(Nv).is_squarefree(): continue
    E = EllipticCurve([0,0,0,-Nv**2,0])
    if not E.is_global_minimal_model(): bad.append(Nv)
print("бесквадратные N <= 2000, для которых y^2=x^3-N^2x НЕ глобально минимальна:", bad[:20], "всего", len(bad))
