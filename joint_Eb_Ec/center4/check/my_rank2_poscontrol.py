# Положительный контроль детектора связей: трёхпростая решётка (5,13,17), показатели <= 2 — должна найтись тройка Саллоуса
# (центр 1105, u = 264, 576, 1104: в нормировке Im w = u/1105).
from fractions import Fraction
def gmul(a, b): return (a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0])
def gpow(a, k):
    r = (1, 0)
    for _ in range(k): r = gmul(r, a)
    return r
pis = {5: (1, 2), 13: (3, 2), 17: (1, 4)}
E = 2; classes = {}
for r in range(-E, E+1):
    for s in range(-E, E+1):
        for t in range(-E, E+1):
            if r == s == t == 0: continue
            rho = (1, 0)
            for (p, k) in ((5, r), (13, s), (17, t)):
                g = pis[p] if k >= 0 else (pis[p][0], -pis[p][1])
                rho = gmul(rho, gpow(g, abs(k)))
            N = rho[0]**2 + rho[1]**2; sq = gmul(rho, rho)
            for eps in range(2):
                v = Fraction(abs(sq[1] if eps == 0 else sq[0]), N)
                if v == 0 or v == 1: continue
                classes.setdefault((eps, abs(r) % 2, abs(s) % 2, abs(t) % 2), set()).add(v)
hits = []
for cls, S in classes.items():
    L = sorted(S)
    for i in range(len(L)):
        for j in range(i+1, len(L)):
            x, y = L[i], L[j]
            for c in (x+y, y-x, 2*y-x, 2*x-y, x+2*y, y+2*x, y-2*x, (x+y)/2, (y-x)/2):
                if c > 0 and c != x and c != y and c in S: hits.append((cls, x, y, c))
target = {Fraction(264, 1105), Fraction(576, 1105), Fraction(1104, 1105)}
print('hits:', len(hits), '; тройка Саллоуса найдена:', any({x, y, c} == target for _, x, y, c in hits))
