# -*- coding: utf-8 -*-
# Независимая перепроверка базовых фактов о (15,8), на которые опирается вывод.
from sage.all import *
import json
m, n = 15, 8
s = QQ(m*m + n*n)/2
print("s =", s, " m^2+n^2 =", m*m+n*n, " квадрат?", is_square(m*m+n*n))

# 1. E_B: Y^2 = X(X-(m^2-n^2)^2)(X-(m^2+n^2)^2) — универсальный фактор Codex
EB = EllipticCurve([0, -((m*m-n*n)**2 + (m*m+n*n)**2), 0, (m*m-n*n)**2*(m*m+n*n)**2, 0])
print("E_B ainvs:", EB.ainvs())
print("E_B torsion:", EB.torsion_subgroup().invariants())
print("E_B pari ellrank:", EB.pari_curve().ellrank())
print("E_B eclib rank_bound:", EB.rank_bound())
try:
    print("E_B analytic_rank:", EB.analytic_rank())
except Exception as e:
    print("analytic_rank err", e)

# 2. Девять клеток как функции t; независимый вывод всех рациональных вырождений
R = PolynomialRing(QQ, 't'); t = R.gen()
F0 = m*m + n*n*t**2
F4 = s*(1 + t**2)
F8 = n*n + m*m*t**2
cells = [F0, (m*t+n)**2, F4 - 2*m*n*t,
         (m*t-n)**2, F4, (m+n*t)**2,
         F4 + 2*m*n*t, (m-n*t)**2, F8]
# контроль магичности
rows = [cells[0:3], cells[3:6], cells[6:9]]
sums = [sum(r) for r in rows] + [sum(rows[i][j] for i in range(3)) for j in range(3)] \
       + [cells[0]+cells[4]+cells[8], cells[2]+cells[4]+cells[6]]
print("магический при всех t:", all(x == 3*F4 for x in sums))

deg_t = set()
why = {}
for i, c in enumerate(cells):
    for r in c.roots(QQ, multiplicities=False):
        deg_t.add(r); why.setdefault(r, []).append("клетка %d = 0" % i)
for i in range(9):
    for j in range(i+1, 9):
        d = cells[i] - cells[j]
        if d == 0:
            why.setdefault('тождество', []).append("клетки %d,%d равны тождественно" % (i, j))
            continue
        for r in d.roots(QQ, multiplicities=False):
            deg_t.add(r); why.setdefault(r, []).append("клетки %d=%d" % (i, j))
print("степени разностей <= 2:",
      all((cells[i]-cells[j]).degree() <= 2 for i in range(9) for j in range(i+1,9) if cells[i] != cells[j]))
print("ВСЕ рациональные вырожденные t:", sorted(deg_t))
print("число:", len(deg_t))
# какие из них дают девять рациональных квадратов
for r in sorted(deg_t):
    vals = [c(r) for c in cells]
    if all(v.is_square() for v in vals):
        print("   t =", r, "-> девять квадратов:", vals, " различных:", len(set(vals)))

# 3. бесконечность
print("на бесконечности ведущие коэффициенты:", [c.leading_coefficient() for c in cells],
      " квадратность центра:", (F4.leading_coefficient()).is_square())

# 4. контроль t=1 — пять форм
L = F4 - 2*m*n*t; U = F4 + 2*m*n*t
print("t=1 пять форм:", [F0(1), F4(1), F8(1), L(1), U(1)],
      "все квадраты:", all(QQ(x).is_square() for x in [F0(1), F4(1), F8(1), L(1), U(1)]))
