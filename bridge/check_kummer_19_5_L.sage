#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Строгость численного неравенства L'(E,1) != 0 для (19,5).
import sys
m, n = 19, 5
s = QQ(m**2+n**2)/2; b = s*m**2*n**2
e1, e2, e3 = -b, -s*m**4, -s*n**4
E = EllipticCurve([0, -(e1+e2+e3), 0, e1*e2+e1*e3+e2*e3, -e1*e2*e3])
print("N =", E.conductor(), " w =", E.root_number())
print("ЛОГИКА: L'(E,1) != 0  =>  ord_{s=1} L(E,s) <= 1  =>  (модулярность + Гросс-Загир")
print("+ Колывагин)  rank E(Q) = ord <= 1.  Знак функционального уравнения для этого НЕ нужен.")
print()
for k in [10**5, 5*10**5, 2*10**6]:
    v, err = E.lseries().deriv_at1(k)
    print("k=%-9d L'(E,1) ~ %-18s граница ошибки %-14s |L'|>err ? %s"
          % (k, v, err, bool(abs(v) > err)))
    sys.stdout.flush()
print()
print("L(E,1) =", E.lseries().at1(10**6))
print("analytic_rank(pari) =", E.analytic_rank(algorithm='pari'))
