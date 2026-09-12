#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Фаза 3: строгость численного неравенства L'(E,1) != 0 и контроль по формуле BSD.
import sys
m, n = 19, 5
s = QQ(m**2+n**2)/2; b = s*m**2*n**2
e1, e2, e3 = -b, -s*m**4, -s*n**4
E = EllipticCurve([0, -(e1+e2+e3), 0, e1*e2+e1*e3+e2*e3, -e1*e2*e3])
Em = E.minimal_model()
G = E(QQ(-28561963657)/1369, QQ(-2089086742828800)/50653)
print("минимальная модель:", Em.ainvs())
print("N =", E.conductor(), "  w =", E.root_number())

print()
print("--- J.1 оценка погрешности L'(E,1) при росте числа членов ---")
for k in [10**5, 3*10**5, 10**6, 3*10**6]:
    try:
        v, err = E.lseries().deriv_at1(k)
        print("k=%-9d  L'(E,1) ~ %s   граница погрешности %s   нуль исключён: %s"
              % (k, v, err, bool(abs(v) > err)))
    except Exception as ex:
        print("k=%d ОШИБКА: %s" % (k, ex))
    sys.stdout.flush()

print()
print("--- J.2 PARI lfun: значение и производные в s=1 ---")
L = pari(E).lfuncreate()
for d in range(0, 4):
    try:
        val = pari.lfun(L, 1, d)
        print("L^(%d)(E,1) =" % d, val)
    except Exception as ex:
        print("d=%d ОШИБКА: %s" % (d, ex))
sys.stdout.flush()

print()
print("--- J.3 контроль по формуле BSD (ранг 1) ---")
Om = E.period_lattice().omega()
ncomp = 2 if E.discriminant() > 0 else 1
Reg = G.height()
tam = E.tamagawa_product()
tor = E.torsion_order()
Lp = E.lseries().deriv_at1(10**6)[0]
print("Omega (Sage omega(), уже с числом компонент) =", Om)
print("число вещественных компонент =", ncomp)
print("Reg = h(G) =", Reg, "  (G насыщен, индекс 1)")
print("prod c_p =", tam, "   |E(Q)_tors| =", tor)
sha = Lp * tor**2 / (Om * Reg * tam)
print("Ш_an = L'(1)*|T|^2/(Omega*Reg*prod c_p) =", sha)
print("округление:", round(sha), "  отклонение:", abs(sha - round(sha)))
print("ожидание при dim Ш[2]=2 и отсутствии прочего: Ш_an = 4")
print("(если Ш_an = 4.000..., то ранг 1, G — генератор, и Ш[2] порядка 4 —")
print(" то есть требуемый класс (1,193,193) действительно сидит в Ш[2])")

print()
print("--- J.4 попытка улучшить границу eclib большими параметрами ---")
sys.stdout.flush()
try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    ec = mwrank_EllipticCurve(list(Em.ainvs()))
    ec.set_verbose(0)
    ec.two_descent(second_limit=18, n_aux=-1, second_descent=True)
    print("eclib (минимальная модель, second_limit=18): rank =", ec.rank(),
          " certain =", ec.certain(), " rank_bound =", ec.rank_bound(),
          " selmer_rank =", ec.selmer_rank())
except Exception as ex:
    print("eclib ОШИБКА:", ex)
print("вывод: eclib делает только 2-спуск -> его верхняя граница 3 непреодолима")
print("без 4-спуска или спаривания Касселса. Это НЕ независимое подтверждение rank<=1.")
