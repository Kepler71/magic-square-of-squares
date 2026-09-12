# -*- coding: utf-8 -*-
# Раунд 5: НЕЗАВИСИМОЕ (не-спусковое) подтверждение rank E(Q) = 1
# через Гросс–Загира–Колывагина: w(E) = -1 и L'(E,1) != 0.
# В раунде 2 я взял слишком мало членов ряда (100) и оценка ошибки
# оказалась бесполезной (7625 при значении 7.93). Берём нужное число.

import sys
from sage.all import *

S = QQ(113); M = ZZ(15); N = ZZ(1); B = S*M**2*N**2
e = [QQ(-B), QQ(-S*M**4), QQ(-S*N**4)]
Rx = PolynomialRing(QQ,'x'); x = Rx.gen()
f = (x-e[0])*(x-e[1])*(x-e[2]); c = f.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]]).minimal_model()
Ncond = E.conductor()
print("  кондуктор N = %s,  sqrt(N) ~ %.1f" % (Ncond, float(sqrt(Ncond))))
print("  корневое число w(E) = %s" % E.root_number())
print("  w = -1  =>  L(E,1) = 0 точно (функциональное уравнение)")
sys.stdout.flush()

for k in [30000, 100000, 300000]:
    try:
        val, err = E.lseries().deriv_at1(k)
        print("  k=%-8d  L'(E,1) ~ %s   оценка ошибки <= %s   ненулевость: %s"
              % (k, val, err, bool(abs(val) > 10*abs(err))))
        sys.stdout.flush()
    except Exception as ex:
        print("  k=%d: ошибка %s" % (k, ex))

print("\n  Если L'(E,1) != 0 с надёжной оценкой ошибки, то аналитический ранг = 1,")
print("  и по теореме Колывагина (с Гросс–Загиром) rank E(Q) = 1, Sha конечна.")
print("  Это НЕ использует второй 2-спуск и спаривание Касселса–Тейта.")
print("  СТАТУС: численно (оценка ошибки Кремоны, реализация Sage), но метод независим.")
print("ГОТОВО")
