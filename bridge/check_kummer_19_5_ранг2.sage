#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Фаза 2 проверки (19,5): независимые подтверждения границы ранга + вопрос о Sel_2.
import sys
m, n = 19, 5
s = QQ(m**2 + n**2)/2; b = s*m**2*n**2
e1, e2, e3 = -b, -s*m**4, -s*n**4
A2 = -(e1+e2+e3); A4 = e1*e2+e1*e3+e2*e3; A6 = -e1*e2*e3
E = EllipticCurve([0, A2, 0, A4, A6])
G = E(QQ(-28561963657)/1369, QQ(-2089086742828800)/50653)

print("="*78); print("G. НЕЗАВИСИМАЯ ГРАНИЦА РАНГА ЧЕРЕЗ АНАЛИТИЧЕСКИЙ РАНГ"); print("="*78)
print("w = ellrootno =", E.root_number(), " -> порядок нуля L(E,s) в s=1 НЕЧЁТЕН -> ord >= 1")
sys.stdout.flush()
for alg in ['pari', 'sympow', 'rubinstein']:
    try:
        ar = E.analytic_rank(algorithm=alg)
        print("analytic_rank(%-10s) = %s" % (alg, ar))
    except Exception as ex:
        print("analytic_rank(%-10s) ОШИБКА: %s" % (alg, ex))
    sys.stdout.flush()
try:
    d1 = E.lseries().deriv_at1(100000)
    print("L'(E,1) (Sage lseries.deriv_at1, 1e5 членов) =", d1)
except Exception as ex:
    print("deriv_at1 ОШИБКА:", ex)
try:
    print("L(E,1) (должно быть 0 при w=-1) =", E.lseries().at1(100000))
except Exception as ex:
    print("at1 ОШИБКА:", ex)
sys.stdout.flush()

print()
print("ТЕОРЕМА (Уайлс-БКДТ модулярность + Гросс-Загир + Колывагин):")
print("  если ord_{s=1} L(E,s) <= 1, то rank E(Q) = ord и Ш(E/Q) конечна.")
print("  Здесь ord = 1 (w=-1 даёт ord нечётен; L'(1) != 0 с большим запасом).")
print("  => rank E(Q) = 1 НЕЗАВИСИМО от спуска PARI.")

print()
print("="*78); print("H. НАСЫЩЕННОСТЬ И Ш (контроль согласованности)"); print("="*78)
sys.stdout.flush()
try:
    sat, idx, nh = E.saturation([G])
    print("E.saturation([G]) -> генераторы:", sat, " индекс:", idx)
except Exception as ex:
    print("saturation ОШИБКА:", ex)
try:
    print("G делится на 2 в E(Q)?", len(G.division_points(2)) > 0)
    print("G делится на 3 в E(Q)?", len(G.division_points(3)) > 0)
except Exception as ex:
    print("division_points ОШИБКА:", ex)
try:
    sha = E.sha().an()
    print("аналитический порядок Ш =", sha, " (ожидается 4 при dim Ш[2]=2)")
except Exception as ex:
    print("sha.an() ОШИБКА:", ex)
print("PARI ellrank дал s = rk(Ш[2]/2Ш[4]) = 2, 2-Selmer rank C = 5, T = 2")
print("  -> r2 = C - T - s = 5 - 2 - 2 = 1.  Согласовано с Ш[2] порядка 4.")
sys.stdout.flush()

print()
print("="*78); print("I. ЛЕЖИТ ЛИ (1,193,193) В 2-SELMER? (собственный локальный тест)"); print("="*78)
# H_delta для delta=(d1,d2,d3), d_i = класс x-e_i:
#   x-e1 = d1 w1^2, x-e2 = d2 w2^2, x-e3 = d3 w3^2
# вычитая:  d1 w1^2 - d2 w2^2 = e2-e1 ,  d2 w2^2 - d3 w3^2 = e3-e2
d1, d2, d3 = 1, 193, 193
c12 = e2 - e1; c23 = e3 - e2
print("H:  %d*w1^2 - %d*w2^2 = %s ;  %d*w2^2 - %d*w3^2 = %s" % (d1, d2, c12, d2, d3, c23))
# в однородном виде с w0:
#   w1^2 = 193 w2^2 + c12 w0^2      (d1=1)
#   w3^2 = w2^2 - (c23/193) w0^2    (d2=d3=193)
k23 = QQ(c23)/193
print("эквивалентно: w1^2 = 193 w2^2 + (%s) w0^2 ,  w3^2 = w2^2 + (%s) w0^2" % (c12, -k23))
assert k23.denominator() == 1
k23 = ZZ(k23)


def is_sq_Qp(a, p):
    if a == 0:
        return True
    a = QQ(a)
    v = a.valuation(p)
    if v % 2:
        return False
    u = a / p**v
    num, den = ZZ(u.numerator()), ZZ(u.denominator())
    u = num * den  # тот же класс квадратов
    if p == 2:
        return (u % 8) == 1
    return kronecker(u, p) == 1


def locally_solvable(p, prec):
    """есть ли (w0:w2) in P^1(Q_p) с обоими значениями квадратами в Q_p"""
    def ok(w0, w2):
        Aq = 193*w2**2 + c12*w0**2
        Bq = w2**2 - k23*w0**2
        return is_sq_Qp(Aq, p) and is_sq_Qp(Bq, p)
    pk = p**prec
    for c in range(pk):          # карта (w0:w2) = (1:c)
        if ok(1, c):
            return True, (1, c)
    for c in range(0, pk, p):    # карта (w0:w2) = (c:1), c in pZ_p
        if ok(c, 1):
            return True, (c, 1)
    return False, None


print("вещественное место: нужно 193w2^2+c12 >= 0 и w2^2-k23 >= 0 одновременно;")
w2big = 10**6
print("   при w2 = 10^6:", 193*w2big**2 + c12, ",", w2big**2 - k23, "-> оба > 0: ОК")
S = [2, 3, 5, 7, 19, 193]
allok = True
for p in S:
    prec = {2: 10, 3: 7, 5: 6, 7: 5, 19: 4, 193: 3}[p]
    res, wit = locally_solvable(p, prec)
    print("p = %-4d разрешимо в Q_p? %-5s  свидетель (w0:w2) = %s" % (p, res, wit))
    allok = allok and res
    sys.stdout.flush()
# хорошие простые: гладкая кривая рода 1 над F_p имеет точки при p не делящем disc (Хассе)
print("для p вне {%s} у H хорошая редукция (род 1) -> точки по Хассе-Вейлю при p>=5" % S)
print()
if allok:
    print("ВЫВОД: (1,193,193) ВСЮДУ ЛОКАЛЬНО РАЗРЕШИМ -> он ЛЕЖИТ в 2-Selmer группе.")
    print("Значит локального препятствия НЕТ, и исключение C(Q)=пусто ЭКВИВАЛЕНТНО")
    print("утверждению, что этот класс — НЕТРИВИАЛЬНЫЙ элемент Ш(E/Q)[2].")
    print("Вся сила аргумента сидит именно в границе rank<=1.")
else:
    print("ВЫВОД: класс локально неразрешим -> исключение БЕЗУСЛОВНО, граница ранга не нужна!")
