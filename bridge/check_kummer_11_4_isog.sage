# Часть 6: НЕЗАВИСИМАЯ верхняя граница ранга.
# Идея: ранг постоянен в классе изогении. Если у какой-то изогенной кривой
# eclib/mwrank даёт 2-Selmer ранг, дающий границу 1, получаем АЛГЕБРАИЧЕСКУЮ
# границу НЕ от PARI ellrank (т.е. независимо от вычисления спаривания Касселса в PARI).
def hdr(t):
    print("\n" + "="*78); print(t); print("="*78)

M, N = ZZ(11), ZZ(4)
S = QQ(M^2+N^2)/2; B = S*M^2*N^2
Rx.<XX> = PolynomialRing(QQ)
cub = (XX+B)*(XX+S*M^4)*(XX+S*N^4); c = cub.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]]).minimal_model()
print("Emin =", E, " проводник", factor(E.conductor()))

hdr("I1. КЛАСС ИЗОГЕНИИ: ранг одинаков у всех; ищем лучшую eclib-границу")
IC = E.isogeny_class()
print("размер класса изогении:", len(IC.curves))
print("матрица изогений:")
print(IC.matrix())
best = None
for i, C in enumerate(IC.curves):
    try:
        mc = C.mwrank_curve()
        sr = mc.selmer_rank(); rb = mc.rank_bound(); cert = mc.certain(); rk = mc.rank()
        tor2 = len([P for P in C.torsion_subgroup().points() if P != C(0) and 2*P == C(0)])
        print("  кривая %d: %s" % (i, C.ainvs()))
        print("     eclib: selmer_rank=%s rank_bound=%s rank=%s certain=%s | 2-кручение(точек порядка2)=%d"
              % (sr, rb, rk, cert, tor2))
        if best is None or rb < best[0]:
            best = (rb, i, cert)
    except Exception as ex:
        print("  кривая %d FAILED: %s" % (i, ex))
print("\nЛУЧШАЯ eclib-граница по классу изогении:", best)

hdr("I2. PARI ellrank на КАЖДОЙ кривой класса (проверка согласованности)")
from sage.libs.pari import pari
for i, C in enumerate(IC.curves):
    try:
        ep = pari.ellinit(list(C.a_invariants()))
        print("  кривая %d: ellrank =" % i, pari.ellrank(ep, 1))
    except Exception as ex:
        print("  кривая %d FAILED: %s" % (i, ex))

hdr("I3. two_descent (eclib) с ростом second_limit — пытаемся получить границу 1")
for sl in [8, 12, 15, 18]:
    try:
        C2 = EllipticCurve(E.a_invariants())
        alarm(500)
        C2.two_descent(second_limit=sl, verbose=False)
        cancel_alarm()
        mc = C2.mwrank_curve()
        print("  second_limit=%d -> eclib rank_bound=%s certain=%s selmer_rank=%s"
              % (sl, mc.rank_bound(), mc.certain(), mc.selmer_rank()))
    except Exception as ex:
        try: cancel_alarm()
        except: pass
        print("  second_limit=%d FAILED: %s" % (sl, ex))

hdr("I4. L-функция: аналитический ранг и НЕЗАВИСИМОСТЬ от 2-спуска")
print("  root_number =", E.root_number(), " => L(E,1)=0, аналитический ранг НЕЧЁТЕН")
try:
    L = E.lseries().dokchitser(100)
    print("  L(E,1)  =", L(1))
    print("  L'(E,1) =", L.derivative(1,1))
    print("  L''(E,1)=", L.derivative(1,2))
    print("  => если L'(E,1) != 0, то аналитический ранг = 1")
except Exception as ex:
    print("  Dokchitser FAILED:", ex)
print("  Гросс-Загье + Колывагин: аналитический ранг 1 => алгебраический ранг 1 и Sha конечна.")
print("  ЭТО ТЕОРЕМА, не гипотеза. Независима от 2-спуска PARI.")

hdr("I5. rank(only_use_mwrank=False) — путь через L-функцию")
try:
    alarm(900)
    print("  E.rank(only_use_mwrank=False) =", E.rank(only_use_mwrank=False))
    cancel_alarm()
except Exception as ex:
    try: cancel_alarm()
    except: pass
    print("  FAILED:", ex)
print("ГОТОВО")
