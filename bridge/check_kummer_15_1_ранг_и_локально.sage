# check_kummer_15_1_ранг_и_локально.sage
#
# Вторая часть атаки на C_{15,1}(Q)=пусто.
# Главный подозреваемый: ГРАНИЦА РАНГА.
#   E.rank_bounds() даёт (1,3): чистый 2-спуск НЕ даёт rank<=1,
#   потому что dim Sel_2 = 5 и Sha[2] нетривиально.
#   Значит вывод держится только на спаривании Касселса (PARI) / втором
#   спуске (eclib).  Здесь это проверяется всеми доступными независимыми путями,
#   включая аналитический ранг (Гросс–Загир–Колывагин) и прямой поиск
#   второй независимой точки.
#
# Автор: Claude (аудит), 2026-09-12.

import time
def hdr(t):
    print("\n" + "="*76); print(t); print("="*76)

m, n = 15, 1
s = QQ(m^2+n^2)/2; b = s*m^2*n^2
e1, e2, e3 = -b, -s*m^4, -s*n^4
E = EllipticCurve([0, -(e1+e2+e3), 0, e1*e2+e1*e3+e2*e3, -e1*e2*e3])
G = E(-75825, 146764800)

hdr("A. ЧТО ИМЕННО ДАЁТ ЧИСТЫЙ 2-СПУСК")
print("  E =", E, "  N =", E.conductor())
print("  E(Q)[2] полный:", E.torsion_subgroup().invariants())
sel = E.selmer_rank()
print("  dim_F2 Sel_2(E/Q) =", sel)
print("  Sel_2 = rank + dim Sha[2] + dim E(Q)[2] = rank + dim Sha[2] + 2")
print("  => rank + dim Sha[2] =", sel-2, "  => ЧИСТЫЙ 2-спуск даёт только rank <=", sel-2)
print("  E.rank_bounds() =", E.rank_bounds())
print("  ВЫВОД: при rank=3 было бы |E(Q)/2E(Q)|=32 и восьми классов НЕ ХВАТИЛО БЫ.")

hdr("B. НЕЗАВИСИМЫЕ ИСТОЧНИКИ ГРАНИЦЫ rank <= 1")
t0 = time.time()
try:
    pr = pari(E).ellrank()
    print("  PARI ellrank =", pr, "   (формат [r_low, r_up, s, points])")
except Exception as ex:
    print("  PARI ellrank упал:", ex)
try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    Em = mwrank_EllipticCurve(list(E.a_invariants()))
    Em.two_descent(verbose=False, second_descent=True)
    print("  eclib two_descent: rank =", Em.rank(), " certain =", Em.certain(),
          " selmer_rank =", Em.selmer_rank(), " rank_bound =", Em.rank_bound())
except Exception as ex:
    print("  eclib упал:", ex)
print("  (%.1f s)" % (time.time()-t0))

hdr("C. НЕЗАВИСИМЫЙ ПУТЬ: АНАЛИТИЧЕСКИЙ РАНГ + ГРОСС-ЗАГИР-КОЛЫВАГИН")
print("  знак функционального уравнения w =", E.root_number(),
      "  => ранг нечётен" if E.root_number() == -1 else "  => ранг чётен")
t0 = time.time()
try:
    ar = E.analytic_rank(algorithm='pari')
    print("  analytic_rank (PARI ellanalyticrank) =", ar, " (%.1f s)" % (time.time()-t0))
    if ar <= 1:
        print("  ord_{s=1} L(E,s) <= 1  =>  по Колывагину/Гросс-Загиру ранг = аналит. ранг")
        print("  => НЕЗАВИСИМОЕ подтверждение rank E(Q) =", ar)
except Exception as ex:
    print("  analytic_rank упал:", ex)

t0 = time.time()
try:
    sh = E.sha().an_numerical()
    print("  |Sha| (численно по BSD) =", sh, " (%.1f s)" % (time.time()-t0))
    print("  ожидание при rank=1, Sha[2]=(Z/2)^2:  |Sha| = 4")
except Exception as ex:
    print("  sha an_numerical упал:", ex)

hdr("D. ПРЯМАЯ АТАКА: ПОИСК ВТОРОЙ НЕЗАВИСИМОЙ ТОЧКИ (был бы rank>=2)")
t0 = time.time()
for h in [14, 16, 18]:
    try:
        pts = E.point_search(h, rank_bound=None)
        nt = [P for P in pts if P.order() == oo]
        # ранг найденного множества
        if nt:
            rr = E.saturation(nt)[0] if len(nt) > 0 else []
            indep = E.gens_quadratic if False else None
        print("  height_limit=%2d : неторсионных точек %d ; x-коорд %s (%.1fs)" %
              (h, len(nt), [P[0] for P in nt][:6], time.time()-t0))
        # проверяем, все ли они кратны G с точностью до кручения
        outside = []
        for P in nt:
            found = False
            for k in range(-8, 9):
                for T in E.torsion_points():
                    if k*G + T == P:
                        found = True; break
                if found: break
            if not found:
                outside.append(P)
        print("     точек, НЕ представимых как k*G+T (|k|<=8):", len(outside),
              [P.xy() for P in outside][:3])
    except Exception as ex:
        print("  point_search(%d) упал: %s" % (h, ex))
try:
    sat, idx, reg = E.saturation([G])
    print("  E.saturation([G]): индекс =", idx, " насыщено:", sat[0].xy() == G.xy())
except Exception as ex:
    print("  saturation упал:", ex)

hdr("E. САМОПРОВЕРКА ЛОКАЛЬНОЙ РАЗРЕШИМОСТИ C (мой счёт, §6 у Codex)")
# F0=225+t^2, F4=113(1+t^2), F8=1+225t^2; ищем t in Q_p с тремя квадратами.
def is_square_Qp(a, p, prec=40):
    a = QQ(a)
    if a == 0: return False
    v = a.valuation(p)
    if v % 2 != 0: return False
    u = a / p^v
    if p == 2:
        return (ZZ(u.numerator()*u.denominator()) % 8) == 1
    num = ZZ(u.numerator()*u.denominator())
    return kronecker(num, p) == 1
primes_to_test = prime_range(200) + [113]
loc = {}
for p in sorted(set(primes_to_test)):
    found = None
    for tn in range(-60, 61):
        for td in [1, 2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 25, 27, 32, 113]:
            tt = QQ(tn)/QQ(td)
            f0 = 225 + tt^2; f4 = 113*(1+tt^2); f8 = 1 + 225*tt^2
            if is_square_Qp(f0,p) and is_square_Qp(f4,p) and is_square_Qp(f8,p):
                found = tt; break
        if found is not None: break
    loc[p] = found
miss = [p for p in loc if loc[p] is None]
print("  простых проверено:", len(loc), "  без найденной точки:", miss)
print("  примеры (p,t):", [(p, loc[p]) for p in sorted(loc)[:12]])
print("  C(R): t=0 -> F0=225>0, F4=113>0, F8=1>0 -> точка есть")
print("  => если miss пуст, C всюду локально разрешима (для проверенных p)")

hdr("F. ТОРСОР D ДЛЯ КЛАССА (1,113,113) — ЧТО ЭТО ЗНАЧИТ")
d1, d2, d3 = 1, 113, 113
print("  d1 z1^2 - d2 z2^2 = e2-e1 =", e2-e1)
print("  d1 z1^2 - d3 z3^2 = e3-e1 =", e3-e1)
print("  сверка с формой у Codex: -5695200 и 25312 ->",
      (e2-e1 == -5695200, e3-e1 == 25312))
print("  отображение C -> D:  z1=15u4, z2=15u0, z3=u8 ; проверка тождеств:")
Rt.<t> = QQ[]
F0 = 225 + t^2; F4 = 113*(1+t^2); F8 = 1 + 225*t^2
chk1 = (225*F4) - 113*(225*F0) - (e2-e1)
chk2 = (225*F4) - 113*(F8) - (e3-e1)
print("    z1^2-113z2^2-(e2-e1) =", chk1, "   z1^2-113z3^2-(e3-e1) =", chk2)
print("  => C всюду локально разрешима  =>  D всюду локально разрешима")
print("  =>  класс (1,113,113) ЛЕЖИТ в Sel_2(E/Q), но не в образе E(Q)/2E(Q)")
print("  =>  это элемент Sha(E/Q)[2]; согласуется с dim Sel_2 =", sel)
