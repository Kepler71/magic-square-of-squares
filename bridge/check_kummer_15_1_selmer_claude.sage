# -*- coding: utf-8 -*-
# Атака на РАНГОВОЕ звено заявления Codex про (15,1).
#
# Вся конструкция держится на одном CAS-звене: rank E(Q) <= 1.
# Здесь я:
#   (1) независимо получаю границы ранга из нескольких источников;
#   (2) сам проверяю всюду-локальную разрешимость торсора D для класса (1,113,113):
#       если он НЕ разрешим локально где-то, то вывод верен БЕЗ границы ранга;
#       если разрешим — ранговое звено действительно нагружено, надо сказать честно;
#   (3) ищу вторую независимую точку на E (если найду — заявление рухнет).
import sys, time
T0 = time.time()

m = Integer(15); n = Integer(1)
s = (m*m+n*n)/2; b = s*m*m*n*n
e1 = -b; e2 = -s*m**4; e3 = -s*n**4
E = EllipticCurve([0, -(e1+e2+e3), 0, e1*e2+e1*e3+e2*e3, -e1*e2*e3])
Em = E.minimal_model()
print("E  =", E)
print("Em =", Em.ainvs())
print("проводник:", E.conductor().factor())
sys.stdout.flush()

print()
print("--- 2. ТОРСОР ДЛЯ КЛАССА (1,113,113): ЛОКАЛЬНАЯ РАЗРЕШИМОСТЬ ---")
# x-e1 = z1^2, x-e2 = 113 z2^2, x-e3 = 113 z3^2
#   z1^2 - 113 z2^2 = e2-e1 = -5695200
#   z1^2 - 113 z3^2 = e3-e1 =  25312
A12 = e2 - e1; A13 = e3 - e1
print("  z1^2 - 113 z2^2 = %s = %s" % (A12, factor(A12)))
print("  z1^2 - 113 z3^2 = %s = %s"  % (A13, factor(A13)))
print("  совпадает с D у Codex:", (A12 == -5695200 and A13 == 25312))
# 113 | A12, A13  =>  113 | z1;  z1 = 113 w:
print("  113 делит обе правые части => z1 = 113 w, система сводится к")
print("     z2^2 - 113 w^2 = 50400,   z3^2 - 113 w^2 = -224")
assert -A12/113 == 50400 and -A13/113 == -224

def solv(p, wmax=400, jmax=6):
    """ищем w в Q_p (w = a/p^j, a целое) такое, что 50400+113w^2 и -224+113w^2 квадраты"""
    K = Qp(p, prec=40)
    for j in range(0, jmax+1):
        for a in range(0, wmax+1):
            w = QQ(a)/QQ(p**j)
            v2 = 50400 + 113*w*w
            v3 = -224 + 113*w*w
            if v2 == 0 or v3 == 0:
                continue
            if K(v2).is_square() and K(v3).is_square():
                return w
    return None

for p in [2, 3, 5, 7, 113, 11, 13]:
    w = solv(p)
    print("  p=%-4s w =" % p, w, "  локально разрешимо:", w is not None)
    sys.stdout.flush()
# вещественное место
print("  R: при w -> большое оба 50400+113w^2>0 и -224+113w^2>0 => разрешимо (например w=2:",
      RR(50400+113*4), RR(-224+113*4), ")")
sys.stdout.flush()

print()
print("--- 1. ГРАНИЦЫ РАНГА, НЕЗАВИСИМЫЕ ИСТОЧНИКИ ---")
sys.stdout.flush()
for nm, fn in [("rank_bounds", lambda: E.rank_bounds()),
               ("rank()", lambda: E.rank()),
               ("analytic_rank", lambda: E.analytic_rank()),
               ("root_number", lambda: E.root_number()),
               ("selmer_rank", lambda: E.selmer_rank()),
               ("sha().an_numerical", lambda: E.sha().an_numerical())]:
    try:
        print("  Sage %-20s =" % nm, fn())
    except Exception as ex:
        print("  Sage %-20s : ИСКЛЮЧЕНИЕ %s" % (nm, ex))
    sys.stdout.flush()
try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    mw = mwrank_EllipticCurve(list(Em.ainvs()))
    mw.two_descent(verbose=False)
    print("  mwrank rank=%s rank_bound=%s selmer_rank=%s certain=%s"
          % (mw.rank(), mw.rank_bound(), mw.selmer_rank(), mw.certain()))
except Exception as ex:
    print("  mwrank: ИСКЛЮЧЕНИЕ", ex)
sys.stdout.flush()
try:
    import subprocess
    code = "E=ellinit(%s); print(ellrank(E));" % str(list(Em.ainvs())).replace(" ", "")
    r = subprocess.run(["gp", "-q"], input=code, capture_output=True, text=True, timeout=1800)
    print("  PARI ellrank:", r.stdout.strip().replace("\n", " | "), r.stderr.strip()[:200])
except Exception as ex:
    print("  PARI: ИСКЛЮЧЕНИЕ", ex)
sys.stdout.flush()

print()
print("--- 3. ПОИСК ВТОРОЙ НЕЗАВИСИМОЙ ТОЧКИ (попытка сломать rank<=1) ---")
sys.stdout.flush()
try:
    pts = E.point_search(13)
    nt = [P for P in pts if P.order() == Infinity]
    print("  point_search(13): всего", len(pts), " неторсионных", len(nt))
    for P in nt:
        print("   ", P)
    if nt:
        sat, idx, reg = E.saturation(nt)
        print("  насыщенная подгруппа:", sat, " индекс", idx, " регулятор", reg)
        print("  РАНГ найденной подгруппы =", len(sat))
except Exception as ex:
    print("  point_search: ИСКЛЮЧЕНИЕ", ex)

print()
print("время: %.1f c" % (time.time()-T0))
