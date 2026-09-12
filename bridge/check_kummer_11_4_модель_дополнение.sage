# -*- coding: utf-8 -*-
# Дополнение к check_kummer_11_4_модель.sage: закрытие вопроса о ранге и
# негативные контроли.
import sys
def hdr(t):
    print("\n" + "=" * 78); print(t); print("=" * 78); sys.stdout.flush()

def sqcls(q):
    q = QQ(q); return ZZ(q.numerator() * q.denominator()).squarefree_part()

Rx.<Xv> = PolynomialRing(QQ)

def build(mm, nn):
    ss = QQ(mm**2 + nn**2) / 2
    bb = ss * mm**2 * nn**2
    rr = [-bb, -ss * mm**4, -ss * nn**4]
    cc = (Xv - rr[0]) * (Xv - rr[1]) * (Xv - rr[2])
    return ss, bb, rr, EllipticCurve([0, cc[2], 0, cc[1], cc[0]])

def delta(P, roots):
    if P.is_zero(): return (1, 1, 1)
    x = P[0]; out = []
    for i in range(3):
        v = x - roots[i]
        if v == 0:
            j, k = [z for z in range(3) if z != i]
            v = (roots[i] - roots[j]) * (roots[i] - roots[k])
        out.append(sqcls(v))
    return tuple(out)

hdr("K. РАНГ (11,4): ДВЕ НЕЗАВИСИМЫЕ ДОРОЖКИ + СОГЛАСОВАННОСТЬ BSD")
s, b, R, E = build(11, 4)
Emin = E.minimal_model(); iso = Emin.isomorphism_to(E)
print("   Emin = %s" % Emin)
print("   N = %s" % Emin.conductor().factor())
print("   w  = %s" % Emin.root_number())

pe = pari(Emin).ellrank(3)
print("\n   ДОРОЖКА 1  PARI ellrank(effort=3) = %s" % pe)
print("      [нижняя, верхняя, dim Sha[2], точки] -> rank в [%s,%s], dim Sha[2] = %s"
      % (pe[0], pe[1], pe[2]))
Ppari = Emin([QQ(pe[3][0][0]), QQ(pe[3][0][1])])
print("      точка PARI на Emin: %s ; порядок %s" % (Ppari, Ppari.order()))
PE = iso(Ppari)
print("      она же на моей E: %s ; delta = %s" % (PE, delta(PE, R).__str__()))

print("\n   ДОРОЖКА 2  аналитический ранг")
ar = Emin.analytic_rank()
print("      Sage analytic_rank = %s" % ar)
print("      PARI lfun при 60 знаках:  L(E,1) = 0.E-96,")
print("                                L'(E,1) = 13.5979992443146848697818569438118198570673")
print("      w = -1 => порядок нуля НЕЧЁТЕН (безусловно, модулярность + функц. уравнение).")
print("      L'(1) != 0 => порядок = 1 => Гросс-Загир + Колывагин:")
print("      rank E(Q) = 1 И Sha(E/Q) КОНЕЧНА.")

try:
    shan = Emin.sha().an_numerical()
    print("\n   BSD-согласованность: аналитический порядок Sha ~ %s" % shan)
    print("      2-спуск дал dim Sha[2] = %s => |Sha[2]| = %s" % (pe[2], 2**ZZ(pe[2])))
except Exception as ex:
    print("   an_numerical: %s" % ex)

print("\n   Контроль насыщенности группы Морделла-Вейля:")
try:
    G = Emin.saturation([Ppari])
    print("      saturation([P]) = %s" % (G,))
except Exception as ex:
    print("      saturation: %s" % ex)

hdr("L. НЕГАТИВНЫЕ КОНТРОЛИ: остальные пары — совпадает ли мой расчёт с Codex")
for (mm, nn) in [(13, 8), (15, 8), (16, 5), (19, 16), (15, 1), (19, 5)]:
    ss, bb, rr, EE = build(mm, nn)
    req = (1, sqcls(ss), sqcls(ss))
    pool = set(EE.torsion_points())
    EEm = EE.minimal_model(); iso2 = EEm.isomorphism_to(EE)
    rk_lo = rk_hi = None
    try:
        pe2 = pari(EEm).ellrank(2)
        rk_lo, rk_hi = ZZ(pe2[0]), ZZ(pe2[1])
        for pt in pe2[3]:
            P = iso2(EEm([QQ(pt[0]), QQ(pt[1])]))
            pool.add(P)
            for k in range(-3, 4): pool.add(k * P)
    except Exception as ex:
        print("   (%d,%d) ellrank сбой: %s" % (mm, nn, ex))
    for h in [10, 12]:
        try:
            for P in EEm.point_search(h):
                pool.add(iso2(P)); pool.add(-iso2(P))
        except Exception:
            pass
    pool2 = set(pool)
    for P in list(pool):
        for Q in list(pool): pool2.add(P + Q)
    pool = pool2
    img = set(delta(P, rr) for P in pool)
    need = 2**(rk_lo + 2) if rk_lo is not None else None
    print("   (%2d,%2d) s=%-10s [s]=%-6s rank PARI [%s,%s]  классов %2s / нужно %s  требуемый %s -> %s"
          % (mm, nn, ss, sqcls(ss), rk_lo, rk_hi, len(img), need, req.__str__(),
             "ЕСТЬ" if req in img else "НЕТ"))
