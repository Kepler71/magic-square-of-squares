#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# ПРОВЕРКА (не подтверждение — попытка сломать) заявленного исключения G1 (m,n)=(19,5).
# Claude, 2026-09-12.  Всё пересчитывается с нуля, данные Codex используются только как мишень.
#
# Что проверяется:
#  A. арифметика семейства: s, b, e_i, тождества X-e_i, требуемый класс (1,s,s);
#  B. структура кручения E(Q) (важно для счёта 2^(r+2));
#  C. верхняя граница ранга ДВУМЯ независимыми средствами: eclib (mwrank) и PARI ellrank;
#  D. 8 классов Куммера: принадлежность точек E(Q), различность, групповая структура,
#     отсутствие требуемого класса;
#  E. лежит ли требуемый класс в 2-Selmer (то есть нужна ли вообще граница ранга);
#  F. ветви, потерянные при переходе C -> E: u_i = 0, бесконечность, особые точки.

import sys

m, n = 19, 5
print("=" * 78)
print("A. АРИФМЕТИКА СЕМЕЙСТВА (m,n) = (%d,%d)" % (m, n))
print("=" * 78)
assert gcd(m, n) == 1, "m,n не взаимно просты"
s = QQ(m**2 + n**2) / 2
b = s * m**2 * n**2
print("s = (m^2+n^2)/2 =", s, "   (Codex: 193)")
print("b = s*m^2*n^2   =", b, "   (Codex: 1741825)")
e1 = -b
e2 = -s * m**4
e3 = -s * n**4
print("e1,e2,e3 =", (e1, e2, e3))
print("Codex e  =", (-1741825, -25151953, -120625))
assert (s, b) == (193, 1741825)
assert (e1, e2, e3) == (-1741825, -25151953, -120625)
assert len({e1, e2, e3}) == 3

R.<t> = QQ[]
F0 = m**2 + n**2 * t**2
F4 = s * (1 + t**2)
F8 = n**2 + m**2 * t**2
print("F0 =", F0, "  F4 =", F4, "  F8 =", F8)

# тождества отображения C -> E : X = b t^2
X = b * t**2
id1 = (X - e1) - (m * n)**2 * F4
id2 = (X - e2) - s * m**2 * F0
id3 = (X - e3) - s * n**2 * F8
print("X-e1 - (mn)^2*F4 =", id1)
print("X-e2 - s*m^2*F0  =", id2)
print("X-e3 - s*n^2*F8  =", id3)
assert id1 == 0 and id2 == 0 and id3 == 0, "ТОЖДЕСТВА ОТОБРАЖЕНИЯ НЕВЕРНЫ"
print("тождества X-e_i выполняются тождественно по t  -> ПОДТВЕРЖДЕНО")


def sqclass(x):
    """представитель класса в Q*/Q*^2: бесквадратное целое."""
    x = QQ(x)
    assert x != 0
    return ZZ(x.numerator() * x.denominator()).squarefree_part()


req = (sqclass(1), sqclass(s), sqclass(s))
print("требуемый класс delta = (1, s, s) =", req, "   (Codex: (1,193,193))")
assert req == (1, 193, 193)
print("произведение координат требуемого класса (должно быть квадратом):",
      sqclass(req[0] * req[1] * req[2]))

print()
print("=" * 78)
print("B. КРИВАЯ E И ЕЁ КРУЧЕНИЕ")
print("=" * 78)
A2 = -(e1 + e2 + e3)
A4 = e1 * e2 + e1 * e3 + e2 * e3
A6 = -e1 * e2 * e3
E = EllipticCurve([0, A2, 0, A4, A6])
print("E:", E)
print("disc  =", E.discriminant())
print("j     =", E.j_invariant())
N = E.conductor()
print("conductor N =", N, " = ", factor(N))
T = E.torsion_subgroup()
print("torsion subgroup:", T, "  order =", T.order())
print("E(Q)[2] =", E.torsion_subgroup().invariants(), " -> dim_F2 E(Q)[2] =",
      len([d for d in T.invariants() if d % 2 == 0]))
t2 = len([d for d in T.invariants() if d % 2 == 0])
print("|T/2T| = |T[2]| =", 2**t2)
print("ВЫВОД: |E(Q)/2E(Q)| = 2^r * |T[2]| = 2^(r+%d)" % t2)
print("(замечание: |T/2T|=|T[2]| для любой конечной абелевой T, поэтому даже")
print(" кручение Z/2xZ/4 или Z/2xZ/8 НЕ изменило бы счёт 2^(r+2))")
print("root number w =", E.root_number())

print()
print("=" * 78)
print("D. ВОСЕМЬ КЛАССОВ КУММЕРА — ПЕРЕСЧЁТ С НУЛЯ")
print("=" * 78)
GX = QQ(-28561963657) / 1369
GY = QQ(-2089086742828800) / 50653
print("точка Codex G = (%s, %s)" % (GX, GY))
lhs = GY**2
rhs = (GX - e1) * (GX - e2) * (GX - e3)
print("V^2 - prod(X-e_i) =", lhs - rhs)
assert lhs == rhs, "ТОЧКА G НЕ ЛЕЖИТ НА E"
G = E(GX, GY)
print("G на E -> ПОДТВЕРЖДЕНО;  порядок G:", "бесконечный" if G.order() == oo else G.order())
print("высота Нерона-Тэйта h(G) =", G.height())

T1 = E(e1, 0)
T2 = E(e2, 0)
T3 = E(e3, 0)
roots = [e1, e2, e3]


def delta(P):
    """полный 2-спуск: E(Q) -> (Q*/Q*^2)^3."""
    if P == E(0):
        return (1, 1, 1)
    x = P[0]
    out = []
    for i in range(3):
        d = x - roots[i]
        if d == 0:
            j, k = [q for q in range(3) if q != i]
            d = (roots[i] - roots[j]) * (roots[i] - roots[k])
        out.append(sqclass(d))
    return tuple(out)


pts = [("O", E(0)), ("T1", T1), ("T2", T2), ("T1+T2", T1 + T2),
       ("G", G), ("G+T1", G + T1), ("G+T2", G + T2), ("G+T1+T2", G + T1 + T2)]
codex_claim = [(1, 1, 1), (-1, 4053, -4053), (-4053, 386, -42), (4053, 42, 386),
               (-3, 386, -1158), (3, 42, 14), (1351, 1, 1351), (-1351, 4053, -3)]
mine = []
print("%-9s %-28s %s" % ("точка", "мой delta", "delta Codex"))
for (nm, P), cc in zip(pts, codex_claim):
    d = delta(P)
    mine.append(d)
    flag = "СОВПАДАЕТ" if d == cc else "!!! РАСХОЖДЕНИЕ !!!"
    print("%-9s %-28s %-28s %s" % (nm, d, cc, flag))
    # контроль: произведение координат должно быть квадратом
    assert sqclass(d[0] * d[1] * d[2]) == 1, "произведение координат не квадрат: " + nm

print("все 8 различны?", len(set(mine)) == 8, " (найдено различных: %d)" % len(set(mine)))
assert len(set(mine)) == 8


def mul(u, v):
    return tuple(sqclass(u[i] * v[i]) for i in range(3))


closed = all(mul(u, v) in set(mine) for u in mine for v in mine)
print("множество замкнуто относительно умножения (это подгруппа)?", closed)
assert closed, "предъявленные классы НЕ образуют группу — противоречие с гомоморфностью"
print("требуемый класс (1,193,193) среди восьми?", req in set(mine))
assert req not in set(mine)
print("ПОДТВЕРЖДЕНО: образ delta (при rank<=1) есть группа порядка 8, (1,193,193) в неё не входит")

print()
print("=" * 78)
print("C. ВЕРХНЯЯ ГРАНИЦА РАНГА — ДВА НЕЗАВИСИМЫХ СРЕДСТВА")
print("=" * 78)
sys.stdout.flush()

print("--- C.1 eclib / mwrank ---")
try:
    from sage.libs.eclib.interface import mwrank_EllipticCurve
    ec = mwrank_EllipticCurve(list(E.ainvs()))
    ec.set_verbose(0)
    lo = ec.rank()
    print("mwrank rank()      =", lo)
    print("mwrank certain()   =", ec.certain())
    print("mwrank rank_bound()=", ec.rank_bound())
    print("mwrank selmer_rank()=", ec.selmer_rank())
except Exception as ex:
    print("eclib ОШИБКА:", ex)
try:
    print("Sage E.rank_bounds() =", E.rank_bounds())
except Exception as ex:
    print("rank_bounds ОШИБКА:", ex)
try:
    print("Sage E.selmer_rank() =", E.selmer_rank())
except Exception as ex:
    print("selmer_rank ОШИБКА:", ex)
sys.stdout.flush()

print("--- C.2 PARI ellrank (r1 <= rank <= r2, s = rk Sha[2]/2Sha[4]) ---")
pE = pari(E)
for eff in [0, 1, 2, 4]:
    try:
        r = pE.ellrank(eff)
        print("effort=%d -> [r1,r2,s] = [%s,%s,%s], найдено точек: %d"
              % (eff, r[0], r[1], r[2], len(r[3])))
    except Exception as ex:
        print("effort=%d ОШИБКА: %s" % (eff, ex))
    sys.stdout.flush()

print()
print("=" * 78)
print("E. ЛЕЖИТ ЛИ ТРЕБУЕМЫЙ КЛАСС (1,193,193) В 2-SELMER?")
print("=" * 78)
print("если НЕТ — исключение безусловно и граница ранга не нужна;")
print("если ДА  — всё держится ровно на rank<=1, и (1,193,193) есть элемент Sha[2].")
# 2-накрытие для delta=(d1,d2,d3): d1 z1^2 - d2 z2^2 = e2-e1, d1 z1^2 - d3 z3^2 = e3-e1
d1, d2, d3 = req
print("2-накрытие H_delta:  %d*z1^2 - %d*z2^2 = %s ,  %d*z1^2 - %d*z3^2 = %s"
      % (d1, d2, e2 - e1, d1, d3, e3 - e1))
badp = set([2]) | set(ZZ(2 * m * n * (m**2 - n**2) * (m**2 + n**2)).prime_factors()) \
    | set(ZZ(E.discriminant().numerator()).prime_factors())
print("плохие простые для локального анализа:", sorted(badp))

print()
print("=" * 78)
print("F. ПОТЕРЯННЫЕ ВЕТВИ ПРИ ПЕРЕХОДЕ C -> E")
print("=" * 78)
print("F.1 нули u_i при конечном рациональном t:")
for nmF, F in [("F0", F0), ("F4", F4), ("F8", F8)]:
    rr = F.roots(QQ)
    print("   %s рациональные корни: %s ; знак при вещественном t: %s"
          % (nmF, rr, "всегда > 0" if F(0) > 0 and F.discriminant() < 0 else "проверить"))
    assert rr == [], "у %s есть рациональный корень — ветвь u=0 не исключена" % nmF
print("   -> ни один u_i не обращается в 0 в рациональной точке: ветвь потеряна НЕ БЫЛА")
print("F.2 гладкость аффинной модели: якобиан по (u0,u4,u8) = diag(2u0,2u4,2u8);")
print("    особая точка требует u_i=0 И F_i'(t)=0 одновременно.")
for nmF, F in [("F0", F0), ("F4", F4), ("F8", F8)]:
    g = gcd(F, F.derivative())
    print("   gcd(%s, %s') = %s  -> кратных корней нет: %s" % (nmF, nmF, g, g.degree() == 0))
    assert g.degree() == 0
print("F.3 точки над t = бесконечность: z=1/t, U_i=u_i/t даёт U4^2 = s.")
print("    s =", s, " квадрат в Q?", QQ(s).is_square())
assert not QQ(s).is_square()
print("    -> рациональных точек над бесконечностью нет")
print("F.4 образ точки C всегда НЕторсионный на E: V = b*u0*u4*u8 != 0, X != e_i,")
print("    значит delta(P) = [X-e1,X-e2,X-e3] без регуляризации. Подтверждено F.1.")

print()
print("=" * 78)
print("ИТОГ БЛОКОВ A,B,D,F — все проверки пройдены. Остаётся блок C (граница ранга).")
print("=" * 78)
