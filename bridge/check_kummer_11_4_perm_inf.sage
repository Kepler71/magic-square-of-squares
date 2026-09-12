# Часть 4: F) инвариантность относительно перестановки корней,
#          G) бесконечность в проективных координатах,
#          H) положительный контроль (15,8),
#          I) прямой поиск точек на C_{11,4}.
def hdr(t):
    print("\n" + "="*78); print(t); print("="*78)
def sqfree(q):
    q = QQ(q); a=q.numerator(); d=q.denominator()
    return ZZ(a*d).squarefree_part()

M, N = ZZ(11), ZZ(4)
S = QQ(M^2+N^2)/2; B = S*M^2*N^2
ROOTS = [QQ(-B), QQ(-S*M^4), QQ(-S*N^4)]
Rx.<XX> = PolynomialRing(QQ)
cub = (XX-ROOTS[0])*(XX-ROOTS[1])*(XX-ROOTS[2]); c = cub.coefficients(sparse=False)
E = EllipticCurve([0, c[2], 0, c[1], c[0]])
Gs = E.gens(proof=False)
T2 = [P for P in E.torsion_subgroup().points() if P != E(0) and 2*P == E(0)]
gens_all = list(Gs) + T2
from itertools import product as iproduct, permutations
grp = []
for co in iproduct([0,1], repeat=len(gens_all)):
    P = E(0)
    for cc,gg in zip(co,gens_all):
        if cc: P = P+gg
    grp.append(P)
grp = list(dict.fromkeys(grp))
print("точек в представителях E(Q)/2E(Q):", len(grp))

def delta(P, roots):
    if P.is_zero(): return (1,1,1)
    x = QQ(P[0]); out=[]
    for i in range(3):
        d = x - roots[i]
        if d == 0:
            j,k = [q for q in range(3) if q != i]
            d = (roots[i]-roots[j])*(roots[i]-roots[k])
        out.append(sqfree(d))
    return tuple(out)

hdr("F. ПЕРЕСТАНОВКА КОРНЕЙ: цель и образ должны переставляться СОГЛАСОВАННО")
print("порядок по определению: e1=-b, e2=-s*m^4, e3=-s*n^4 =", ROOTS)
print("требуемый класс по тождествам в этом порядке: (1, s, s) =", (1, sqfree(S), sqfree(S)))
print()
base = [1, sqfree(S), sqfree(S)]
for perm in permutations([0,1,2]):
    rts = [ROOTS[i] for i in perm]
    req_p = tuple(base[i] for i in perm)
    ims = set(delta(P, rts) for P in grp)
    print("  perm %s: корни(e1',e2',e3') = %s" % (str(perm), rts))
    print("       цель = %s   |образ| = %d   ЦЕЛЬ В ОБРАЗЕ: %s"
          % (str(req_p), len(ims), req_p in ims))

print()
print("ВАЖНО: 'почти совпадение' — в образе есть класс (274,1,274),")
print("       который является ПЕРЕСТАНОВКОЙ цели (1,274,274).")
print("       Проверяем, что это не артефакт порядка: цель и образ переставляются ВМЕСТЕ,")
print("       поэтому ответ 'цель в образе' одинаков для всех 6 перестановок (см. выше).")

hdr("F2. Из какой точки берётся (274,1,274)? и почему это НЕ цель")
for P in grp:
    d = delta(P, ROOTS)
    if sorted([abs(x) for x in d]) == sorted([1,274,274]):
        print("  точка", P, " delta =", d)
        if not P.is_zero():
            x = QQ(P[0])
            print("    x - e1 =", x-ROOTS[0], " класс", sqfree(x-ROOTS[0]))
            print("    x - e2 =", x-ROOTS[1], " класс", sqfree(x-ROOTS[1]))
            print("    x - e3 =", x-ROOTS[2], " класс", sqfree(x-ROOTS[2]))
            print("    для точки C нужно x-e1 = (mn)^2*u4^2 > 0 и КВАДРАТ; здесь класс =", sqfree(x-ROOTS[0]))
            print("    также нужно X/b = t^2 квадрат: X/b =", x/B, " квадрат:", QQ(x/B).is_square())

hdr("G. БЕСКОНЕЧНОСТЬ — проективно/локально, без потери ветвей")
Rz.<z> = PolynomialRing(QQ)
G0 = M^2*z^2 + N^2; G4 = S*(z^2+1); G8 = N^2*z^2 + M^2
Fq = Rz.fraction_field()
print("t = 1/z,  U_i = u_i*z  =>  U_i^2 = z^2 * F_i(1/z):")
print("  z^2*F0(1/z) =", Fq(M^2 + N^2/z^2)*z^2, " == G0 :", Fq(M^2+N^2/z^2)*z^2 == Fq(G0))
print("  z^2*F4(1/z) =", Fq(S*(1+1/z^2))*z^2, " == G4 :", Fq(S*(1+1/z^2))*z^2 == Fq(G4))
print("  z^2*F8(1/z) =", Fq(N^2 + M^2/z^2)*z^2, " == G8 :", Fq(N^2+M^2/z^2)*z^2 == Fq(G8))
print()
print("значения при z=0:  U0^2 =", G0(0), "  U4^2 =", G4(0), "  U8^2 =", G8(0))
print("  G0(0) = n^2 =", N^2, " квадрат в Q:", QQ(G0(0)).is_square(), " => sqrt(G0) в Q((z)) (Гензель)")
print("  G8(0) = m^2 =", M^2, " квадрат в Q:", QQ(G8(0)).is_square(), " => sqrt(G8) в Q((z)) (Гензель)")
print("  G4(0) = s   =", G4(0), " квадрат в Q:", QQ(G4(0)).is_square(),
      " => sqrt(G4) в Q((z)) ТОЛЬКО если s квадрат")
print()
print("ВЫВОД: поле вычетов всех мест над z=0 содержит Q(sqrt(s)); s =", S, "не квадрат =>")
print("       все места над t=inf имеют степень 2 => рациональных точек над бесконечностью НЕТ.")
print()
print("ПЕРЕКРЁСТНАЯ ПРОВЕРКА через H: Y^2 = F0*F4*F8 (секстика по t)")
Rt.<tt> = PolynomialRing(QQ)
H6 = (M^2+N^2*tt^2)*(S*(1+tt^2))*(N^2+M^2*tt^2)
print("  H6 =", H6)
print("  deg =", H6.degree(), " старший коэфф =", H6.leading_coefficient(), " = b:", H6.leading_coefficient()==B)
print("  b =", B, " квадрат в Q:", QQ(B).is_square(), " класс:", sqfree(B))
print("  => у гладкой модели H нет рациональных точек на бесконечности; C->H => и у C нет.")
print()
print("ПЕРЕКРЁСТНАЯ ПРОВЕРКА 2: сколько точек над t=inf на гладкой модели C (над Qbar)?")
print("  F0,F8 имеют старшие коэффициенты n^2,m^2 (квадраты) => каждое sqrt расщепляется: 2*2=4")
print("  F4 имеет старший коэффициент s (не квадрат) => не расщепляется над Q")
print("  всего 8 точек над t=inf в Qbar, объединённых в 4 пары, сопряжённые над Q(sqrt(s)).")
print()
print("ОБНУЛЕНИЕ u_i при конечном t? F0,F4,F8 > 0 для всех вещественных t:")
for nm, pol in [("F0", M^2+N^2*tt^2), ("F4", S*(1+tt^2)), ("F8", N^2+M^2*tt^2)]:
    print("   %s = %s  вещественных корней: %d  дискриминант: %s"
          % (nm, pol, len(pol.roots(RR)), pol.discriminant()))
print("  => u0,u4,u8 != 0 => V = b*u0*u4*u8 != 0 => образ НЕ является 2-кручением;")
print("     delta берётся напрямую, подмена нулевой координаты НЕ нужна. Множитель не теряется.")
print()
print("ПРОВЕРКА X != e_i: X-e1=(mn)^2 u4^2>0, X-e2=s m^2 u0^2>0, X-e3=s n^2 u8^2>0 (s>0:", S>0, ")")

hdr("H. ПОЛОЖИТЕЛЬНЫЙ КОНТРОЛЬ (15,8): точка C известна (t=1) — метод обязан её увидеть")
def build(mm,nn):
    ss = QQ(mm^2+nn^2)/2; bb = ss*mm^2*nn^2
    rr = [QQ(-bb), QQ(-ss*mm^4), QQ(-ss*nn^4)]
    cu = (XX-rr[0])*(XX-rr[1])*(XX-rr[2]); cc = cu.coefficients(sparse=False)
    return ss, bb, rr, EllipticCurve([0,cc[2],0,cc[1],cc[0]])
for (mm,nn,t0) in [(15,8,QQ(1))]:
    ss,bb,rr,EE = build(mm,nn)
    f0 = mm^2+nn^2*t0^2; f4 = ss*(1+t0^2); f8 = nn^2+mm^2*t0^2
    print("(m,n)=(%d,%d) s=%s b=%s t=%s" % (mm,nn,ss,bb,t0))
    print("  F0=%s кв:%s | F4=%s кв:%s | F8=%s кв:%s" %
          (f0,f0.is_square(),f4,f4.is_square(),f8,f8.is_square()))
    u0=QQ(f0).sqrt(); u4=QQ(f4).sqrt(); u8=QQ(f8).sqrt()
    Xp = bb*t0^2; Vp = bb*u0*u4*u8
    P = EE(Xp,Vp)
    print("  точка на E:", P, " ЛЕЖИТ на E: True")
    dd = tuple(sqfree(Xp-rr[i]) for i in range(3))
    print("  delta =", dd, "  цель (1,s,s) =", (1,sqfree(ss),sqfree(ss)),
          "  СОВПАЛО:", dd == (1,sqfree(ss),sqfree(ss)))
    print("  X/b = t^2 =", Xp/bb, " квадрат:", QQ(Xp/bb).is_square())

hdr("I. ПРЯМОЙ ПОИСК точек на C_{11,4} (санити-контроль)")
LIM = 500
found=[]
for q in range(1, LIM+1):
    for p in range(0, LIM+1):
        if gcd(p,q)!=1: continue
        tv = QQ(p)/q
        a0 = M^2+N^2*tv^2
        if not a0.is_square(): continue
        a4 = S*(1+tv^2)
        if not a4.is_square(): continue
        a8 = N^2+M^2*tv^2
        if a8.is_square(): found.append(tv)
print("t=p/q, 0<=p<=%d, 1<=q<=%d: найдено %d точек: %s" % (LIM,LIM,len(found),found))
print("(отсутствие находки — НЕ доказательство; это контроль согласованности)")
print("ГОТОВО")
