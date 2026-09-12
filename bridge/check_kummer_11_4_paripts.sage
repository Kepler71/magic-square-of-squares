# Часть 7: КРИТИЧЕСКАЯ ПРОВЕРКА. PARI ellrank вернул ДРУГИЕ точки, чем Sage gens().
# Если хоть одна из них даёт класс (1,274,274) — исключение Codex немедленно РУШИТСЯ.
def hdr(t):
    print("\n"+"="*78); print(t); print("="*78)
def sqfree(q):
    q=QQ(q); a=q.numerator(); d=q.denominator()
    return ZZ(a*d).squarefree_part()

M,N = ZZ(11),ZZ(4)
S = QQ(M^2+N^2)/2; B = S*M^2*N^2
ROOTS = [QQ(-B), QQ(-S*M^4), QQ(-S*N^4)]
Rx.<XX> = PolynomialRing(QQ)
cub=(XX-ROOTS[0])*(XX-ROOTS[1])*(XX-ROOTS[2]); c=cub.coefficients(sparse=False)
E = EllipticCurve([0,c[2],0,c[1],c[0]])
Emin = E.minimal_model()
phi = E.isomorphism_to(Emin)
print("E    =", E)
print("Emin =", Emin)
print("изоморфизм E->Emin:", phi)

# корни Emin = образы корней E
Rm.<YY> = PolynomialRing(QQ)
fmin = YY^3 + Emin.a2()*YY^2 + Emin.a4()*YY + Emin.a6()
rmin = sorted([r for r,_ in fmin.roots(QQ)])
print("корни Emin:", rmin)
# соответствие: 2-кручение E -> Emin
corr = []
for i,e in enumerate(ROOTS):
    P = E(e,0); Q = phi(P)
    corr.append(QQ(Q[0]))
print("образы корней e1,e2,e3 на Emin (В ТОМ ЖЕ ПОРЯДКЕ):", corr)
print("это перестановка корней Emin:", sorted(corr)==rmin)

u = phi.u if hasattr(phi,'u') else None
print("масштаб u изоморфизма:", u, " u^2 квадрат => квадратные классы СОХРАНЯЮТСЯ")

def delta_min(P):
    if P.is_zero(): return (1,1,1)
    x = QQ(P[0]); out=[]
    for i in range(3):
        d = x - corr[i]
        if d==0:
            j,k=[q for q in range(3) if q!=i]
            d=(corr[i]-corr[j])*(corr[i]-corr[k])
        out.append(sqfree(d))
    return tuple(out)

REQ = (1, sqfree(S), sqfree(S))
print("требуемый класс:", REQ)

hdr("проверка: delta сохраняется при изоморфизме E -> Emin")
def delta_E(P):
    if P.is_zero(): return (1,1,1)
    x=QQ(P[0]); out=[]
    for i in range(3):
        d=x-ROOTS[i]
        if d==0:
            j,k=[q for q in range(3) if q!=i]
            d=(ROOTS[i]-ROOTS[j])*(ROOTS[i]-ROOTS[k])
        out.append(sqfree(d))
    return tuple(out)
Gs = E.gens(proof=False)
T2 = [P for P in E.torsion_subgroup().points() if P!=E(0) and 2*P==E(0)]
for P in [E(0)]+T2+list(Gs)+[Gs[0]+T2[0], Gs[0]+T2[1], 2*Gs[0], 3*Gs[0]]:
    print("  P=%-32s delta_E=%-24s delta_Emin=%-24s совпало:%s"
          % (P, delta_E(P), delta_min(phi(P)), delta_E(P)==delta_min(phi(P))))

hdr("ТОЧКИ, НАЙДЕННЫЕ САМИМ PARI ellrank — их квадратные классы")
pari_pts = [(QQ(121840950)/169, QQ(1814346676800)/2197),
            (QQ(6609976), QQ(16178953076))]
for (xx,yy) in pari_pts:
    try:
        P = Emin(xx,yy)
        d = delta_min(P)
        print("  PARI-точка (%s, %s)" % (xx,yy))
        print("     ЛЕЖИТ на Emin: True")
        print("     delta =", d, "   == требуемый %s : %s" % (str(REQ), d==REQ))
        # выразить через известный генератор
        Pb = Emin.point((xx,yy))
        print("     высота:", Pb.height())
    except Exception as ex:
        print("  PARI-точка (%s,%s) FAILED: %s" % (xx,yy,ex))

hdr("СВЯЗЬ PARI-точек с генератором Sage: независимы ли они?")
Gm = phi(Gs[0])
print("  phi(G) =", Gm, " высота:", Gm.height())
pts = [Emin(p[0],p[1]) for p in pari_pts] + [Gm]
try:
    sat, idx, reg = Emin.saturation(pts)
    print("  saturation([PARI-точки, phi(G)]) -> базис размера %d: %s" % (len(sat), sat))
    print("  => все точки ЗАВИСИМЫ, ранг подгруппы =", len(sat))
except Exception as ex:
    print("  saturation FAILED:", ex)
try:
    print("  матрица высот (регулятор):")
    print(Emin.height_pairing_matrix(pts))
    print("  ранг матрицы высот:", Emin.height_pairing_matrix(pts).rank(),
          " (если 1 — все точки кратны одному генератору)")
except Exception as ex:
    print("  height_pairing_matrix FAILED:", ex)

hdr("ИСЧЕРПЫВАЮЩИЙ ПОИСК точек E с delta=(1,274,274) среди всех найденных")
allpts = set()
for h in [10,14,18]:
    try:
        alarm(400)
        for P in Emin.point_search(h, rank_bound=4):
            allpts.add(P)
        cancel_alarm()
    except Exception as ex:
        try: cancel_alarm()
        except: pass
        print("  point_search(%d): %s" % (h,ex))
allpts |= set(Emin(p[0],p[1]) for p in pari_pts)
allpts.add(Gm)
print("  всего различных найденных точек:", len(allpts))
hits = [P for P in allpts if delta_min(P)==REQ]
print("  точек с требуемым классом:", len(hits), hits)
classes = sorted(set(delta_min(P) for P in allpts))
print("  все встреченные классы (%d):" % len(classes))
for cl in classes: print("     ", cl)
print("  требуемый класс среди них:", REQ in classes)
print("ГОТОВО")
