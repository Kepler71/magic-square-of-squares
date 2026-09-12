# Контроль перечисления прообразов: сверяем два пути для (7,1) и печатаем все исключения.
load('/home/kep/magicKube/corners/corner_search.sage')
Rp = PolynomialRing(QQ, 'p'); pv = Rp.gen(); RT = PolynomialRing(QQ, 'T')
COEF = [(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]
def cells_G1(m, n):
    m, n = QQ(m), QQ(n)
    P = m*pv + n; Q = m - n*pv; R = m*pv - n; S = m + n*pv
    c = (m^2 + n^2)*(pv^2 + 1)/2
    M = matrix(QQ, [[COEF[1][1], COEF[1][2]], [COEF[3][1], COEF[3][2]]]); Mi = M.inverse()
    ab = Mi*vector(Rp, [Rp(P^2 - c), Rp(R^2 - c)])
    a, b = Rp(ab[0]), Rp(ab[1])
    return [Rp(cc*c + ca*a + cb*b) for cc, ca, cb in COEF]
m, n = 7, 1
C9 = cells_G1(m, n)
print("клетки G1 при (7,1):")
for i, c in enumerate(C9): print(f"  c{i} = {c}")
print("\nсверка с формулами Astra: F0 = m²+n²t², F4 = (m²+n²)(t²+1)/2, F8 = n²+m²t²")
print("  c0 == m²+n²p²:", C9[0] == Rp(m^2 + n^2*pv^2), "   c4 == (m²+n²)(p²+1)/2:", C9[4] == Rp(QQ(m^2+n^2)*(pv^2+1)/2), "   c8 == n²+m²p²:", C9[8] == Rp(n^2 + m^2*pv^2))
q = RT((C9[0]*C9[8]).list())
print("\nквартика F0·F8 =", q)
qc = QuarticCurve(q, QQ(0)); E = qc.E.minimal_model()
print("E =", E.ainvs(), " кручение", E.torsion_order())
iso = E.isomorphism_to(qc.E)          # E -> qc.E (обратная сторона: .inverse() в этой версии Sage нет)
ts = {QQ(0)}
for P_ in E.torsion_points():
    try:
        XY = qc.from_E(iso(P_))
        print("   точка", P_, "→", ("t = %s" % (XY[0] + qc.t0)) if XY else "вырождена (None)")
        if XY: ts.add(XY[0] + qc.t0)
    except Exception as ex:
        print("   точка", P_, "→ ИСКЛЮЧЕНИЕ", type(ex).__name__, str(ex)[:60])
print("итог t:", sorted(ts))
print("\nпроверка каждого t на три квадрата:")
for t1 in sorted(ts):
    v = [C9[0](t1), C9[4](t1), C9[8](t1)]
    print(f"   t={t1}: {v}  квадраты: {[x.is_square() for x in v]}  все девять клеток: {[c(t1) for c in C9]}")
