# G1 (две рёберные пары полные), тройка свободных клеток {0,4,8}; параметр (m,n) = пара сечения
# (тогда при p = 0 центр — квадрат, есть базовая рациональная точка).
# Якобиан кривой рода 5 расщепляется; если хоть один фактор v² = F_i·F_j имеет ранг 0, то точек конечное число
# и они перечисляются явно — безусловное закрытие без Чабо. Маршрут предложен Astra (INTEGRALS_BRIDGE).
import sys
load('/home/kep/magicKube/corners/corner_search.sage')     # QuarticCurve (Коннелл)
Rp = PolynomialRing(QQ, 'p'); pv = Rp.gen()
RT = PolynomialRing(QQ, 'T')
COEF = [(1,1,0),(1,-1,-1),(1,0,1),(1,-1,1),(1,0,0),(1,1,-1),(1,0,-1),(1,1,1),(1,-1,0)]

def cells_G1(m, n):
    m, n = QQ(m), QQ(n)
    P = m*pv + n; Q = m - n*pv; R = m*pv - n; S = m + n*pv
    c = (m^2 + n^2)*(pv^2 + 1)/2
    M = matrix(QQ, [[COEF[1][1], COEF[1][2]], [COEF[3][1], COEF[3][2]]]); Mi = M.inverse()
    ab = Mi*vector(Rp, [Rp(P^2 - c), Rp(R^2 - c)])
    a, b = Rp(ab[0]), Rp(ab[1])
    cells = [Rp(cc*c + ca*a + cb*b) for cc, ca, cb in COEF]
    assert cells[1] == Rp(P^2) and cells[7] == Rp(Q^2) and cells[3] == Rp(R^2) and cells[5] == Rp(S^2)
    return cells

def rk(E):
    r = pari(E.minimal_model().ainvs()).ellinit().ellrank(); return ZZ(r[0]), ZZ(r[1])

secs = []
for line in open('/home/kep/magicKube/corners/sections_500.txt'):
    parts = line.split()
    if len(parts) >= 3:
        try: secs.append(tuple(int(x) for x in parts[:3]))
        except ValueError: pass
lim = int(sys.argv[1]) if len(sys.argv) > 1 else 25
secs = secs[:lim]
print(f"сечений: {len(secs)}; клетки (0,4,8); ищем факторы ранга 0")
closed = 0; opened = 0; hits = []
for (b, h, n) in secs:
    C9 = cells_G1(b, h)
    F = [C9[0], C9[4], C9[8]]
    try:
        facs = {}
        for (i, j) in [(0,1),(0,2),(1,2)]:
            q = RT((F[i]*F[j]).list())
            qc = QuarticCurve(q, QQ(0))
            facs[(i,j)] = (qc, qc.E.minimal_model())
        ranks = {k_: rk(E) for k_, (qc, E) in facs.items()}
    except Exception as ex:
        print(f"  ({b},{h},{n}): пропуск — {type(ex).__name__}: {str(ex)[:40]}"); continue
    zero = [k_ for k_, (lo, hi) in ranks.items() if hi == 0]
    rs = " ".join(f"{k_}:{ranks[k_][0]}..{ranks[k_][1]}" for k_ in ranks)
    if zero:
        k_ = zero[0]; qc, E = facs[k_]
        ts = {QQ(0)}
        iso = E.isomorphism_to(qc.E)          # E -> qc.E (обратная сторона: .inverse() в этой версии Sage нет)
        for P_ in E.torsion_points():
            try:
                XY = qc.from_E(iso(P_))
                if XY: ts.add(XY[0] + qc.t0)
            except Exception as ex:
                print(f"      ИСКЛЮЧЕНИЕ на точке {P_}: {type(ex).__name__}"); raise
        good = [t1 for t1 in ts if all(Fi(t1) > 0 and Fi(t1).is_square() for Fi in F)]
        nd = []
        for t1 in good:
            vals = [ci(t1) for ci in C9]
            if len(set(vals)) == 9 and min(vals) > 0: nd.append((t1, vals))
        closed += 1
        print(f"  ({b},{h},{n}): {rs}  → ЗАКРЫТО (ранг 0): точек {len(ts)}, три квадрата у {len(good)}, невырожденных {len(nd)}")
        if nd: hits.append((b, h, n, nd))
    else:
        opened += 1
        print(f"  ({b},{h},{n}): {rs}")
print(f"\nзакрыто безусловно: {closed}, осталось открытыми: {opened}")
print("невырожденные попадания (7 квадратов):", hits if hits else "нет")
