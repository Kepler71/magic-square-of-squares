# Полное определение C(Q) для пилота (m,n)=(7,1): фактор B имеет ранг 0, значит точек конечно и их можно перечислить.
load('/home/kep/magicKube/corners/corner_search.sage')
R = PolynomialRing(QQ, 'T'); T = R.gen()
m, n = 7, 1
F0 = m^2 + n^2*T^2; F4 = (m^2+n^2)*(T^2+1)/2; F8 = n^2 + m^2*T^2
qc = QuarticCurve(R(F0*F8), QQ(0)); B = qc.E
print("B =", B.minimal_model().ainvs(), " ранг 0, кручение", B.torsion_order())
pts = B.torsion_points()
print("рациональных точек на B:", len(pts))
ts = set()
for P in pts:
    XY = qc.from_E(P)
    if XY is None: continue
    t1 = XY[0] + qc.t0
    ts.add(t1)
ts |= {QQ(0)}                      # базовая точка t = 0
print("прообразы: значения t на C08(Q):", sorted(ts))
print("\nпроверяем каждое t: квадратны ли F0, F4, F8 по отдельности (это и есть точки кривой рода 5)")
found = []
for t1 in sorted(ts):
    v0, v4, v8 = F0(t1), F4(t1), F8(t1)
    sq = [v0.is_square(), v4.is_square(), v8.is_square()]
    print(f"   t = {t1}:  F0 = {v0} {'□' if sq[0] else '—'},  F4 = {v4} {'□' if sq[1] else '—'},  F8 = {v8} {'□' if sq[2] else '—'}")
    if all(sq): found.append(t1)
print("\nточки кривой рода 5 (все три клетки квадратны):", found)
for t1 in found:
    u0, u4, u8 = sqrt(F0(t1)), sqrt(F4(t1)), sqrt(F8(t1))
    print(f"   t = {t1}: корни клеток (u0,u4,u8) = ({u0}, {u4}, {u8})  — различны: {len({u0,u4,u8})==3}")
