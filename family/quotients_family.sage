# Поиск фактор-кривых ранга 0 дополненного покрытия для всех сечений (обобщение s15/quotients.sage).
# Формы (t = X^2): f1 = t-1, f2 = t+1, f3 = A t - C, f4 = C t - A; класс (d1,d2,d3), d4 = d3.
# |S|=2: C_S (чётная квартика в X) род 1; |S|=3: E+ (кубика в t), E- (квартика t*G); |S|=4: E (квартика в t).
# Если у компоненты верхняя граница ранга (mwrank 2-спуск) = 0 -> ранг 0 -> X^2 в конечном множестве.
import itertools, functools
print = functools.partial(print, flush=True)
Rx = PolynomialRing(QQ, 'X'); X = Rx.gen()
Rt = PolynomialRing(QQ, 't'); t = Rt.gen()
def jac_quartic(g):
    a, b, c, d, e = [g[i] for i in (4, 3, 2, 1, 0)]
    I = 12*a*e - 3*b*d + c^2; J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return EllipticCurve([-27*I, -27*J]).minimal_model()
def jac_cubic(g):
    a3, a2, a1, a0 = [g[i] for i in (3, 2, 1, 0)]
    return EllipticCurve([0, a2, 0, a1*a3, a0*a3^2]).minimal_model()
secs = [(17,7,13,[(5,13,65),(6,26,39)]), (7,1,5,[(3,5,15),(2,10,5)]), (23,7,17,[(1,34,34),(15,17,255)]),
        (31,17,25,[(7,1,7),(3,2,6)]), (41,1,29,[(21,29,609),(10,58,145)]), (47,23,37,[(35,37,1295),(6,74,111)]),
        (49,31,41,[(1,41,41),(5,82,410)]), (73,17,53,[(5,53,265),(14,106,371)]), (71,49,61,[(11,61,671),(30,122,915)]),
        (89,23,65,[(7,130,910),(33,65,2145)]), (79,47,65,[(7,65,455),(2,130,65)])]
hits = []
for (b, h, n, classes) in secs:
    A, C = ZZ((h^2 + n^2)/2), ZZ((b^2 + n^2)/2)
    Ft = [t - 1, t + 1, A*t - C, C*t - A]
    summary = []
    for cls in classes:
        d = (cls[0], cls[1], cls[2], cls[2])
        for kk in (2, 3, 4):
            for S in itertools.combinations(range(4), kk):
                Dd = prod(d[i] for i in S); G = prod(Ft[i] for i in S)
                if kk == 2: comps = [("C_S", jac_quartic(Rx(Dd * G(X^2))))]
                elif kk == 3: comps = [("E+", jac_cubic(Rt(Dd*G))), ("E-", jac_quartic(Rt(Dd*t*G)))]
                else: comps = [("E", jac_quartic(Rt(Dd*G)))]
                for nm, E in comps:
                    try:
                        lo, hi = E.rank_bounds()
                    except Exception:
                        lo, hi = None, None
                    summary.append(hi)
                    if hi == 0:
                        lab = "{" + ",".join(str(i+1) for i in S) + "}"
                        hits.append(((b, h, n), cls, lab, nm, E.ainvs()))
                        print(f"  *** RANK 0: section ({b},{h},{n}) class {cls} S={lab} {nm}: {E.ainvs()}  torsion {E.torsion_subgroup().invariants()}")
    his = [x for x in summary if x is not None]
    print(f"({b},{h},{n}): {len(summary)} components, upper bounds: min {min(his) if his else None}, histogram {sorted(set((x, his.count(x)) for x in his))}")
print("\nRANK-0 COMPONENTS:", hits if hits else "none")
