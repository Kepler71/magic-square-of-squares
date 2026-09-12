# G1, клетки 0,4,8: для каждого (m,n) строим три фактора; если хоть один ранга 0 — перечисляем все точки кривой рода 5.
import sys
load('/home/kep/magicKube/corners/corner_search.sage')
R = PolynomialRing(QQ, 'T'); T = R.gen()
def rk(C):
    r = pari(C.minimal_model().ainvs()).ellinit().ellrank(); return ZZ(r[0]), ZZ(r[1])
MAX = int(sys.argv[1]) if len(sys.argv) > 1 else 12
closed = []; open_ = []; hits = []
for m in range(2, MAX+1):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        F0 = m^2 + n^2*T^2; F4 = QQ(m^2+n^2)*(T^2+1)/2; F8 = n^2 + m^2*T^2
        facs = {}
        try:
            for nm, q in (('A', F0*F4), ('B', F0*F8), ('C48', F4*F8)):
                qc = QuarticCurve(R(q), QQ(0)); facs[nm] = (qc, qc.E.minimal_model())
            ranks = {nm: rk(E) for nm, (qc, E) in facs.items()}
        except Exception as ex:
            print(f"({m},{n}): ошибка {type(ex).__name__}"); continue
        zero = [nm for nm, (lo, hi) in ranks.items() if hi == 0]
        line = f"({m},{n}): ранги " + ", ".join(f"{nm} {ranks[nm][0]}..{ranks[nm][1]}" for nm in ranks)
        if zero:
            nm = zero[0]; qc, E = facs[nm]
            ts = {QQ(0)}
            for P in E.torsion_points():
                try:
                    XY = qc.from_E(qc.E.isomorphism_to(E).inverse()(P) if E is not qc.E else P)
                except Exception:
                    XY = None
                if XY: ts.add(XY[0] + qc.t0)
            good = [t1 for t1 in ts if F0(t1).is_square() and F4(t1).is_square() and F8(t1).is_square()]
            nondeg = []
            for t1 in good:
                u = (sqrt(F0(t1)), sqrt(F4(t1)), sqrt(F8(t1)))
                if len(set(u)) == 3 and all(x != 0 for x in u): nondeg.append((t1, u))
            closed.append((m, n, nm, len(ts), good))
            if nondeg: hits.append((m, n, nondeg))
            print(line + f"  → ЗАКРЫТО через {nm} (ранг 0): точек t {sorted(ts)}, с тремя квадратами {sorted(good)}")
        else:
            open_.append((m, n, ranks)); print(line)
print(f"\nвсего пар: {len(closed)+len(open_)}, закрыто безусловно: {len(closed)}, осталось: {len(open_)}")
print("невырожденных попаданий (7 квадратов):", hits if hits else "нет")
