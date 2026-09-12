# -*- coding: utf-8 -*-
# G1 (m,n)=(19,16): поиск генератора E(Q) (rank 1 по Codex/PARI).
# Подход M0: mwrank/eclib с большим параметром поиска + point_search.
# Всё критическое — точной рациональной арифметикой.

import sys, time

m = 19; n = 16
s = QQ(m^2 + n^2)/2          # 617/2
b = s*m^2*n^2                 # 28510336
e1 = -b
e2 = -s*m^4
e3 = -s*n^4

print("m,n =", m, n)
print("s   =", s)
print("b   =", b, " (integer? %s)" % (b in ZZ))
print("e1,e2,e3 =", e1, e2, e3)

# E : V^2 = (X-e1)(X-e2)(X-e3)
R.<X> = QQ[]
cub = (X-e1)*(X-e2)*(X-e3)
print("cubic =", cub)
a2 = cub[2]; a4 = cub[1]; a6 = cub[0]
E0 = EllipticCurve([0, a2, 0, a4, a6])
print("E0 =", E0)
print("E0 disc =", E0.discriminant().factor())

Emin = E0.minimal_model()
print("Emin =", Emin)
print("conductor =", Emin.conductor().factor())
iso  = E0.isomorphism_to(Emin)     # E0 -> Emin
isoi = Emin.isomorphism_to(E0)     # Emin -> E0
print("iso E0->Emin :", iso)
print("iso Emin->E0 :", isoi)

# ---- проверка, что изоморфизм — квадратный масштаб u^2 (иначе классы delta поменяются)
# Weierstrass iso: (x,y) -> (u^2 x + r, u^3 y + ...)
u_iso = isoi.u   # для Emin->E0
print("u(Emin->E0) =", u_iso, "  u^2 =", u_iso^2)

T = Emin.torsion_subgroup()
print("torsion =", T.invariants(), " order =", T.order())
for P in Emin.torsion_points():
    print("   tors pt:", P)

sys.stdout.flush()

# ---------------- rank bounds -----------------
t0 = time.time()
try:
    rl = Emin.rank(only_use_mwrank=False, proof=False)
    print("rank(sage, proof=False) =", rl)
except Exception as ex:
    print("rank() exception:", ex)
sys.stdout.flush()
print("rank_bounds =", Emin.rank_bounds())
print("analytic_rank =", Emin.analytic_rank())
print("time so far", time.time()-t0)
sys.stdout.flush()

# ---------------- ПОИСК ТОЧКИ: эскалация высоты -----------------
found = []

def report(P, tag):
    if P.is_zero():
        return False
    if P.has_finite_order():
        return False
    print("### НАЙДЕНА неторсионная точка [%s]: %s" % (tag, P))
    print("###   height =", P.height())
    found.append(P)
    return True

# (1) point_search с растущей высотой
for h in [10, 14, 18, 22, 26, 30, 34, 38]:
    t1 = time.time()
    try:
        pts = Emin.point_search(h, rank_bound=1)
    except Exception as ex:
        print("point_search(%s) exception: %s" % (h, ex)); sys.stdout.flush(); continue
    dt = time.time()-t1
    print("point_search(h=%d): %d pts, %.1f s" % (h, len(pts), dt))
    sys.stdout.flush()
    ok = False
    for P in pts:
        if report(P, "point_search h=%d" % h):
            ok = True
    if ok:
        break

print("=== итог поиска: найдено неторсионных:", len(found))
for P in found:
    print("   ", P, " h=", P.height())
