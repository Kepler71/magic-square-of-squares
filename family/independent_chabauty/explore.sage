load('/home/kep/magicKube/family/independent_chabauty/setup.sage')
import sys
b, h, n, D = [ZZ(a) for a in sys.argv[1:5]]
H = int(sys.argv[5]) if len(sys.argv) > 5 else 40
S = Section(b, h, n, D)
print(f"(b,h,n)=({b},{h},{n})  A={S.A} C={S.C}  beta'={S.bet}  k=Q(sqrt({S.dk}))  sb={S.sb}")
print(f"D={D}  c0={S.c0}  roots={S.roots}")
E = S.E
print(f"E = {E.ainvs()}")
print(f"disc = {E.discriminant().norm().factor()}   j = {E.j_invariant()}")
print(f"torsion: {E.torsion_subgroup().invariants()}  ; points: {[P for P in E.torsion_points()]}")
# символическая проверка cz
Rt = PolynomialRing(QQ, 'T'); T = Rt.gen()
f = T*(T^2 - 1)*(S.A*T - S.C)*(S.C*T - S.A)
Rz = PolynomialRing(QQ, 'z'); z = Rz.gen()
Fz = Rz(((1 - z)^6 * f((1 + z)/(1 - z))).simplify_full()) if False else Rz((1-z)^6 * f.subs(T=(1+z)/(1-z)))
print("control cz:", Fz == S.cz * z*(z^2-1)*(z^2-S.bet^2))
# Рациональные u с k-точкой обязаны иметь N(rhs) in Q^2, т.е. u^2-4beta' in Q^2,
# т.е. u = z + beta'/z с z in Q.  Поэтому перебираем z.
us = set()
for q in range(1, H + 1):
    for p0 in range(-H, H + 1):
        if p0 == 0 or gcd(p0, q) != 1: continue
        z0 = QQ(p0)/q
        us.add(z0 + S.bet/z0)
found = []
for u in sorted(us):
    v = S.rhs(S.k(u)*S.scale)
    if v == 0 or v.is_square():
        found.append(u)
print(f"rational u with a k-point, |num|,den <= {H}: {found}")
for u in found:
    P = S.point_from_u(u)
    ts = S.t_of_u(u)
    print(f"   u={u}: t={ts}  X=sqrt(t)->{[tt.sqrt() if QQ(tt).is_square() else None for tt in ts]}  P={None if P is None else P[0]}")
