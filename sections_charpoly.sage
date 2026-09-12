# Галуа-орбиты целых сечений через матрицу умножения в Q[a..g]/I (dim 224).
# charpoly(M_l) = prod_{точки} (T - l(pt))^{кратность}; неприводимые множители = орбиты,
# если l разделяет точки (проверяем двумя разными l: разбиение должно совпасть).
import time
Rs.<s> = QQ[]
P = 9*s^4 - 120*s^3 + 664*s^2 - 480*s - 1392
Q = 9*s^4 - 120*s^3 + 760*s^2 - 480*s - 2160
A4, A6 = -432*P, 3456*(s-6)*(3*s-2)*Q
B.<a,b,c,d,e,f,g> = PolynomialRing(QQ, order='degrevlex')
Bs.<t> = B[]
X = a*t^2 + b*t + c; Y = d*t^3 + e*t^2 + f*t + g
I = B.ideal((Y^2 - (X^3 + Bs(A4)(t)*X + Bs(A6)(t))).coefficients())
t0 = time.time(); G = I.groebner_basis(); NB = I.normal_basis()
print("dim Q[..]/I =", len(NB), f"({time.time()-t0:.1f}s)", flush=True)
idx = {m: i for i, m in enumerate(NB)}
def mult_matrix(l):
    M = matrix(QQ, len(NB), len(NB))
    for j, m in enumerate(NB):
        r = (l*m).reduce(G)
        for coef, mon in zip(r.coefficients(), r.monomials()):
            M[idx[mon], j] = coef
    return M
res = {}
for name, l in [("l1", a + 3*b - 2*c + 5*d + 7*e - 11*f + 13*g), ("l2", 2*a - b + 7*c - 3*d + e + 17*f - 5*g), ("a", a), ("c", c)]:
    t0 = time.time()
    cp = mult_matrix(l).charpoly()
    fl = cp.factor()
    res[name] = fl
    print(f"\n[{name}] charpoly factored ({time.time()-t0:.1f}s): orbit (deg, mult):", sorted((h.degree(), m) for h, m in fl), flush=True)
    print(f"      distinct points: {sum(h.degree() for h, m in fl)},  total with mult: {sum(h.degree()*m for h, m in fl)}")
save(res, 'sections_charpoly.sobj')
# поля определения орбит (по l1)
for h, m in sorted(res["l1"], key=lambda hm: hm[0].degree()):
    if h.degree() == 1:
        print(f"orbit size 1, mult {m}: rational")
    else:
        K = NumberField(h, 'z')
        print(f"orbit size {h.degree()}, mult {m}: disc {K.discriminant().factor()}", flush=True)
        sub = [F.degree() for F, _, _ in K.subfields() if F.degree() > 1 and F.degree() < h.degree()]
        quad = [F for F, _, _ in K.subfields() if F.degree() == 2]
        print(f"      proper subfield degrees: {sorted(sub)}; quadratic subfields: {[F.discriminant() for F in quad]}", flush=True)
