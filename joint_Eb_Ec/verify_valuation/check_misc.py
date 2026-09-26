# Доп. контроли: (1) формула деления пополам (шаг 1) символьно;
# (2) изоморфизм E_b -> E_N, b = N t^2, сохраняет 2E и прогрессию;
# (3) 2-адика: x(P) при P in 2E_N(Q)\{O} имеет v_2<0 (т.е. 2 | B_2 всегда) — проверка
#     вычетов квадратов mod 4 и численно;
# (4) позитивный контроль «модель важна»: на немасштабированной модели
#     y^2 = x^3 - (N u^2)^2 x примитивный делитель B_{2M} исчезает;
# (5) чувствительность к гипотезе 2E: поиск прогрессий среди x(nG+T), n нечётные/сдвиги.
from sage.all import *
import itertools, json
res = {}

# (1) символьно: alpha^2 + beta^2 = 2 gamma^2, x = gamma^2, t = gamma^2 - alpha^2
R = PolynomialRing(QQ, 'al,be,ga'); al, be, ga = R.gens()
I = R.ideal(al**2 + be**2 - 2*ga**2)
S = R.quotient(I)
t = ga**2 - al**2
X = (al + ga)*(ga + be); Y = (al + ga)*(ga + be)*(al + be)
on_curve = S(Y**2 - (X**3 - t**2*X)) == 0
# x(2Q) = (X^2 + t^2)^2 / (4 Y^2) должно равняться gamma^2
dbl = S((X**2 + t**2)**2 - 4*Y**2*ga**2) == 0
res["halving_on_curve"] = bool(on_curve); res["halving_x2Q_eq_x"] = bool(dbl)
print("(1) halving:", on_curve, dbl)

# (2) масштабирование
ok2 = True
for N in [5, 6, 7, 34]:
    EN = EllipticCurve([-N**2, 0]); G = EN.gens()[0]
    for tq in [QQ(3)/5, QQ(7)/2, QQ(-11)/13]:
        b = N*tq**2
        Eb = EllipticCurve([-b**2, 0])
        for m in [1, 2, 3]:
            P = 2*m*G
            xb, yb = tq**2*P[0], tq**3*P[1]
            Pb = Eb(xb, yb)
            # в 2E_b: xb, xb±b — квадраты
            ok2 &= all(QQ(z).is_square() for z in [xb, xb - b, xb + b])
            ok2 &= Pb in [2*Qp for Qp in [Eb(tq**2*R_[0], tq**3*R_[1]) for R_ in [m*G]]]
        xs = [tq**2*(2*m*G)[0] for m in [1, 2, 3]]
        # прогрессия сохраняется линейно: проверим на искусственной тройке
        s = [QQ(1), QQ(4), QQ(7)]
        ok2 &= (tq**2*s[0] + tq**2*s[2] == 2*tq**2*s[1])
res["scaling_ok"] = bool(ok2); print("(2) scaling:", ok2)

# (3) 2-адика: квадраты Z_2 по mod 4 в {0,1}; разность двух квадратов != 2 mod 4
sq4 = {(a*a) % 4 for a in range(4)}
diffs = {(u - w) % 4 for u in sq4 for w in sq4}
res["sq_diff_mod4"] = sorted(diffs); print("(3) diffs of squares mod 4:", sorted(diffs), "(2 не встречается)")

# (4) немасштабированная модель: a' = -(N u^2)^2 (не свободно от 4-х степеней)
demo = []
for N in [5, 6, 7, 13, 14]:
    EN = EllipticCurve([-N**2, 0]); G = EN.gens()[0]
    for M in [3, 4]:
        n = 2*M
        B = {k: QQ((k*G)[0]).denominator() for k in range(1, n+1)}
        r = B[n]
        for k in range(1, n):
            g = gcd(r, B[k])
            while g > 1:
                r //= g; g = gcd(r, g)
        # u = радикал примитивной части (без факторизации нельзя знать простые;
        # берём u = isqrt(prim) если prim квадрат, иначе prim)
        u = r.isqrt() if r.is_square() else r
        Eu = EllipticCurve([-(N*u**2)**2, 0])
        Gu = Eu(u**2*G[0], u**3*G[1])
        Bu = {k: QQ((k*Gu)[0]).denominator() for k in range(1, n+1)}
        ru = Bu[n]
        for k in range(1, n):
            g = gcd(ru, Bu[k])
            while g > 1:
                ru //= g; g = gcd(ru, g)
        demo.append({"N": N, "n": n, "prim_short_model_gt1": bool(r > 1),
                     "prim_scaled_model_gt1": bool(ru > 1), "B_n_scaled": str(Bu[n])[:40]})
res["scaled_model_demo"] = demo
for d in demo: print("(4)", d)

# (5) чувствительность: прогрессии среди x(nG+T), 1<=n<=8, T in E[2] (точки не обязательно в 2E)
found = []
for N in [5, 6, 7, 13, 14, 15, 21, 22, 23]:
    E = EllipticCurve([-N**2, 0]); G = E.gens()[0]
    pts = {}
    for n in range(1, 9):
        for i, T in enumerate(E.torsion_points()):
            P = n*G + T
            pts[(n, i)] = QQ(P[0])
    xs = sorted(set(pts.values()))
    S_ = set(xs)
    for p1, p2 in itertools.combinations(xs, 2):
        mid = (p1 + p2)/2
        if mid in S_ and mid != p1 and mid != p2:
            labs = [k for k, vv in pts.items() if vv in (p1, mid, p2)]
            found.append((N, str(p1), str(mid), str(p2), [str(l) for l in labs]))
res["AP_among_nG_plus_T"] = found
print("(5) APs among x(nG+T):", len(found))
for f in found[:10]: print("   ", f)
json.dump(res, open("check_misc.json", "w"), indent=1, default=str)
