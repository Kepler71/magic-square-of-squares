# Сечение s=1/5, оставшиеся классы. Дополненное покрытие (Codex, §2–3):
#   f1 = p^2 - q^2 = d1 r^2,  f2 = p^2 + q^2 = d2 s^2,  f3 = 109p^2 - 229q^2 = d3 w^2,  f4 = 229p^2 - 109q^2 = d3 v^2
# (Z/2)^4-накрытие P^1, 8 точек ветвления, род 17. Jac ~ prod_S Jac(C_S), C_S: D_S y^2 = prod_{i in S} f_i(X,1),
# D_S = prod_{i in S} d_i  (Kani–Rosen). Все f_i чётны по X, t = X^2:
#   |S|=2: C_S род 1 (чётная квартика в X);
#   |S|=3: C_S род 2, бисэллиптична: E+ : D y^2 = G(t),  E- : D W^2 = t G(t);
#   |S|=4: C_S род 3:  E : D y^2 = G(t) (род 1)  и  J2 : D W^2 = t G(t) (род 2, не считаем).
# Эллиптическая компонента ранга 0 => X^2 в конечном множестве => класс закрывается.
import itertools, functools
print = functools.partial(print, flush=True)
Rx.<X> = QQ[]
Rt.<t> = QQ[]
Ft = [t - 1, t + 1, 109*t - 229, 229*t - 109]           # f_i(X,1) = Ft[i](X^2)

def jac_quartic(g):
    # Якобиан y^2 = g(x), deg g <= 4 (Cremona): Y^2 = X^3 - 27 I X - 27 J
    a, b, c, d, e = [g[i] for i in (4, 3, 2, 1, 0)]
    I = 12*a*e - 3*b*d + c^2
    J = 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3
    return EllipticCurve([-27*I, -27*J]).minimal_model()

def jac_cubic(g):
    # y^2 = g(t), deg g = 3; X = a3 t, Y = a3 y
    a3, a2, a1, a0 = [g[i] for i in (3, 2, 1, 0)]
    return EllipticCurve([0, a2, 0, a1*a3, a0*a3^2]).minimal_model()

def npts(D, poly, p):
    # число точек гладкой модели D y^2 = poly(x) над F_p
    Fp = GF(p); P = poly.change_ring(Fp); Dp = Fp(D)
    n = sum(1 if Dp*P(x) == 0 else (2 if (Dp*P(x)).is_square() else 0) for x in Fp)
    if P.degree() % 2 == 1: n += 1
    else: n += 2 if (Dp * P.leading_coefficient()).is_square() else 0
    return n

def rank_str(E):
    try:
        return str(E.rank(proof=True))
    except Exception:
        lo, hi = E.rank_bounds()
        return f"[{lo},{hi}] an={E.analytic_rank()}"

classes = [(5, 13, 65), (6, 26, 39), (-5, 13, -65), (-6, 26, -39)]
good = [p for p in primes(11, 120) if p not in (13, 109, 229)]
summary = {}
for d3t in classes:
    d = (d3t[0], d3t[1], d3t[2], d3t[2])          # у f4 тот же коэффициент d3
    print(f"\n================ class (d1,d2,d3) = {d3t}")
    zero = []
    for k in (2, 3, 4):
        for S in itertools.combinations(range(4), k):
            D = prod(d[i] for i in S); G = prod(Ft[i] for i in S)
            PX = Rx(G(X^2))
            if k == 2:
                ells = [("C_S", jac_quartic(Rx(D * PX)))]
            elif k == 3:
                ells = [("E+", jac_cubic(Rt(D * G))), ("E-", jac_quartic(Rt(D * t * G)))]
            else:
                ells = [("E", jac_quartic(Rt(D * G)))]
            parts = []
            for nm, E in ells:
                rs = rank_str(E)
                if rs == "0": zero.append((S, nm, E))
                parts.append(f"{nm}: rank {rs} [{E.conductor().factor()}]")
            if k <= 3:   # контроль разложения по числу точек
                bad = [p for p in good if npts(D, PX, p) != p + 1 - sum(E.ap(p) for _, E in ells)]
                parts.append("ctrl OK" if not bad else f"ctrl FAIL {bad}")
            lab = "{" + ",".join(str(i+1) for i in S) + "}"
            print(f"  S={lab:10s} D={D:7d} g={k-1}: " + " | ".join(parts))
    summary[d3t] = zero
print("\n==== rank-0 components:")
for c, z in summary.items():
    print(c, [("{" + ",".join(str(i+1) for i in S) + "}", nm, E.ainvs()) for S, nm, E in z])
