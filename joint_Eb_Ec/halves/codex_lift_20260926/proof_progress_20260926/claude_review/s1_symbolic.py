# -*- coding: utf-8 -*-
"""s1: символьные проверки (Sage 10.9) формул для половинок и их образов при 2-спуске / phi-спуске.
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s1_symbolic.py
"""
from sage.all import QQ, PolynomialRing, factor, EllipticCurve, Integer
import itertools, json

out = {}
R = PolynomialRing(QQ, ['al', 'ga', 'be', 'm', 'n'])
al, ga, be, m, n = R.gens()

# ---------- (a) тождества на тройке alpha^2 = x-t, gamma^2 = x, beta^2 = x+t по модулю alpha^2+beta^2-2gamma^2
I = R.ideal([al**2 + be**2 - 2 * ga**2])
t = ga**2 - al**2
A, B, C = al + ga, ga + be, al + be
u = A * B
v = A * B * C
checks = {
    "u-t == A*C": (u - t - A * C) in I,
    "u+t == B*C": (u + t - B * C) in I,
    "v^2 == u^3 - t^2 u": (v**2 - (u**3 - t**2 * u)) in I,
    # x(2Q) = ((u^2+t^2)/(2v))^2 = gamma^2  <=>  (u^2+t^2)^2 - 4 v^2 gamma^2 in I  (при v != 0)
    "x(2Q) == gamma^2": ((u**2 + t**2)**2 - 4 * v**2 * ga**2) in I,
    # также y(2Q)^2 = x(x-t)(x+t) с x = gamma^2 -- следует из того, что 2Q на кривой; проверим x(2Q)-t, x(2Q)+t
    "x(2Q)-t == ((u^2-2tu-t^2)/(2v))^2": ((u**2 - 2 * t * u - t**2)**2 - 4 * v**2 * (ga**2 - t)) in I,
    "x(2Q)+t == ((u^2+2tu-t^2)/(2v))^2": ((u**2 + 2 * t * u - t**2)**2 - 4 * v**2 * (ga**2 + t)) in I,
}
out["a_identities"] = {k: bool(val) for k, val in checks.items()}

# ---------- (b) рациональная параметризация тройки: (alpha, gamma, beta) = s*(m^2-2mn-n^2, m^2+n^2, m^2+2mn-n^2)
sub = {al: m**2 - 2 * m * n - n**2, ga: m**2 + n**2, be: m**2 + 2 * m * n - n**2}
par = {
    "AP": (sub[al]**2 + sub[be]**2 - 2 * sub[ga]**2) == 0,
    "t": str(factor(t.subs(sub))),
    "A": str(factor(A.subs(sub))), "B": str(factor(B.subs(sub))), "C": str(factor(C.subs(sub))),
    "delta_x=AB": str(factor(u.subs(sub))),
    "delta_x-t=AC": str(factor((A * C).subs(sub))),
    "delta_x+t=BC": str(factor((B * C).subs(sub))),
}
out["b_parametrization"] = par

# ---------- (c) смена знаков корней = сдвиг половинки на 2-кручение: delta(Q') = delta(Q) * delta(T)
# delta(T) для E_t: y^2 = x(x-t)(x+t), координаты (x, x-t, x+t):
#   T0=(0,0): (-1, -t, t);  Tp=(t,0): (t, 2, 2t);  Tm=(-t,0): (-t, -2t, 2)
tpar = t.subs(sub)
deltaT = {"O": (1, 1, 1), "T0": (-1, -tpar, tpar), "Tp": (tpar, 2, 2 * tpar), "Tm": (-tpar, -2 * tpar, 2)}


def is_rational_square(f):
    """f -- ненулевой многочлен из Q[m,n]; квадрат в Q(m,n)?"""
    F = factor(R(f))
    return all(e % 2 == 0 for _, e in F) and QQ(F.unit()).is_square()


flip_table = {}
base = [sub[al], sub[ga], sub[be]]
old = (u.subs(sub), (A * C).subs(sub), (B * C).subs(sub))
for signs in itertools.product((1, -1), repeat=3):
    a2, g2, b2 = [sg * x for sg, x in zip(signs, base)]
    A2, B2, C2 = a2 + g2, g2 + b2, a2 + b2
    new = (A2 * B2, A2 * C2, B2 * C2)
    found = [name for name, dT in deltaT.items()
             if all(is_rational_square(new[k] * old[k] * dT[k]) for k in range(3))]
    flip_table[str(signs)] = found
out["c_sign_flips_to_torsion"] = flip_table

# ---------- (d) phi-спуск: для изогении с ядром <(e,0)> отображение x -> x - e (коэффициент спуска двойственной изогении)
# численная проверка на E_34: для точек G1, G2 и половинок -- delta-координаты совпадают с x-e mod квадраты,
# и для P = 2Q все три x-e -- квадраты.
E = EllipticCurve([-34**2, 0])
gens = E.gens()
res = []
for G in gens:
    P = 2 * G
    xs = [P[0] - e for e in (0, 34, -34)]
    res.append({"G": str(G), "x(2G)-e squares": [bool(QQ(x).is_square()) for x in xs]})
    # половинка по формуле через корни
    x = P[0]
    alq, gaq, beq = (x - 34).sqrt(), x.sqrt(), (x + 34).sqrt()
    uq = (alq + gaq) * (gaq + beq)
    vq = uq * (alq + beq)
    Q = E(uq, vq)
    res[-1]["2Q==±P"] = bool(2 * Q == P or 2 * Q == -P)
    res[-1]["Q - G in E[2]"] = bool((Q - G).order() in (1, 2) or (Q + G).order() in (1, 2))
out["d_numeric_E34"] = res

print(json.dumps(out, ensure_ascii=False, indent=1))
json.dump(out, open("s1_symbolic.json", "w"), ensure_ascii=False, indent=1)
