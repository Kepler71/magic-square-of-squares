# -*- coding: utf-8 -*-
"""s5: контроль на семилинейных квадратах (lambda != 0) соотношения (R0), обязательного при lambda = 0.

Семилинейный квадрат в форме Codex = три точки P-, P0, P+ из 2E_b(Q) с x = a-c, a+lambda, a+c
(любые три точки 2E_b(Q) с различными x дают такой квадрат; lambda = 0 <=> x в прогрессии).
Половинки Q_j берём по формуле через положительные корни строк.

Соотношение (R0) [NOTE.md, §3, следствие T3]: для нечётного простого q | c, q ∤ b (в целой нормировке корней)
    delta_q(Q_0) * delta_q(Q_{+-1})  in  delta_q(E_b[2])          -- при lambda = 0 доказано;
Контрольное (R+-): delta_q(Q_+) * delta_q(Q_-) in delta_q(E_b[2]) -- верно для ЛЮБОГО семилинейного квадрата
(x_+ - x_- = 2c), должно выполняться всегда.
Тест невырожден, только если delta_q(E_b[2]) -- собственная подгруппа неразветвлённых классов (есть 4-кручение mod q).
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s5_seven_line_control.py
"""
from sage.all import EllipticCurve, QQ, ZZ, lcm, factor, Integer
import itertools, json, sys, random, time
from fractions import Fraction
sys.path.insert(0, '.')
from grid import vp_int, legendre

random.seed(925)
from sage.all import prime_range
SMALLP = [int(q) for q in prime_range(3, 20000)]  # пробное деление: простые q < 2*10^4


def lclass(n, q):
    """класс целого n != 0 в Q_q^*/квадраты: (v mod 2, символ единичной части)."""
    n = int(n)
    v = vp_int(abs(n), q)
    u = n // q ** v
    return (v % 2, legendre(u, q))


def mul(x, y):
    return ((x[0] + y[0]) % 2, x[1] * y[1])


def delta_T(t, q):
    t = int(t)
    cl = lambda z: lclass(z, q)
    return [((0, 1), (0, 1)), (cl(-1), cl(-t)), (cl(t), cl(2)), (cl(-t), cl(-2 * t))]


def half_class(al, ga, be, q):
    A, B, C = al + ga, ga + be, al + be
    return (lclass(A * B, q), lclass(A * C, q))


def rational_sqrt(x):
    x = QQ(x)
    return x.sqrt() if x.is_square() else None


def seven_line_from_points(b, xs):
    """xs = (x_-, x_0, x_+) -- абсциссы точек 2E_b. Возвращает целую нормировку корней строк и b, c."""
    rows = []
    for x in xs:
        r = [rational_sqrt(x - b), rational_sqrt(x), rational_sqrt(x + b)]
        if any(v is None for v in r):
            return None
        rows.append(r)
    D = lcm([v.denominator() for r in rows for v in r])
    for scale in (D, 2 * D):
        rr = [[ZZ(v * scale) for v in r] for r in rows]
        bi = QQ(b) * scale ** 2
        ci = (QQ(xs[2]) - QQ(xs[0])) / 2 * scale ** 2
        li = (QQ(xs[1]) - (QQ(xs[2]) + QQ(xs[0])) / 2) * scale ** 2
        if bi in ZZ and ci in ZZ and li in ZZ:
            return rr, ZZ(bi), ZZ(ci), ZZ(li)
    return None


def test_square(rr, bi, ci, li):
    """Проверка (R0) и (R+-) во всех нечётных q | c, q ∤ b."""
    res = []
    for q in SMALLP:
        if ci % q:
            continue
        e = vp_int(abs(int(ci)), q)
        if bi % q == 0:
            continue
        H = [half_class(abs(r[0]), abs(r[1]), abs(r[2]), q) for r in rr]  # строки -, 0, +
        TT = delta_T(bi, q)
        tors = set(TT)
        vacuous = (len(tors) == 4)
        rel = lambda X, Y: (mul(X[0], Y[0]), mul(X[1], Y[1])) in tors
        res.append({"q": q, "v_q(c)": int(e), "q|lambda": bool(li % q == 0), "|delta_q(E[2])|": len(tors),
                    "R0-": rel(H[1], H[0]), "R0+": rel(H[1], H[2]), "R+-": rel(H[2], H[0])})
    return res


out = {"curves": {}, "summary": {}}
t0 = time.time()
BLIST = [5, 6, 7, 13, 14, 15, 21, 22, 23, 29, 30, 31, 34, 37, 38, 39, 41, 46, 47, 65, 210]
tot = {"tests": 0, "nonvacuous": 0, "R0_fail_nonvac": 0, "R0_checks_nonvac": 0, "Rpm_fail": 0,
       "R0_fail_when_q_divides_lambda": 0, "squares": 0}
for b in BLIST:
    E = EllipticCurve([-b ** 2, 0])
    try:
        G = E.gens()
    except Exception as ex:
        out["curves"][b] = {"error": str(ex)}
        continue
    r = len(G)
    pts = {}
    rng = range(-2, 3) if r >= 2 else range(-4, 5)
    for coeffs in itertools.product(rng, repeat=r):
        if all(c == 0 for c in coeffs):
            continue
        R = sum((c * g for c, g in zip(coeffs, G)), E(0))
        P = 2 * R
        if P.is_zero():
            continue
        x = P[0]
        if x not in pts and len(str(x)) < 400:
            pts[x] = P
    xs = sorted(pts.keys(), key=lambda z: len(str(z)))[:14]
    nsq, ntests = 0, 0
    recs = []
    for trip in itertools.permutations(xs, 3):
        if trip[0] >= trip[2]:
            continue  # x_- < x_+ (ориентация); середина -- любая
        sl = seven_line_from_points(b, trip)
        if sl is None:
            continue
        rr, bi, ci, li = sl
        if li == 0:
            print("!!! lambda = 0 найдено", b, trip)
            continue
        res = test_square(rr, bi, ci, li)
        nsq += 1
        for t in res:
            tot["tests"] += 1
            if not t["|delta_q(E[2])|"] == 4:
                tot["nonvacuous"] += 1
                tot["R0_checks_nonvac"] += 2
                nf = (not t["R0-"]) + (not t["R0+"])
                tot["R0_fail_nonvac"] += nf
                if t["q|lambda"]:
                    tot["R0_fail_when_q_divides_lambda"] += nf
            tot["Rpm_fail"] += (not t["R+-"])
        if len(recs) < 5:
            recs.append({"x": [str(v) for v in trip], "b_int": int(bi), "c_int": int(ci), "lambda_int": int(li), "tests": res})
    tot["squares"] += nsq
    out["curves"][b] = {"rank": r, "gens": [str(g) for g in G], "points_2E_used": len(xs),
                        "seven_line_squares": nsq, "examples": recs}
    print(b, "rank", r, "squares", nsq, dict(tot), "%.0fs" % (time.time() - t0), flush=True)

# Саллоус (a,b,c,lambda) = (4420, -3360, 1056, 8349): строки x = a-c, a+lambda, a+c
a, b, c, lam = 4420, -3360, 1056, 8349
sl = seven_line_from_points(b, (a - c, a + lam, a + c))
rr, bi, ci, li = sl
out["sallows"] = {"b_int": int(bi), "c_int": int(ci), "lambda_int": int(li), "tests": test_square(rr, bi, ci, li)}
print("Саллоус:", out["sallows"])
out["summary"] = tot
json.dump(out, open("s5_seven_line_control.json", "w"), ensure_ascii=False, indent=1, default=str)
print("ИТОГ", tot)
