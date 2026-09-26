# -*- coding: utf-8 -*-
"""s6: контроль теоремы D -- для ЛЮБЫХ невырожденных b, c в каждом месте v существует точка
X_{b,c}(Q_v) (центр a = t^2, |t|_v велико), в которой все 24 дельта-координаты половинок
(8 линий x 3 координаты) -- квадраты в Q_v. Проверка p-адикой Sage (p = 2 включительно) и R.
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s6_trivial_adelic_point.py
"""
from sage.all import Qp, QQ, RealField, prime_divisors, ceil
import random, json, itertools

random.seed(26)
IDX = (-1, 0, 1)
LINES = []
for j in IDX:
    LINES.append(((-1, j), (0, j), (1, j)))
for i in IDX:
    LINES.append(((i, -1), (i, 0), (i, 1)))
LINES.append(((-1, -1), (0, 0), (1, 1)))
LINES.append(((-1, 1), (0, 0), (1, -1)))


def check_place(p, b, c):
    b, c = QQ(b), QQ(c)
    vb = b.valuation(p)
    vc = c.valuation(p)
    K = max(1, int(ceil((8 - min(vb, vc)) / 2)))
    F = Qp(p, prec=120)
    t = F(p) ** (-K)
    a = t ** 2
    r = {}
    for i in IDX:
        for j in IDX:
            cell = a + i * F(b) + j * F(c)
            s = cell.sqrt()
            r[(i, j)] = s if (s / t - 1).valuation() > 0 else -s
    ok = True
    for (P1, P2, P3) in LINES:
        A, B, C = r[P1] + r[P2], r[P2] + r[P3], r[P1] + r[P3]
        for val in (A * B, A * C, B * C):
            ok &= bool(val.is_square())
    return ok, K


def check_real(b, c):
    R = RealField(200)
    t = R(10) ** 6 * (abs(R(b)) + abs(R(c)))
    a = t ** 2
    r = {(i, j): (a + i * R(b) + j * R(c)).sqrt() for i in IDX for j in IDX}
    return all(v > 0 for (P1, P2, P3) in LINES
               for v in ((r[P1] + r[P2]) * (r[P2] + r[P3]), (r[P1] + r[P2]) * (r[P1] + r[P3]),
                         (r[P2] + r[P3]) * (r[P1] + r[P3])))


out = []
pairs = [(34, 3400), (-3360, 1056), (24, 120), (QQ(7) / 4, QQ(-15) / 8)]
for _ in range(20):
    b = QQ(random.randint(-10 ** 6, 10 ** 6)) / random.randint(1, 50)
    c = QQ(random.randint(-10 ** 6, 10 ** 6)) / random.randint(1, 50)
    if 0 in (b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c):
        continue
    pairs.append((b, c))
tot, good = 0, 0
for b, c in pairs:
    b, c = QQ(b), QQ(c)
    ps = set([2, 3, 5, 7, 11])
    for x in (b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c):
        ps |= set(prime_divisors(x.numerator())) | set(prime_divisors(x.denominator()))
    ps = sorted(ps)[:25]
    res = {}
    for p in ps:
        ok, K = check_place(p, b, c)
        res[int(p)] = ok
        tot += 1
        good += ok
    rok = check_real(b, c)
    tot += 1
    good += rok
    out.append({"b": str(b), "c": str(c), "places_ok": res, "real_ok": rok})
    print(b, c, "все места:", all(res.values()) and rok, flush=True)
print("ИТОГ: мест проверено", tot, "из них все 24 координаты -- квадраты:", good)
json.dump({"pairs": out, "total_places": tot, "all_squares": good}, open("s6_trivial_adelic_point.json", "w"),
          ensure_ascii=False, indent=1)
