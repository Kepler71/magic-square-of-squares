# -*- coding: utf-8 -*-
"""s12: глобальный контроль леммы 2 (NOTE.md §7) над Q в Sage.
Лемма 2: если x-t, x, x+t -- квадраты в Q_2 (x != 0), то t/x in 8Z_2 и для P = (x, y) in E_t:
    P in 4E_t(Q_2)  <=>  v_2(t/x) >= 5;
класс канонической половинки (корни = 1 mod 4 после масштаба) в координатах (x, x-t):
    v_2(t/x)=3: (5, 3) при (t/x)/8 = 1 mod 4, (5, 7) при (t/x)/8 = 3 mod 4;  v_2(t/x)=4: (1, 5);  >= 5: (1, 1).
Проверки:
 (1) тройки квадратов x(P)-N, x(P), x(P)+N из P = 2R, R in E_N(Q) (N конгруэнтные): P in 4E_N(Q) (Sage division_points)
     => v_2(N/x) >= 5; таблица по v_2;
 (2) настоящая половинка Q из положительных рациональных корней: её 2-адический класс = kappa_2 * delta_2(T(метки)),
     метки -- знаки корней по модулю 4 относительно канонических;
 (3) строки семилинейных квадратов из s5 (тройки точек 2E_b) -- то же.
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s12_lemma2_global.py
"""
from sage.all import EllipticCurve, QQ, ZZ
import itertools, json
from collections import Counter


def v2(q):
    return QQ(q).valuation(2)


def cls2(q):
    """класс рационального q в Q_2^*/Q_2^*2: (v mod 2, единичная часть mod 8)."""
    q = QQ(q)
    v = q.valuation(2)
    u = q / QQ(2) ** v
    return (int(v % 2), int((u.numerator() * pow(int(u.denominator()), -1, 8)) % 8))


def cmul(x, y):
    return ((x[0] + y[0]) % 2, (x[1] * y[1]) % 8)


def kappa2(s):
    k = v2(s)
    if k == 3:
        sig = s / 8
        m4 = int((sig.numerator() * pow(int(sig.denominator()), -1, 4)) % 4)
        return ((0, 5), (0, 3 if m4 == 1 else 7))
    if k == 4:
        return ((0, 1), (0, 5))
    return ((0, 1), (0, 1))


def delta2_T(t, lab):
    one = (0, 1)
    return {(0, 0): (one, one), (1, 0): (cls2(t), cls2(2)), (0, 1): (cls2(-t), cls2(-2 * t)),
            (1, 1): (cls2(-1), cls2(-t))}[lab]


def unit_mod4(q):
    q = QQ(q)
    u = q / QQ(2) ** q.valuation(2)
    return int((u.numerator() * pow(int(u.denominator()), -1, 4)) % 4)


def check_ap(x, t):
    """x-t, x, x+t -- рациональные квадраты. Возвращает (v2(t/x), класс половинки == предсказание?)."""
    x, t = QQ(x), QQ(t)
    al, ga, be = (x - t).sqrt(), x.sqrt(), (x + t).sqrt()
    s = t / x
    k = v2(s)
    # метки: знак корня относительно канонического (= 1 mod 4 после деления на общую степень 2)
    eps = [1 if unit_mod4(r) == 1 else -1 for r in (al, ga, be)]
    lab = (int(eps[0] != eps[1]), int(eps[2] != eps[1]))
    A, B, C = al + ga, ga + be, al + be
    got = (cls2(A * B), cls2(A * C))
    kap = kappa2(s)
    dT = delta2_T(t, lab)
    pred = (cmul(kap[0], dT[0]), cmul(kap[1], dT[1]))
    return k, got == pred


out = {"ap": Counter(), "in4E_by_v2": Counter(), "total_by_v2": Counter(), "class_formula_ok": 0,
       "class_formula_total": 0, "lemma_violations": 0}
for N in [5, 6, 7, 13, 14, 15, 21, 22, 23, 29, 30, 31, 34, 37, 38, 39, 41, 46, 47, 65, 210, 1254]:
    E = EllipticCurve([-N ** 2, 0])
    G = E.gens()
    rng = range(-3, 4) if len(G) >= 2 else range(-8, 9)
    if len(G) >= 3:
        rng = range(-2, 3)
    seen = set()
    for co in itertools.product(rng, repeat=len(G)):
        R = sum((c * g for c, g in zip(co, G)), E(0))
        if R.is_zero():
            continue
        P = 2 * R
        x = P[0]
        if x in seen or len(str(x)) > 120:
            continue
        seen.add(x)
        k, ok = check_ap(x, N)
        in4 = len(P.division_points(4)) > 0
        out["total_by_v2"][int(k)] += 1
        out["in4E_by_v2"][int(k)] += in4
        out["class_formula_total"] += 1
        out["class_formula_ok"] += ok
        if in4 and k < 5:
            out["lemma_violations"] += 1
    print(N, "rank", len(G), "точек", len(seen), dict(out["total_by_v2"]), "в 4E:", dict(out["in4E_by_v2"]),
          "формула класса:", out["class_formula_ok"], "/", out["class_formula_total"], flush=True)

# (3) строки семилинейных квадратов (s5): x = a-c, a+lambda, a+c на E_b
s5 = json.load(open("s5_seven_line_control.json"))
rows_checked, rows_ok = 0, 0
for bkey, cur in s5["curves"].items():
    for ex in cur.get("examples", []):
        b = QQ(int(bkey))
        for xs in ex["x"]:
            k, ok = check_ap(QQ(xs), b)
            rows_checked += 1
            rows_ok += ok
out["seven_line_rows"] = {"checked": rows_checked, "class_formula_ok": rows_ok}
# Саллоус
a, b, c, lam = 4420, -3360, 1056, 8349
sal = []
for x in (a - c, a + lam, a + c):
    k, ok = check_ap(x, b)
    E = EllipticCurve([-b ** 2, 0])
    y = (QQ(x) * (x - b) * (x + b)).sqrt()
    P = E(x, y)
    sal.append({"x": x, "v2(b/x)": int(k), "formula_ok": ok, "P_in_4E": len(P.division_points(4)) > 0})
out["sallows_rows"] = sal
print("семилинейные строки:", out["seven_line_rows"], "Саллоус:", sal)
print("ИТОГ: нарушений леммы 2 (P in 4E при v2 < 5):", out["lemma_violations"])
out = {k: (dict(v) if isinstance(v, Counter) else v) for k, v in out.items()}
json.dump(out, open("s12_lemma2_global.json", "w"), ensure_ascii=False, indent=1, default=str)
