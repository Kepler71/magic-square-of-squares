# -*- coding: utf-8 -*-
"""s8: контроль предложения S (NOTE.md, §6): архимедов размер против меток половинок.
Для примитивной целой АП квадратов alpha^2, gamma^2, beta^2 (разность t > 0, НОД корней 1) и нечётных p | t:
  метка s_alpha(p) = [p | gamma+alpha], s_beta(p) = [p | beta+gamma];
  Pi1 = prod_{s_alpha=0} p^{v_p(t)} <= gamma - alpha < sqrt(t);
  Pi2 = prod_{s_beta=0}  p^{v_p(t)} <= beta - gamma < sqrt(t)/2;
  Pi3 = prod_{s_alpha=s_beta} p^{v_p(t)} <= beta - alpha < sqrt(2t).
Следствие: если t_odd := prod_{p нечётно, v_p(t) нечётно} p^{v_p(t)} >= sqrt(t)/2, то P = (gamma^2, ...) не лежит в 4E_t(Q).
Проверка: АП из точек 2E_N(Q) (N конгруэнтные) -> примитивная целая нормировка; сверка неравенств и
(для следствия) прямая проверка P in 4E через Sage (P.division_points(4)).
Запуск: env DOT_SAGE=/tmp/claude_halves_sage python3 s8_size_labels.py
"""
from sage.all import EllipticCurve, QQ, ZZ, lcm, gcd, factor, sqrt, RR
import itertools, json

out = {"aps": 0, "ineq_fail": 0, "corollary_applicable": 0, "corollary_fail": 0, "in4E_total": 0, "examples": []}
for N in [5, 6, 7, 13, 14, 15, 21, 22, 23, 30, 34, 41, 65, 210]:
    E = EllipticCurve([-N ** 2, 0])
    G = E.gens()
    rng = range(-3, 4) if len(G) >= 2 else range(-6, 7)
    seen = set()
    for co in itertools.product(rng, repeat=len(G)):
        R = sum((c * g for c, g in zip(co, G)), E(0))
        if R.is_zero():
            continue
        P = 2 * R
        x = P[0]
        if x in seen or len(str(x)) > 70:
            continue
        seen.add(x)
        roots = [(x - N).sqrt(), x.sqrt(), (x + N).sqrt()]
        D = lcm([r.denominator() for r in roots])
        al, ga, be = [abs(ZZ(r * D)) for r in roots]
        g = gcd([al, ga, be])
        al, ga, be = al // g, ga // g, be // g
        t = ga ** 2 - al ** 2
        assert be ** 2 - ga ** 2 == t and t > 0
        Pi1 = Pi2 = Pi3 = 1
        todd = 1
        for p, e in factor(t):
            if p == 2:
                continue
            sa = int((ga + al) % p == 0)
            sb = int((be + ga) % p == 0)
            if sa == 0:
                Pi1 *= p ** e
            if sb == 0:
                Pi2 *= p ** e
            if sa == sb:
                Pi3 *= p ** e
            if e % 2:
                todd *= p ** e
        ok = (Pi1 ** 2 < t) and (4 * Pi2 ** 2 < t) and (Pi3 ** 2 < 2 * t)
        out["aps"] += 1
        out["ineq_fail"] += (not ok)
        # следствие: t_odd >= sqrt(t)/2  =>  P not in 4E_t
        Et = EllipticCurve([-t ** 2, 0])
        Pt = Et(ga ** 2, al * ga * be)
        in4 = len(Pt.division_points(4)) > 0
        out["in4E_total"] += in4
        if 4 * todd ** 2 >= t:
            out["corollary_applicable"] += 1
            out["corollary_fail"] += in4
        if len(out["examples"]) < 12:
            out["examples"].append({"N": N, "t": int(t), "roots": [int(al), int(ga), int(be)],
                                    "Pi1/sqrt(t)": float(Pi1 / RR(t).sqrt()), "Pi2/sqrt(t)": float(Pi2 / RR(t).sqrt()),
                                    "Pi3/sqrt(2t)": float(Pi3 / RR(2 * t).sqrt()), "P_in_4E": in4})
    print(N, {k: v for k, v in out.items() if k != "examples"}, flush=True)
json.dump(out, open("s8_size_labels.json", "w"), ensure_ascii=False, indent=1)
