# Скептическая проверка шага примитивного делителя (Claude §2 шаг 5; Codex §3).
# Запуск: env DOT_SAGE=/tmp/claude_vval_sage python3 check_valuation.py
# Все проверки точные (рациональная арифметика), без факторизации: примитивная
# часть B_n получается выбрасыванием общих простых с B_k, k<n, через НОД.
from sage.all import *
import json, sys, time, itertools

K = int(sys.argv[1]) if len(sys.argv) > 1 else 12      # m = 1..K, индексы 2m <= 2K
out = {"K": K, "curves": []}
t0 = time.time()

def strip(r, B):
    g = gcd(r, B)
    while g > 1:
        r //= g
        g = gcd(r, g)
    return r

def v(q, p):
    q = QQ(q)
    return Infinity if q == 0 else q.valuation(p)

def analyse(E, N, G, label, K):
    rec = {"N": int(N), "G": label, "x(G)": str(G[0])}
    xs = {}
    Q = E(0)
    for n in range(1, 2*K+1):
        Q = Q + G
        xs[n] = QQ(Q[0])
    B = {n: xs[n].denominator() for n in xs}
    # (a) чётность: все B_n — точные квадраты
    rec["all_B_square"] = all(B[n].is_square() for n in B)
    # 2 | B_2 (тогда 2 никогда не примитивный делитель B_{2M}, M>=2)
    rec["v2_B2"] = int(B[2].valuation(2))
    rec["v2_x2m_negative"] = all(v(xs[2*m], 2) < 0 for m in range(1, K+1))
    # (b) примитивная часть B_{2m}
    prim = {}
    for m in range(1, K+1):
        n = 2*m
        r = B[n]
        for k in range(1, n):
            r = strip(r, B[k])
            if r == 1:
                break
        prim[m] = r
    rec["prim_exists_2m_ge4"] = all(prim[m] > 1 for m in range(2, K+1))
    rec["prim_missing_m"] = [m for m in range(1, K+1) if prim[m] == 1]
    rec["two_is_primitive_m"] = [m for m in range(1, K+1) if prim[m] % 2 == 0]
    # (c) никакие три различных x(2mG) не в прогрессии (все расстановки)
    X = {m: xs[2*m] for m in range(1, K+1)}
    ap = []
    for a_, b_, c_ in itertools.combinations(range(1, K+1), 3):
        for (o1, mid, o2) in ((a_, b_, c_), (a_, c_, b_), (b_, a_, c_)):
            if X[o1] + X[o2] == 2*X[mid]:
                ap.append((o1, mid, o2))
    rec["AP_found"] = ap
    # (d) сам механизм: для каждой тройки m1<m2<m3 примитивная часть r = prim[m3]
    # (все простые r — примитивные) делит знаменатель дефекта X[o1]+X[o2]-2X[mid]
    # во всех трёх расстановках, т.е. нормировка дефекта < 0 по этим простым.
    bad = 0; checked = 0; mid_cases = 0
    for m1, m2, m3 in itertools.combinations(range(1, K+1), 3):
        r = prim[m3]
        if r == 1:
            continue
        # простые r не делят B_{2m1}, B_{2m2}
        assert gcd(r, B[2*m1]) == 1 and gcd(r, B[2*m2]) == 1
        for (o1, mid, o2) in ((m1, m2, m3), (m1, m3, m2), (m2, m1, m3)):
            d = (X[o1] + X[o2] - 2*X[mid])
            checked += 1
            if mid == m3:
                mid_cases += 1
            # каждый простой делитель r делит знаменатель дефекта:
            # проверяем без факторизации: strip(r, den(d)) должен дать 1
            if d == 0 or strip(r, d.denominator()) != 1:
                bad += 1
    rec["mechanism_checked"] = checked
    rec["mechanism_mid_cases"] = mid_cases
    rec["mechanism_bad"] = bad
    # то же для сдвига подъёма G+T (T in E[2]): чётные члены не меняются
    return rec, xs, B

curves_r1 = []
curves_r2 = []
for N in range(1, 400):
    if not Integer(N).is_squarefree():
        continue
    E = EllipticCurve([-N**2, 0])
    assert E.is_minimal()
    try:
        r = E.rank(only_use_mwrank=True, proof=True)
    except Exception as e:
        continue
    if r == 1 and len(curves_r1) < 22:
        curves_r1.append(N)
    elif r >= 2 and len(curves_r2) < 8:
        curves_r2.append(N)
    if len(curves_r1) >= 22 and len(curves_r2) >= 8:
        break
out["rank1_N"] = [int(n) for n in curves_r1]
out["rank_ge2_N"] = [int(n) for n in curves_r2]
print("rank1:", curves_r1, "rank>=2:", curves_r2, flush=True)

summary = {"curves": 0, "lines": 0, "all_B_square": True, "v2_B2_ge2": True,
           "v2_x2m_neg": True, "prim_all": True, "two_primitive_any_m_ge2": False,
           "AP_any": False, "mech_checked": 0, "mech_mid": 0, "mech_bad": 0,
           "torsion_4": True}

def lines_for(E, gens):
    T = [t for t in E.torsion_points()]
    L = []
    if len(gens) == 1:
        g = gens[0]
        L = [(g, "g1")] + [(g + t, "g1+T%s" % (i,)) for i, t in enumerate(T) if not t.is_zero()]
    else:
        g1, g2 = gens[0], gens[1]
        combos = [(1,0),(0,1),(1,1),(1,-1),(2,1),(1,2),(3,-1),(2,-3)]
        for (i, j) in combos:
            L.append((i*g1 + j*g2, "%d*g1%+d*g2" % (i, j)))
        L.append((g1 + g2 + T[1] if not T[1].is_zero() else g1+g2+T[2], "g1+g2+T"))
    return L

for N in curves_r1 + curves_r2:
    E = EllipticCurve([-N**2, 0])
    if E.torsion_order() != 4:
        summary["torsion_4"] = False
    gens = E.gens(proof=True)
    summary["curves"] += 1
    evens_by_lift = {}
    for G, lab in lines_for(E, gens):
        rec, xs, B = analyse(E, N, G, lab, K)
        out["curves"].append(rec)
        summary["lines"] += 1
        summary["all_B_square"] &= rec["all_B_square"]
        summary["v2_B2_ge2"] &= (rec["v2_B2"] >= 2)
        summary["v2_x2m_neg"] &= rec["v2_x2m_negative"]
        summary["prim_all"] &= rec["prim_exists_2m_ge4"]
        summary["two_primitive_any_m_ge2"] = summary.get("two_primitive_any_m_ge2", False) or any(m >= 2 for m in rec["two_is_primitive_m"])
        summary["AP_any"] |= bool(rec["AP_found"])
        summary["mech_checked"] += rec["mechanism_checked"]
        summary["mech_mid"] += rec["mechanism_mid_cases"]
        summary["mech_bad"] += rec["mechanism_bad"]
        if len(gens) == 1:
            evens_by_lift[lab] = tuple(xs[2*m] for m in range(1, K+1))
        print(N, lab, "Bsq", rec["all_B_square"], "v2B2", rec["v2_B2"],
              "prim", rec["prim_exists_2m_ge4"], "miss", rec["prim_missing_m"],
              "2prim", rec["two_is_primitive_m"], "AP", rec["AP_found"],
              "mech", rec["mechanism_checked"], "bad", rec["mechanism_bad"],
              "t=%.0fs" % (time.time()-t0), flush=True)
    if evens_by_lift:
        # подъём G и G+T дают одинаковые x(2mG)
        vals = list(evens_by_lift.values())
        summary.setdefault("lift_independent", True)
        summary["lift_independent"] &= all(vv == vals[0] for vv in vals)

out["summary"] = summary
print(json.dumps(summary, indent=1))
with open("check_valuation.json", "w") as f:
    json.dump(out, f, indent=1, default=str)
