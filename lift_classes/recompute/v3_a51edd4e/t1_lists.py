# -*- coding: utf-8 -*-
"""t1: независимый пересчёт допустимых классов d_T, d_L по ФОРМУЛИРОВКЕ теоремы Codex.
Код Codex не импортируется; из его пакета читается только таблица class_lists.json (данные) для сверки.

  D_-(r,s) = {d>0 бесквадр.: p|d => p|(r-s), p=1 mod 4; d=1 mod 24}   (кандидаты d_T)
  D_+(r,s) = то же с r+s                                               (кандидаты d_L)

Реализация A: своё пробное деление + перебор подмножеств простых.
Реализация B (s<=SMAX_B): Sage divisors(m), moebius(d)!=0, factor(d), d%24==1 — другой маршрут.
Дополнительно (помечено как ГИПОТЕЗА до проверки в t3): уточнённые списки D'_,
  где простой p с чётной v_p(r-s) допускается только при p=1 mod 8.

Наклоны: 1 <= r < s <= SMAX, gcd(r,s)=1 (пара {r,s} и r<->s дают тот же набор списков;
для -r/s списки D_- и D_+ меняются местами, доля "оба тривиальны" от знака не зависит).
Все циклы ограничены: s <= SMAX.
"""
import json, sys, time, itertools
from math import gcd, prod

OUT = "t1_lists"
SMAX_LIST = [300, 600, 1000, 2000]
SMAX = max(SMAX_LIST)
SMAX_B = 300

def factor_td(n):
    """Пробное деление, n>0. Ограничено: d*d<=n."""
    assert n > 0
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f

def admissible_A(m, refined=False):
    """Список допустимых d для m=r-s или r+s (m != 0)."""
    f = factor_td(abs(m))
    P = [p for p, e in f.items() if p % 4 == 1 and (not refined or e % 2 == 1 or p % 8 == 1)]
    out = []
    for k in range(len(P) + 1):
        for S in itertools.combinations(P, k):
            d = prod(S)
            if d % 24 == 1:
                out.append(d)
    return sorted(out)

def main():
    t0 = time.time()
    try:
        cl = json.load(open("/home/kep/magicKube/lift_classes/codex/proof_progress_20260926/lift/class_lists.json"))
    except Exception as e:
        cl = None
        print("class_lists.json Codex не прочитан:", e)
    res = {}
    # --- три наклона
    for (r, s) in [(126, 451), (73, 362), (265, 298)]:
        rec = {}
        for sign in (1, -1):
            rr = sign * r
            DT, DL = admissible_A(rr - s), admissible_A(rr + s)
            rec["%d/%d" % (rr, s)] = {"r-s": rr - s, "fac(r-s)": factor_td(abs(rr - s)),
                                     "r+s": rr + s, "fac(r+s)": factor_td(abs(rr + s)),
                                     "D_T": DT, "D_L": DL,
                                     "refined_D_T": admissible_A(rr - s, True),
                                     "refined_D_L": admissible_A(rr + s, True)}
        key = "%d/%d" % (r, s)
        if cl is not None and key in cl:
            rec["codex_T"] = cl[key]["T_classes"]; rec["codex_L"] = cl[key]["L_classes"]
            rec["match_codex"] = (rec[key]["D_T"] == cl[key]["T_classes"] and rec[key]["D_L"] == cl[key]["L_classes"])
        res[key] = rec
        print(key, json.dumps(rec, default=str), flush=True)

    # --- массовый перебор
    stats = {S: dict(slopes=0, both_trivial=0, only_T=0, only_L=0, both_nontriv=0,
                     refined_both_trivial=0, refined_differs=0) for S in SMAX_LIST}
    viol_gcd = 0; viol_73 = 0; mismatch_AB = 0; checked_B = 0
    maxpairs = (0, None)
    first_refined_diff = []
    sage_ok = True
    try:
        import os
        os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
        from sage.all import divisors, moebius, factor as sfactor
    except Exception as e:
        sage_ok = False
        print("Sage недоступен:", e)
    listing300 = {}
    for s in range(2, SMAX + 1):
        for r in range(1, s):
            if gcd(r, s) != 1:
                continue
            DT, DL = admissible_A(r - s), admissible_A(r + s)
            rDT, rDL = admissible_A(r - s, True), admissible_A(r + s, True)
            # следствия теоремы, которые должны выполняться автоматически
            for a in DT:
                for b in DL:
                    if gcd(a, b) != 1:
                        viol_gcd += 1
            if 0 < abs(r - s) < 73 and DT != [1]:
                viol_73 += 1
            if 0 < abs(r + s) < 73 and DL != [1]:
                viol_73 += 1
            npairs = len(DT) * len(DL)
            if npairs > maxpairs[0]:
                maxpairs = (npairs, (r, s, DT, DL))
            if (rDT, rDL) != (DT, DL) and len(first_refined_diff) < 15:
                first_refined_diff.append((r, s, DT, DL, rDT, rDL))
            for S in SMAX_LIST:
                if s <= S:
                    st = stats[S]
                    st["slopes"] += 1
                    tT, tL = DT == [1], DL == [1]
                    if tT and tL: st["both_trivial"] += 1
                    elif tL: st["only_T"] += 1
                    elif tT: st["only_L"] += 1
                    else: st["both_nontriv"] += 1
                    if rDT == [1] and rDL == [1]: st["refined_both_trivial"] += 1
                    if (rDT, rDL) != (DT, DL): st["refined_differs"] += 1
            if s <= SMAX_B:
                if not (DT == [1] and DL == [1]):
                    listing300["%d/%d" % (r, s)] = [DT, DL]
                if sage_ok:
                    for m, D in ((r - s, DT), (r + s, DL)):
                        B = sorted(int(d) for d in divisors(abs(m))
                                   if moebius(d) != 0 and d % 24 == 1
                                   and all(p % 4 == 1 for p, _ in sfactor(d)))
                        checked_B += 1
                        if B != D:
                            mismatch_AB += 1
                            if mismatch_AB < 5:
                                print("РАСХОЖДЕНИЕ A/B", r, s, m, D, B)
        if s % 200 == 0 or s == 300:
            print("  s=%d  t=%.1fs" % (s, time.time() - t0), flush=True)
    for S in SMAX_LIST:
        st = stats[S]
        st["frac_both_trivial"] = st["both_trivial"] / st["slopes"]
        st["frac_refined_both_trivial"] = st["refined_both_trivial"] / st["slopes"]
        print("s<=%d:" % S, st, flush=True)
    summary = dict(stats=stats, viol_gcd=viol_gcd, viol_73=viol_73, checked_B=checked_B,
                   mismatch_AB=mismatch_AB, sage_ok=sage_ok, maxpairs=maxpairs,
                   first_refined_diff=first_refined_diff, time_s=round(time.time() - t0, 1))
    print("нарушений gcd:", viol_gcd, " нарушений следствия 73:", viol_73)
    print("сверка A/B (Sage) при s<=%d: списков %d, расхождений %d" % (SMAX_B, checked_B, mismatch_AB))
    print("макс. число пар (d_T,d_L):", maxpairs)
    print("первые наклоны, где уточнённые списки отличаются:", first_refined_diff[:8])
    json.dump(dict(three=res, summary=summary, nontrivial_s_le_300=listing300),
              open(OUT + ".json", "w"), indent=1, default=str)
    print("время %.1f с" % (time.time() - t0))

if __name__ == "__main__":
    main()
