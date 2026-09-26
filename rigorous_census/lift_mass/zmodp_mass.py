#!/usr/bin/env python3
"""zmodp_mass.py -- массовая проверка утверждения REVIEW_LIFT_CLASSES_2026-09-26 §5:
  у всех 2760 наклонов 0<r<s<=300 с нетривиальными списками Codex (D_-(r-s) для d_T, D_+(r+s) для d_L)
  после фильтра по простым p = 3 mod 4, 7<=p<=400, остаётся только (d_T,d_L) = (1,1).

Фильтр: пара (d_T,d_L) из D_- x D_+ выживает при p, если
  (класс d_T в Q_p, класс d_L в Q_p) лежит в образе {([T]_p,[L]_p)} по всем z in Q_p
  с ненулевыми клетками и восемью квадратами-произведениями в Q_p.
При p = 3 mod 4 полюса нет (доказано в NOTE.md, проверено в zmodp_controls.py), поэтому образ =
образ по z in Z_p, который local_image_Zp считает ТОЧНО прямым вычислением классов девяти клеток.
Фильтр совместный (по паре); маргинальный (по d_T и d_L отдельно) считается для сравнения.
Запуск: python3 zmodp_mass.py > zmodp_mass.log  (все циклы конечны: s<=300, p<=400, KMAX в движке)
"""
import json
import sys
import time
from math import gcd

from zmodp_common import (codex_list, local_image_Zp, primes_3mod4, tl_of, unit_code)

SMAX = 300
P = primes_3mod4(7, 400)


def log(*a):
    print(*a, flush=True)


def main():
    t0 = time.time()
    log("P (%d простых):" % len(P), P)
    slopes = [(r, s) for s in range(2, SMAX + 1) for r in range(1, s) if gcd(r, s) == 1]
    log("наклонов 0<r<s<=%d, gcd=1: %d" % (SMAX, len(slopes)))
    nontriv = []
    for (r, s) in slopes:
        Dm = codex_list(s - r)
        Dp = codex_list(r + s)
        if len(Dm) > 1 or len(Dp) > 1:
            nontriv.append((r, s, Dm, Dp))
    log("с нетривиальными списками: %d" % len(nontriv))

    records = []
    p7_insufficient = []
    discrepancies = []
    marginal_discrepancies = []
    odd_val_vectors = 0          # векторы образа с нечётной оценкой где-либо (теорема §2: должно быть 0)
    trivial_missing = 0          # (0,0) обязан лежать в каждом образе (z = 0 mod p)
    tl_union = {p: set() for p in P}
    first_kill_hist = {}
    needed_prefix_max = 0
    nodes = {}
    for n_done, (r, s, Dm, Dp) in enumerate(nontriv, 1):
        TL = {}
        for p in P:
            im = local_image_Zp(r, s, p, nodes)
            for v in im:
                if any(c & 1 for c in v):
                    odd_val_vectors += 1
            TL[p] = set(tl_of(v) for v in im)
            if (0, 0) not in TL[p]:
                trivial_missing += 1
            tl_union[p] |= TL[p]
        pairs = [(a, b) for a in Dm for b in Dp]
        pair_info = []
        survivors = []
        survivors_p7 = []
        prefix_needed = 0  # индекс в P, после которого все нетривиальные пары убиты
        for (a, b) in pairs:
            killers = [p for p in P if (unit_code(a, p), unit_code(b, p)) not in TL[p]]
            if not killers:
                survivors.append((a, b))
            if 7 not in killers:
                survivors_p7.append((a, b))
            if (a, b) != (1, 1):
                pair_info.append({"dT": a, "dL": b, "killers": killers})
                if killers:
                    fk = killers[0]
                    first_kill_hist[fk] = first_kill_hist.get(fk, 0) + 1
                    prefix_needed = max(prefix_needed, P.index(fk))
        needed_prefix_max = max(needed_prefix_max, prefix_needed)
        # маргинальный фильтр
        mT = [a for a in Dm if all(any(t == unit_code(a, p) for (t, l) in TL[p]) for p in P)]
        mL = [b for b in Dp if all(any(l == unit_code(b, p) for (t, l) in TL[p]) for p in P)]
        if survivors != [(1, 1)]:
            discrepancies.append((r, s, survivors))
        if mT != [1] or mL != [1]:
            marginal_discrepancies.append((r, s, mT, mL))
        if survivors_p7 != [(1, 1)]:
            p7_insufficient.append((r, s))
        nontriv_TL = {p: sorted(TL[p]) for p in P if TL[p] != {(0, 0)}}
        records.append({"r": r, "s": s, "D_minus_dT": Dm, "D_plus_dL": Dp,
                        "pairs": pair_info, "survivors": survivors,
                        "p7_survivors": survivors_p7,
                        "last_needed_prime": P[prefix_needed] if pair_info else None,
                        "nontrivial_TL_images": {str(p): v for p, v in nontriv_TL.items()}})
        if n_done % 250 == 0 or n_done == len(nontriv):
            log("  [%6.1f c] обработано %d/%d; расхождений %d; p=7 недостаточно: %d"
                % (time.time() - t0, n_done, len(nontriv), len(discrepancies), len(p7_insufficient)))

    log("")
    log("ИТОГ")
    log("  наклонов всего: %d (утверждение: 27397)" % len(slopes))
    log("  с нетривиальными списками: %d (утверждение: 2760)" % len(nontriv))
    log("  совместный фильтр: наклонов, где выжило что-то кроме (1,1): %d" % len(discrepancies))
    log("  маргинальный фильтр: наклонов, где выжило что-то кроме 1/1: %d" % len(marginal_discrepancies))
    log("  одного p=7 недостаточно (совместно): %d (утверждение: 504)" % len(p7_insufficient))
    log("  векторов образа с нечётной оценкой: %d (теорема §2: 0)" % odd_val_vectors)
    log("  образов без (0,0): %d (должно быть 0)" % trivial_missing)
    log("  наибольшее простое, реально нужное хоть одному наклону: %d" % P[needed_prefix_max])
    log("  первое убивающее простое (гистограмма по нетривиальным парам):",
        dict(sorted(first_kill_hist.items())))
    log("  объединение образов ([T],[L]) по всем 2760 наклонам (коды 0=квадрат, 2=невычет):")
    for p in P:
        log("    p=%3d: %s" % (p, sorted(tl_union[p])))
    log("  узлов перебора: %d" % nodes.get("nodes", 0))
    if discrepancies:
        log("  РАСХОЖДЕНИЯ (совместный):", discrepancies[:50])
    if marginal_discrepancies:
        log("  маргинальные расхождения (первые 30):", marginal_discrepancies[:30])
    npairs = sum(len(rec["pairs"]) for rec in records)
    log("  нетривиальных пар (d_T,d_L) всего: %d" % npairs)
    log("  время: %.1f c" % (time.time() - t0))
    with open("zmodp_mass.json", "w") as f:
        json.dump({"P": P, "n_slopes": len(slopes), "n_nontrivial": len(nontriv),
                   "discrepancies": discrepancies,
                   "marginal_discrepancies": marginal_discrepancies,
                   "p7_insufficient": p7_insufficient,
                   "odd_val_vectors": odd_val_vectors,
                   "first_kill_hist": first_kill_hist,
                   "tl_union": {str(p): sorted(v) for p, v in tl_union.items()},
                   "records": records}, f, ensure_ascii=False, indent=0)


if __name__ == "__main__":
    sys.exit(main())
