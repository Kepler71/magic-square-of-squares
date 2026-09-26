#!/usr/bin/env python3
"""zres_mass.py -- массовая проверка REVIEW_LIFT_CLASSES_2026-09-26 §5 прямым вычислением классов клеток.

Утверждение: у всех 2760 наклонов 0<r<s<=300 (gcd=1) с нетривиальными списками Codex
(D_-(r,s) для d_T, D_+(r,s) для d_L) после фильтра по p = 3 mod 4, 7<=p<=400, остаётся только (1,1).

Фильтр (совместный, по паре): (d_T,d_L) выживает при p, если ([d_T]_p,[d_L]_p) лежит в ТОЧНОМ образе
image_Zp(p,r,s) из zres_core (прямое вычисление классов девяти клеток по дереву z mod p^k).
Почему достаточно z in Z_p: при p = 3 mod 4 полюса нет (NOTE.md §2; контроль -- zres_controls.py, блок D).
Классы [d]_p вычисляются по самому d (class_of_int), без предположения p ∤ d.
Ранней остановки нет: считаются все 39 простых для всех 2760 наклонов.
Все циклы конечны: s<=300, p<=400, глубина дерева <= KMAX.
Запуск: python3 zres_mass.py  (лог -- zres_mass.log, данные -- zres_mass.json)
"""
import json
import sys
import time
from math import gcd

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from zres_core import (chi_table, class_of_int, codex_lists, image_Zp, primes_upto)

SMAX = 300
P = [p for p in primes_upto(400) if p % 4 == 3 and p >= 7]
OUT = __file__.rsplit("/", 1)[0]
logf = open(OUT + "/zres_mass.log", "w")


def log(*a):
    msg = " ".join(str(x) for x in a)
    print(msg, flush=True)
    logf.write(msg + "\n")
    logf.flush()


def main():
    t0 = time.time()
    chis = {p: chi_table(p) for p in P}
    log("zres_mass.py, запуск", time.strftime("%Y-%m-%d %H:%M:%S"))
    log("простые p=3 mod 4, 7<=p<=400: %d шт., %s ... %s" % (len(P), P[:6], P[-3:]))
    slopes = [(r, s) for s in range(2, SMAX + 1) for r in range(1, s) if gcd(r, s) == 1]
    lists = {sl: codex_lists(*sl) for sl in slopes}
    nontriv = [sl for sl in slopes if len(lists[sl][0]) > 1 or len(lists[sl][1]) > 1]
    log("наклонов 0<r<s<=%d взаимно простых: %d; с нетривиальными списками: %d" % (SMAX, len(slopes), len(nontriv)))
    npairs = sum(len(lists[sl][0]) * len(lists[sl][1]) - 1 for sl in nontriv)
    log("нетривиальных пар (d_T,d_L) всего: %d" % npairs)
    nT = sum(1 for sl in nontriv if len(lists[sl][0]) > 1)
    nL = sum(1 for sl in nontriv if len(lists[sl][1]) > 1)
    nB = sum(1 for sl in nontriv if len(lists[sl][0]) > 1 and len(lists[sl][1]) > 1)
    log("из них: |D_-|>1 у %d, |D_+|>1 у %d, оба у %d" % (nT, nL, nB))

    # --- контроль p=7 для ВСЕХ 27397 наклонов: точный образ {(0,0)} и ни одной нетривиальной клетки
    bad7 = []
    for sl in slopes:
        img, info = image_Zp(7, sl[0], sl[1], chis[7])
        if img != {(0, 0)} or info["nontriv_cell"]:
            bad7.append(sl)
    log("p=7, все %d наклонов: образ != {(0,0)} или нетривиальная клетка -- у %d наклонов" % (len(slopes), len(bad7)))

    records = []
    val_bit_seen = 0          # классы с нечётной оценкой в образах (по §2 их быть не должно)
    total_nodes = 0
    total_refined = 0
    maxdepth = 0
    survivors_all = []
    need_more_than_7 = 0
    first_kill_hist = {}
    marg_diff = 0
    img_size_hist = {}
    for idx, (r, s) in enumerate(nontriv):
        DT, DL = lists[(r, s)]
        pairs = [(a, b) for a in DT for b in DL if (a, b) != (1, 1)]
        excl = {pr: [] for pr in pairs}
        margT = {a: [] for a in DT}
        margL = {b: [] for b in DL}
        imgs = {}
        for p in P:
            img, info = image_Zp(p, r, s, chis[p])
            total_nodes += info["nodes"]
            total_refined += info["refined"]
            maxdepth = max(maxdepth, info["maxdepth"])
            img_size_hist[len(img)] = img_size_hist.get(len(img), 0) + 1
            if any((cT | cL) & 1 for cT, cL in img):
                val_bit_seen += 1
            imgs[p] = img
            cls = {}
            for d in set(DT) | set(DL):
                cls[d] = class_of_int(d, p, chis[p])
            projT = {cT for cT, _ in img}
            projL = {cL for _, cL in img}
            for pr in pairs:
                if (cls[pr[0]], cls[pr[1]]) not in img:
                    excl[pr].append(p)
            for a in DT:
                if cls[a] not in projT:
                    margT[a].append(p)
            for b in DL:
                if cls[b] not in projL:
                    margL[b].append(p)
        surv = [pr for pr in pairs if not excl[pr]]
        surv7 = [pr for pr in pairs if 7 not in excl[pr]]
        if surv7:
            need_more_than_7 += 1
        # наименьшее P0 такое, что простые <= P0 исключают все пары
        fk = max((min(excl[pr]) if excl[pr] else 10 ** 9) for pr in pairs)
        first_kill_hist[fk] = first_kill_hist.get(fk, 0) + 1
        margsurv = [(a, b) for a in DT for b in DL if (a, b) != (1, 1) and not margT[a] and not margL[b]]
        if sorted(margsurv) != sorted(surv):
            marg_diff += 1
        if surv:
            survivors_all.append({"r": r, "s": s, "surv": surv})
        records.append({"r": r, "s": s, "DT": DT, "DL": DL,
                        "excl": {"%d,%d" % pr: excl[pr] for pr in pairs},
                        "surv": surv, "surv_after_p7": surv7, "P0": fk,
                        "marg_surv": margsurv})
        if (idx + 1) % 250 == 0 or idx + 1 == len(nontriv):
            log("  прогресс: %d/%d наклонов, выживших нетривиальных пар пока: %d, %.1f c"
                % (idx + 1, len(nontriv), sum(len(x["surv"]) for x in survivors_all), time.time() - t0))

    log("")
    log("ИТОГ")
    log("  наклонов с нетривиальными списками: %d, нетривиальных пар: %d" % (len(nontriv), npairs))
    log("  наклонов, где после ВСЕХ p остаётся нетривиальная пара: %d" % len(survivors_all))
    for x in survivors_all[:50]:
        log("    РАСХОЖДЕНИЕ:", x)
    log("  наклонов, где одного p=7 недостаточно: %d" % need_more_than_7)
    log("  P0 = наибольший из 'первых исключающих' простых по парам наклона; гистограмма P0:")
    for k in sorted(first_kill_hist):
        log("    P0=%s: %d" % (k if k < 10 ** 9 else "не исключено", first_kill_hist[k]))
    log("  совместный и маргинальный фильтр дают разные множества выживших у %d наклонов" % marg_diff)
    log("  образы с классом нечётной оценки (не должно быть по §2): %d" % val_bit_seen)
    log("  гистограмма размеров образа по (наклон,p): %s" % dict(sorted(img_size_hist.items())))
    log("  узлов дерева всего %d, из них делений %d, максимальная глубина %d" % (total_nodes, total_refined, maxdepth))
    # сколько простых исключает каждую пару (устойчивость вывода)
    nex = {}
    for rec in records:
        for k, v in rec["excl"].items():
            nex[len(v)] = nex.get(len(v), 0) + 1
    log("  число исключающих простых у пары (из %d) -> число пар: %s" % (len(P), dict(sorted(nex.items()))))
    for (r, s) in [(126, 451), (73, 362), (265, 298)]:
        rec = [x for x in records if x["r"] == r and x["s"] == s]
        if rec:
            rec = rec[0]
            log("  %d/%d: D_-=%s D_+=%s; исключающие простые: %s" % (r, s, rec["DT"], rec["DL"],
                {k: v[:8] for k, v in rec["excl"].items()}))
        else:
            log("  %d/%d: списки тривиальны (в 2760 не входит)" % (r, s))
    log("  время %.1f c" % (time.time() - t0))
    json.dump({"P": P, "n_slopes": len(slopes), "n_nontriv": len(nontriv), "n_pairs": npairs,
               "bad7": bad7, "survivors": survivors_all, "need_more_than_7": need_more_than_7,
               "records": records}, open(OUT + "/zres_mass.json", "w"))


if __name__ == "__main__":
    main()
