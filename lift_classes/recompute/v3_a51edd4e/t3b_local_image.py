# -*- coding: utf-8 -*-
"""t3b: выборочный локальный ОБРАЗ пары классов ([T],[L]) в (Q_p^*/Q_p^{*2})^2 на точках восьми произведений,
и проверка: лежат ли локальные классы каждого кандидата (d_T,d_L) из списков Codex в этом образе.
Выборка (не полный перебор) => отсутствие класса в образе -- НАБЛЮДЕНИЕ, не доказательство.
Использует функции t3_local.py (свой код).
"""
import json, sys, random, time
from math import gcd
from t3_local import cell_classes, eight_ok, cls_int, prime_factors, samples_a, T_CELLS, L_CELLS
from t1_lists import admissible_A

def local_image(r, s, p, NA, rng, EMIN=-7, EMAX=6):
    img = {}
    for e in range(EMIN, EMAX + 1):
        for a in samples_a(p, r, s, rng, NA):
            C = cell_classes(r, s, p, e, a)
            if C is None or not eight_ok(C): continue
            cT = C[T_CELLS[0]] ^ C[T_CELLS[1]] ^ C[T_CELLS[2]]
            cL = C[L_CELLS[0]] ^ C[L_CELLS[1]] ^ C[L_CELLS[2]]
            img[(cT, cL)] = img.get((cT, cL), 0) + 1
    return img

def main():
    t0 = time.time()
    rng = random.Random(7)
    slopes = [(126, 451), (-126, 451), (73, 362), (265, 298), (12, 85), (2, 3), (4, 9), (8, 9)]
    out = {}
    for (r, s) in slopes:
        DT, DL = admissible_A(r - s), admissible_A(r + s)
        primes = sorted(set(p for p in range(2, 60) if all(p % q for q in range(2, p)))
                        | set(prime_factors(r)) | set(prime_factors(s))
                        | set(prime_factors(r - s)) | set(prime_factors(r + s)))
        imgs = {p: local_image(r, s, p, 6000, rng) for p in primes}
        res = {}
        for dT in DT:
            for dL in DL:
                miss = []
                for p in primes:
                    key = (cls_int(dT, p), cls_int(dL, p))
                    if key not in imgs[p]:
                        miss.append((p, key))
                res["(%d,%d)" % (dT, dL)] = miss
        out["%d/%d" % (r, s)] = dict(candidates=res,
                                     image={p: sorted(list(map(list, imgs[p].keys()))) for p in primes})
        print("%d/%d  [%.1fs]" % (r, s, time.time() - t0), flush=True)
        for k, v in res.items():
            print("   кандидат %s: не найден локально при %s" % (k, v if v else "— (все места ок)"), flush=True)
        for p in primes:
            if len(imgs[p]) < 16 and p > 3:
                print("     образ при p=%d: %s" % (p, sorted(imgs[p].keys())))
    json.dump(out, open("t3b_local_image.json", "w"), indent=1, default=str)
    print("время %.1f" % (time.time() - t0))

if __name__ == "__main__":
    main()
