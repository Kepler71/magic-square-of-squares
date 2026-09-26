# -*- coding: utf-8 -*-
"""t3e: контроль универсальной леммы p=7 на реальных наклонах и новая граница вместо 73.
(1) Для каждого из 48 классов (r0,s0) mod 7 берём реальный взаимно простой наклон (r,s) с такими вычетами
    (и отдельно наклоны с 7 | r-s, 7 | r+s, 49 | r+-s), перебираем z = 7^e a, a -- единицы mod 7^4, e in [-3,4],
    ТОЧНО проверяем восемь произведений в Q_7 и то, что T, L -- квадраты в Q_7.
(2) Наименьший d>1: бесквадратный, простые = 1 mod 4, d = 1 mod 24, И квадрат mod 7.
Циклы ограничены явно.
"""
import json, time
from math import gcd
from t3_local import cell_classes, eight_ok, T_CELLS, L_CELLS
t0 = time.time()
reps = {}
for s in range(2, 400):
    for r in range(-s + 1, s):
        if r == 0 or gcd(r, s) != 1: continue
        key = (r % 7, s % 7)
        if key not in reps: reps[key] = (r, s)
extra = [(r, s) for (r, s) in [(1, 8), (1, 13), (1, 48), (3, 46), (2, 51), (5, 54), (1, 97), (24, 25)] if gcd(r, s) == 1]
slopes = sorted(set(reps.values())) + extra
tested = passed = bad = 0
cls_counts = {}
for (r, s) in slopes:
    for e in range(-3, 5):
        for a in range(1, 7 ** 4):
            if a % 7 == 0: continue
            for sg in (1, -1):
                C = cell_classes(r, s, 7, e, sg * a)
                if C is None: continue
                tested += 1
                if not eight_ok(C): continue
                passed += 1
                cT = C[T_CELLS[0]] ^ C[T_CELLS[1]] ^ C[T_CELLS[2]]
                cL = C[L_CELLS[0]] ^ C[L_CELLS[1]] ^ C[L_CELLS[2]]
                cls_counts[(cT, cL)] = cls_counts.get((cT, cL), 0) + 1
                if cT or cL: bad += 1
print("(1) наклонов %d (классов mod 7: %d), тестов %d, прошло %d, нетривиальный [T]/[L] в Q_7: %d, классы %s, t=%.1fs"
      % (len(slopes), len(reps), tested, passed, bad, cls_counts, time.time() - t0))

def ok_codex(d):
    n = d; p = 2; seen = set()
    while p * p <= n:
        while n % p == 0:
            if p in seen or p % 4 != 1: return False
            seen.add(p); n //= p
        p += 1
    if n > 1:
        if n in seen or n % 4 != 1: return False
    return d % 24 == 1
QR7 = {1, 2, 4}
first = [d for d in range(2, 5000) if ok_codex(d)][:20]
first7 = [d for d in range(2, 5000) if ok_codex(d) and d % 7 in QR7][:20]
print("(2) первые допустимые по Codex:", first)
print("    из них квадраты mod 7:", first7)
json.dump(dict(slopes=len(slopes), residue_classes=len(reps), tested=tested, passed=passed, nontrivial=bad,
               class_counts={str(k): v for k, v in cls_counts.items()}, first_codex=first, first_with_QR7=first7),
          open("t3e_p7_universal.json", "w"), indent=1)
