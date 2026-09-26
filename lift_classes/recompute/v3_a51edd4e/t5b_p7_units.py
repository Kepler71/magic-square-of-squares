# -*- coding: utf-8 -*-
"""t5b: прицельный контроль леммы p=7 на z -- 7-адических ЕДИНИЦАХ (нетривиальная ветвь леммы).
Для случайных наклонов берём z0 in F_7^*, разрешённые леммой (x=r z0, y=s z0, x+-y не равны +-3 mod 7, c(x)=c(y),
c(x+y)c(x-y)=c(x)), поднимаем z = z0 + 7 t (t случайно mod 7^7), либо z = корень клетки + 7^k * u (прицельно в
высокую оценку клетки). Прошедшие 8 условий (Sage, точно) проверяем: T, L и все 9 клеток -- квадраты в Q_7.
Контроль: доля прошедших среди z0, ЗАПРЕЩЁННЫХ леммой, должна быть 0.
"""
import os, random, time, json
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
from math import gcd
from sage.all import QQ
t0 = time.time()
IDX = (-1, 0, 1); CELLS = [(i, j) for i in IDX for j in IDX]
LINES = ([[(i, j) for j in IDX] for i in IDX] + [[(i, j) for i in IDX] for j in IDX]
         + [[(-1, -1), (0, 0), (1, 1)], [(-1, 1), (0, 0), (1, -1)]])
QR = {1, 2, 4}
def c(w):
    w %= 7
    a, b = (1 + w) % 7, (1 - w) % 7
    ca = 0 if a == 0 else (1 if a in QR else -1); cb = 0 if b == 0 else (1 if b in QR else -1)
    if ca and cb and ca != cb: return None
    return ca if ca else cb
def allowed(r, s, z0):
    x, y = r * z0, s * z0
    cs = [c(x), c(y), c(x + y), c(x - y)]
    if None in cs: return False
    return cs[0] == cs[1] and cs[2] * cs[3] == cs[0]
rng = random.Random(7007)
stat = dict(allowed_tested=0, allowed_passed=0, bad=0, forbidden_tested=0, forbidden_passed=0,
            passed_with_cell_div7=0)
for it in range(40000):
    s = rng.randrange(2, 2000); r = rng.randrange(-s + 1, s)
    if r == 0 or gcd(r, s) != 1 or abs(r) == s: continue
    z0 = rng.randrange(1, 7)
    ok = allowed(r, s, z0)
    if rng.random() < 0.5:
        z = QQ(z0 + 7 * rng.randrange(0, 7 ** 7))
    else:
        # прицельно: сделать некоторую клетку 1 + lam z делящейся на высокую степень 7
        lam = rng.choice([r, s, r + s, r - s, -r, -s, -(r + s), -(r - s)])
        if lam % 7 == 0 or (lam * z0 + 1) % 7 != 0: 
            z = QQ(z0 + 7 * rng.randrange(0, 7 ** 7))
        else:
            k = rng.randrange(2, 7); M = 7 ** k
            zr = (-pow(lam, -1, M)) % M
            z = QQ(zr + M * rng.randrange(0, 7 ** 4))
    if z == 0 or z.valuation(7) != 0: continue
    f = {(i, j): 1 + (i * r + j * s) * z for (i, j) in CELLS}
    if any(v == 0 for v in f.values()): continue
    passed = all(QQ(f[l[0]] * f[l[1]] * f[l[2]]).is_padic_square(7) for l in LINES)
    if ok:
        stat["allowed_tested"] += 1
        if passed:
            stat["allowed_passed"] += 1
            if any(v.valuation(7) > 0 for v in f.values()): stat["passed_with_cell_div7"] += 1
            T = f[(-1, 0)] * f[(1, 1)] * f[(0, -1)]; L = f[(-1, 0)] * f[(1, -1)] * f[(0, 1)]
            if not (T.is_padic_square(7) and L.is_padic_square(7) and all(v.is_padic_square(7) for v in f.values())):
                stat["bad"] += 1
    else:
        stat["forbidden_tested"] += 1
        if passed: stat["forbidden_passed"] += 1
    if it % 10000 == 9999: print(" it=%d %s t=%.1fs" % (it + 1, stat, time.time() - t0), flush=True)
print("ИТОГ:", stat, "t=%.1fs" % (time.time() - t0))
json.dump(stat, open("/home/kep/magicKube/lift_classes/recompute/v3_a51edd4e/t5b_p7_units.json", "w"), indent=1)
