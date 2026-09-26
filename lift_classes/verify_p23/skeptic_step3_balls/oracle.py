# -*- coding: utf-8 -*-
"""
Контроль реализации balls.py независимым оракулом Sage (QQ(x).is_padic_square(p)).
1) class_code: тривиальность <=> квадрат в Q_p (Sage), мультипликативность.
2) mask-таблицы: для случайных X из шара X0 + p^N Z класс X лежит в маске шара;
   и каждая классовая метка маски реально достигается (полнота => ослабление точное поклеточно).
3) Точная рациональная проверка утверждения Шага 3 ДРУГОЙ дорогой (Sage, без шаров):
   случайные и прицельные (b,c) = (r z, s z), r,s взаимно просты; восемь произведений
   квадраты в Q_p  =>  v_2(z) >= 3 (p=2), v_3(z) >= 1 (p=3); и тогда все 9 клеток квадраты в Q_p.
Все циклы с явными границами.
"""
import sys, random, time, json
from sage.all import QQ, ZZ, gcd, valuation
sys.path.insert(0, '.')
from balls import class_code, residue_mask_table, CELLS, LINES

random.seed(20260926)
t0 = time.time()
out = {}

# 1) class_code против Sage
bad1 = 0
for p in (2, 3):
    for _ in range(20000):                               # граница
        x = random.choice([-1, 1]) * random.randint(1, 10 ** 6) * p ** random.randint(0, 7)
        triv = class_code(x, p) == 0
        if triv != QQ(x).is_padic_square(p):
            bad1 += 1
        y = random.choice([-1, 1]) * random.randint(1, 10 ** 6)
        if class_code(x * y, p) != class_code(x, p) ^ class_code(y, p):
            bad1 += 1
print(f"[1] class_code vs Sage / мультипликативность: ошибок {bad1}", flush=True)
out['class_code_errors'] = bad1

# 2) маски шаров: корректность и полнота
bad2 = 0; incomplete = 0
for p, N in ((2, 1), (2, 2), (2, 3), (2, 5), (3, 1), (3, 2), (3, 4)):
    tab = residue_mask_table(p, N)
    pN = p ** N
    for X0 in range(pN):                                  # граница
        seen = set()
        for t in range(-400, 401):                        # граница
            X = X0 + pN * t
            if X == 0:
                continue
            c = class_code(X, p)
            seen.add(c)
            if not (int(tab[X0]) >> c) & 1:
                bad2 += 1
        mask_set = {c for c in range(16) if (int(tab[X0]) >> c) & 1}
        if seen != mask_set:
            incomplete += 1
print(f"[2] маски: нарушений включения {bad2}, шаров с недостигнутым классом маски {incomplete}", flush=True)
out['mask_errors'] = bad2
out['mask_incomplete_in_window'] = incomplete


# 3) точная проверка утверждения другой дорогой
def eight_square(b, c, p):
    X = [1 + i * b + j * c for (i, j) in CELLS]
    if any(x == 0 for x in X):
        return None
    return all((X[a] * X[bb] * X[cc]).is_padic_square(p) for a, bb, cc in LINES), X


stats = {}
for p, thr in ((2, 3), (3, 1)):
    tested = passed = viol = not_all9 = 0
    by_v = {}
    def check(r, s, z):
        global tested, passed, viol, not_all9
        b, c = r * z, s * z
        res = eight_square(b, c, p)
        if res is None:
            return
        ok, X = res
        tested += 1
        if ok:
            passed += 1
            vz = valuation(z, p)
            by_v[vz] = by_v.get(vz, 0) + 1
            if vz < thr:
                viol += 1
                print("   НАРУШЕНИЕ:", p, r, s, z, flush=True)
            if not all(x.is_padic_square(p) for x in X):
                not_all9 += 1
    # (a) случайные
    for _ in range(60000):                                 # граница
        r = random.randint(-60, 60); s = random.randint(-60, 60)
        if r == 0 or s == 0 or abs(r) == abs(s) or gcd(r, s) != 1:
            continue
        e = random.randint(-3, 5)
        z = QQ(random.randint(-500, 500)) / random.randint(1, 500) * QQ(p) ** e
        if z == 0:
            continue
        check(r, s, z)
    # (b) прицельные: клетка 1 + (ir+js) z почти ноль p-адически (z = -1/lam + p^k t)
    for _ in range(60000):                                 # граница
        r = random.randint(-60, 60); s = random.randint(-60, 60)
        if r == 0 or s == 0 or abs(r) == abs(s) or gcd(r, s) != 1:
            continue
        i, j = random.choice([(i, j) for (i, j) in CELLS if (i, j) != (0, 0)])
        lam = i * r + j * s
        if lam == 0:
            continue
        k = random.randint(1, 12)
        z = QQ(-1) / lam + QQ(p) ** k * QQ(random.randint(-50, 50)) / random.randint(1, 30)
        if z == 0:
            continue
        check(r, s, z)
    # (c) прицельные: z с v_p(z) ровно 0,1,2 и случайной единичной частью (плотная выборка)
    for _ in range(60000):                                 # граница
        r = random.randint(-200, 200); s = random.randint(-200, 200)
        if r == 0 or s == 0 or abs(r) == abs(s) or gcd(r, s) != 1:
            continue
        e = random.randint(-2, thr + 1)
        u = random.randint(1, 10 ** 5)
        if u % p == 0:
            continue
        z = QQ(p) ** e * u / random.choice([1, 5, 7, 11, 13, 25, 35, 49])
        check(r, s, z)
    stats[p] = dict(tested=tested, eight_square_Qp=passed, violations=viol,
                    passed_but_not_all9=not_all9, passed_by_vz=by_v)
    print(f"[3] p={p}: проверено {tested}, 8 квадратов в Q_p: {passed} (по v(z): {dict(sorted(by_v.items()))}), "
          f"нарушений вывода {viol}, не все 9 квадраты {not_all9}, {time.time() - t0:.0f}s", flush=True)
out['exact'] = {str(k): {kk: (vv if not isinstance(vv, dict) else {str(a): b for a, b in vv.items()})
                         for kk, vv in v.items()} for k, v in stats.items()}
json.dump(out, open('oracle.json', 'w'), indent=1, ensure_ascii=False)
print(f"Готово за {time.time() - t0:.0f}s")
