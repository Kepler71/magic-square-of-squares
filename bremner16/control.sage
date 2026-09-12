# Контроли движка: (1) совпадение F1 с corners/bremner_type.sage, (2) поиск примера Бремнера–Саллоуса.
import functools, sys
print = functools.partial(print, flush=True)
load('/home/kep/magicKube/bremner16/engine.sage')

SALLOWS = [373^2, 289^2, 565^2, 360721, 425^2, 23^2, 205^2, 527^2, 222121]

print("=== 0. Проверка параметризации (3) Бремнера ===")
lines = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
a_, b_, c_ = QQ(7), QQ(-3), QQ(11)
vv = [cc*c_ + ca*a_ + cb*b_ for cc, ca, cb in COEF]
print("  все 8 линий равны:", len(set(sum(vv[i] for i in L) for L in lines)) == 1)
print("  квадрат Саллоуса магический:", len(set(sum(SALLOWS[i] for i in L) for L in lines)) == 1,
      " квадратные клетки:", [i for i in range(9) if ZZ(SALLOWS[i]).is_square()])

print("\n=== 1. F1 против corners/bremner_type.sage (b,h)=(31,17) ===")
# bremner_type: cell1 = b^2 (верх-середина), cell7 = h^2 (низ-середина), cell5 = f^2 (середина справа)
cells, pair, j, free = family_cells('F1', 31, 17)
nm = {0: 'TL', 2: 'TR', 3: 'D', 6: 'BL', 8: 'BR'}
print("  C =", cells[4][0], " пара =", pair, " третья =", j, " свободные =", free)
for k in free:
    print(f"    {nm[k]} = {cells[k][0]} + ({cells[k][1]})*t^2")
print("  ожидалось из bremner_type: TL=(h^2+f^2)/2, TR=(b^2+2h^2-f^2)/2, BL=(b^2+f^2)/2, "
      "BR=(2b^2+h^2-f^2)/2, D=b^2+h^2-f^2")
b, h = QQ(31), QQ(17)
exp = {0: (h^2/2, QQ(1)/2), 2: ((b^2+2*h^2)/2, -QQ(1)/2), 6: ((b^2)/2, QQ(1)/2),
       8: ((2*b^2+h^2)/2, -QQ(1)/2), 3: (b^2+h^2, QQ(-1))}
print("  совпадение:", all(cells[k] == exp[k] for k in free))

print("\n=== 2. Пример Бремнера–Саллоуса как точка семейства F1, (P,Q)=(17,31), t=23/17 ===")
cells, pair, j, free = family_cells('F1', 17, 31)
n, dist, pos, v, fl = score9(cells, QQ(23)/17)
print("  клетки:", v)
print("  квадраты:", fl, " всего:", n, " различны:", dist, " все>0:", pos)
print("  совпадает с (1) Бремнера с точностью до множителя 289:",
      [x*289 for x in v] == SALLOWS)

print("\n=== 3. Вырожденные t и доступные базовые пары для (P,Q)=(17,31) ===")
C = cells[4][0]
ts = degenerate_ts(cells, pair, j, C)
print("  C =", C, " квадрат:", C.is_square(), "  вырожденные t:", ts)
for t0 in ts:
    sq = [NAMES[k] for k in free if (cells[k][0] + cells[k][1]*t0^2) > 0
          and (cells[k][0] + cells[k][1]*t0^2).is_square()]
    print(f"    t0={t0}: квадратны свободные клетки {sq}")

print("\n=== 4. Поиск: находит ли движок пример Бремнера–Саллоуса? ===")
base_cnt = 3 + int(C.is_square())
found = []
for c1, c2 in itertools.combinations(free, 2):
    ok = False
    for t0 in ts:
        val = cells[c1][0] + cells[c1][1]*t0^2
        if val <= 0 or not val.is_square():
            continue
        try:
            line, hits = run_base(cells, free, base_cnt, c1, c2, t0, K=1500, nprimes=60, hcap=1500)
        except Exception as ex:
            print(f"  {NAMES[c1]}/{NAMES[c2]}@{t0}: ERROR {type(ex).__name__}: {str(ex)[:70]}")
            continue
        print("  " + line)
        for hh in hits:
            print("     HIT", hh[0], "квадратов, t =", hh[2], " положительны:", hh[1], " клетки:", hh[3])
            found.append(hh[2])
        ok = True
        break
    if not ok:
        print(f"  {NAMES[c1]}/{NAMES[c2]}: базовая кривая не построена")
print("\n  Пример Саллоуса (t = +-23/17) найден:", any(abs(t) == QQ(23)/17 for t in found))
