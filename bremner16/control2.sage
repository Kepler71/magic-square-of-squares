# Контроли движка-2 (F- и G-семейства).
import functools, sys, itertools
print = functools.partial(print, flush=True)
sys.argv = ['x']
load('/home/kep/magicKube/bremner16/engine2.sage')

LINES = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]

print("=== A. G-семейства: магичность и автоматические квадраты ===")
for gt in GFAM:
    cells, auto, free = cells_G(gt, 5, 2)
    ok = len(set(sum((cells[i] for i in L), Rp(0)) for L in LINES)) == 1
    sq = [k for k in auto]
    print(f"  {gt}: пары {GFAM[gt]}  магический: {ok}  авто-квадраты: {sq}  свободные: {free}"
          f"  степени клеток: {[cells[k].degree() for k in range(9)]}")

print("\n=== B. F-семейства через движок-2: контроль Бремнера–Саллоуса (F1, (P,Q)=(17,31), t=23/17) ===")
cells, auto, free = cells_F('F1', 17, 31)
n, dist, pos, v, fl = score9(cells, QQ(23)/17)
print("  клетки:", v, "\n  квадратов:", n, " различны:", dist, " положительны:", pos)

print("\n=== C. G3: пример Бремнера–Саллоуса как точка (m,n)=(61,-23), p=7/6 ===")
cells, auto, free = cells_G('G3', 61, -23)
n, dist, pos, v, fl = score9(cells, QQ(7)/6)
print("  клетки:", v)
print("  x36:", [x*36 for x in v])
print("  квадратов:", n, " различны:", dist, " положительны:", pos, " флаги:", fl)
print("  авто:", auto, " свободные:", free, " центр квадрат:", v[4].is_square())

print("\n=== D. Поиск в G3 при (m,n)=(61,-23): находит ли движок пример? ===")
ts = degenerate_ts(cells)
print("  вырожденные p:", ts)
base_cnt = 4
found = []
for c1, c2 in itertools.combinations(free + [4], 2):
    ok = False
    for t0 in ts:
        for (x, y) in ((c1, c2), (c2, c1)):
            val = cells[x](t0)
            if val <= 0 or not val.is_square():
                continue
            try:
                line, hits = run_base(cells, free + [4], base_cnt, x, y, t0, K=1200, nprimes=50, hcap=1200)
            except Exception as ex:
                continue
            print("  " + line)
            for hh in hits:
                print("     HIT", hh[0], "квадратов, p =", hh[2], " клетки:", hh[3])
                found.append(hh[2])
            ok = True
            break
        if ok:
            break
    if not ok:
        print(f"  {NAMES[c1]}/{NAMES[c2]}: кривая не построена")
print("\n  Пример Саллоуса (p = 7/6) найден:", QQ(7)/6 in found or QQ(-7)/6 in found)
