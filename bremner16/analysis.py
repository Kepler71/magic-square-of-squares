#!/usr/bin/env python3
# Структурный разбор 16 конфигураций и их привязка к (тип семейства, базовая пара) нашего движка.
from itertools import combinations
from configs import d8, orbit, classify, FIG1, decode, CELLS

G = d8()
classes = classify()
REPS = [frozenset(rep) for rep, orb in classes]
ORBS = [orb for rep, orb in classes]


def cls_of(S):
    S = frozenset(S)
    for k, orb in enumerate(ORBS, 1):
        if S in orb:
            return k
    return None


LAB = {}
for name, s in FIG1.items():
    sols, mode = decode(name, s)
    good = [x for x in sols if len(x[1]) == 6]
    LAB[cls_of(good[0][1])] = name

LINES = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]
PAIRS = [(0, 8), (1, 7), (2, 6), (3, 5)]
FAMILIES = {'F1': ((1, 7), 5), 'F2': ((1, 7), 0), 'F3': ((0, 8), 2), 'F4': ((0, 8), 1)}
SALLOWS_SQ = frozenset([0, 1, 2, 4, 5, 6, 7])       # квадратные клетки примера (1) Бремнера

# --- какие (семейство, базовая пара [, «везучая» клетка]) обслуживают каждый класс ---
serve_center = {k: [] for k in range(1, 17)}   # 6 гарантированных (центр в конфигурации, C = квадрат)
serve_nocent = {k: [] for k in range(1, 17)}   # 5 гарантированных + 1 везение
for ft, ((i, ip), j) in FAMILIES.items():
    free = [k for k in range(9) if k not in (i, ip, j, 4)]
    for c1, c2 in combinations(free, 2):
        S = frozenset({i, ip, j, 4, c1, c2})
        serve_center[cls_of(S)].append((ft, (c1, c2)))
        for x in free:
            if x in (c1, c2):
                continue
            S2 = frozenset({i, ip, j, c1, c2, x})
            if len(S2) == 6:
                serve_nocent[cls_of(S2)].append((ft, (c1, c2), x))

print("Конфигурация | Бремнер | клетки | дополнение | центр | полных пар | полных линий | "
      "в примере (1)? | обслуживают")
for k, (rep, orb) in enumerate(classes, 1):
    S = set(rep)
    comp = sorted(set(CELLS) - S)
    npairs = sum(1 for p in PAIRS if set(p) <= S)
    nlines = sum(1 for L in LINES if set(L) <= S)
    known = any(frozenset(x) <= SALLOWS_SQ for x in orb)
    sv = serve_center[k] if 4 in S else serve_nocent[k]
    fams = sorted(set(x[0] for x in sv))
    print(f"  {k:2d} | {LAB[k]:5s} | {sorted(S)} | {comp} | {'да' if 4 in S else 'нет':3s} | "
          f"{npairs} | {nlines} | {'ДА' if known else 'нет':3s} | {len(sv)} вариантов, семейства {fams}")

print()
print("Детально — по одному варианту на класс (гарантированные квадраты подчёркнуты):")
for k, (rep, orb) in enumerate(classes, 1):
    S = set(rep)
    sv = serve_center[k] if 4 in S else serve_nocent[k]
    if not sv:
        print(f"  {k:2d} ({LAB[k]}): НЕТ обслуживающего варианта!")
        continue
    v = sv[0]
    ft = v[0]; (i, ip), j = FAMILIES[ft]
    if 4 in S:
        print(f"  {k:2d} ({LAB[k]:5s}): {ft} пара{{{i},{ip}}}=P^2,Q^2  третья c{j}=t^2  центр c4=C=□  "
              f"базовая пара {{c{v[1][0]},c{v[1][1]}}}  -> 6 гарантированных; "
              f"7-я клетка из {sorted(set(CELLS)-S)}")
    else:
        print(f"  {k:2d} ({LAB[k]:5s}): {ft} пара{{{i},{ip}}}=P^2,Q^2  третья c{j}=t^2  "
              f"базовая пара {{c{v[1][0]},c{v[1][1]}}}  + «везение» c{v[2]} -> 5 гарантированных + 1; "
              f"7-я клетка из {sorted(set(CELLS)-S)}")

print()
print("Проверка полноты: каждый класс обслуживается хотя бы одним вариантом:",
      all((serve_center[k] if 4 in set(REPS[k-1]) else serve_nocent[k]) for k in range(1, 17)))

# --- шесть-подмножества примера (1) ---
print()
print("Классы шестиэлементных подмножеств известного примера (1) Бремнера–Саллоуса:")
for sub in combinations(sorted(SALLOWS_SQ), 6):
    print(f"  {list(sub)} -> класс {cls_of(sub)} ({LAB[cls_of(sub)]})")
