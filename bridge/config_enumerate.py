#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ТОЧНОЕ ПЕРЕЧИСЛЕНИЕ конфигураций квадратных клеток магического 3x3.

Считаем ПЕРЕБОРОМ, не рассуждением.  Чистый Python, без Sage.

Разметка клеток (построчно, индексы 0..8):

      0 1 2        a b c
      3 4 5   =    d e f
      6 7 8        g h i

Противоположные через центр пары (ассоциативность магического 3x3):
    УГЛОВЫЕ  : (0,8) = (a,i),  (2,6) = (c,g)
    РЁБЕРНЫЕ : (1,7) = (b,h),  (3,5) = (d,f)
Центр 4 = e ни в какую пару не входит.

Что считается:
  1. число ПОЛНЫХ пар (обе клетки в подмножестве);
  2. ТИПЫ полных пар (углов./рёбер./смешанно);
  3. входит ли центр;
  4. классы эквивалентности относительно D8 (8 элементов), представители;
  5. проверка ЛЕММЫ о двух полных парах перебором;
  6. сверка: 6 клеток <-> 16 конфигураций Бремнера, 7 клеток <-> 8 конфигураций Бойера;
  7. контроль Бёрнсайда на число орбит;
  8. контроль: переводит ли D8 рёберные пары в угловые (проверяется перебором, не рассуждением);
  9. привязка известного примера Бремнера-Саллоуса.

Запуск: python3 config_enumerate.py
"""

from itertools import combinations

CELLS = tuple(range(9))
CENTER = 4
CORNER_PAIRS = ((0, 8), (2, 6))
EDGE_PAIRS = ((1, 7), (3, 5))
ALL_PAIRS = CORNER_PAIRS + EDGE_PAIRS
PAIR_TYPE = {(0, 8): "C", (2, 6): "C", (1, 7): "E", (3, 5): "E"}
PAIR_NAME = {(0, 8): "(a,i)", (2, 6): "(c,g)", (1, 7): "(b,h)", (3, 5): "(d,f)"}
LETTER = "abcdefghi"


# ----------------------------------------------------------------------------- D8
def build_d8():
    """8 симметрий квадрата как перестановки индексов 0..8."""
    def apply(k, refl):
        g = []
        for i in CELLS:
            r, c = divmod(i, 3)
            if refl:
                r, c = c, r
            for _ in range(k):
                r, c = c, 2 - r          # поворот на 90 градусов
            g.append(3 * r + c)
        return tuple(g)
    gens = {apply(k, refl) for k in range(4) for refl in (False, True)}
    return sorted(gens)


D8 = build_d8()
assert len(D8) == 8, len(D8)
# группа: замкнутость и обратимость — проверяем перебором, не рассуждением
_S = set(D8)
for _g in D8:
    assert sorted(_g) == list(CELLS)
    for _h in D8:
        assert tuple(_g[_h[i]] for i in CELLS) in _S


def img(g, S):
    return frozenset(g[i] for i in S)


def orbit(S):
    return frozenset(img(g, S) for g in D8)


def stab_size(S):
    return sum(1 for g in D8 if img(g, S) == frozenset(S))


# ------------------------------------------------- контроль: D8 на четырёх парах
def pair_action_check():
    """Перебором: куда D8 переводит каждую из четырёх пар."""
    rows = []
    edge_to_corner = False
    corner_to_edge = False
    for p in ALL_PAIRS:
        images = set()
        for g in D8:
            q = tuple(sorted(img(g, p)))
            images.add(q)
        types = {PAIR_TYPE[q] for q in images}
        if PAIR_TYPE[p] == "E" and "C" in types:
            edge_to_corner = True
        if PAIR_TYPE[p] == "C" and "E" in types:
            corner_to_edge = True
        rows.append((p, sorted(images)))
    return rows, edge_to_corner, corner_to_edge


# --------------------------------------------------------------- анализ подмножества
def analyze(S):
    S = frozenset(S)
    full = [p for p in ALL_PAIRS if p[0] in S and p[1] in S]
    nc = sum(1 for p in full if PAIR_TYPE[p] == "C")
    ne = sum(1 for p in full if PAIR_TYPE[p] == "E")
    if nc and ne:
        kind = "смешанно"
    elif nc:
        kind = "только угловые"
    elif ne:
        kind = "только рёберные"
    else:
        kind = "нет полных пар"
    return {
        "set": S,
        "full": full,
        "n_full": len(full),
        "n_corner": nc,
        "n_edge": ne,
        "kind": kind,
        "center": CENTER in S,
        "combo": (nc, ne),
    }


def schema(S):
    r = []
    for row in range(3):
        r.append("".join("X" if 3 * row + c in S else "." for c in range(3)))
    return " / ".join(r)


def cellnames(S):
    return "".join(LETTER[i] for i in sorted(S))


def classes_of_size(k):
    """Все C(9,k) подмножеств, разбитые на орбиты D8.  Представитель — лексминимум."""
    seen, out = set(), []
    for S in combinations(CELLS, k):
        S = frozenset(S)
        if S in seen:
            continue
        orb = orbit(S)
        seen |= orb
        rep = min(tuple(sorted(x)) for x in orb)
        out.append({"rep": rep, "orbit": orb, "orbit_size": len(orb)})
    out.sort(key=lambda d: d["rep"])
    for i, d in enumerate(out, 1):
        d["idx"] = i
        d.update(analyze(d["rep"]))
    return out


def burnside(k):
    """Контроль числа орбит леммой Бёрнсайда — независимый счёт."""
    tot = 0
    per_g = []
    for g in D8:
        fix = sum(1 for S in combinations(CELLS, k) if img(g, S) == frozenset(S))
        per_g.append(fix)
        tot += fix
    assert tot % 8 == 0
    return tot // 8, per_g


# --------------------------------------------------------------------------- печать
def line(ch="-", n=100):
    print(ch * n)


def print_table(k, cls):
    print()
    line("=")
    print(f"РАЗМЕР {k}:  всего подмножеств C(9,{k}) = "
          f"{len(list(combinations(CELLS, k)))},  классов D8 = {len(cls)}")
    line("=")
    hdr = (f"{'#':>2} | {'представитель':<12} | {'схема':<13} | {'орбита':>6} | "
           f"{'полн.пар':>8} | {'типы (C,E)':<10} | {'какие пары':<26} | {'центр':<5} | тип")
    print(hdr)
    line()
    for d in cls:
        pairs = ",".join(PAIR_NAME[p] for p in d["full"]) or "—"
        print(f"{d['idx']:>2} | {cellnames(d['rep']):<12} | {schema(d['rep']):<13} | "
              f"{d['orbit_size']:>6} | {d['n_full']:>8} | "
              f"{'(%d,%d)' % (d['n_corner'], d['n_edge']):<10} | {pairs:<26} | "
              f"{'да' if d['center'] else 'нет':<5} | {d['kind']}")
    nb, per_g = burnside(k)
    ok = (nb == len(cls))
    print()
    print(f"  контроль Бёрнсайда: неподвижных по 8 элементам = {per_g}, "
          f"сумма/8 = {nb}  -> {'СОВПАЛО' if ok else 'РАСХОЖДЕНИЕ!'}")
    assert ok
    tot_orb = sum(d["orbit_size"] for d in cls)
    print(f"  контроль суммы орбит: {tot_orb} = C(9,{k}) = "
          f"{len(list(combinations(CELLS, k)))} -> "
          f"{'СОВПАЛО' if tot_orb == len(list(combinations(CELLS, k))) else 'РАСХОЖДЕНИЕ!'}")
    assert tot_orb == len(list(combinations(CELLS, k)))


def combo_breakdown(k, cls):
    """Разбиение классов по комбинации типов полных пар."""
    print()
    print(f"  --- размер {k}: комбинации типов полных пар (классы D8) ---")
    buckets = {}
    for d in cls:
        buckets.setdefault(d["combo"], []).append(d)
    names = {
        (0, 0): "нет полных пар",
        (0, 1): "{рёберная}",
        (1, 0): "{угловая}",
        (0, 2): "{рёберная, рёберная}",
        (1, 1): "{угловая, рёберная}  (смешанно)",
        (2, 0): "{угловая, угловая}",
        (1, 2): "три пары: 1 угловая + 2 рёберные",
        (2, 1): "три пары: 2 угловые + 1 рёберная",
        (2, 2): "четыре пары (все)",
        (0, 3): "невозможно (рёберных пар всего 2)",
    }
    for combo in sorted(buckets):
        ds = buckets[combo]
        nm = names.get(combo, str(combo))
        reps = ", ".join(cellnames(d["rep"]) for d in ds)
        n_sets = sum(d["orbit_size"] for d in ds)
        print(f"    (C={combo[0]}, E={combo[1]})  {nm:<36} : классов {len(ds):>2}, "
              f"подмножеств {n_sets:>3}  | представители: {reps}")
    for combo in [(0, 2), (1, 1), (2, 0), (1, 2), (2, 1), (2, 2)]:
        if combo not in buckets:
            print(f"    (C={combo[0]}, E={combo[1]})  {names.get(combo, ''):<36} : "
                  f"классов  0, подмножеств   0  | НЕВОЗМОЖНО")


# ------------------------------------------------------------------------------ main
def main():
    print("ТОЧНОЕ ПЕРЕЧИСЛЕНИЕ КОНФИГУРАЦИЙ (перебор, не рассуждение)")
    print("клетки:  0 1 2 / 3 4 5 / 6 7 8  =  a b c / d e f / g h i")
    print("пары: УГЛОВЫЕ (a,i)=(0,8), (c,g)=(2,6);  РЁБЕРНЫЕ (b,h)=(1,7), (d,f)=(3,5);  центр e=4")

    # 0. действие D8 на парах
    line("=")
    print("0. ДЕЙСТВИЕ D8 НА ЧЕТЫРЁХ ПАРАХ (перебор всех 8 симметрий)")
    line("=")
    rows, e2c, c2e = pair_action_check()
    for p, images in rows:
        print(f"  пара {PAIR_NAME[p]} = {p}  тип {PAIR_TYPE[p]}  ->  образы: "
              f"{', '.join(PAIR_NAME[q] for q in images)}")
    print(f"  рёберная -> угловая когда-нибудь? {'ДА' if e2c else 'НЕТ'}")
    print(f"  угловая -> рёберная когда-нибудь? {'ДА' if c2e else 'НЕТ'}")
    print("  ВЫВОД (проверено перебором): D8 сохраняет тип пары; крест в икс не переводится.")

    all_cls = {}
    for k in (6, 7, 8, 9):
        cls = classes_of_size(k)
        all_cls[k] = cls
        print_table(k, cls)
        combo_breakdown(k, cls)

    # 1. ЛЕММА
    line("=")
    print("1. ПРОВЕРКА ЛЕММЫ О ДВУХ ПОЛНЫХ ПАРАХ — ПОЛНЫЙ ПЕРЕБОР")
    line("=")
    counter = []
    for k in (7, 8, 9):
        for S in combinations(CELLS, k):
            a = analyze(S)
            if a["n_full"] < 2:
                counter.append((k, S, a["n_full"]))
    print(f"  проверено подмножеств размера >= 7: "
          f"{sum(len(list(combinations(CELLS,k))) for k in (7,8,9))}")
    if counter:
        print("  !!! КОНТРПРИМЕРЫ НАЙДЕНЫ — ЛЕММА ЛОЖНА !!!")
        for k, S, n in counter:
            print(f"      размер {k}: {cellnames(S)} — полных пар {n}")
    else:
        print("  КОНТРПРИМЕРОВ НЕТ: у каждого подмножества размера >= 7 полных пар >= 2.")
        print("  ЛЕММА ПОДТВЕРЖДЕНА ПЕРЕБОРОМ.")
    # точное распределение
    for k in (6, 7, 8, 9):
        dist = {}
        for S in combinations(CELLS, k):
            dist[analyze(S)["n_full"]] = dist.get(analyze(S)["n_full"], 0) + 1
        print(f"  размер {k}: распределение числа полных пар "
              f"{ {n: dist[n] for n in sorted(dist)} }  (мин = {min(dist)})")
    # усиление: минимум по размеру 6
    m6 = min(analyze(S)["n_full"] for S in combinations(CELLS, 6))
    print(f"  ДОБАВКА: размер 6 — минимум полных пар = {m6} "
          f"(то есть 'две пары' для шести клеток НЕ гарантированы).")

    # 2. сверка с Бремнером / Бойером
    line("=")
    print("2. СВЕРКА С ЛИТЕРАТУРОЙ")
    line("=")
    print(f"  размер 6: классов D8 = {len(all_cls[6])};  Бремнер I–XVI = 16  -> "
          f"{'СОВПАДАЕТ' if len(all_cls[6]) == 16 else 'РАСХОЖДЕНИЕ'}")
    print(f"  размер 7: классов D8 = {len(all_cls[7])};  Бойер 7.I–7.VIII = 8  -> "
          f"{'СОВПАДАЕТ' if len(all_cls[7]) == 8 else 'РАСХОЖДЕНИЕ'}")
    print(f"  размер 8: классов D8 = {len(all_cls[8])};  Бойер 8.I–8.III = 3  -> "
          f"{'СОВПАДАЕТ' if len(all_cls[8]) == 3 else 'РАСХОЖДЕНИЕ'}")
    print(f"  размер 9: классов D8 = {len(all_cls[9])} (тривиально 1)")

    # 3. покрытие G-семействами и теоремой об отношениях
    line("=")
    print("3. ПОКРЫТИЕ G-СЕМЕЙСТВАМИ (G1 = две рёберные, G2 = две угловые, G3 = рёберная+угловая)")
    print("   и теоремой об исключении отношений (она СФОРМУЛИРОВАНА ТОЛЬКО ДЛЯ РЁБЕРНЫХ ПАР)")
    line("=")
    for k in (7, 6):
        print(f"  --- размер {k} ---")
        for d in all_cls[k]:
            fams = []
            if d["n_edge"] >= 2:
                fams.append("G1")
            if d["n_corner"] >= 2:
                fams.append("G2")
            if d["n_corner"] >= 1 and d["n_edge"] >= 1:
                fams.append("G3")
            cov = "ДА" if d["n_edge"] >= 1 else "НЕТ (нет ни одной полной рёберной пары)"
            print(f"    класс {d['idx']:>2} {cellnames(d['rep']):<12} {schema(d['rep']):<13} "
                  f"пар (C,E)=({d['n_corner']},{d['n_edge']})  "
                  f"G-семейства: {','.join(fams) if fams else '—':<10} "
                  f"| есть рёберная пара: {cov}")
        n_no_edge = sum(1 for d in all_cls[k] if d["n_edge"] == 0)
        n_g1 = sum(1 for d in all_cls[k] if d["n_edge"] >= 2)
        n_g2 = sum(1 for d in all_cls[k] if d["n_corner"] >= 2)
        n_g3 = sum(1 for d in all_cls[k] if d["n_corner"] >= 1 and d["n_edge"] >= 1)
        print(f"    ИТОГО размер {k}: классов {len(all_cls[k])}; "
              f"покрыто G1 = {n_g1}, G2 = {n_g2}, G3 = {n_g3}; "
              f"без единой полной рёберной пары = {n_no_edge}")

    # 3bis. привязка к нумерации Бремнера и к расщеплению Бойера
    line("=")
    print("3bis. ПРИВЯЗКА К ЛИТЕРАТУРНЫМ НУМЕРАЦИЯМ")
    line("=")
    # Нумерация классов размера 6 (лексминимум представителя) совпадает с таблицей §1
    # BREMNER16_2026-09-12.md, где привязка I..XVI к классам сделана расшифровкой Fig.1.
    BREMNER = {1: "V", 2: "I", 3: "X", 4: "VIII", 5: "IX", 6: "XII", 7: "XVI", 8: "XV",
               9: "VII", 10: "XIV", 11: "III", 12: "VI", 13: "XI", 14: "II",
               15: "IV", 16: "XIII"}
    print("  размер 6 -> Бремнер (метки взяты из BREMNER16_2026-09-12.md §1, там привязка сделана")
    print("  расшифровкой Fig.1; здесь сверяю только число полных пар — независимый контроль):")
    B16_PAIRS = {1: 1, 2: 1, 3: 1, 4: 1, 5: 2, 6: 2, 7: 2, 8: 2, 9: 2, 10: 2,
                 11: 3, 12: 2, 13: 2, 14: 3, 15: 1, 16: 1}
    bad = []
    for d in all_cls[6]:
        ok = (d["n_full"] == B16_PAIRS[d["idx"]])
        if not ok:
            bad.append(d["idx"])
        print(f"    класс {d['idx']:>2} = Бремнер {BREMNER[d['idx']]:<5} "
              f"{cellnames(d['rep']):<8} полных пар {d['n_full']} "
              f"(в BREMNER16: {B16_PAIRS[d['idx']]}) {'ok' if ok else 'РАСХОЖДЕНИЕ'}")
    print(f"  столбец 'полных пар' сверен со старой таблицей: "
          f"{'все 16 совпали' if not bad else 'РАСХОЖДЕНИЯ в ' + str(bad)}")
    assert not bad

    print()
    print("  размер 7 -> Бойер: он делит восемь конфигураций на 7.I–7.VI (КВАДРАТНЫЙ центр)")
    print("  и 7.VII–7.VIII (неквадратный центр). Проверяю расщепление перебором:")
    with_c = [d["idx"] for d in all_cls[7] if d["center"]]
    without_c = [d["idx"] for d in all_cls[7] if not d["center"]]
    print(f"    с центром   : классы {with_c}  -> {len(with_c)} штук (у Бойера 6)")
    print(f"    без центра  : классы {without_c}  -> {len(without_c)} штук (у Бойера 2)")
    print(f"    расщепление 6+2: "
          f"{'СОВПАДАЕТ с 7.I–7.VI / 7.VII–7.VIII' if (len(with_c), len(without_c)) == (6, 2) else 'РАСХОЖДЕНИЕ'}")

    print()
    print("  'песочные часы' (hourglass, три линии через центр = дополнение есть ПОЛНАЯ пара):")
    for p in ALL_PAIRS:
        S = frozenset(CELLS) - frozenset(p)
        for d in all_cls[7]:
            if S in d["orbit"]:
                print(f"    дыра = пара {PAIR_NAME[p]} ({PAIR_TYPE[p]}): подмножество "
                      f"{cellnames(S)} -> класс {d['idx']} "
                      f"(C,E)=({d['n_corner']},{d['n_edge']})")

    # 4. известный пример Бремнера-Саллоуса
    line("=")
    print("4. ПРИВЯЗКА ИЗВЕСТНОГО ПРИМЕРА (Бремнер–Саллоус)")
    line("=")
    known = [373**2, 289**2, 565**2,
             360721, 425**2, 23**2,
             205**2, 527**2, 222121]
    from math import isqrt
    sq = [i for i, v in enumerate(known) if isqrt(v) ** 2 == v]
    print(f"  квадратные клетки: {cellnames(sq)} = {sorted(sq)};  неквадратные: "
          f"{cellnames(set(CELLS) - set(sq))}")
    # магичность — контроль
    s = sum(known[0:3])
    lines_ = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
    print(f"  контроль магичности: все 8 линий = {s}? "
          f"{all(sum(known[i] for i in L) == s for L in lines_)}")
    a = analyze(sq)
    print(f"  полных пар: {a['n_full']}  типы (C,E) = ({a['n_corner']},{a['n_edge']})  "
          f"пары: {', '.join(PAIR_NAME[p] for p in a['full'])}")
    for d in all_cls[7]:
        if frozenset(sq) in d["orbit"]:
            print(f"  класс D8 размера 7: #{d['idx']}, представитель {cellnames(d['rep'])} "
                  f"({schema(d['rep'])}), орбита {d['orbit_size']}")
    # какие шестёрки Бремнера он реализует
    real6 = set()
    for six in combinations(sq, 6):
        for d in all_cls[6]:
            if frozenset(six) in d["orbit"]:
                real6.add(d["idx"])
    print(f"  реализуемые им классы шести клеток: {sorted(real6)} "
          f"(всего {len(real6)} из 16)")

    # 5. семейство G1 из задания: F0,F4,F8
    line("=")
    print("5. СЕМЕЙСТВО G1 ИЗ ЗАДАНИЯ (u0^2=F0, u4^2=F4, u8^2=F8) — какие клетки квадратны")
    line("=")
    g1_auto = {1, 7, 3, 5}                       # две рёберные пары — автоматические квадраты
    g1_curve = {0, 4, 8}                         # три условия кривой C
    S = g1_auto | g1_curve
    a = analyze(S)
    print(f"  автоматические (две рёберные пары): {cellnames(g1_auto)}")
    print(f"  кривая C навешивает: {cellnames(g1_curve)}")
    print(f"  итого квадратных клеток: {len(S)} = {cellnames(S)}; "
          f"неквадратные: {cellnames(set(CELLS) - S)}")
    print(f"  полных пар {a['n_full']}, типы (C,E) = ({a['n_corner']},{a['n_edge']}): "
          f"{', '.join(PAIR_NAME[p] for p in a['full'])}")
    for d in all_cls[7]:
        if frozenset(S) in d["orbit"]:
            print(f"  класс D8 размера 7: #{d['idx']}, представитель {cellnames(d['rep'])}")
    print("  ВЫВОД: G1 = ДВЕ РЁБЕРНЫЕ пары (крест). Угловых пар в его основании нет.")

    # 6. полный список 36 семиклеточных подмножеств по классам
    line("=")
    print("6. ВСЕ 36 СЕМИКЛЕТОЧНЫХ ПОДМНОЖЕСТВ ПОИМЁННО (по классам, через дополнение)")
    line("=")
    for d in all_cls[7]:
        items = []
        for S in sorted(d["orbit"], key=lambda x: tuple(sorted(x))):
            comp = sorted(set(CELLS) - set(S))
            items.append(f"{cellnames(S)}(дыры {cellnames(comp)})")
        print(f"  класс {d['idx']} [{d['orbit_size']} шт, (C,E)=({d['n_corner']},{d['n_edge']}), "
              f"центр {'да' if d['center'] else 'нет'}]:")
        for j in range(0, len(items), 4):
            print("      " + "  ".join(items[j:j + 4]))

    # 7. итоговая сводка покрытия теоремой
    line("=")
    print("7. ИТОГ: ПОКРЫТИЕ ТЕОРЕМОЙ ОБ ИСКЛЮЧЕНИИ ОТНОШЕНИЙ (только рёберные пары)")
    line("=")
    print("  Девять клеток (полное решение): 4 полные пары = 2 угловые + 2 рёберные.")
    print("  Теорема ограничивает 2 из 4 пар (рёберные). Угловые (a,i) и (c,g) НЕ ограничены ничем.")
    print()
    print("  Семь клеток (приз Бойера), 8 классов D8:")
    cov = [d["idx"] for d in all_cls[7] if d["n_edge"] >= 1]
    nocov = [d["idx"] for d in all_cls[7] if d["n_edge"] == 0]
    g1cov = [d["idx"] for d in all_cls[7] if d["n_edge"] >= 2]
    print(f"    есть хотя бы одна полная РЁБЕРНАЯ пара (рёберный аргумент применим): "
          f"классы {cov} — {len(cov)} из 8")
    print(f"    полных рёберных пар НЕТ вообще (рёберный аргумент неприменим в принципе): "
          f"классы {nocov} — {len(nocov)} из 8")
    print(f"    ДВЕ полные рёберные пары (основание семейства G1, где вся арифметика проекта): "
          f"классы {g1cov} — {len(g1cov)} из 8")
    for i in nocov:
        d = all_cls[7][i - 1]
        print(f"    класс {i} = {cellnames(d['rep'])} ({schema(d['rep'])}), дыры "
              f"{cellnames(set(CELLS) - set(d['rep']))} — только угловые пары.")

    # 8. ШАГИ ПАР: какие тройки пар связаны соотношением x3 = x1 + x2 (перебор, не алгебра)
    line("=")
    print("8. ШАГИ ПАР И СООТНОШЕНИЕ x3 = x1 + x2 (перебор по случайным магическим квадратам)")
    line("=")
    import random
    rnd = random.Random(20260912)

    def build(e, x, y):
        """Общий магический 3x3 с центром e: два свободных параметра x, y."""
        A, B, C, D = y, -(x + y), x, x - y      # шаги пар (a,i),(b,h),(c,g),(d,f)
        q = [0] * 9
        q[0], q[8] = e + A, e - A               # a, i
        q[1], q[7] = e + B, e - B               # b, h
        q[2], q[6] = e + C, e - C               # c, g
        q[3], q[5] = e + D, e - D               # d, f
        q[4] = e
        return q, {(0, 8): A, (1, 7): B, (2, 6): C, (3, 5): D}

    lines_ = [(0,1,2),(3,4,5),(6,7,8),(0,3,6),(1,4,7),(2,5,8),(0,4,8),(2,4,6)]
    samples = []
    for _ in range(400):
        e = rnd.randint(10, 10**6)
        x = rnd.randint(-10**5, 10**5)
        y = rnd.randint(-10**5, 10**5)
        q, st = build(e, x, y)
        assert all(sum(q[i] for i in L) == 3 * e for L in lines_), "не магический!"
        samples.append(st)
    print(f"  построено {len(samples)} случайных магических квадратов, магичность всех 8 линий — ok")
    print("  (общий вид: шаги пар (a,i)=y, (b,h)=-(x+y), (c,g)=x, (d,f)=x-y; x, y свободны)")
    print()
    print("  Для каждой ТРОЙКИ пар: есть ли тождество e1*s1 + e2*s2 + e3*s3 = 0 с e_i = ±1")
    print("  (то есть 'один шаг = сумма двух других' — условие Бьюэлла/Моргенштерна x3 = x1+x2):")
    from itertools import product as iproduct
    for trip in combinations(ALL_PAIRS, 3):
        found = None
        for eps in iproduct((1, -1), repeat=3):
            if all(sum(e_ * st[p] for e_, p in zip(eps, trip)) == 0 for st in samples):
                found = eps
                break
        types = "".join(PAIR_TYPE[p] for p in trip)
        nm = ", ".join(PAIR_NAME[p] for p in trip)
        # какой класс размера 7 отвечает этой тройке полных пар
        holes = tuple(sorted(set(CELLS) - {i for p in trip for i in p} - {CENTER}))
        S = frozenset(CELLS) - frozenset(holes)
        kls = [d["idx"] for d in all_cls[7] if S in d["orbit"]]
        # соотношение с коэффициентом 2, если ±1 не нашлось
        found2 = None
        if not found:
            for eps in iproduct((1, -1), repeat=3):
                for k2 in range(3):
                    co = [e_ * (2 if j == k2 else 1) for j, e_ in enumerate(eps)]
                    if all(sum(c * st[p] for c, p in zip(co, trip)) == 0 for st in samples):
                        found2 = tuple(co)
                        break
                if found2:
                    break
        rel = (f"ЕСТЬ ±1: {found}" if found
               else (f"НЕТ ±1; есть с коэффициентом 2: {found2}" if found2 else "нет никакого"))
        print(f"    {types}  {nm:<28} класс 7-клеточный {kls}  ->  {rel}")
    print()
    print("  ЧТЕНИЕ: 'x3 = x1 + x2' (Бьюэлл, Моргенштерн) выполняется РОВНО для троек,")
    print("  содержащих ОБЕ угловые пары. Это класс 7 = abceghi = XXX / .X. / XXX —")
    print("  форма буквальных 'песочных часов'. Класс 8 = abdefhi = XX. / XXX / .XX")
    print("  (крест + угловая пара) тоже даёт три линии через центр, но соотношение там")
    print("  с коэффициентом 2, а не x3 = x1 + x2.")


if __name__ == "__main__":
    main()
