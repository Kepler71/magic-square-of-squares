#!/usr/bin/env python3
# 1) Независимое перечисление всех C(9,6)=84 шестиэлементных подмножеств клеток 3x3 с точностью до D8.
# 2) Расшифровка Fig.1 из Bremner, Acta Arith. 99 (2001): pdftotext теряет верхние индексы, поэтому
#    токен "492" — это либо 492, либо 49^2.  Перебираем 2^9 прочтений, оставляем те, что дают
#    магический квадрат ровно с 6 квадратными клетками.  Так конфигурации I..XVI привязываются к нашим классам.
from itertools import product
from math import isqrt

CELLS = list(range(9))          # 0..8, построчно


def d8():
    """8 симметрий квадрата как перестановки индексов 0..8 (позиция i -> позиция g[i])."""
    def rc(i):
        return divmod(i, 3)
    gens = []
    for k in range(4):
        for refl in (False, True):
            g = []
            for i in CELLS:
                r, c = rc(i)
                if refl:
                    r, c = c, r
                for _ in range(k):          # поворот на 90 градусов
                    r, c = c, 2 - r
                g.append(3 * r + c)
            gens.append(tuple(g))
    return sorted(set(gens))


G = d8()
assert len(G) == 8


def orbit(S):
    return frozenset(frozenset(g[i] for i in S) for g in G)


def classify():
    seen, classes = set(), []
    for S in product(*[[0, 1]] * 9):
        if sum(S) != 6:
            continue
        Sf = frozenset(i for i in CELLS if S[i])
        if Sf in seen:
            continue
        orb = orbit(Sf)
        seen |= orb
        classes.append((min(sorted(tuple(sorted(x))) for x in orb), orb))
    classes.sort(key=lambda t: t[0])
    return classes


# ---- Fig. 1: токены как их выдал pdftotext (построчно) ----
FIG1 = {
 'I':    "492 1432 1552 1932 1252 -10999 852 10801 28849",
 'II':   "1945 12 372 232 1105 412 292 472 265",
 'III':  "5412 4212 492 -132839 157441 447721 5592 3712 1492",
 'IV':   "93961 1912 43801 892 2412 3292 2692 79681 1492",
 'V':    "1153 1057 313 12 292 412 372 252 232",
 'VI':   "3385 472 7081 892 652 232 372 792 5065",
 'VII':  "889 697 172 52 252 352 312 553 192",
 'VIII': "1561 312 12 -719 292 492 412 721 112",
 'IX':   "2713 673 352 72 1537 552 432 492 192",
 'X':    "3001 -1679 612 492 412 312 -359 712 192",
 'XI':   "10585 -1679 1132 972 852 712 412 1272 3865",
 'XII':  "22009 1192 9265 492 15145 1672 1452 1272 912",
 'XIII': "313 232 412 472 292 -527 12 1153 372",
 'XIV':  "52 1561 172 889 252 192 312 -311 352",
 'XV':   "265 12 132 72 145 241 112 172 52",
 'XVI':  "372 5089 672 6769 3649 232 532 472 772",
}


def readings(tok):
    """Возможные значения токена: само число и (если оканчивается на '2') квадрат отсечённого числа."""
    out = [(int(tok), False)]
    if tok.endswith('2') and len(tok.lstrip('-')) > 1:
        base = tok[:-1]
        if base not in ('', '-'):
            out.append((int(base) ** 2, True))
    return out


def is_magic(v):
    S = sum(v[:3])
    lines = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]
    return all(sum(v[i] for i in L) == S for L in lines)


def issq(n):
    return n > 0 and isqrt(n) ** 2 == n


def decode(name, s, allow_wild=True):
    toks = s.split()
    sols = []
    for combo in product(*[readings(t) for t in toks]):
        v = [c[0] for c in combo]
        if not is_magic(v):
            continue
        sq = frozenset(i for i in range(9) if issq(v[i]))
        sols.append((tuple(v), sq, tuple(c[1] for c in combo)))
    if sols or not allow_wild:
        return sols, None
    # ни одно прочтение не магическое: считаем, что один токен искажён при извлечении текста.
    for bad in range(9):
        for combo in product(*[readings(t) if i != bad else [(None, None)]
                               for i, t in enumerate(toks)]):
            v = [c[0] for c in combo]
            known = [i for i in range(9) if v[i] is not None]
            # восстановить недостающую клетку из магического условия
            lines = [(0, 1, 2), (3, 4, 5), (6, 7, 8), (0, 3, 6), (1, 4, 7), (2, 5, 8), (0, 4, 8), (2, 4, 6)]
            L = next((L for L in lines if bad in L and all(x in known for x in L if x != bad)), None)
            if L is None:
                continue
            S = None
            for L2 in lines:
                if all(x in known for x in L2):
                    S = sum(v[x] for x in L2)
                    break
            if S is None:
                continue
            v[bad] = S - sum(v[x] for x in L if x != bad)
            if not is_magic(v):
                continue
            sq = frozenset(i for i in range(9) if issq(v[i]))
            sols.append((tuple(v), sq, None, bad))
    return sols, 'wild'


if __name__ == '__main__':
    classes = classify()
    print("Число классов D8 шестиэлементных подмножеств:", len(classes))
    reps = {}
    for k, (rep, orb) in enumerate(classes, 1):
        reps[frozenset(rep)] = k
        print(f"  класс {k:2d}: клетки {sorted(rep)}  дополнение {sorted(set(CELLS)-set(rep))}"
              f"  |орбита|={len(orb)}")
    print()
    print("Расшифровка Fig. 1:")
    lab2cls = {}
    for name, s in FIG1.items():
        sols, mode = decode(name, s)
        good = [x for x in sols if len(x[1]) == 6]
        sqsets = set(x[1] for x in good)
        if len(sqsets) == 1:
            sq = next(iter(sqsets))
            cls = next(k for r, k in reps.items() if frozenset(sq) in orbit(r) or sq == r
                       or frozenset(sq) in orbit(r))
            # надёжнее: найти класс по орбите
            cls = None
            for r, k in reps.items():
                if frozenset(sq) in orbit(r):
                    cls = k
                    break
            lab2cls[name] = (cls, tuple(sorted(sq)), good[0][0], mode)
            print(f"  {name:5s} -> класс {cls:2d}, квадратные клетки {sorted(sq)}"
                  f"{'  (одна клетка восстановлена)' if mode else ''}   {good[0][0]}")
        else:
            print(f"  {name:5s} -> НЕОДНОЗНАЧНО/НЕТ: {len(good)} решений, множества {sqsets}")
    print()
    got = sorted(v[0] for v in lab2cls.values())
    print("Классы, покрытые Fig.1:", got, "  все 16 различны:", len(set(got)) == 16)
