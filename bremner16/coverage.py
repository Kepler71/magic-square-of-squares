#!/usr/bin/env python3
# Покрытие по конфигурациям: для каждого класса — сколько базовых кривых реально построено
# в логах (по строкам "  G1 c2/c3@..." / "  F1 c0/c2@...") и сколько пропущено.
import sys, re
from itertools import combinations
from configs import classify, CELLS

classes = classify()


def cls_of(S):
    S = frozenset(S)
    for k, (rep, orb) in enumerate(classes, 1):
        if S in orb:
            return k


FFAM = {'F1': ((1, 7), 5), 'F2': ((1, 7), 0), 'F3': ((0, 8), 2), 'F4': ((0, 8), 1)}
GFAM = {'G1': ((1, 7), (3, 5)), 'G2': ((0, 8), (2, 6)), 'G3': ((1, 7), (0, 8))}
LAB = {1: 'V', 2: 'I', 3: 'X', 4: 'VIII', 5: 'IX', 6: 'XII', 7: 'XVI', 8: 'XV', 9: 'VII',
       10: 'XIV', 11: 'III', 12: 'VI', 13: 'XI', 14: 'II', 15: 'IV', 16: 'XIII'}


def serving(fam, c1, c2, csq):
    """Класс конфигурации, все шесть клеток которой гарантированы этой (семейство, база).
       Для F при неквадратном центре гарантированных только 5 — возвращаем список из трёх классов
       (5 гарантированных + одна из трёх оставшихся свободных клеток)."""
    if fam in GFAM:
        (i, ip), (k, kp) = GFAM[fam]
        return [cls_of({i, ip, k, kp, c1, c2})]
    (i, ip), j = FFAM[fam]
    if csq:
        return [cls_of({i, ip, j, 4, c1, c2})]
    free = [x for x in range(9) if x not in (i, ip, j, 4)]
    return [cls_of({i, ip, j, c1, c2, x}) for x in free if x not in (c1, c2)]


built = {k: 0 for k in range(1, 17)}
missed = {k: 0 for k in range(1, 17)}
pat_ok = re.compile(r'^  ([FG]\d) c(\d)/c(\d)@')
pat_no = re.compile(r'^  ([FG]\d) c(\d)/c(\d): (?:кривая не построена|нет рацион)')
csq = True
for fn in sys.argv[1:]:
    for line in open(fn):
        if line.startswith('=='):
            csq = ('SQ' in line) or ('m,n' in line)
        m = pat_ok.match(line) or pat_no.match(line)
        if not m:
            continue
        fam, a, b = m.group(1), int(m.group(2)), int(m.group(3))
        tgt = built if pat_ok.match(line) else missed
        for k in serving(fam, a, b, csq):
            if k:
                tgt[k] += 1

print("| класс | Бремнер | построенных базовых кривых | пропущено (нет точки) |")
print("|---|---|---|---|")
for k in range(1, 17):
    print(f"| {k} | {LAB[k]} | {built[k]} | {missed[k]} |")
