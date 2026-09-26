# -*- coding: utf-8 -*-
"""
Полный перебор вычетов mod p^k (k = 1..KFLAT) + адаптивное доуточнение неопределённых узлов.
Запуск:  python3 p23_enum.py p M KFLAT KMAX
  p     -- 2 или 3
  M     -- максимальный порядок полюса min(v(b),v(c)) = -m, m = 0..M
  KFLAT -- до какого k делать ПЛОСКИЙ полный перебор (по условию задачи 8)
  KMAX  -- граница глубины адаптивного дерева (явная)
Все циклы ограничены; бюджет узлов дерева NODE_BUDGET.
"""
import sys, time, json
import numpy as np
from p23_masks import (CELLS, LINES, CENTER, NONC, GBITS, THR, mask_of_residue,
                       build_solution_list, selftest_SU, brute_all_G8)

NODE_BUDGET = 5_000_000
t0 = time.time()
p, M, KFLAT, KMAX = (int(a) for a in sys.argv[1:5])
t = THR[p]
G = 1 << GBITS[p]
FULL = (1 << G) - 1

# --- решатель: список решений 8 линий (независимо от S,U) ---
sols = build_solution_list(p)
ok_su, nsol, nsu = selftest_SU(p, sols)
print(f'p={p}: решений линий {nsol}, образ S,U {nsu}, совпадают: {ok_su}', flush=True)
if p == 3:
    bf = brute_all_G8(p)
    print(f'  контроль полным перебором G^8 без отсечения: {len(bf)} решений, совпадает: {bf == set(sols)}', flush=True)
    assert bf == set(sols)
assert ok_su
SOLS = np.array(sols, dtype=np.int64)          # (nsol, 9)


def vec_masks(N0, k, m):
    """векторная маска для массива N0 (значения mod p^k)"""
    q = p ** k
    N0 = N0 % q
    zero = N0 == 0
    w = np.zeros_like(N0)
    Nw = N0.copy()
    for _ in range(k):                     # явная граница
        d = (Nw % p == 0) & (~zero) & (w < k)
        w = w + d
        Nw = np.where(d, Nw // p, Nw)
    u = Nw
    if p == 2:
        prec = k - w
        c = (w & 1) | (((u & 3) == 3).astype(np.int64) << 1) | ((((u & 7) == 3) | ((u & 7) == 5)).astype(np.int64) << 2)
        c0 = (w & 1) | (((u & 3) == 3).astype(np.int64) << 1)
        mask = np.where(prec >= 3, np.int64(1) << c,
                        np.where(prec == 2, (np.int64(1) << c0) | (np.int64(1) << (c0 | 4)),
                                 np.where((w & 1) == 0, 0x55, 0xAA)))
        mask = np.where(zero, FULL, mask)
        if m & 1:
            mask = ((mask & 0x55) << 1) | ((mask & 0xAA) >> 1)
    else:
        c = (w & 1) | (((u % 3) == 2).astype(np.int64) << 1)
        mask = np.where(zero, FULL, np.int64(1) << c)
        if m & 1:
            mask = ((mask & 0x5) << 1) | ((mask & 0xA) >> 1)
    return mask


def vec_classify(Barr, Carr, k, m):
    """Barr, Carr -- массивы одинаковой длины. Возвращает (consistent, alldet, masks)."""
    q = p ** k
    pm = pow(p, m, q) if m < 64 else 0
    masks = []
    for (i, j) in CELLS:
        masks.append(vec_masks((pm + i * Barr + j * Carr) % q, k, m))
    cons = np.zeros(Barr.shape, dtype=bool)
    for s in SOLS:
        ok = np.ones(Barr.shape, dtype=bool)
        for n in NONC:
            ok &= ((masks[n] >> int(s[n])) & 1).astype(bool)
        cons |= ok
    alldet = np.ones(Barr.shape, dtype=bool)
    for n in NONC:
        mk = masks[n]
        alldet &= (mk & (mk - 1)) == 0
    return cons, alldet, masks


def scalar_classify(k, m, B0, C0):
    q = p ** k
    pm = pow(p, m, q)
    mk = [mask_of_residue((pm + i * B0 + j * C0) % q, k, p, m) for (i, j) in CELLS]
    cons = any(all(mk[n] >> s[n] & 1 for n in NONC) for s in sols)
    alldet = all((mk[n] & (mk[n] - 1)) == 0 for n in NONC)
    return cons, alldet, mk


def in_target(k, m, B0, C0):
    """узел целиком в целевой области min v(b),v(c) >= t"""
    if m > 0:
        return False
    if k < t:
        return False
    tt = p ** t
    return B0 % tt == 0 and C0 % tt == 0


def meets_target(k, m, B0, C0):
    if m > 0:
        return False
    kk = min(k, t)
    return B0 % p ** kk == 0 and C0 % p ** kk == 0


report = {'p': p, 'M': M, 'KFLAT': KFLAT, 'KMAX': KMAX, 'flat': {}, 'tree': {}}
CE = []
BUG = []
for m in range(M + 1):
    undet_final = []
    for k in range(1, KFLAT + 1):
        q = p ** k
        cnt = {'target': 0, 'fail': 0, 'CE': 0, 'undet': 0, 'target_touching_fail_BUG': 0}
        Carr_all = np.arange(q, dtype=np.int64)
        for B0 in range(q):                               # явная граница
            Carr = Carr_all
            if m > 0 and B0 % p == 0:
                Carr = Carr_all[Carr_all % p != 0]
            if Carr.size == 0:
                continue
            Barr = np.full(Carr.shape, B0, dtype=np.int64)
            cons, alldet, _ = vec_classify(Barr, Carr, k, m)
            if m == 0 and k >= t and B0 % p ** t == 0:
                tgt = (Carr % p ** t) == 0
            else:
                tgt = np.zeros(Carr.shape, dtype=bool)
            # узлы, пересекающие целевую область, не могут быть ПРОВАЛ (там все клетки -- квадраты)
            if m == 0:
                kk = min(k, t)
                touch = (B0 % p ** kk == 0) & ((Carr % p ** kk) == 0)
                bad = touch & ~cons
                if bad.any():
                    cnt['target_touching_fail_BUG'] += int(bad.sum())
                    BUG.append((m, k, B0, int(Carr[bad][0])))
            cnt['target'] += int(tgt.sum())
            rest = ~tgt
            cnt['fail'] += int((rest & ~cons).sum())
            ce = rest & cons & alldet
            if ce.any():
                cnt['CE'] += int(ce.sum())
                for C0 in Carr[ce][:5]:
                    CE.append((m, k, B0, int(C0)))
            un = rest & cons & ~alldet
            cnt['undet'] += int(un.sum())
            if k == KFLAT and un.any():
                undet_final.extend((B0, int(C0)) for C0 in Carr[un])
        report['flat'][f'm={m},k={k}'] = cnt
        print(f'p={p} m={m} k={k}: {cnt}  [{time.time()-t0:.1f}s]', flush=True)
    # --- адаптивное дерево от неопределённых узлов уровня KFLAT ---
    stack = [(KFLAT, B0, C0) for (B0, C0) in undet_final]
    nodes = 0
    maxdepth = KFLAT
    tree_cnt = {'target': 0, 'fail': 0, 'CE': 0, 'unresolved_at_KMAX': 0}
    unresolved = []
    while stack:
        if nodes >= NODE_BUDGET:
            print('  БЮДЖЕТ УЗЛОВ ИСЧЕРПАН', flush=True)
            tree_cnt['budget_exhausted'] = len(stack)
            break
        k, B0, C0 = stack.pop()
        q = p ** k
        for a in range(p):
            for a2 in range(p):
                k1 = k + 1
                B1, C1 = B0 + a * q, C0 + a2 * q
                nodes += 1
                maxdepth = max(maxdepth, k1)
                if in_target(k1, m, B1, C1):
                    tree_cnt['target'] += 1
                    continue
                cons, alldet, _ = scalar_classify(k1, m, B1, C1)
                if not cons:
                    if meets_target(k1, m, B1, C1):
                        BUG.append(('tree', m, k1, B1, C1))
                    tree_cnt['fail'] += 1
                elif alldet:
                    tree_cnt['CE'] += 1
                    CE.append(('tree', m, k1, B1, C1))
                elif k1 >= KMAX:
                    tree_cnt['unresolved_at_KMAX'] += 1
                    if len(unresolved) < 20:
                        unresolved.append((k1, B1, C1))
                else:
                    stack.append((k1, B1, C1))
    tree_cnt['nodes'] = nodes
    tree_cnt['maxdepth'] = maxdepth
    tree_cnt['start_undet'] = len(undet_final)
    tree_cnt['unresolved_examples'] = unresolved
    report['tree'][f'm={m}'] = tree_cnt
    print(f'  дерево p={p} m={m}: {tree_cnt}  [{time.time()-t0:.1f}s]', flush=True)

report['CE'] = CE[:50]
report['BUG'] = BUG[:50]
report['n_CE'] = len(CE)
report['n_BUG'] = len(BUG)
print('КОНТРПРИМЕРЫ:', CE[:10], ' ОШИБКИ МАСОК:', BUG[:10], flush=True)
with open(f'enum_p{p}_M{M}_k{KFLAT}.json', 'w') as f:
    json.dump(report, f, indent=1, ensure_ascii=False)
print(f'готово за {time.time()-t0:.1f}s', flush=True)
