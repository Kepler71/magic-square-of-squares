# -*- coding: utf-8 -*-
"""
ПОЛНЫЙ перебор шаров (b, c) mod p^k для p = 2, 3, k = 1..KMAX, порядок полюса m = 0..MMAX(k).
b = B/p^m, c = C/p^m, B, C пробегают ВСЕ вычеты mod p^(k+m); при m >= 1 берём только
примитивные (B, C) (не оба делятся на p), т.е. полюс ровно порядка m.

Классы кодируются битами (реализация независима от p23lib.py):
  p=2: (чётность v, [u = 3 mod 4], [u = +-3 mod 8]) in F_2^3   (характеры (Z/8)^*)
  p=3: (чётность v, [u = 2 mod 3])                    in F_2^2
Линия определена, если определены все три клетки; её класс = XOR кодов.
Категории шара:
  naive: EXCL (есть определённая нетривиальная линия) / PASS (всё определено, всё тривиально) / UND
  pairs: то же после подстановки класса противоположной клетки вместо неопределённого
         (законно: 4 центральные линии дают равенство классов противоположных клеток);
         конфликт определённых классов в паре = EXCL.
Выход: JSON-таблица + проверки утверждений. Время каждого (p,m,k) печатается.
"""
import sys, json, time
import numpy as np

CELLS = [(i, j) for i in (-1, 0, 1) for j in (-1, 0, 1)]
IDX = {ij: n for n, ij in enumerate(CELLS)}
LINES = ([[IDX[(i, j)] for j in (-1, 0, 1)] for i in (-1, 0, 1)] +
         [[IDX[(i, j)] for i in (-1, 0, 1)] for j in (-1, 0, 1)] +
         [[IDX[(t, t)] for t in (-1, 0, 1)], [IDX[(t, -t)] for t in (-1, 0, 1)]])
OPP = [IDX[(-i, -j)] for (i, j) in CELLS]
NEED = {2: 3, 3: 1}
THR = {2: 3, 3: 1}


def classes(X, m, p, K):
    """X: int64-массив значений p^m * x mod p^K. Возвращает (det, code)."""
    PK = p ** K
    X = np.mod(X, PK)
    nz = X != 0
    Y = X.copy()
    v = np.zeros_like(X)
    for _ in range(K):                      # явная граница
        div = nz & (Y % p == 0)
        if not div.any():
            break
        Y = np.where(div, Y // p, Y)
        v += div
    det = nz & ((K - v) >= NEED[p])
    par = np.mod(v - m, 2)
    if p == 2:
        u = np.mod(Y, 8)
        e1 = (np.mod(u, 4) == 3).astype(np.int64)
        e2 = ((u == 3) | (u == 5)).astype(np.int64)
        code = par * 4 + e1 * 2 + e2
    else:
        u = np.mod(Y, 3)
        code = par * 2 + (u == 2).astype(np.int64)
    code = np.where(det, code, 0)
    return det, code


def categorize(det, code):
    """det, code: списки из 9 массивов. Возвращает массив категорий 0=EXCL,1=PASS,2=UND."""
    n = det[0].shape[0]
    excl = np.zeros(n, bool)
    allline = np.ones(n, bool)
    for L in LINES:
        d = det[L[0]] & det[L[1]] & det[L[2]]
        c = code[L[0]] ^ code[L[1]] ^ code[L[2]]
        excl |= d & (c != 0)
        allline &= d
    alldet = np.ones(n, bool)
    for t in range(9):
        alldet &= det[t]
    cat = np.full(n, 2, np.int8)
    cat[allline & alldet & ~excl] = 1
    cat[excl] = 0
    return cat


def run(p, m, k):
    K = k + m
    PK = p ** K
    pm = p ** m
    Cs = np.arange(PK, dtype=np.int64)
    vC = np.zeros_like(Cs)
    Yc = Cs.copy()
    nzc = Cs != 0
    for _ in range(K):
        div = nzc & (Yc % p == 0)
        Yc = np.where(div, Yc // p, Yc); vC += div
    vC = np.where(nzc, vC, K)
    tally = {}   # (method, cat, below) -> count
    t0 = time.time()
    for B in range(PK):                                  # явная граница: p^K итераций
        vB = K if B == 0 else next(t for t in range(K + 1) if B % p ** (t + 1) != 0)
        minv = np.minimum(vC, vB)
        if m >= 1:
            sel = minv == 0                              # полюс ровно порядка m
            if not sel.any():
                continue
            C = Cs[sel]; mv = minv[sel] - m
        else:
            C = Cs; mv = minv
        below = mv < THR[p]
        det, code = [], []
        for (i, j) in CELLS:
            d, c = classes(pm + i * B + j * C, m, p, K)
            det.append(d); code.append(c)
        cat_n = categorize(det, code)
        # метод пар
        conflict = np.zeros(C.shape[0], bool)
        det2, code2 = [], []
        for t in range(9):
            o = OPP[t]
            conflict |= det[t] & det[o] & (code[t] != code[o])
            det2.append(det[t] | det[o])
            code2.append(np.where(det[t], code[t], code[o]))
        cat_p = categorize(det2, code2)
        cat_p[conflict] = 0
        for meth, cat in (('naive', cat_n), ('pairs', cat_p)):
            for cv in (0, 1, 2):
                for bl in (True, False):
                    cnt = int(np.count_nonzero((cat == cv) & (below == bl)))
                    if cnt:
                        key = (meth, ('EXCL', 'PASS', 'UND')[cv], 'below' if bl else 'atleast')
                        tally[key] = tally.get(key, 0) + cnt
    return tally, time.time() - t0


PLAN = {
    # p: список (k, m_max)
    2: [(k, 3) for k in range(1, 9)] + [(k, 8) for k in (1, 2, 3, 4)],
    3: [(k, 2) for k in range(1, 9)] + [(k, 6) for k in (1, 2, 3)],
}

if __name__ == '__main__':
    only = int(sys.argv[1]) if len(sys.argv) > 1 else None
    out = {}
    viol = []
    seen = set()
    for p in (2, 3):
        if only and p != only:
            continue
        for (k, mmax) in PLAN[p]:
            for m in range(0, mmax + 1):
                if (p, k, m) in seen:
                    continue
                seen.add((p, k, m))
                if p ** (2 * (k + m)) > 5 * 10 ** 7:          # явный предел размера
                    continue
                tally, dt = run(p, m, k)
                key = f'p={p},k={k},m={m}'
                out[key] = {f'{a}/{b}/{c}': n for (a, b, c), n in sorted(tally.items())}
                # утверждения
                npb = tally.get(('naive', 'PASS', 'below'), 0)
                if npb:
                    viol.append((key, 'naive PASS below threshold', npb))
                # при k >= THR все шары с minv >= THR должны быть PASS (все клетки = 1 mod 8 / mod 3)
                if k >= THR[p]:
                    for cat in ('EXCL', 'UND'):
                        n_ = tally.get(('naive', cat, 'atleast'), 0)
                        if n_:
                            viol.append((key, f'naive {cat} at/above threshold', n_))
                ppb = tally.get(('pairs', 'PASS', 'below'), 0)
                pub = tally.get(('pairs', 'UND', 'below'), 0)
                print(f'{key}: {out[key]}  [{dt:.1f}s]', flush=True)
                out[key]['_pairs_not_excluded_below'] = ppb + pub
    out['_violations'] = viol
    print('VIOLATIONS:', viol, flush=True)
    fn = f'enum_full_p{only}.json' if only else 'enum_full.json'
    with open(fn, 'w') as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
