# -*- coding: utf-8 -*-
"""
Контроль классификатора:
 (1) вторая реализация (p23lib, чистый Python) против numpy-перебора: те же счётчики при малых k;
     плюс монотонность: потомки EXCL/PASS-шара на уровне k+1 имеют ту же категорию.
 (2) Sage: точная p-адическая квадратность (QQ(x).is_padic_square(p)) восьми произведений
     на случайных и прицельных рациональных (b, c); сравнение с категорией шара, и главное:
     ни одной рациональной точки с min v_p(b,c) < порога и всеми 8 квадратами в Q_p.
Запуск: env DOT_SAGE=/tmp/dotsage_verify_p23 python3 validate.py
"""
import json, random, time, sys
from fractions import Fraction as Fr
import p23lib as L
import enum_full as E
import numpy as np

random.seed(20260926)
out = {}

# ---------- (1) две реализации и монотонность ----------
def cats_python(p, m, k):
    K = k + m
    PK = p ** K
    tally = {}
    catmap = {}
    for B in range(PK):
        for C in range(PK):
            if m >= 1 and B % p == 0 and C % p == 0:
                continue
            cls = L.ball_cells(B, C, m, p, K)
            cn = L.classify_naive(cls, p)
            cp = L.classify_pairs(cls, p)
            mv = L.minval_ball(B, C, m, p, K)
            bl = 'below' if mv < L.THR[p] else 'atleast'
            for meth, c in (('naive', cn), ('pairs', cp)):
                tally[(meth, c, bl)] = tally.get((meth, c, bl), 0) + 1
            catmap[(B, C)] = (cn, cp)
    return tally, catmap

mismatch = []
mono_bad = 0
mono_checked = 0
for p, plan in ((2, [(k, m) for k in range(1, 6) for m in range(0, 3)]),
                (3, [(k, m) for k in range(1, 4) for m in range(0, 2)])):
    for (k, m) in plan:
        tp, cmap = cats_python(p, m, k)
        tn, _ = E.run(p, m, k)
        if tp != tn:
            mismatch.append((p, k, m, str(tp), str(tn)))
        # монотонность: сравнить с уровнем k+1 (тот же m)
        if p ** (2 * (k + 1 + m)) <= 3 * 10 ** 5:
            K1 = k + 1 + m
            PK = p ** (k + m)
            for B1 in range(p ** K1):
                for C1 in range(p ** K1):
                    if m >= 1 and B1 % p == 0 and C1 % p == 0:
                        continue
                    parent = cmap[(B1 % PK, C1 % PK)]
                    cls = L.ball_cells(B1, C1, m, p, K1)
                    child = (L.classify_naive(cls, p), L.classify_pairs(cls, p))
                    mono_checked += 1
                    for a, b in zip(parent, child):
                        if a in ('EXCL', 'PASS') and a != b:
                            mono_bad += 1
        print(f'impl-compare p={p} k={k} m={m} ok={tp == tn}', flush=True)
out['impl_mismatch'] = mismatch
out['monotonicity'] = {'checked_children': mono_checked, 'bad': mono_bad}

# ---------- (2) Sage ----------
from sage.all import QQ, Integer

def vfr(x, p):
    x = QQ(x)
    return x.valuation(p) if x != 0 else 10 ** 9

def eight_products(b, c):
    cell = {(i, j): 1 + i * b + j * c for i in (-1, 0, 1) for j in (-1, 0, 1)}
    prods = []
    for line in L.LINES:
        pr = QQ(1)
        for n in line:
            pr *= cell[L.CELLS[n]]
        prods.append(pr)
    return cell, prods

def ball_of(b, c, p, k):
    """(B, C, m) шара уровня k, содержащего рациональные (b, c)."""
    m = max(0, -min(vfr(b, p), vfr(c, p)))
    K = k + m
    PK = p ** K
    def red(x):
        y = QQ(x) * p ** m
        num, den = int(y.numerator()), int(y.denominator())
        return (num * pow(den, -1, PK)) % PK
    return red(b), red(c), m, K

def rand_q(p, lo, hi):
    """случайное рациональное с v_p в [lo, hi] и 'случайной' единичной частью."""
    v = random.randint(lo, hi)
    num = random.randint(1, 10 ** 6) * random.choice((-1, 1))
    while num % p == 0:
        num //= p
    den = random.randint(1, 10 ** 4)
    while den % p == 0:
        den //= p
    return QQ(num) / den * QQ(p) ** v

stats = {}
bad_rational = []
consistency_bad = []
NS = 60000
t0 = time.time()
for p in (2, 3):
    st = {'tested': 0, 'all8_square': 0, 'all8_square_below': 0,
          'targeted': 0, 'ball_PASS_but_not_all8': 0, 'ball_EXCL_but_all8': 0}
    for it in range(NS):                       # явная граница
        mode = it % 4
        if mode == 0:        # общие
            b = rand_q(p, -3, 5); c = rand_q(p, -3, 5)
        elif mode == 1:      # около особых прямых: b = +-1 + p^t w, c любое малое/целое
            t = random.randint(1, 12)
            b = random.choice((1, -1)) + rand_q(p, t, t)
            c = rand_q(p, 0, 6)
        elif mode == 2:      # около b + c = +-1 и b - c = +-1
            t = random.randint(1, 12)
            c = rand_q(p, 0, 4)
            b = random.choice((1, -1)) - random.choice((1, -1)) * c + rand_q(p, t, t)
        else:                # точки с малым min v, но 'почти' проходящие: b, c = единицы*p^e, e = thr-1
            e = L.THR[p] - 1
            b = rand_q(p, e, e + 2); c = rand_q(p, e, e + 2)
            if random.random() < 0.5:
                c = b * QQ(random.choice((1, -1))) + rand_q(p, e + 1, e + 6)
        if b == 0 or c == 0 or b == c or b == -c:
            continue
        cell, prods = eight_products(b, c)
        if any(v == 0 for v in cell.values()):
            continue
        st['tested'] += 1
        if mode:
            st['targeted'] += 1
        ok8 = all(QQ(x).is_padic_square(p) for x in prods)
        mv = min(vfr(b, p), vfr(c, p))
        if ok8:
            st['all8_square'] += 1
            if mv < L.THR[p]:
                st['all8_square_below'] += 1
                bad_rational.append((p, str(b), str(c)))
            # подъём: все 9 клеток квадраты в Q_p
            if not all(QQ(x).is_padic_square(p) for x in cell.values()):
                bad_rational.append((p, 'no-lift', str(b), str(c)))
        # согласие с категорией шара уровня k = 6 (p=2) / 4 (p=3)
        kk = 6 if p == 2 else 4
        B, C, m, K = ball_of(b, c, p, kk)
        cls = L.ball_cells(B, C, m, p, K)
        cat = L.classify_naive(cls, p)
        if cat == 'PASS' and not ok8:
            st['ball_PASS_but_not_all8'] += 1; consistency_bad.append((p, str(b), str(c), cat))
        if cat == 'EXCL' and ok8:
            st['ball_EXCL_but_all8'] += 1; consistency_bad.append((p, str(b), str(c), cat))
        # и класс каждой клетки из шара = точный класс (где определён)
        for n, (i, j) in enumerate(L.CELLS):
            cl = cls[n]
            if cl is None:
                continue
            x = QQ(cell[(i, j)])
            v = x.valuation(p)
            u = x / QQ(p) ** v
            um = (int(u.numerator()) * pow(int(u.denominator()), -1, L.UMOD[p])) % L.UMOD[p]
            if cl != (v % 2, um):
                consistency_bad.append((p, 'cellclass', str(b), str(c), n, cl, (v % 2, um)))
        if it % 20000 == 0:
            print(f'sage p={p} it={it} {st} [{time.time() - t0:.0f}s]', flush=True)
    stats[p] = st
out['sage_random'] = stats
out['bad_rational'] = bad_rational[:50]
out['consistency_bad'] = consistency_bad[:50]
print(json.dumps(out, indent=1, default=str), flush=True)
with open('validate.json', 'w') as f:
    json.dump(out, f, indent=1, default=str)
