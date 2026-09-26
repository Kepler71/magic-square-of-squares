# -*- coding: utf-8 -*-
"""Независимая перепроверка (проверяющий, 26.09.2026) утверждений 3, 4, 5, 7, 12 записки center4/NOTE.md.
Не использует circle_lib.py и sympy. Другая параметризация:
  a - u = alpha^2, a + u = beta^2  <=>  a = x^2 + y^2, u = 2xy, x = (beta+alpha)/2, y = (beta-alpha)/2, x > y > 0.
Тогда Omega_u = alpha*beta + i u = (x + i y)^2 = z^2, w_u = z / conj(z).
Нормирования: для p = 1 (mod 4) берём p-адический корень iota из -1 (гензелевский подъём), гомоморфизм
phi: Z[i] -> Z/p^K, i -> iota; v_pi(z) = v_p(x + iota*y), v_{pi-bar}(z) = v_p(x - iota*y).
Проверки: лемма о двойном минимуме (через gcd(alpha, beta)), ранг матрицы v_pi(w_u) (свой алгоритм),
неравенство D по классам ориентаций, шаги доказательства D (Z' = sum n tau T != 0 и phi_p(Z') = 0 mod p^{2e}),
отрицательный контроль. Запуск: python3 my_configs.py A BLOCK SEED"""
import sys, json, time, random
import numpy as np
from math import isqrt, gcd
from fractions import Fraction
from itertools import product

A = int(sys.argv[1]); BLOCK = int(sys.argv[2]) if len(sys.argv) > 2 else 10**8
SEED = int(sys.argv[3]) if len(sys.argv) > 3 else 12345

# --- простые до sqrt(A) для пробного деления
M = isqrt(A) + 2
sv = bytearray([1]) * (M + 1); sv[0] = sv[1] = 0
for i in range(2, isqrt(M) + 1):
    if sv[i]:
        sv[i * i::i] = bytearray(len(sv[i * i::i]))
PR = [i for i in range(2, M + 1) if sv[i]]

def factor(n):
    f = {}
    for p in PR:
        if p * p > n: break
        while n % p == 0:
            f[p] = f.get(p, 0) + 1; n //= p
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def vp(n, p):
    if n == 0: return 10**9
    k = 0
    while n % p == 0:
        n //= p; k += 1
    return k

def iota(p, K):
    """корень из -1 по модулю p^K (p = 1 mod 4), гензелевский подъём"""
    for g in range(2, p):
        t = pow(g, (p - 1) // 4, p)
        if t * t % p == p - 1: break
    m = p
    for _ in range(K - 1):
        m2 = m * p
        # t <- t - (t^2+1)/(2t) mod m2
        t = (t - (t * t + 1) * pow(2 * t, -1, m2)) % m2
        m = m2
    assert (t * t + 1) % (p ** K) == 0
    return t

def rank_int(rows):
    R = [[Fraction(v) for v in r] for r in rows]
    if not R or not R[0]: return 0
    nr, nc = len(R), len(R[0]); rk = 0; col = 0
    for col in range(nc):
        piv = None
        for r in range(rk, nr):
            if R[r][col] != 0: piv = r; break
        if piv is None: continue
        R[rk], R[piv] = R[piv], R[rk]
        for r in range(nr):
            if r != rk and R[r][col] != 0:
                f = R[r][col] / R[rk][col]
                R[r] = [R[r][c] - f * R[rk][c] for c in range(nc)]
        rk += 1
        if rk == nr: break
    return rk

# все связи sum n u = 0 с мультимножеством |n| = {1,1,1} или {1,1,2}, с точностью до общего знака
NVECS = []
for n in product((-2, -1, 1, 2), repeat=3):
    if sorted(abs(v) for v in n) in ([1, 1, 1], [1, 1, 2]) and n[0] > 0:
        NVECS.append(n)

def relations(tr):
    return [n for n in NVECS if sum(a * b for a, b in zip(n, tr)) == 0]

def analyse(a, tr, rels):
    """tr — три положительных шага (отсортированы). Возвращает словарь проверок."""
    out = dict(a=a, us=list(tr), rels=[list(n) for n in rels])
    xy = []
    for u in tr:
        al, be = isqrt(a - u), isqrt(a + u)
        assert al > 0 and al * al == a - u and be * be == a + u
        assert (al + be) % 2 == 0
        x, y = (be + al) // 2, (be - al) // 2
        assert x * x + y * y == a and 2 * x * y == u
        xy.append((x, y, al, be))
    f = factor(a)
    out['factor'] = {str(p): e for p, e in f.items()}
    # лемма: для нечётных p | a — min v_p(gcd(alpha,beta)) достигается >= 2 раз
    lemma_ok = True
    for p in f:
        if p == 2: continue
        ts = [vp(gcd(al, be), p) for (_, _, al, be) in xy]
        if ts.count(min(ts)) < 2: lemma_ok = False
    out['lemma_ok'] = lemma_ok
    split = sorted(p for p in f if p % 4 == 1)
    out['nsplit'] = len(split)
    cols, prof = [], {}
    for p in split:
        e = f[p]; K = 2 * e + 2; m = p ** K
        io = iota(p, K)
        vs = []
        for (x, y, _, _) in xy:
            r = min(vp((x + io * y) % m, p), K); rb = min(vp((x - io * y) % m, p), K)
            assert r + rb == e, (a, p, r, rb, e)
            vs.append((r, rb))
        prof[p] = (e, io, vs)
        cols.append([r - rb for (r, rb) in vs])
    Mrows = [[cols[j][i] for j in range(len(cols))] for i in range(3)]
    out['rank'] = rank_int(Mrows) if cols else 0
    # классы: все три линии p-чистые, sigma-вектор с точностью до знака
    classes = {}
    for p in split:
        e, io, vs = prof[p]
        if all(min(r, rb) == 0 for r, rb in vs):
            sig = tuple(1 if r > rb else -1 for r, rb in vs)
            if sig[0] < 0: sig = tuple(-s for s in sig)
            classes.setdefault(sig, []).append(p)
    out['classes'] = []
    D_ok = True; steps_ok = True; maxr = Fraction(0)
    for n in rels:
        NL = sum(abs(v) for v in n)
        for sig, ps in classes.items():
            Pi = 1
            for p in ps: Pi *= p ** f[p]
            r = Fraction(Pi * Pi, NL * a)
            if r > 1: D_ok = False
            maxr = max(maxr, r)
            # шаги доказательства: Z' = sum n_u tau_u T_u, T_u in {z^2, conj(z)^2}, phi_p(T_u) — единица
            Zr, Zi = 0, 0
            iotas = {}
            for p in ps:
                e, io, vs = prof[p]
                m = p ** (2 * e + 2)
                # ориентируем iota так, чтобы у первой линии v_p(x + iota y) = e (sigma=+1)
                x1, y1 = xy[0][0], xy[0][1]
                if min(vp((x1 + io * y1) % m, p), 2 * e + 2) != e: io = (-io) % m
                iotas[p] = (io, m, e)
            for (x, y, _, _), nu in zip(xy, n):
                taus = set()
                for p, (io, m, e) in iotas.items():
                    taus.add(1 if vp((x + io * y) % m, p) == 0 else -1)
                if len(taus) != 1: steps_ok = False; continue
                tau = taus.pop()
                yy = y if tau == 1 else -y
                Tr, Ti = x * x - yy * yy, 2 * x * yy
                Zr += nu * tau * Tr; Zi += nu * tau * Ti
            if Zr == 0 and Zi == 0: steps_ok = False
            for p, (io, m, e) in iotas.items():
                if (Zr + io * Zi) % (p ** (2 * e)) != 0: steps_ok = False
                if (Zr - io * Zi) % (p ** (2 * e)) != 0: steps_ok = False
            if Zr * Zr + Zi * Zi > (NL * a) ** 2: steps_ok = False
            out['classes'].append(dict(NL=NL, sig=list(sig), ps=ps, ratio=str(r)))
    out['D_ok'] = D_ok; out['steps_ok'] = steps_ok; out['maxr'] = maxr
    return out

def D_ok_norel(a, tr, NL=4):
    """для тройки без связи: нарушилось бы неравенство D?"""
    f = factor(a)
    classes = {}
    for p in sorted(q for q in f if q % 4 == 1):
        e = f[p]; K = 2 * e + 2; m = p ** K; io = iota(p, K)
        vs = []
        for u in tr:
            al, be = isqrt(a - u), isqrt(a + u); x, y = (be + al) // 2, (be - al) // 2
            vs.append((min(vp((x + io * y) % m, p), K), min(vp((x - io * y) % m, p), K)))
        if all(min(r, rb) == 0 for r, rb in vs):
            sig = tuple(1 if r > rb else -1 for r, rb in vs)
            if sig[0] < 0: sig = tuple(-s for s in sig)
            classes.setdefault(sig, []).append(p)
    for sig, ps in classes.items():
        Pi = 1
        for p in ps: Pi *= p ** f[p]
        if Pi * Pi > NL * a: return False
    return True

t0 = time.time()
res = dict(A=A, configs=0, config_rel_pairs=0, by_rel={}, square_centers=[], lemma_fail=[], rank_hist={},
           nsplit_hist={}, D_fail=[], steps_fail=[], max_ratio='0', max_ratio_example=None, multi_rel=[])
maxr = Fraction(0)
rng = random.Random(SEED)
neg_pool = []
pool_seen = 0
lo = 1
while lo <= A:
    hi = min(A, lo + BLOCK - 1)
    As, Us = [], []
    for x in range(2, isqrt(hi) + 1):
        x2 = x * x
        ymin2 = lo - x2
        ymin = 1 if ymin2 <= 1 else isqrt(ymin2 - 1) + 1
        if hi - x2 < 1: continue
        ymax = min(x - 1, isqrt(hi - x2))
        if ymin > ymax: continue
        y = np.arange(ymin, ymax + 1, dtype=np.int64)
        As.append(x2 + y * y); Us.append(2 * x * y)
    a_arr = np.concatenate(As); u_arr = np.concatenate(Us); del As, Us
    o = np.lexsort((u_arr, a_arr)); a_arr = a_arr[o]; u_arr = u_arr[o]; del o
    brk = np.flatnonzero(np.diff(a_arr)) + 1
    st = np.concatenate(([0], brk)); en = np.concatenate((brk, [len(a_arr)]))
    big = np.flatnonzero(en - st >= 3)
    for gi in big:
        s0, e0 = st[gi], en[gi]
        a = int(a_arr[s0]); L = [int(v) for v in u_arr[s0:e0]]; S = set(L)
        assert len(S) == len(L)
        found = set()
        for i in range(len(L)):
            p_ = L[i]
            for j in range(i + 1, len(L)):
                q_ = L[j]
                for c in (p_ + q_, q_ - p_, 2 * q_ - p_, 2 * p_ - q_, p_ + 2 * q_, q_ + 2 * p_, q_ - 2 * p_):
                    if c > 0 and c in S and c != p_ and c != q_:
                        found.add(tuple(sorted((p_, q_, c))))
                for c2 in (p_ + q_, q_ - p_):
                    if c2 % 2 == 0 and (c2 // 2) in S and c2 // 2 not in (p_, q_):
                        found.add(tuple(sorted((p_, q_, c2 // 2))))
        for tr in found:
            rels = relations(tr)
            if not rels: continue
            res['configs'] += 1; res['config_rel_pairs'] += len(rels)
            if len(rels) > 1: res['multi_rel'].append((a, tr, rels))
            for n in rels:
                k = str(tuple(sorted(abs(v) for v in n)))
                res['by_rel'][k] = res['by_rel'].get(k, 0) + 1
            if isqrt(a) ** 2 == a: res['square_centers'].append((a, tr))
            o_ = analyse(a, tr, rels)
            if not o_['lemma_ok']: res['lemma_fail'].append((a, tr))
            res['rank_hist'][o_['rank']] = res['rank_hist'].get(o_['rank'], 0) + 1
            res['nsplit_hist'][o_['nsplit']] = res['nsplit_hist'].get(o_['nsplit'], 0) + 1
            if not o_['D_ok']: res['D_fail'].append((a, tr))
            if not o_['steps_ok']: res['steps_fail'].append((a, tr))
            if o_['maxr'] > maxr:
                maxr = o_['maxr']; res['max_ratio'] = str(maxr)
                res['max_ratio_example'] = dict(a=a, us=list(tr), factor=o_['factor'], classes=o_['classes'])
        # пул для отрицательного контроля (резервуарная выборка групп)
        pool_seen += 1
        if len(neg_pool) < 20000:
            neg_pool.append((a, L))
        else:
            k = rng.randrange(0, pool_seen)
            if k < 20000: neg_pool[k] = (a, L)
    print(f'блок [{lo},{hi}] пар {len(a_arr)} конф {res["configs"]} ранги {res["rank_hist"]} {time.time()-t0:.0f}s', flush=True)
    del a_arr, u_arr
    lo = hi + 1

# отрицательный контроль: случайные тройки без связи
neg = dict(triples=0, would_fail=0)
for a, L in neg_pool[:20000]:
    tr = tuple(sorted(rng.sample(L, 3)))
    if relations(tr): continue
    neg['triples'] += 1
    if not D_ok_norel(a, tr, 4): neg['would_fail'] += 1
res['negative_control'] = neg
res['time_s'] = round(time.time() - t0, 1)
res['multi_rel'] = res['multi_rel'][:20]
json.dump(res, open(f'/home/kep/magicKube/joint_Eb_Ec/center4/check/my_configs_A{A}.json', 'w'), indent=1, default=str)
print(json.dumps({k: v for k, v in res.items() if k not in ('lemma_fail', 'D_fail', 'steps_fail')}, default=str)[:3000])
