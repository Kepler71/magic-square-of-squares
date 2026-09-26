# -*- coding: utf-8 -*-
"""t4: независимая сверка 34-точечного сертификата над F_10007 (claude_review/independence_certificate.json).
Код Codex/Claude НЕ импортируется; из JSON берутся только данные (a,b,c, корни, маски).
Определение столбцов восстановлено по описанию column_order:
  для каждой из 8 линий (строки j=-1,0,1 при меняющемся i; столбцы i=-1,0,1; диагональ; антидиагональ)
  с тройкой корней (alpha,gamma,beta) = (r[-1], r[0], r[+1]) вдоль линии: A=alpha+gamma, B=gamma+beta, C=alpha+beta,
  столбцы A*B (=x половинки) и A*C (=x-t); затем b,c,b+c,b-c,b+2c,b-2c,2b+c,2b-c; 9 корней в лекс. порядке (i,j); константа.
Символ Лежандра -- своим алгоритмом Якоби (взаимность), сверка с Sage kronecker. Ранг -- своя элиминация по столбцам
и Sage matrix(GF(2)).rank(); проверяем также det != 0 и тавтологию BC = AB*AC (mod квадраты) в каждой точке.
Дополнительно: свежая выборка своим кодом при p = 10009 (p = 1 mod 4, корни через Sage) и p = 10007 -- ранг 34.
Все циклы ограничены.
"""
import os, json, random, time
os.environ.setdefault("DOT_SAGE", "/tmp/sage_a51edd4e")
t0 = time.time()
CERT = "/home/kep/magicKube/lift_classes/codex/proof_progress_20260926/claude_review/independence_certificate.json"
data = json.load(open(CERT))
IDX = (-1, 0, 1)
LEX = [(i, j) for i in IDX for j in IDX]

def jacobi(a, n):
    assert n > 0 and n % 2 == 1
    a %= n; res = 1
    for _ in range(10 ** 6):          # ограниченный цикл
        if a == 0:
            return 0 if n != 1 else res
        while a % 2 == 0:
            a //= 2
            if n % 8 in (3, 5): res = -res
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3: res = -res
        a %= n
        if n == 1: return res
    raise RuntimeError("jacobi: не сошлось")

def is_prime_td(n):
    if n < 2: return False
    d = 2
    while d * d <= n:
        if n % d == 0: return False
        d += 1
    return True

def lines(R):
    L = []
    for j in IDX: L.append((R[(-1, j)], R[(0, j)], R[(1, j)]))
    for i in IDX: L.append((R[(i, -1)], R[(i, 0)], R[(i, 1)]))
    L.append((R[(-1, -1)], R[(0, 0)], R[(1, 1)]))
    L.append((R[(-1, 1)], R[(0, 0)], R[(1, -1)]))
    return L

def values(a, b, c, R):
    v = []
    for (al, ga, be) in lines(R):
        A, B, C = al + ga, ga + be, al + be
        v += [A * B, A * C]
    v += [b, c, b + c, b - c, b + 2 * c, b - 2 * c, 2 * b + c, 2 * b - c]
    v += [R[k] for k in LEX]
    return v

def bits_of(vals, p):
    out = []
    for x in vals:
        L = jacobi(x % p, p)
        if L == 0: return None
        out.append(0 if L == 1 else 1)
    return out + [1]

def rank_cols(rows, ncols):
    """своя элиминация: ранг по столбцам (транспонируем)"""
    cols = []
    for k in range(ncols):
        m = 0
        for i, row in enumerate(rows):
            m |= row[k] << i
        cols.append(m)
    piv = {}
    rk = 0
    for m in cols:
        x = m
        for _ in range(len(rows) + 1):
            if x == 0: break
            h = x.bit_length() - 1
            if h in piv: x ^= piv[h]
            else:
                piv[h] = x; rk += 1; break
    return rk

from sage.all import GF, matrix, kronecker, is_prime, sqrt as ssqrt
p = data["prime"]
res = dict(prime=p, prime_td=is_prime_td(p), prime_sage=bool(is_prime(p)), n_witness=len(data["witnesses"]))
rows = []; bad_roots = bad_zero = bad_mask = bad_kron = bad_taut = 0
pts = set()
for w in data["witnesses"]:
    a, b, c = w["a"], w["b"], w["c"]
    R = dict(zip(LEX, w["roots_lex"]))
    pts.add((a, b, c, tuple(w["roots_lex"])))
    for (i, j) in LEX:
        if (R[(i, j)] ** 2 - (a + i * b + j * c)) % p: bad_roots += 1
        if (a + i * b + j * c) % p == 0: bad_zero += 1
    vals = values(a, b, c, R)
    assert len(vals) == 33
    bt = bits_of(vals, p)
    if bt is None:
        bad_zero += 1; continue
    for x, bb in zip(vals, bt):
        if int(kronecker(x % p, p)) != (1 if bb == 0 else -1): bad_kron += 1
    mask = sum(bb << k for k, bb in enumerate(bt))
    if mask != w["mask"]: bad_mask += 1
    # тавтология: [B*C] = [A*B] + [A*C] для каждой линии
    for (al, ga, be) in lines(R):
        A, B, C = al + ga, ga + be, al + be
        if jacobi(B * C % p, p) != jacobi(A * B % p, p) * jacobi(A * C % p, p): bad_taut += 1
    rows.append(bt)
res.update(distinct_points=len(pts), bad_roots=bad_roots, bad_zero=bad_zero, bad_mask=bad_mask,
           bad_kronecker=bad_kron, bad_tautology=bad_taut)
res["rank_own"] = rank_cols(rows, 34)
M = matrix(GF(2), rows)
res["rank_sage"] = int(M.rank()); res["det_sage"] = int(M.det())
res["rank_33_without_const"] = int(matrix(GF(2), [r[:33] for r in rows]).rank())
print("сертификат:", res, flush=True)

# свежая выборка своим кодом
def fresh(p, npts, seed):
    rng = random.Random(seed)
    F = GF(p)
    rows = []
    tries = 0
    while len(rows) < npts and tries < 200000:
        tries += 1
        a, b, c = (rng.randrange(1, p) for _ in range(3))
        cells = {(i, j): (a + i * b + j * c) % p for (i, j) in LEX}
        if any(jacobi(v, p) != 1 for v in cells.values()): continue
        R = {}
        for k, v in cells.items():
            s0 = int(F(v).sqrt())
            R[k] = s0 if rng.random() < 0.5 else p - s0
        bt = bits_of(values(a, b, c, R), p)
        if bt is None: continue
        rows.append(bt)
    return rows, tries
fr = {}
for (pp, n, sd) in [(10007, 200, 1), (10009, 200, 2), (20021, 200, 3)]:
    assert is_prime(pp)
    rws, tries = fresh(pp, n, sd)
    fr[pp] = dict(p_mod4=pp % 4, points=len(rws), tries=tries, rank=int(matrix(GF(2), rws).rank()),
                  rank_own=rank_cols(rws, 34))
    print("свежая выборка p=%d:" % pp, fr[pp], "t=%.1fs" % (time.time() - t0), flush=True)
res["fresh"] = fr
res["time_s"] = round(time.time() - t0, 1)
json.dump(res, open("/home/kep/magicKube/lift_classes/recompute/v3_a51edd4e/t4_certificate.json", "w"), indent=1)
print("время %.1f с" % (time.time() - t0))
