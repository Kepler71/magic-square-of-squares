"""Сертификат по классам вычетов для Шагов 3-4 (p=2,3).
Точка (b,c) в Q_p^2 <-> примитивная тройка (t,beta,gamma) в Z_p^3, b=beta/t, c=gamma/t, t!=0.
Карты: A: t=1 (b,c в Z_p);  B: v(t)>=1, beta=1;  C: v(t)>=1, v(beta)>=1, gamma=1.  Покрывают все Q_p^2.
Класс клетки f_ij = [t+i*beta+j*gamma] - [t]  в Q_p*/Q_p*^2.
Для каждого класса вычетов (t,beta,gamma) mod p^K строим МНОЖЕСТВО возможных классов каждой клетки
(все ненулевые p-адические числа с этим вычетом) и проверяем, существует ли пара (S,U) в G^2,
совместимая с матрицей шага 1 (а значит — с квадратностью 8 произведений).
Утверждение Шагов 3-4: совместимость возможна лишь при min(v(b),v(c)) >= m0 (m0=3 при p=2, 1 при p=3),
и тогда только (S,U)=(0,0).  Ноль «выживших» при K=K0 влечёт ноль при всех K>=K0 (монотонность)."""
import sys, time, json
import numpy as np

p = int(sys.argv[1]); KMAX = int(sys.argv[2])
LINESET = sys.argv[3] if len(sys.argv) > 3 else "all8"
m0 = 3 if p == 2 else 1
d = 3 if p == 2 else 2          # dim_F2 Q_p*/Q_p*^2
NG = 1 << d
FULL = (1 << NG) - 1

def cls_code(vbit, u):
    """u -- единица по модулю 8 (p=2) или 3 (p=3); код класса."""
    if p == 2:
        a = (u % 4 == 3).astype(np.int64); bb = ((u % 8 == 3) | (u % 8 == 5)).astype(np.int64)
        return vbit | (a << 1) | (bb << 2)
    return vbit | ((u % 3 == 2).astype(np.int64) << 1)

def classmask(x, K):
    """x: массив вычетов mod p^K (0..p^K-1). Возвращает битовую маску возможных классов."""
    P = p**K
    x = np.mod(x, P).astype(np.int64)
    out = np.full(x.shape, FULL, dtype=np.int64)
    nz = x != 0
    v = np.zeros(x.shape, dtype=np.int64); u = x.copy()
    for _ in range(K):                          # явная граница цикла
        m = nz & (u % p == 0)
        if not m.any(): break
        u = np.where(m, u // p, u); v = np.where(m, v + 1, v)
    prec = K - v                                # число известных p-адических цифр единицы
    vbit = v & 1
    if p == 2:
        # prec>=3: единица известна mod 8 -> один класс
        m3 = nz & (prec >= 3)
        out = np.where(m3, 1 << cls_code(vbit, u % 8), out)
        # prec==2: известен u mod 4 -> 2 класса (u или u+4 mod 8)
        m2 = nz & (prec == 2)
        u4 = u % 4
        out = np.where(m2, (1 << cls_code(vbit, u4)) | (1 << cls_code(vbit, u4 + 4)), out)
        # prec==1: u нечётно, всё прочее неизвестно -> 4 класса с данной чётностью оценки
        m1 = nz & (prec == 1)
        allu = np.zeros(x.shape, dtype=np.int64)
        for uu in (1, 3, 5, 7): allu |= (1 << cls_code(vbit, np.full(x.shape, uu)))
        out = np.where(m1, allu, out)
    else:
        out = np.where(nz, 1 << cls_code(vbit, u % 3), out)
    return out, np.where(nz, v, K)             # v (или K, если вычет 0: тогда v>=K)

# перенос маски на класс tau:  {c ^ tau}
SHIFT = np.zeros((NG, 1 << NG), dtype=np.int64)
for tau in range(NG):
    for mask in range(1 << NG):
        SHIFT[tau, mask] = sum(1 << (c ^ tau) for c in range(NG) if mask >> c & 1)

def idx(i, j): return (i+1)*3 + (j+1)
NONC = [k for k in range(9) if k != 4]
# решение системы линий: для all8 -- M(S,U); для других наборов -- явный базис над F2
LINES = {}
for i in (-1, 0, 1): LINES[f"row{i:+d}"] = [idx(i, j) for j in (-1, 0, 1)]
for j in (-1, 0, 1): LINES[f"col{j:+d}"] = [idx(i, j) for i in (-1, 0, 1)]
LINES["diag"] = [idx(-1, -1), 4, idx(1, 1)]; LINES["anti"] = [idx(-1, 1), 4, idx(1, -1)]
if LINESET == "all8": names = list(LINES)
elif LINESET == "central4": names = ["row+0", "col+0", "diag", "anti"]
else: names = [n for n in LINES if n != LINESET.replace("drop_", "")]
from itertools import product
sols = []
for bits in product((0, 1), repeat=8):
    x = [0]*9
    for k, bt in zip(NONC, bits): x[k] = bt
    if all(sum(x[k] for k in LINES[n]) % 2 == 0 for n in names): sols.append(x)
# базис пространства решений (гауссом по F2)
basis = []
def reduce(v, B):
    v = v[:]
    for bvec, piv in B:
        if v[piv]: v = [a ^ b for a, b in zip(v, bvec)]
    return v
for s in sols:
    r = reduce(s, basis)
    if any(r): basis.append((r, r.index(1)))
basis = [b for b, _ in basis]
delta = len(basis)
combos = list(product(range(NG), repeat=delta))   # |G|^delta, явная граница
print(f"p={p} lines={LINESET} ({len(names)}) solution dim over F2 = {delta}, combos={len(combos)}", flush=True)
# для all8 проверим, что базис даёт M(S,U): координаты S=x[idx(1,1)], U=x[idx(1,0)]

def run_block(K, chart, t, be, ga):
    n = t.size
    tmask, vt = classmask(t, K)
    gm = {}
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if (i, j) != (0, 0): gm[idx(i, j)] = classmask(t + i*be + j*ga, K)[0]
    _, vb = classmask(be, K); _, vg = classmask(ga, K)
    if chart == "A":
        viol = np.minimum(vb, vg) < m0          # при K >= m0 известно точно
    else:
        viol = np.ones(n, dtype=bool)            # полюс: min(v(b),v(c)) = -v(t) < 0
    anycons = np.zeros(n, dtype=bool); nontriv = np.zeros(n, dtype=bool); trivok = np.zeros(n, dtype=bool)
    for tau in range(NG):
        tin = (tmask >> tau) & 1 == 1
        if not tin.any(): continue
        cm = {k: SHIFT[tau][gm[k]] for k in gm}
        for g in combos:
            x = [0]*9
            for gk, bv in zip(g, basis):
                if gk: x = [xi ^ (gk if bi else 0) for xi, bi in zip(x, bv)]
            ok = tin.copy()
            for k in NONC:
                ok &= (cm[k] >> x[k]) & 1 == 1
            anycons |= ok
            if any(x): nontriv |= ok
            else: trivok |= ok
    surv = viol & anycons
    ex = [[int(t[e]), int(be[e]), int(ga[e])] for e in np.nonzero(surv)[0][:5]]
    return dict(points=int(n), viol_region=int(viol.sum()), survivors_in_viol=int(surv.sum()),
                good_region=int((~viol).sum()), good_nontriv_consistent=int((~viol & nontriv).sum()),
                good_trivial_consistent=int((~viol & trivok).sum())), ex

def run_chart(K, chart):
    P = p**K
    ar = np.arange(P, dtype=np.int64); arp = np.arange(0, P, p, dtype=np.int64); one = np.ones(1, dtype=np.int64)
    T, Bt, Gm = {"A": (one, ar, ar), "B": (arp, one, ar), "C": (arp, arp, one)}[chart]
    # внешняя ось -- первая неединичная, куски по ~1.5e6 точек
    inner = Gm.size if chart != "C" else Bt.size
    outer_len = (Bt.size if chart == "A" else T.size)
    step = max(1, 1500000 // max(1, inner))
    tot = None; exs = []
    for s0 in range(0, outer_len, step):            # явная граница
        if chart == "A":
            tt, bb, gg = np.meshgrid(T, Bt[s0:s0+step], Gm, indexing="ij")
        else:
            tt, bb, gg = np.meshgrid(T[s0:s0+step], Bt, Gm, indexing="ij")
        r, ex = run_block(K, chart, tt.ravel(), bb.ravel(), gg.ravel())
        tot = r if tot is None else {k: tot[k] + r[k] for k in tot}
        exs += ex
    if exs: tot["examples(t,beta,gamma)"] = exs[:5]
    return tot

out = {"p": p, "lines": LINESET, "dim": delta, "by_K": {}}
for K in range(1, KMAX + 1):
    t0 = time.time()
    out["by_K"][K] = {ch: run_chart(K, ch) for ch in ("A", "B", "C")}
    tot = sum(v["survivors_in_viol"] for v in out["by_K"][K].values())
    print(f"K={K} ({time.time()-t0:.1f}s) survivors_total={tot} ", json.dumps(out["by_K"][K]), flush=True)
json.dump(out, open(f"classset_cert_p{p}_{LINESET}.json", "w"), indent=1)
