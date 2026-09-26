# (б) Независимая проверка:
#  (1) gcd(D(P),D(Q)) | gcd(D(dG),D(dH)),  d = det;
#  (3) gcd(D(P1),D(P2),D(P3)) | gcd(D(Delta G), D(Delta H)), Delta = НОД 2x2-миноров (групповая часть, без AP);
#  INDEX: p≡3 mod 8, p∤2N, все P_j в ядре редукции, p ∤ gcd(D(A),D(B)), A=2G,B=2H  =>  gcd(d,(p+1)/4)>1,
#         а точнее |red(ZA+ZB)| делит d; плюс #E_N(F_p)=p+1 при p≡3 mod 4.
# Все циклы с явными границами. Запуск: env DOT_SAGE=/tmp/claude_verify_bridge_sage python3 b_lattice.py
from sage.all import *
import json, itertools, time, sys
from pathlib import Path
OUT = Path(__file__).resolve().parent
src = Path('/home/kep/magicKube/joint_Eb_Ec/lattice2d/rank2_N_le_700.json')
curves = json.loads(src.read_text())
NCURVES = 6          # явная граница
BOX_PAIR = 3         # |m|,|n| <= 3 для леммы (1)
BOX_IDX = 5          # |m|,|n| <= 5 для INDEX
PMAX = 400           # простые p ≡ 3 mod 8 до PMAX
t0 = time.time()

def D(P):
    return ZZ(0) if P.is_zero() else ZZ(P[0].denominator()).sqrt(extend=False)

def vecs(B):
    return [(m, n) for m in range(-B, B+1) for n in range(-B, B+1) if m > 0 or (m == 0 and n > 0)]

summary = dict(curves=[], pair_checks=0, pair_fail=0, triple_checks=0, triple_fail=0,
               index_instances=0, index_nonvacuous=0, index_fail=0, kernel_mismatch=0, count_fail=0)
chosen = [c for c in curves if all(QQ(x).denominator() == 1 for g in c['gens'] for x in g)][:NCURVES]
for ci, c in enumerate(chosen):
    N = ZZ(c['N']); E = EllipticCurve([-N**2, 0])
    G0, H0 = [E(QQ(g[0]), QQ(g[1])) for g in c['gens']]
    # --- лемма (1), G,H = 2G0,2H0 (как у Codex) и G,H = G0,H0 (лемма не требует 2E) ---
    for tag, (G, H) in (('2E', (2*G0, 2*H0)), ('E', (G0, H0))):
        V = vecs(BOX_PAIR)
        pts = {v: v[0]*G + v[1]*H for v in V}
        Dv = {v: D(pts[v]) for v in V}
        cache = {}
        def capf(dd):
            dd = abs(dd)
            if dd not in cache:
                cache[dd] = gcd(D(dd*G), D(dd*H))
            return cache[dd]
        for v, w in itertools.combinations(V, 2):
            det = v[0]*w[1] - v[1]*w[0]
            if det == 0: continue
            summary['pair_checks'] += 1
            if capf(det) % gcd(Dv[v], Dv[w]) != 0:
                summary['pair_fail'] += 1
                print('PAIR FAIL', N, tag, v, w, flush=True)
        # --- (3) групповая часть: тройки ---
        for tri in itertools.combinations(V, 3):
            dets = [tri[i][0]*tri[j][1] - tri[i][1]*tri[j][0] for i, j in ((0,1),(0,2),(1,2))]
            Delta = gcd(dets)
            if Delta == 0: continue
            summary['triple_checks'] += 1
            g3 = gcd([Dv[t] for t in tri])
            if capf(Delta) % g3 != 0:
                summary['triple_fail'] += 1
                print('TRIPLE FAIL', N, tag, tri, flush=True)
    print(f'[{time.time()-t0:.0f}s] N={N}: pair/triple done', flush=True)
    # --- INDEX-лемма ---
    A, B = 2*G0, 2*H0
    V = vecs(BOX_IDX)
    pts = {v: v[0]*A + v[1]*B for v in V}
    Dv = {v: D(pts[v]) for v in V}
    DA, DB = D(A), D(B)
    for p in prime_range(3, PMAX):
        if p % 8 != 3 or (2*N) % p == 0: continue
        Ep = E.change_ring(GF(p))
        if Ep.cardinality() != p + 1:
            summary['count_fail'] += 1; print('COUNT FAIL', N, p, flush=True)
        a_, b_ = Ep(A), Ep(B)
        # ядро редукции K_p = {(m,n): m a_ + n b_ = O}; прямо сверяем с делимостью D
        inK = {v: (v[0]*a_ + v[1]*b_).is_zero() for v in V}
        byD = {v: Dv[v] % p == 0 for v in V}
        if inK != byD:
            summary['kernel_mismatch'] += 1; print('KERNEL MISMATCH', N, p, flush=True)
        # порядок образа Gamma
        h = len(set((m*a_ + n*b_) for m in range(0, a_.order()) for n in range(0, b_.order())))
        assert ((p + 1)//4) % h == 0
        Vp = [v for v in V if byD[v]]
        MAXT = 3000  # явная граница на число троек на (кривая, p)
        cnt = 0
        for tri in itertools.combinations(Vp, 3):
            if cnt >= MAXT: break
            dets = [tri[i][0]*tri[j][1] - tri[i][1]*tri[j][0] for i, j in ((0,1),(0,2),(1,2))]
            d = gcd(dets)
            if d == 0: continue
            cnt += 1
            summary['index_instances'] += 1
            ok = (d % h == 0)
            if gcd(DA, DB) % p != 0:
                summary['index_nonvacuous'] += 1
                ok = ok and gcd(d, (p+1)//4) > 1
            if not ok:
                summary['index_fail'] += 1; print('INDEX FAIL', N, p, tri, d, h, flush=True)
    summary['curves'].append(int(N))
    print(f'[{time.time()-t0:.0f}s] N={N}: index done; running totals {summary}', flush=True)

(OUT / 'b_lattice.json').write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
assert summary['pair_fail'] == summary['triple_fail'] == summary['index_fail'] == summary['kernel_mismatch'] == summary['count_fail'] == 0
print('ALL (b) CHECKS PASSED')
