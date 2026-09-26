# (д) Дополнительные попытки опровержения записки Codex «восемь линий».
#  1. Локальный критерий для AP квадратов по ВСЕМ оценкам (полюса и нули), все p, включая 2:
#     перебор примитивных AP квадратов (m^2-2mn-n^2)^2, (m^2+n^2)^2, (m^2+2mn-n^2)^2 и рациональных масштабов.
#  2. Подгрупповое свойство фильтрации {v_p(D) >= k} на модели y^2=x^3-N^2x при p=2 (N чётно) и p | N.
#  3. Лемма (1) с G,H, содержащими кручение, и на кривой ранга 3 (N=1254).
# Все циклы с явными границами. Запуск: env DOT_SAGE=/tmp/claude_verify_bridge_sage python3 e_extra.py
from sage.all import *
import json, itertools, time, random
from pathlib import Path
OUT = Path(__file__).resolve().parent
t0 = time.time(); res = {}
random.seed(20260926)

def allowed(X, mid, p):
    """True если тройка X с серединой mid проходит локальный критерий Codex в p (по всем оценкам)."""
    v = [x.valuation(p) for x in X]; mu = min(v); at = [i for i in range(3) if v[i] == mu]
    if len(at) == 3: return True
    if len(at) == 1: return False
    if p == 2: return False
    need = 2 if mid in at else -1
    return kronecker(need, p) == 1

# --- 1 ---
MMAX = 40; NSCALE = 6
cnt = 0; bad = []; eq_fail = []
pattern = {}
for m in range(2, MMAX+1):
    for n in range(1, m):
        if gcd(m, n) != 1: continue
        base = [QQ((m*m-2*m*n-n*n)**2), QQ((m*m+n*n)**2), QQ((m*m+2*m*n-n*n)**2)]
        if 0 in base: continue
        assert base[0] + base[2] == 2*base[1]
        for _ in range(NSCALE):
            u = random.randint(1, 10**4); w = random.randint(1, 10**4)
            X = [x*QQ(u)**2/QQ(w)**2 for x in base]
            order = random.sample(range(3), 3)       # перемешиваем, чтобы середина была не всегда 1
            Y = [X[i] for i in order]; mid = order.index(1)
            ps = sorted(set(sum([prime_divisors(y.numerator()) + prime_divisors(y.denominator()) for y in Y], [2])))
            for p in ps:
                cnt += 1
                if not allowed(Y, mid, p): bad.append((m, n, u, w, int(p)))
                v = [y.valuation(p) for y in Y]
                if (p == 2 or p % 8 == 3) and len(set(v)) != 1: eq_fail.append((m, n, u, w, int(p)))
                if len(set(v)) == 2 and p != 2:
                    mu = min(v); hi = [i for i in range(3) if v[i] != mu][0]
                    key = (int(p % 8), 'mid_high' if hi == mid else 'end_high')
                    pattern[key] = pattern.get(key, 0) + 1
res['ap_criterion'] = dict(prime_checks=cnt, false_rejections=len(bad), eq_fail_2_or_3mod8=len(eq_fail),
                           two_level_patterns={'%d:%s' % k: c for k, c in sorted(pattern.items())})
assert not bad and not eq_fail
# ожидание: только (1,*), (5,mid_high), (7,end_high)
assert all((k[0] == 1) or (k[0] == 5 and k[1] == 'mid_high') or (k[0] == 7 and k[1] == 'end_high') for k in pattern)
print(f'[{time.time()-t0:.0f}s] 1: {res["ap_criterion"]}', flush=True)

# --- 2 ---
def Dp(P, p):
    return 10**9 if P.is_zero() else max(0, -P[0].valuation(p)) // 2
sub_checks = 0; sub_fail = 0
for N, gens in [(34, [(-16, 120), (-2, 48)]), (138, None), (154, None), (65, None), (1254, None)]:
    E = EllipticCurve([-N**2, 0])
    if gens is None:
        G = E.gens(proof=False)
    else:
        G = [E(*g) for g in gens]
    G = list(G)[:3]
    T = [E(0, 0), E(N, 0), E(-N, 0)]
    BOX = 3 if len(G) == 2 else 2
    pts = []
    for coeffs in itertools.product(range(-BOX, BOX+1), repeat=len(G)):
        for t in [E(0)] + T:
            P = sum((c*g for c, g in zip(coeffs, G)), E(0)) + t
            if not P.is_zero(): pts.append(P)
    pts = pts[:400]   # явная граница
    primes = sorted(set([2] + prime_divisors(N) + sum([prime_divisors(P[0].denominator()) for P in pts[:60]], [])))[:12]
    for p in primes:
        for P, Q in itertools.islice(itertools.combinations(pts, 2), 6000):
            k = min(Dp(P, p), Dp(Q, p))
            if k == 0: continue
            sub_checks += 1
            for R in (P + Q, P - Q):
                if Dp(R, p) < k: sub_fail += 1; print('SUBGROUP FAIL', N, p, P, Q, flush=True)
    print(f'[{time.time()-t0:.0f}s] 2: N={N} rank-gens={len(G)} checks so far {sub_checks}', flush=True)
res['filtration_subgroup'] = dict(checks=sub_checks, fails=sub_fail)
assert sub_fail == 0

# --- 3: лемма (1) с кручением в G,H и на кривой ранга 3 ---
def D(P): return ZZ(0) if P.is_zero() else ZZ(P[0].denominator()).sqrt()
l1 = 0; l1f = 0
for N in (34, 1254):
    E = EllipticCurve([-N**2, 0]); g = list(E.gens(proof=False))
    tors = E(N, 0)
    pairs = [(g[0] + tors, g[1])] + ([(g[0] + g[2], g[1] - g[2]), (g[2], g[0] + E(0, 0))] if len(g) >= 3 else [])
    for G, H in pairs:
        V = [(m, n) for m in range(-3, 4) for n in range(-3, 4) if (m, n) > (0, 0)]
        Dv = {v: D(v[0]*G + v[1]*H) for v in V}
        for v, w in itertools.combinations(V, 2):
            det = v[0]*w[1] - v[1]*w[0]
            if det == 0: continue
            l1 += 1
            if gcd(D(det*G), D(det*H)) % gcd(Dv[v], Dv[w]) != 0: l1f += 1; print('L1 FAIL', N, v, w, flush=True)
    print(f'[{time.time()-t0:.0f}s] 3: N={N} rank={len(g)}', flush=True)
res['lemma1_torsion_rank3'] = dict(checks=l1, fails=l1f)
assert l1f == 0
(OUT / 'e_extra.json').write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str)); print('ALL (e) CHECKS PASSED')
