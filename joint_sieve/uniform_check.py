# Проверка леммы равномерности: для множителей ранга 1 с образующей P_T(0) (точка z=0) и кручением (ℤ/2)²
# таблицы решета при ℓ зависят только от k mod ℓ. Берём k = 204/247 и k' ≡ k (mod 5·7·11·13) с тем же порядком клеток
# (r > s/2), те же T-образцы (по меткам s, r, d = s−r, u = s+r), решето только по ℓ ∈ {5,7,11,13}: выжившие должны совпасть.
# Отрицательный контроль: наклон k'' ≢ k — выжившие, как правило, отличаются.
from jsieve import *
from run_slope import degenerate_points
import sys

LAB = ['s', '-s', 'r', '-r', 'd', '-d', 'u', '-u']


def labels(r, s):
    return {'s': s, '-s': -s, 'r': r, '-r': -r, 'd': s - r, '-d': r - s, 'u': s + r, '-u': -s - r}


def uniform_factor(r, s, pat):
    lab = labels(r, s); T = sorted(lab[p] for p in pat)
    f = Factor(T)
    if not (f.rank_proved and f.rank == 1 and sorted(f.tinv) == [2, 2]): return None
    P0 = f.point(0, 1); cls = f.decompose(P0)
    if abs(cls[0]) != 1: return None
    # равномерный базис: образующая P0 (E(ℚ) = ℤ P0 ⊕ (ℤ/2)²), кручение — t_λ для двух наибольших по метке λ ∈ T \ {c0}
    f.gens = [P0]
    others = [lam for lam in f.T if lam != f.c0]
    tg = []
    for lam in others[:2]:
        w = 1 / (QQ(-1) / lam + QQ(1) / f.c0); X = f.L * w; tg.append(f.E(X, 0))
    f.tgens = tg; f.tinv = [2, 2]
    f.tdict = {}
    for co in itertools.product(range(2), range(2)):
        P = f.E(0)
        for c, g in zip(co, tg): P = P + c * g
        f.tdict[P] = co
    assert len(f.tdict) == 4
    f.lcache = {}
    return f


def run(r, s, pats, primes, N=4):
    cells = slope_cells(r, s)
    F = []
    for pat in pats:
        f = uniform_factor(r, s, pat)
        if f is None: return None, pat
        F.append(f)
    JS = JointSieve(r, s, cells, F, verbose=False)
    # порядок клеток в cellset — по меткам, чтобы знаковые векторы сопоставлялись одинаково
    lab = labels(r, s); JS.cellset = [lab[p] for p in LAB]; JS.cidx = {lam: i for i, lam in enumerate(JS.cellset)}
    N0 = {i: 4 for i in range(len(F))}
    good = [l for l in primes if all(l not in f.bad for f in F)]
    cand, idx, cols, Nn, log = JS.run(list(range(len(F))), N0, good, lifts=[])
    S = set(map(tuple, cand.tolist()))
    return S, good


if __name__ == '__main__':
    r, s = 204, 247
    # T-образцы: тройки без пары ±λ, ранга 1 на 204/247 (подбираем автоматически)
    lab = labels(r, s)
    pats = []
    for pat in itertools.combinations(LAB, 3):
        if any(p.lstrip('-') == q.lstrip('-') for p, q in itertools.combinations(pat, 2)): continue
        f = uniform_factor(r, s, pat)
        if f is not None: pats.append(pat)
        if len(pats) >= 6: break
    print('образцы:', pats, flush=True)
    Fb = [uniform_factor(r, s, p) for p in pats]; badp = set().union(*[f.bad for f in Fb])
    primes = [int(l) for l in prime_range(5, 200) if l not in badp][:3]
    M = prod(primes)
    print('простые', primes, 'M =', M, flush=True)
    S0, good0 = run(r, s, pats, primes)
    print(f'{r}/{s}: простые {good0}, выживших {len(S0)}', flush=True)
    # ищем k' ≡ k (mod M), r' > s'/2, с теми же свойствами множителей
    found = []
    for s2 in range(50, 4000):
        for r2 in range(s2 // 2 + 1, s2):
            if gcd(r2, s2) != 1 or (r2, s2) == (r, s): continue
            if (r2 * s - r * s2) % M != 0: continue
            if any(x % l == 0 for x in (s2, r2, s2 - r2, s2 + r2) for l in primes): continue
            S2, good2 = run(r2, s2, pats, primes)
            if S2 is None: print(f'  {r2}/{s2}: образец {good2} не подходит (ранг/кручение/индекс)'); continue
            print(f'  {r2}/{s2} ≡ {r}/{s} (mod {M}): простые {good2}, выживших {len(S2)}, совпадение множеств: {S2 == S0}', flush=True)
            found.append((r2, s2, S2 == S0))
            if len(found) >= 3: break
        if len(found) >= 3: break
    # отрицательный контроль: несравнимые наклоны
    neg = 0
    for s2 in range(60, 400):
        for r2 in range(s2 // 2 + 1, s2):
            if gcd(r2, s2) != 1 or (r2 * s - r * s2) % primes[0] == 0: continue
            if any(x % l == 0 for x in (s2, r2, s2 - r2, s2 + r2) for l in primes): continue
            S2, good2 = run(r2, s2, pats, primes)
            if S2 is None: continue
            print(f'  контроль {r2}/{s2} ≢ {r}/{s} (mod {primes[0]}): выживших {len(S2)}, совпадение: {S2 == S0}', flush=True)
            neg += 1
            if neg >= 3: break
        if neg >= 3: break
