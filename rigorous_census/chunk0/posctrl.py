# Claude, 26.09.2026. Положительный контроль всей машины rig_dem.run_T.
# Ослабленная задача: S' = T ∪ −T, «решение» — z, при котором ∏_{λ∈T}(1+λz) и ∏_{λ∈T}(1−λz) — квадраты.
# (Для z0 = 1 и |λ| <= 60 таких троек нет — проверено; поэтому z0 = p/q, p,q <= 12.)
# Доказательство теоремы 7.1 использует только рациональность P_T(z) и P_T(−z), поэтому для неё
# верно то же: все такие z лежат в Z_A ∪ Z_B. Берём целые тройки T и z0, для которых z0 — такое решение,
# E_T ранга [1,1] с точкой, P_T(z0) бесконечного порядка; машина обязана выдать ±z0 среди кандидатов.
from sage.all import *
import sys, json, time, itertools, multiprocessing as mp
sys.path.insert(0, '/home/kep/magicKube/rigorous_census/chunk0')
from rig_dem import run_T, is_sq, model3, clean

NMAX = 60          # |λ| <= NMAX
NTEST = 40         # сколько троек ранга 1 проверить


def weak(T):
    return lambda z: is_sq(prod(1 + l * z for l in T)) and is_sq(prod(1 - l * z for l in T))


def sqf(n):
    n = ZZ(n)
    return (1 if n > 0 else -1) * n.abs().squarefree_part()


def pick():
    """тройки целых T (|λ| <= NMAX) и z = p/q (1 <= p,q <= 12), при которых q∏(q±λp) — квадраты"""
    out = []
    vals = [v for v in range(-NMAX, NMAX + 1) if v != 0]
    for q in range(1, 13):
        for p in range(1, 13):
            if gcd(p, q) != 1:
                continue
            cls = {}
            info = {}
            for l in vals:
                v = q * (q + l * p); w = q * (q - l * p)
                if v == 0 or w == 0:
                    continue
                info[l] = (sqf(v), sqf(w))
                cls.setdefault(info[l], []).append(l)
            ls = sorted(info)
            for i, a in enumerate(ls):
                for b in ls[i + 1:]:
                    # нужно sqf(v_a v_b v_c * q^3 ...) — считаем для q∏(q+λp)/q^2: класс v_a v_b v_c / q^2 ~ v_a v_b v_c
                    key = (sqf(info[a][0] * info[b][0]), sqf(info[a][1] * info[b][1]))
                    # q*prod(q+λp) = (q(q+ap))(q(q+bp))(q(q+cp)) / q^2 -> класс = класс произведения трёх v
                    for c in cls.get(key, []):
                        if c > b:
                            out.append(([a, b, c], QQ(p) / q))
    return out


def one(arg):
    T, z0 = arg
    Tq = [QQ(t) for t in T]
    E1, L, g = model3(Tq)
    y0 = prod(1 + t * z0 for t in Tq)
    if not y0.is_square() or not prod(1 - t * z0 for t in Tq).is_square():
        return None
    P1 = E1([L * z0, L * sqrt(y0)])       # P_T(z0): X = L z, Y = L y
    if P1.order() != oo:
        return None
    Em = E1.minimal_model()
    lo, hi = [int(v) for v in pari(Em).ellrank()[:2]]
    if (lo, hi) != (1, 1):
        return None
    S = sorted(set(Tq + [-t for t in Tq]))
    o = run_T(0, 0, T, S=S, okfun=weak(Tq))
    targets = (str(z0), str(-z0))
    found = any(z in targets for z in o.get('caseA_nondeg', []) + o.get('caseB_nondeg', []))
    o['z0'] = str(z0)
    o['posctrl_found_z0'] = bool(found)
    o['hits_z0'] = [h for h in o.get('caseA_hits', []) if h[0] in targets]
    o['hP_z0_sage'] = float(Em(E1.isomorphism_to(Em)(P1)).height())
    for k in ('cert', 'E1', 'Emin', 'caseB_roots_z'):
        o.pop(k, None)
    return o


if __name__ == '__main__':
    t0 = time.time()
    cands = pick()
    print('пар (T, z0) с обоими квадратами:', len(cands), flush=True)
    res = []
    with mp.get_context('fork').Pool(3) as pool:
        for o in pool.imap(one, cands[:3000]):
            if o is None:
                continue
            res.append(o)
            print(f"[{len(res)} {time.time()-t0:.0f}s] T={o['T']} z0={o['z0']} M0={o.get('M0_rig')} "
                  f"найдено ±z0: {o['posctrl_found_z0']} (z,n,t): {o['hits_z0']} caseB: {o.get('caseB_nondeg')} "
                  f"hP0={o.get('hP0_sage', 0):.3f} h(P_z0)={o['hP_z0_sage']:.3f}", flush=True)
            if len(res) >= NTEST:
                pool.terminate()
                break
    json.dump(clean(res), open('/home/kep/magicKube/rigorous_census/chunk0/posctrl.json', 'w'),
              ensure_ascii=False, indent=1)
    nf = sum(o['posctrl_found_z0'] for o in res)
    print(f'итог: проверено {len(res)}, ±z0 найдено в {nf}', flush=True)
